# Generic imports
import datetime
import json
import re
from collections import OrderedDict

from django.contrib.auth.models import User
from django.db.models import CharField, Count, F, Func, Prefetch, Value

# Local imports
import core.core_config
import core.models
import core.utils.commercial_kits
import core.utils.common
import core.utils.protocols
import wetlab.api.serializers


def project_table_fields(projects):
    """Gathers projects fields info to be displayed in the record project form

    Parameters
    ----------
    projects
        List of project names. Pe. ["None", "Relecov"]
    samples
        List of sample dict. Pe. [{"sample_name": "2133", "sample_project": "Relecov",...},
                                  {"sample_name": "2134", "sample_project": "Relecov",...}]
        "sample_project" key is mandatory

    Returns
    -------
        project_fields = [
            {
                "project": < SampleProject "Relecov" > , (obj)
                "project_fields": [list of projectField objects],
                "project_samples": [list of samples belonging to this project dict]
            },
            {
                "project": < SampleProject "Mepram" >,
                "project_fields": [list of projectField objects],
                "project_samples": [list of samples belonging to this project dict]
            }
        ]
    """
    projects_fields = []

    projects = core.models.SampleProjects.objects.prefetch_related(
        Prefetch(
            "sample_project_fields",
            queryset=core.models.SampleProjectsFields.objects.prefetch_related(
                "opt_value_prop"
            )
            .exclude(sample_project_field_used=False)
            .order_by("sample_project_field_order"),
            to_attr="project_fields_options",
        )
    ).filter(sample_project_name__in=projects)

    projects_fields = wetlab.api.serializers.ProjectsSerializer(
        projects, many=True
    ).data

    return projects_fields


def sample_table_fields(app_name):
    """Gathers fields info to include in record sample form. It grabs verbose_name from
        Samples model when available and options for dropdown lists.

    Parameters
    ----------
    app_name
        app name. p.e wetlab or drylab

    Returns
    -------
        field_info
        {
            "fields": {
                "sample_name": "Sample Name",
                "sample_collection_date": "Sample collection date"
            },
            "species": ["human", "rat"],
            "lab_request": ["isciii", "hrycl"],
            "sample_type": ["fastq", "blood"],
            "sample_project": ["None","Relecov"],
        }
    """
    # get the choices to be included in the form
    fields_info = {}
    fields_info["fields"] = {}
    for field in core.models.Samples._meta.get_fields():
        try:
            fields_info["fields"][field.name] = field.verbose_name
        except Exception:
            fields_info["fields"][field.name] = field.name

    fields_info["species"] = get_species()
    fields_info["lab_request"] = get_lab_requested()
    fields_info["sample_type"] = list(
        core.models.SampleType.objects.filter(apps_name=app_name).values_list(
            "sample_type", flat=True
        )
    )
    fields_info["sample_project"] = get_defined_sample_projects(app_name)
    fields_info["sample_project"].insert(0, "None")

    return fields_info


def save_recorded_samples(samples_data, req_user, app_name):
    """This function saves samples using record_samples form/view

    Parameters
    ----------
    samples_data
        accepts a list of dicts obtained from the jspreadsheet in the form
        eg.
    req_user
        user recording the samples
    app_name
        app_name. (wetlab, drylab, core, etc.)

    Returns
    -------
        sample sample_data as a list of dicts, including fields like:
        - unique_sample_id
        - patient_core obj if exists
        - sample project obj if exists
        - sample_state
        - complete_date - if only recorded and sample_state is completed
        - success: True / False if recording failed
        - error: if success False the actual error
    """
    for sample in samples_data:
        # Fill fields
        sample["user"] = req_user
        sample["app_name"] = app_name

        sample["sample_code_id"] = req_user + "_" + str(sample["sample_name"])

        # Set unique ID
        if not core.models.Samples.objects.exclude(
            unique_sample_id__isnull=True
        ).exists():
            sample["unique_sample_id"] = "AAA-0001"
        else:
            last_unique_value = (
                core.models.Samples.objects.exclude(unique_sample_id__isnull=True)
                .last()
                .unique_sample_id
            )
            sample["unique_sample_id"] = increase_unique_value(last_unique_value)

        # Check if patient code  already exists on database,
        # If not if will be created giving a sequencial dummy value
        if sample["patient_core"] != "":
            patient_obj = check_patient_code_exists(sample["patient_core"])
            if patient_obj is False:
                # Define the new patient only Patient code is defined
                patient_obj = create_empty_patient(sample["patient_core"])
        else:
            patient_obj = None

        sample["patient_core"] = patient_obj

        # Check if sample project exist and generate de appropiate objetct
        if sample["sample_project"] == "None" or sample["sample_project"] == "":
            sample["sample_project"] = None
        else:
            sample["sample_project"] = core.models.SampleProjects.objects.get(
                sample_project_name__exact=sample["sample_project"]
            )

        # If only recorded set sample to completed state
        if sample["only_recorded"] and sample["sample_project"] is None:
            sample["sample_state"] = "Completed"
            sample["completed_date"] = datetime.datetime.now()
        # If no sample project data needed set to defined
        elif sample["sample_project"] is None:
            sample["sample_state"] = "Defined"
        else:
            sample["sample_state"] = "Pre-Defined"

        try:
            core.models.Samples.objects.create_sample(sample)
            sample["success"] = True
        except Exception as e:
            sample["success"] = False
            sample["error"] = e

    return samples_data


def validate_sample_data(
    sample_data: json,
    req_user: str,
    app_name: str,
    repeat_allowed: bool = False,
    allow_user__repeat: bool = False,
) -> tuple:
    """sample data validation

    Args:
        sample_data (json): sample data formatted in json obtained from jspreadsheet
        req_user (str): requested user name
        app_name (str): application name (wetlab, drylab, core, etc.)
        repeat_allowed (bool, optional): allow or not that sample can be repeated. Defaults to False.
        allow_user__repeat (bool, optional): if allowed to be repeated check if same request user is allowed have repeated sample names. Defaults to False.

    Returns:
        tuple: returns result as boolean and a dictionary with the following keys
        validation
            list with validation info for each sample with format:
            [{"Sample Name": "test 01",
              "Validate": True,
              "Validation error": None
              },
              {"Sample Name": "test 02",
              "Validate": False,
              "Validation error": "Mandatory field is missing."
              }]
    """
    validation = []
    sample_name_list = []
    line = 0
    result = True
    not_allowed_sample_names = {}
    # collect the sample names already in the database in case that the user
    #  is not allowed to repeat
    if not repeat_allowed:
        if not allow_user__repeat:
            if core.models.Samples.objects.filter(
                sample_user__username__exact=req_user
            ).exists():
                existing_sample_name_list = list(
                    core.models.Samples.objects.filter(
                        sample_user__username__iexact=req_user
                    ).values_list("sample_name", flat=True)
                )
        else:
            existing_sample_name_list = list(
                core.models.Samples.objects.all().values_list("sample_name", flat=True)
            )
        # convert to dict to speed up the search
        not_allowed_sample_names = dict.fromkeys(existing_sample_name_list, 0)
    for sample in sample_data:
        line += 1
        sample_dict = {}

        # Check if sample name is not empty
        if sample["sample_name"] == "":
            error_cause = core.core_config.ERROR_EMPTY_SAMPLE_NAME.copy()
            error_cause.insert(1, str(line))
            sample_dict["Sample name"] = "Empty_" + str(line)
            sample_dict["Validate"] = False
            sample_dict["Validation error"] = " ".join(error_cause)
            result = False
            validation.append(sample_dict)
            continue

        if repeat_allowed or (
            sample["sample_name"] not in sample_name_list
            and sample["sample_name"] not in not_allowed_sample_names
        ):
            sample_name_list.append(sample["sample_name"])
        else:
            error_cause = core.core_config.ERROR_REPEATED_SAMPLE_NAME.copy()
            error_cause.insert(1, str(line))
            sample_dict["Sample name"] = (
                sample["sample_name"] + "_repeated_" + str(line)
            )
            sample_dict["Validate"] = False
            sample_dict["Validation error"] = " ".join(error_cause)
            result = False
            validation.append(sample_dict)
            continue

        # Create not error dictionary
        sample_dict["Sample name"] = sample["sample_name"]
        sample_dict["Validate"] = True
        sample_dict["Validation error"] = []

        # Check only recorded format
        if (
            sample["only_recorded"] != ""
            and sample["only_recorded"] is not True
            and sample["only_recorded"] is not False
        ):
            sample_dict["Validate"] = False
            sample_dict["Validation error"].append(
                "".join(core.core_config.ERROR_ONLY_RECORDED_FIELD)
            )

        # Check if sample already in the DB
        if core.models.Samples.objects.filter(
            sample_name__exact=sample["sample_name"],
            sample_user__username__exact=req_user,
        ).exists():
            error_cause = core.core_config.ERROR_SAMPLE_ALREADY_DEFINED.copy()
            error_cause.insert(1, sample["sample_name"])
            sample_dict["Validate"] = False
            sample_dict["Validation error"].append(" ".join(error_cause))

        # Check sample type filled (mandatory), check if its in the DB, check sample type's mandatory fields
        if sample["sample_type"] == "":
            sample_dict["Validate"] = False
            sample_dict["Validation error"].append(
                "".join(core.core_config.ERROR_EMPTY_SAMPLE_TYPE)
            )
        elif defined_sample_type(sample["sample_type"], app_name):
            sample_dict["Validate"] = False
            sample_dict["Validation error"].append(
                defined_sample_type(sample["sample_type"], app_name)
            )
        else:
            # Check mandatory fields
            if check_mandatory_fields_included(sample, sample["sample_type"], app_name):
                sample_dict["Validate"] = False
                sample_dict["Validation error"].append(
                    check_mandatory_fields_included(
                        sample, sample["sample_type"], app_name
                    )
                )

        # Check if laboratory is in the DB
        if sample["lab_request"] != "" and defined_lab_request(
            sample["lab_request"], app_name
        ):
            sample_dict["Validate"] = False
            sample_dict["Validation error"].append(
                defined_lab_request(sample["lab_request"], app_name)
            )

        # Check if species is in the DB
        if sample["species"] != "" and defined_species(sample["species"], app_name):
            sample_dict["Validate"] = False
            sample_dict["Validation error"].append(
                defined_species(sample["species"], app_name)
            )

        if not sample_dict["Validation error"]:
            sample_dict["Validation error"] = ""
        else:
            sample_dict["Validation error"] = ". ".join(
                map(str, sample_dict["Validation error"])
            )
            result = False
        validation.append(sample_dict)

    return validation, result


def validate_project_data(project_data, project_name, sample_validation=False):
    """Sample data validation

    Parameters
    ----------
    project_data
        project data formatted in json obtained from jspreadsheet

    Returns
    -------
        validation
            list with validation info for each sample with format:
            [{"Sample Name": "test 01",
              "Validate": True,
              "Validation error": None
              },
              {"Sample Name": "test 02",
              "Validate": False,
              "Validation error": "Mandatory field is missing."
              }]
    """
    validation = []
    result = True
    for sample in project_data:
        sample_dict = {}
        if sample_validation:
            for dict in sample_validation:
                if dict["Sample name"] == sample["sample_name"]:
                    sample_dict["Sample name"] = dict["Sample name"]
                    sample_dict["Project name"] = project_name
                    sample_dict["Validate"] = dict["Validate"]
                    sample_dict["Validation error"] = []
                    if dict["Validation error"] != "":
                        sample_dict["Validation error"].append(dict["Validation error"])
        else:
            sample_dict["Sample name"] = sample["sample_name"]
            sample_dict["Project name"] = project_name
            sample_dict["Validate"] = True
            sample_dict["Validation error"] = []

        s_proj_obj = core.models.SampleProjects.objects.filter(
            sample_project_name__iexact=project_name
        ).last()
        s_pro_field_objs = core.models.SampleProjectsFields.objects.filter(
            sample_projects_id=s_proj_obj
        )
        for s_pro_field_obj in s_pro_field_objs:
            field_name = s_pro_field_obj.get_field_name()
            field_type = s_pro_field_obj.get_field_type()
            field_options = s_pro_field_obj.get_field_options_list()
            """
            if field_type == "String" and type(sample[field_name]) is not str:
                error_cause = core.core_config.ERROR_PROJECT_FIELD_NOTSTRING.copy()
                error_cause.insert(1, field_name)
                sample_dict["Validate"] = False
                sample_dict["Validation error"].append(" ".join(error_cause))
            """
            if field_type == "Date" and sample[field_name] != "":
                try:
                    datetime.datetime.strptime(sample[field_name], "%Y-%m-%d")
                except Exception:
                    # Unknown string format: string present at position x
                    error_cause = core.core_config.ERROR_PROJECT_FIELD_NODATE.copy()
                    error_cause.insert(1, field_name)
                    sample_dict["Validate"] = False
                    sample_dict["Validation error"].append(" ".join(error_cause))
            elif field_type == "Options List" and sample[field_name] != "":
                if sample[field_name] not in field_options:
                    error_cause = core.core_config.ERROR_PROJECT_FIELD_NOOPTION.copy()
                    error_cause.insert(1, field_name)
                    error_cause.insert(3, ", ".join(field_options))
                    sample_dict["Validate"] = False
                    sample_dict["Validation error"].append(" ".join(error_cause))
        if not sample_dict["Validation error"]:
            sample_dict["Validation error"] = ""
        else:
            sample_dict["Validation error"] = ". ".join(
                map(str, sample_dict["Validation error"])
            )
            result = False
        validation.append(sample_dict)

    return validation, result


def save_project_data(excel_data, project_info):
    """Saves the project form data for each sample

    Parameters
    ----------
    excel_data
        excel form data in json format (list of dicts)
        e.g.
    project_info
        project info according to project_table_fields function.

    Returns
    -------
        Same project info adding keys
         - success: True/False
         - error: error message
    """
    for sample in excel_data:
        sample_id = core.models.Samples.objects.get(
            sample_code_id__exact=sample["sample_code_id"]
        )

        for field in project_info["sample_project_fields"]:
            field_value = {}
            field_value["sample_id"] = sample_id
            field_value["sample_project_field_id"] = (
                core.models.SampleProjectsFields.objects.get(
                    sample_projects_id__exact=core.models.SampleProjects.objects.get(
                        sample_project_name__exact=project_info["sample_project_name"]
                    ),
                    sample_project_field_name__exact=field["sample_project_field_name"],
                )
            )
            field_value["sample_project_field_value"] = sample[
                field["sample_project_field_name"]
            ]
            try:
                core.models.SampleProjectsFieldsValue.objects.create_project_field_value(
                    field_value
                )
                project_info["success"] = True
            except Exception:
                project_info["success"] = False
                project_info["error"] = "Error saving any of the project fields"
        if project_info["success"] and sample["only_recorded"] == "True":
            sample_id.set_state("Completed")
        elif project_info["success"] and sample["only_recorded"] == "False":
            sample_id.set_state("Defined")
        elif project_info["success"] and sample["only_recorded"]:
            sample_id.set_state("Completed")
        elif project_info["success"]:
            sample_id.set_state("Defined")
        else:
            sample_id.set_state("Pre-Defined")
    return project_info


def check_empty_fields(row_data, optional_index):
    """
    Description:
        The function check if row_data contains empty values. If a empty field is defined as optional
        it will be ignored.
    Input:
        row_data:       # data to be checked
        optional_index :   # list with the index values that can be empty
    Return:
        False is all mandatory fields contain data. True if any of them are empty
    """
    for i in range(len(row_data)):
        if row_data[i] == "" and i not in optional_index:
            return True
    return False


def create_empty_patient(p_code_id):
    """
    Description:
        The function create patient in database.
    Input:
        name:       patient name
        surname     patient family name
    Return:
        patient_obj
    """
    patient_data = {}
    patient_data["patientSex"] = "Not Provided"
    patient_data["patientName"] = "Not Provided"
    patient_data["patientSurname"] = "Not Provided"
    patient_data["patientCode"] = p_code_id

    patient_obj = core.models.PatientCore.objects.create_patient(patient_data)

    return patient_obj


def add_molecule_protocol_parameters(data, parameters):
    """Store molecule parameter values asigned to the molecule

    Args:
        data (list): _description_
        parameters (list): _description_
    """
    for row in data:
        molecule_obj = get_molecule_obj_from_id(row["m_ids"])
        prot_obj = molecule_obj.get_protocol_obj()
        molecule_obj = molecule_obj.set_user_lot_kit(
            row["Lot Commercial Kit"], update_usage_kit=True
        )
        for param in parameters:
            molecule_parameter_value = {}
            molecule_parameter_value["molecule_parameter_id"] = (
                core.models.ProtocolParameters.objects.filter(
                    protocol_id=prot_obj, parameter_name__iexact=param
                ).last()
            )
            molecule_parameter_value["molecule_id"] = molecule_obj
            molecule_parameter_value["parameter_value"] = row[param]
            _ = core.models.MoleculeParameterValue.objects.create_molecule_parameter_value(
                molecule_parameter_value
            )

        molecule_obj.set_state("Completed")
        # Update sample state
        sample_obj = molecule_obj.get_sample_obj()
        sample_obj.set_state("Pending for use")
    return


def check_if_molecule_use_defined(app_name):
    """
    Description:    The function check if there are defined the use for molecules

    Input:
        app_name    # application name to assign the right molecule use
    Return:
        True or False #
    """
    if core.models.MoleculeUsedFor.objects.filter(apps_name__exact=app_name).exists():
        return True
    return False


def check_if_sample_project_id_exists(sample_project_id):
    if core.models.SampleProjects.objects.filter(pk__exact=sample_project_id).exists():
        return True
    return False


def check_mandatory_fields_included(data, sample_type, app_name):
    """Check if for the type of sample all mandatory fields are not empty

    Args:
        data (dictionary): contains data recorded by user
        sample_type (string): Type of sample
        app_name (_type_): application name

    Returns:
        _type_: _description_
    """
    missing_mandatory = []
    mandatory_fields = (
        core.models.SampleType.objects.filter(
            sample_type__exact=sample_type, apps_name__exact=app_name
        )
        .last()
        .get_mandatory_values()
    )

    for field in mandatory_fields:
        if field == "only_recorded":
            if not data[field]:
                missing_mandatory.append(field)
        else:
            if data[field] == "":
                missing_mandatory.append(field)

        if len(missing_mandatory) != 0:
            missing_fields = []
            for field in core.models.Samples._meta.get_fields():
                if field.name in missing_mandatory:
                    try:
                        missing_fields.append(field.verbose_name)
                    except Exception:
                        missing_fields.append(field.name)
            error_cause = core.core_config.ERROR_MISSING_MANDATORY.copy()
            error_cause.insert(1, ", ".join(missing_fields))
            return " ".join(error_cause)


def check_patient_code_exists(p_code_id):
    """
    Description:
        The function check if patient name/surname is defined in database.
    Input:
        name:       patient name
        surname     patient family name
    Return:
        False is user is not define or patient_obj is patient exists
    """
    if core.models.PatientCore.objects.filter(patient_code__iexact=p_code_id).exists():
        patient_obj = core.models.PatientCore.objects.filter(
            patient_code__iexact=p_code_id
        ).last()
    else:
        return False

    return patient_obj


def create_new_sample_project(form_data, app_name):
    """
    Description:
        The function create sample project in database.
    Input:
        form_data:  Information collected in the user form
        app_name    application name
    Return:
        patient_obj
    """
    s_project_data = {}
    form_fields = [
        "sampleProyectName",
        "sampleProyectManager",
        "sampleProyectManagerContact",
        "description",
    ]
    db_fields = [
        "sampleProjectName",
        "sampleProjectManager",
        "sampleProjectContact",
        "sampleProjectDescription",
    ]
    for i in range(len(form_fields)):
        s_project_data[db_fields[i]] = form_data[form_fields[i]]
    s_project_data["apps_name"] = app_name
    new_sample_project = core.models.SampleProjects.objects.create_sample_project(
        s_project_data
    )

    new_sample_project_id = new_sample_project.get_id()

    return new_sample_project_id


def create_table_molecule_pending_use(sample_list, app_name):
    """
    Description:
        The function get the type of use that the molecule can have to assign it.
    Input:
        sample_list:  list of samples that are pending the use that the molecuel will have
        app_name : application name
    Return:
        use_type
    """
    use_type = {}
    use_type["data"] = list(
        core.models.MoleculePreparation.objects.filter(
            molecule_used_for=None, sample__in=sample_list
        ).values_list("sample__sample_name", "molecule_code_id", "pk")
    )
    if len(use_type["data"]) > 0:
        if core.models.MoleculeUsedFor.objects.filter(
            apps_name__exact=app_name
        ).exists():
            use_type["types"] = list(
                core.models.MoleculeUsedFor.objects.filter(
                    apps_name__exact=app_name
                ).values_list("used_for", flat=True)
            )

        use_type["heading"] = core.core_config.HEADING_FOR_SELECTING_MOLECULE_USE
    return use_type


def create_table_to_select_molecules(samples_list):
    """
    Description:
        The function return a dictionary with the information to display to user to
        select the molecules.

    Parameters:
    -----------
        samples_list  : sample list object to get the information

    Return:
    -------
        sample_information.
    """
    sample_information = {}
    sample_information["sample_information"] = []
    sample_code_id = []
    for sample in samples_list:
        sample_information["sample_information"].append(
            sample.get_info_in_defined_state()
        )
        sample_code_id.append(sample.get_sample_code())
    sample_information["sample_heading"] = core.core_config.HEADING_FOR_DEFINED_SAMPLES
    # sample_information["sample_code_ids"] = ",".join(sample_code_id)

    return sample_information


def create_table_pending_molecules(molecule_list):
    """The function prepare molecule information to display.

    Parameters:
    -----------
        _owner_molecules:  list of molecules belongs to user

    Return:
    -------
        molecule_data
    """
    molecule_data = {}
    molecule_data["data"] = list(
        molecule_list.values_list(
            "sample__sample_name", "molecule_code_id", "protocol_used__name"
        )
        .annotate(
            extrac_date=Func(
                F("molecule_extraction_date"),
                Value("%Y-%m-%d"),
                function="DATE_FORMAT",
                output_field=CharField(),
            )
        )
        .annotate(pk=F("pk"))
    )
    if len(molecule_data["data"]) > 0:
        molecule_data["molecule_heading"] = (
            core.core_config.HEADING_FOR_PENDING_MOLECULES
        )

    return molecule_data


def define_table_for_sample_project_fields(sample_project_id):
    """
    Description:
        The function return a dictionary with the information to create the table
        for defining the fields used in the sample project
    Input:
        sample_project_id # id  to get sample project information
    Return:
        sample_project_data
    """
    sample_project_data = {}
    sample_project_obj = core.models.SampleProjects.objects.get(
        pk__exact=sample_project_id
    )

    sample_project_data["sample_project_name"] = (
        sample_project_obj.get_sample_project_name()
    )
    sample_project_data["sample_project_id"] = sample_project_id
    sample_project_data["heading"] = core.core_config.HEADING_FOR_SAMPLE_PROJECT_FIELDS
    return sample_project_data


def display_molecule_use(app_name):
    """
    Description:    The function collect the defined molecule use

    Input:
        app_name    # application name to assign the right molecule use
    Return:
        molecule_use_data #
    """
    molecule_use_data = {}
    molecule_use_data["defined_molecule_use"] = []
    if core.models.MoleculeUsedFor.objects.filter(apps_name__exact=app_name).exists():
        molecule_uses = core.models.MoleculeUsedFor.objects.filter(
            apps_name__exact=app_name
        )
        for molecule in molecule_uses:
            massive = molecule.get_massive()
            if massive == "True":
                molecule_use_data["defined_molecule_use"].append(
                    [molecule.get_molecule_use_name(), "YES"]
                )
            else:
                molecule_use_data["defined_molecule_use"].append(
                    [molecule.get_molecule_use_name(), "NO"]
                )
    return molecule_use_data


def display_sample_types(app_name):
    """
    Description:
        The function return a dictionary with the information to define the type of sample
    Input:
        app_name # application name where are sample type are defined
    Return:
        sample_types
    """
    sample_types = {}
    defined_sample_types = []

    if core.models.SampleType.objects.filter(apps_name__exact=app_name).exists():
        s_types = core.models.SampleType.objects.filter(
            apps_name__exact=app_name
        ).order_by("sample_type")
        for s_type in s_types:
            defined_sample_types.append([s_type.get_sample_type_id, s_type.get_name])
        sample_types["defined_sample_types"] = defined_sample_types

    mandatory_value = {}
    for field in core.models.Samples._meta.get_fields():
        try:
            mandatory_value[field.name] = field.verbose_name
        except Exception:
            mandatory_value[field.name] = field.name

    sample_types["mandatory_value"] = mandatory_value
    return sample_types


def get_all_sample_information(sample_id, join_values=False):
    sample_information = {}
    sample_obj = get_sample_obj_from_id(sample_id)
    if not core.models.Samples.objects.filter(pk__exact=sample_id).exists():
        return {"Error": " ".join(core.core_config.ERROR_SAMPLE_NOT_FOUND)}
    sample_information["sample_id"] = sample_id
    sample_information["sample_name"] = sample_obj.get_sample_name()

    sample_information["sample_definition"] = sample_obj.get_info_for_display()
    sample_information["sample_definition_heading"] = (
        core.core_config.HEADING_FOR_SAMPLE_DEFINITION
    )
    if join_values:
        sample_information["sample_definition_join_value"] = list(
            zip(
                sample_information["sample_definition_heading"],
                sample_information["sample_definition"],
            )
        )
    # get the sample project information fields
    sample_project_obj = sample_obj.get_sample_project_obj()
    if sample_project_obj is not None:
        sample_information.update(
            get_sample_project_information(
                sample_project_obj, sample_obj, join_values=True
            )
        )

    # check if molecule information exists for the sample
    if core.models.MoleculePreparation.objects.filter(sample=sample_obj).exists():
        molecules = core.models.MoleculePreparation.objects.filter(sample=sample_obj)
        sample_information["molecule_definition_heading"] = (
            core.core_config.HEADING_FOR_MOLECULE_DEFINITION
        )
        sample_information["molecule_definition"] = []
        sample_information["molecule_parameter_values"] = []
        sample_information["molecule_definition_data"] = []

        for molecule in molecules:
            molecule_definition_data = []
            molecule_definition_data.append(molecule.get_info_for_display())
            protocol_used_obj = molecule.get_protocol_obj()
            if core.models.ProtocolParameters.objects.filter(
                protocol_id=protocol_used_obj
            ).exists():
                parameter_names = core.models.ProtocolParameters.objects.filter(
                    protocol_id=protocol_used_obj
                ).order_by("parameter_order")
                molecule_param_heading = ["Molecule CodeID"]
                mol_param_value = [molecule.get_molecule_code_id()]
                for p_name in parameter_names:
                    molecule_param_heading.append(p_name.get_parameter_name())
                    if core.models.MoleculeParameterValue.objects.filter(
                        molecule_id=molecule
                    ).exists():
                        try:
                            mol_param_value.append(
                                core.models.MoleculeParameterValue.objects.get(
                                    molecule_id=molecule, molecule_parameter_id=p_name
                                ).get_param_value()
                            )
                        except core.models.MoleculeParameterValue.DoesNotExist:
                            # if the parameter was not set at the time the molecule was handeled
                            mol_param_value.append("")

                molecule_definition_data.append(molecule_param_heading)
                molecule_definition_data.append(mol_param_value)
            else:
                # Add empty values in case that there is no paramters defined yet
                molecule_definition_data.append("")
                molecule_definition_data.append("")
            sample_information["molecule_definition_data"].append(
                molecule_definition_data
            )

    return sample_information


def get_defined_samples(register_user):
    """
    Description:
        The function will return the samples created by the user which are in Defined state.
    Input:
        register_user
    Variables:
        laboratories # list containing all laboratory names
    Return:
        samples_list.
    """
    defined_samples = {}
    # sample_ids = []
    defined_samples["heading"] = core.core_config.HEADING_FOR_DEFINED_SAMPLES_STATE
    sample_information = []
    if core.models.Samples.objects.filter(
        sample_state__sample_state_name__exact="Defined",
        sample_user__username__exact=register_user,
    ).exists():
        sample_list = core.models.Samples.objects.filter(
            sample_state__sample_state_name__exact="Defined",
            sample_user__username__exact=register_user,
        ).order_by("generated_at")
        for sample in sample_list:
            sample_information.append(sample.get_info_in_defined_state())
        defined_samples["sample_information"] = sample_information
    return defined_samples


def get_info_to_display_sample_project(sample_project_id):
    """
    Description:
        The function return the information for the requested sample project
    Return:
        info_s_project.
    """
    info_s_project = {}
    if core.models.SampleProjects.objects.filter(pk__exact=sample_project_id).exists():
        sample_project_obj = core.models.SampleProjects.objects.get(
            pk__exact=sample_project_id
        )
        # collect data from project
        info_s_project["sample_project_id"] = sample_project_id
        info_s_project["sample_project_name"] = (
            sample_project_obj.get_sample_project_name()
        )
        info_s_project["main_data"] = list(
            zip(
                core.core_config.SAMPLE_PROJECT_MAIN_DATA,
                sample_project_obj.get_full_info_to_display(),
            )
        )

        if core.models.SampleProjectsFields.objects.filter(
            sample_projects_id=sample_project_obj
        ).exists():
            sample_project_fields = core.models.SampleProjectsFields.objects.filter(
                sample_projects_id=sample_project_obj
            ).order_by("sample_project_field_order")
            s_project_fields_list = []
            for sample_project_field in sample_project_fields:
                s_project_fields_list.append(
                    sample_project_field.get_sample_project_fields_name(
                        include_search=True
                    )
                )
            info_s_project["fields"] = s_project_fields_list
        info_s_project["heading"] = core.core_config.HEADING_FOR_SAMPLE_PROJECT_FIELDS
    else:
        info_s_project["ERROR"] = "Sample project does not exists"
        return info_s_project

    return info_s_project


def get_lab_requested():
    """
    Description:
        The function will return the Sample origin places defined in database.
    Variables:
        lab_requested_places # list containing the place names
    Return:
        lab_requested_places.
    """
    lab_requested_places = []
    if core.models.LabRequest.objects.filter().exists():
        lab_requesteds = core.models.LabRequest.objects.all()

        for lab_requested in lab_requesteds:
            lab_requested_places.append(lab_requested.get_lab_request_code())
    return lab_requested_places


def get_only_recorded_samples_and_dates():
    """
    Description:
        The function api return a list of list of samples which are defined as only recorded,
        project name, sample type, species, recorded date, and sample id
    Return:
        samples_data
    """
    samples_data = []
    if core.models.Samples.objects.filter(only_recorded=True).exists():
        sample_objs = core.models.Samples.objects.filter(only_recorded=True).order_by(
            "generated_at"
        )
        for sample_obj in sample_objs:
            data = [sample_obj.get_sample_name()]
            data.append(sample_obj.get_sample_project())
            data.append(sample_obj.get_sample_type())
            data.append(sample_obj.get_species())
            data.append(sample_obj.get_entry_date())
            data.append(sample_obj.get_sample_id())
            samples_data.append(data)

    return samples_data


def get_parameters_sample_project(sample_project_id):
    """
    Description:
        The function return the parameters definded for the sample project id.
    Input:
        sample_project_id       # id of the sample project
    Return:
        parameters_s_project
    """
    parameters_s_project = {}
    if core.models.SampleProjects.objects.filter(pk__exact=sample_project_id).exists():
        sample_project_obj = core.models.SampleProjects.objects.get(
            pk__exact=sample_project_id
        )
        if core.models.SampleProjectsFields.objects.filter(
            sample_projects_id=sample_project_obj
        ).exists():
            sample_project_fields = core.models.SampleProjectsFields.objects.filter(
                sample_projects_id=sample_project_obj
            ).order_by("sample_project_field_order")
            s_project_fields_list = []
            # parameter_ids = []
            # parameter_names = []
            for sample_project_field in sample_project_fields:
                parameter_data = sample_project_field.get_sample_project_fields_name(
                    include_search=True
                )
                # parameter_names.append(parameter_data[0])
                parameter_data.append(sample_project_field.get_field_id())
                # add empty field to include the new name
                parameter_data.insert(1, "")
                s_project_fields_list.append(parameter_data)
            parameters_s_project["fields"] = s_project_fields_list
        parameters_s_project["heading"] = (
            core.core_config.HEADING_FOR_MODIFY_SAMPLE_PROJECT_FIELDS
        )
        parameters_s_project["sample_project_id"] = sample_project_id
        parameters_s_project["sample_project_name"] = (
            sample_project_obj.get_sample_project_name()
        )
        # parameters_s_project["parameter_names"] = ",".join(parameter_names)
        # parameters_s_project["parameter_ids"] = ",".join(parameter_ids)
    else:
        return "ERROR"

    return parameters_s_project


def get_info_for_defined_sample_projects(app_name):
    """
    Description:
        The function return a list with all defined sample projects.
    Return:
        info_s_projects.
    """
    info_s_projects = []
    if core.models.SampleProjects.objects.filter(apps_name__exact=app_name).exists():
        s_projects = core.models.SampleProjects.objects.filter(
            apps_name__exact=app_name
        )
        for s_project in s_projects:
            s_project_data = s_project.get_info_to_display()
            if core.models.SampleProjectsFields.objects.filter(
                sample_projects_id=s_project
            ).exists():
                s_project_data.append(True)
            else:
                s_project_data.append(False)
            info_s_projects.append(s_project_data)
    return info_s_projects


def get_defined_sample_projects(app_name):
    """
    Description:
        The function will return the samples projects that a sample could has.
    Variables:
        sample_projects # list containing the sample projects defined in database
    Return:
        sample_projects.
    """
    sample_projects = []
    if core.models.SampleProjects.objects.filter(apps_name__exact=app_name).exists():
        s_projects = core.models.SampleProjects.objects.filter(
            apps_name__exact=app_name
        )
        for s_project in s_projects:
            sample_projects.append(s_project.get_sample_project_name())
    return sample_projects


def get_molecule_codeid_from_object(molecule_obj):
    """
    Description:
        The function will return the molecule id form the object class.
    Input:
        molecule_obj : molecule object from where to get the molecule id
    Return:
        molecules_id.
    """
    return molecule_obj.get_molecule_code_id()


def get_molecule_obj_from_id(molecule_id):
    """
    Description:
        The function will return the molecule object that are assigned to the id.
        It returns '' if molecule does not have the id
    Input:
        molecule_id : id number from where to get the molecule
    Return:
        molecules_obj.
    """
    if core.models.MoleculePreparation.objects.filter(pk__exact=molecule_id).exists():
        molecules_obj = core.models.MoleculePreparation.objects.get(
            pk__exact=molecule_id
        )
        return molecules_obj
    else:
        return ""


def get_molecule_objs_from_sample(sample_obj):
    """
    Description:
        The function will return the molecule object that are assigned to the sample.
        It returns '' if no molecule is assigned yet to the sample
    Input:
        sample_obj : sample object from where to get the molecule
    Return:
        molecules_obj.
    """
    if core.models.MoleculePreparation.objects.filter(sample=sample_obj).exists():
        molecules_obj = core.models.MoleculePreparation.objects.filter(
            sample=sample_obj
        )
        return molecules_obj
    else:
        return ""


def get_molecule_objs_in_state(m_state, user=None, friend_list=None):
    """The function returns a list of instance molecules which match the conditions
        molecule instances are filter by the user and the user friend list.
        If no user is set then no filter is done

    Parameters
    ----------
    m_state : string
        state name to get the molecules
    user : user object, optional
        instance of user to filter for their molecules, by default None
    friend_list : Bool
        Boolean value to allow for getting molecule inside the user friend list

    Returns
    -------
    molecule_objs : list of Molecules instances
        contains the objects that were matched on the condition. Empty list if
        not match
    """
    if user:
        if friend_list:
            user_list = core.utils.common.get_friend_list(user)
        else:
            user_list = [user]
        return core.models.MoleculePreparation.objects.filter(
            state__molecule_state_name__exact=m_state,
            molecule_user__in=user_list,
        )
    else:
        return (
            core.models.MoleculePreparation.objects.filter(
                state__molecule_state_name__exact=m_state
            )
            .order_by("molecule_user")
            .order_by("generated_at")
        )


def get_modules_type():
    """
    Description:
        The function will return the list of type of molecules defined.
    Input:
        none
    Variables:
        molecule_type # list containing all type of molecules
    Return:
        molecule_type.
    """
    molecule_type = []
    types = core.models.MoleculeType.objects.all()
    for type in types:
        molecule_type.append(type.get_name())
    return molecule_type


def get_molecule_protocols(apps_name):
    """
    Description:
        The function will return the protocols defined.
    Input:
        apps_name
    Variables:
        protocols # dictionary containing all protocols using "key" as protocol type
    Return:
        protocols.
    """
    protocol_types = []
    protocol_list = []
    protocols = {}
    p_types = core.models.ProtocolType.objects.filter(
        molecule__isnull=False, apps_name__exact=apps_name
    )
    molecule_types = core.models.MoleculeType.objects.filter()
    for molecule in molecule_types:
        protocols[molecule.get_name()] = []
    for p_type in p_types:
        protocol_types.append(p_type.get_name())
    for protocol_type in p_types:
        protocols_in_type = core.models.Protocols.objects.filter(type=protocol_type)
        for p_in_type in protocols_in_type:
            protocol_name = p_in_type.get_name()
            protocols[protocol_type.get_molecule_type()].append(protocol_name)
            protocol_list.append(protocol_name)

    return protocols, protocol_list


def get_molecule_data_and_protocol_parameters(protocol_objs):
    mol_data_parm = {}
    for protocol_obj, mol_ids in protocol_objs.items():
        prot_name = protocol_obj.get_name()
        mol_data_parm[prot_name] = {}
        mol_data_parm[prot_name]["params_type"] = (
            core.utils.protocols.get_protocol_parameters_and_type(protocol_obj)
        )
        mol_data_parm[prot_name][
            "fix_heading"
        ] = core.core_config.HEADING_FOR_MOLECULE_ADDING_PARAMETERS
        mol_data_parm[prot_name]["lot_kit"] = (
            core.utils.commercial_kits.get_lot_commercial_kits(protocol_obj)
        )
        mol_data_parm[prot_name]["param_heading"] = []
        prot_params = core.models.ProtocolParameters.objects.filter(
            protocol_id=protocol_obj, parameter_used=True
        ).order_by("parameter_order")
        for param in prot_params:
            mol_data_parm[prot_name]["param_heading"].append(param.get_parameter_name())
        mol_data_parm[prot_name]["param_heading_in_string"] = ";".join(
            mol_data_parm[prot_name]["param_heading"]
        )
        mol_data_parm[prot_name]["m_data"] = list(
            core.models.MoleculePreparation.objects.filter(pk__in=mol_ids).values_list(
                "pk", "sample__sample_name", "molecule_code_id"
            )
        )

    return mol_data_parm


def get_sample_objs_in_state(s_state, user=None, friend_list=None):
    """Return a list of sample objects that are in the indicate state.
    Sample are filter by the user and the user friend list. If no user
    is set then no filter is done

    Parameters
    ----------
    s_state : string
        state name to get the samples
    user : user object, optional
        instance of user to filter for their samples, by default None
    friend_list : Bool
        Boolean value to allow for getting samples inside the user friend list

    Returns
    -------
    sample_objs : list of Samples instances
        contains the objects that were matched on the condition. Empty list if
        not match
    """
    if user:
        if friend_list:
            user_list = core.utils.common.get_friend_list(user)
        else:
            user_list = [user]
        return core.models.Samples.objects.filter(
            sample_state__sample_state_name__exact=s_state,
            sample_user__in=user_list,
        )
    else:
        return (
            core.models.Samples.objects.filter(
                sample_state__sample_state_name__exact=s_state
            )
            .order_by("sample_user")
            .order_by("sample_entry_date")
        )


def get_sample_obj_from_sample_name(sample_name):
    if core.models.Samples.objects.filter(sample_name__exact=sample_name).exists():
        sample_obj = core.models.Samples.objects.get(sample_name__exact=sample_name)
        return sample_obj
    return


def get_sample_obj_from_id(sample_id):
    """The function will return the class object from id number of the class.

    Parameters
    ----------
    sample_id : string
        id value to get the sample obj

    Returns:
    --------
    sample_obj : Sample instance or None if not match
    """

    if core.models.Samples.objects.filter(pk__exact=sample_id).exists():
        sample_obj = core.models.Samples.objects.get(pk__exact=sample_id)
    else:
        sample_obj = None
    return sample_obj


def get_sample_project_information(sample_project_obj, sample_obj, join_values=False):
    """Get the sample project, fields and value"""
    s_project_info = {}
    s_project_info["sample_project_name"] = sample_project_obj.get_sample_project_name()
    if core.models.SampleProjectsFields.objects.filter(
        sample_projects_id=sample_project_obj
    ).exists():
        s_project_info["sample_project_field_heading"] = []
        s_project_info["sample_project_field_value"] = []
        sample_project_fields = core.models.SampleProjectsFields.objects.filter(
            sample_projects_id=sample_project_obj
        )
        for s_p_field in sample_project_fields:
            s_project_info["sample_project_field_heading"].append(
                s_p_field.get_field_name()
            )
            if core.models.SampleProjectsFieldsValue.objects.filter(
                sample_id=sample_obj, sample_project_field_id=s_p_field
            ).exists():
                field_value = (
                    core.models.SampleProjectsFieldsValue.objects.filter(
                        sample_id=sample_obj, sample_project_field_id=s_p_field
                    )
                    .last()
                    .get_field_value()
                )
                if s_p_field.get_field_type() == "Date":
                    field_value = field_value.replace(" 00:00:00", "")
            else:
                field_value = core.core_config.VALUE_NOT_PROVIDED
            s_project_info["sample_project_field_value"].append(field_value)

    if join_values:
        s_project_info["sample_project_join_value"] = [
            list(x)
            for x in zip(
                s_project_info["sample_project_field_heading"],
                s_project_info["sample_project_field_value"],
            )
        ]
    return s_project_info


def get_sample_states():
    """
    Description:
        The function will return the sample states defined in database.
    Return:
        sample_states.
    """
    sample_states = []
    if core.models.StatesForSample.objects.all().exists():
        states = core.models.StatesForSample.objects.all()
        for state in states:
            sample_states.append(state.get_sample_state())
    return sample_states


def get_species():
    """
    Description:
        The function will return the name of the species defined in database.
    Input:
        none
    Variables:
        species_names # list containing all species names
    Return:
        species_names.
    """
    species_names = []
    if core.models.Species.objects.filter().exists():
        all_species = core.models.Species.objects.all()

        for species in all_species:
            species_names.append(species.get_name())
    return species_names


def get_sample_project_obj_from_id(s_project_id):
    """
    Description:
        The fuction return sample project object from id
    Input:
        s_project_id    # id of the sample project
    Output:
        s_project_obj
    """
    s_project_obj = None
    if core.models.SampleProjects.objects.filter(pk__exact=s_project_id).exists():
        s_project_obj = core.models.SampleProjects.objects.get(pk__exact=s_project_id)
    return s_project_obj


def get_sample_project_field_obj_from_id(sample_project_field_id):
    """
    Description:
        The function will return the sample project field object from id number.
    Input:
        sample_project_field_id
    Return:
        sample_project_field_obj.
    """
    sample_project_field_obj = False
    if core.models.SampleProjectsFields.objects.filter(
        pk__exact=sample_project_field_id
    ).exists():
        sample_project_field_obj = core.models.SampleProjectsFields.objects.get(
            pk__exact=sample_project_field_id
        )
    return sample_project_field_obj


def get_selection_from_excel_data(data, heading, check_field, field_id):
    """_summary_

    Args:
        data (_type_): _description_
        heading (list): _description_
        check_field (string_): _description_
        field_id (_string: _description_

    Returns:
        list: _description_
    """

    selected = []
    selected_row = []
    excel_data = json.loads(data)
    # Convert excel list-list to dictionary with field_names
    excel_json_data = core.utils.common.jspreadsheet_to_dict(heading, excel_data)
    for row in excel_json_data:
        if check_field is not None:
            if row[check_field] is True or row[check_field] != "":
                selected.append(row[field_id])
                selected_row.append(row)
        else:
            selected.append(row[field_id])
            selected_row.append(row)

    return selected, selected_row


def get_table_record_molecule(samples, apps_name):
    """
    Description:    The function get the sample ids to create the molecule table where
            define the type of molecule and the protocol usec from the extracion.
    Input:
        samples     # list of the samples to be include in the table
    Functions:
        get_modules_type         # located at this file
        get_molecule_protocols   # located at this file
    Variables:
        molecule_information # dictionary which collects all info
    Return:
        molecule_information #
    """

    molecule_information = {}
    molecule_information["headings"] = (
        core.core_config.HEADING_FOR_MOLECULE_PROTOCOL_DEFINITION
    )
    sample_code_ids = []
    valid_samples = []
    for sample in samples:
        try:
            if core.models.Samples.objects.filter(pk__exact=int(sample)).exists():
                valid_samples.append(int(sample))
        except Exception:
            continue
        if len(valid_samples) == 0:
            molecule_information["ERROR"] = True
            return molecule_information
    sample_code_id = []
    molecule_information["data"] = []
    for sample in valid_samples:
        sample_obj = get_sample_obj_from_id(sample)
        sample_code_id = sample_obj.get_sample_code()
        sample_code_ids.append(sample_code_id)
        data = [""] * (
            len(core.core_config.HEADING_FOR_MOLECULE_PROTOCOL_DEFINITION) + 1
        )
        data[0] = sample_code_id
        data[1] = sample_obj.get_sample_type()
        data[-1] = sample
        molecule_information["data"].append(data)

    molecule_information["type_of_molecules"] = get_modules_type()
    (
        molecule_information["protocols_dict"],
        molecule_information["protocol_list"],
    ) = get_molecule_protocols(apps_name)
    molecule_information["number_of_samples"] = len(valid_samples)
    molecule_information["table_length"] = len(
        core.core_config.HEADING_FOR_MOLECULE_PROTOCOL_DEFINITION
    )
    molecule_information["protocol_type"] = list(
        molecule_information["protocols_dict"].keys()
    )
    molecule_information["protocol_filter_selection"] = []
    molecule_information["sample_code_ids"] = ",".join(sample_code_ids)
    for key, value in molecule_information["protocols_dict"].items():
        molecule_information["protocol_filter_selection"].append([key, value])
    return molecule_information


def get_type_of_sample_information(sample_type_id):
    """
    Description:
        The function return a dictionary with the information to display the type of sample
        and the optional fields
    Input:
        sample_type_id # id for the type of sample to display
    Return:
        sample_type_data
    """
    mandatory_fields = [
        "patient_core",
        "sample_name",
        "lab_request",
        "sample_type",
        "species",
        "sample_project",
        "sample_entry_date",
        "collection_sample_date",
        "sample_location",
        "only_recorded",
    ]
    sample_type_data = {}
    sample_type_data["mandatory_data"] = []
    if core.models.SampleType.objects.filter(pk__exact=sample_type_id).exists():
        sample_type_obj = core.models.SampleType.objects.get(pk__exact=sample_type_id)
        s_type_mandatory = sample_type_obj.get_mandatory_values()
        sample_type_data["sample_type_name"] = sample_type_obj.get_name()
        avail_fields = core.models.Samples._meta.get_fields()
        for m_field in mandatory_fields:
            for avail_field in avail_fields:
                if m_field == avail_field.name:
                    try:
                        display_name = avail_field.verbose_name
                    except Exception:
                        display_name = avail_field.name

                    if avail_field.name in s_type_mandatory:
                        sample_type_data["mandatory_data"].append(
                            [display_name, "Mandatory"]
                        )
                    else:
                        sample_type_data["mandatory_data"].append(
                            [display_name, "Not required"]
                        )
                    break
    else:
        sample_type_data["ERROR"] = (
            core.core_config.ERROR_TYPE_OF_SAMPLE_ID_DOES_NOT_EXISTS
        )
    return sample_type_data


def group_molecules_by_protocol(molecule_ids):
    protocols = {}
    for molecule_id in molecule_ids:
        mol_obj = get_molecule_obj_from_id(molecule_id)
        mol_protocol_obj = mol_obj.get_protocol_obj()
        if mol_protocol_obj not in protocols:
            protocols[mol_protocol_obj] = []
        protocols[mol_protocol_obj].append(molecule_id)
    return protocols


def increase_unique_value(old_unique_number):
    """
    Description:
        The function will increase the sample unique value.
    Input:
        old_unique_number # contains the last unique value store in Database
    Return:
        The increased number.
    """
    split_value = old_unique_number.split("-")
    number = int(split_value[1]) + 1
    letter = split_value[0]

    if number > 9999:
        number = 1
        index_letter = list(split_value[0])
        if index_letter[2] == "Z":
            if index_letter[1] == "Z":
                index_letter[0] = chr(ord(index_letter[0]) + 1)
                index_letter[1] = "A"
                index_letter[2] = "A"
            else:
                index_letter[1] = chr(ord(index_letter[1]) + 1)
                index_letter[2] = "A"

            index_letter = "".join(index_letter)
        else:
            index_letter[2] = chr(ord(index_letter[2]) + 1)

        letter = "".join(index_letter)

    number_str = str(number)
    number_str = number_str.zfill(4)
    return str(letter + "-" + number_str)


def modify_fields_in_sample_project(form_data):
    """
    Description:    The function get the project field value and check if there is
        some changes. If change then replace the old values by thenew ones
    Input:
        form_data     # form data from user
    Return:
        saved_fields #
    """

    sample_project_id = form_data["sample_project_id"]
    # parameter_ids = form_data["parameter_ids"].split(",")
    # parameter_names = form_data["parameter_names"].split(",")
    sample_project_obj = core.models.SampleProjects.objects.get(
        pk__exact=sample_project_id
    )
    json_data = json.loads(form_data["table_data1"])

    fields = core.core_config.HEADING_FOR_MODIFY_SAMPLE_PROJECT_FIELDS.copy()
    # add samnple parameter id on field values
    fields.append("s_p_field_id")
    excel_json_data = core.utils.common.jspreadsheet_to_dict(fields, json_data)

    saved_fields = {}
    saved_fields["fields"] = []
    saved_fields["heading"] = core.core_config.HEADING_FOR_SAMPLE_PROJECT_FIELDS
    saved_fields["sample_project_name"] = sample_project_obj.get_sample_project_name()

    # Delete existing optionns to add new values. Even they are the same,
    # is easier to remmmve all and created again instaed of infividual checking
    if core.models.SamplesProjectsTableOptions.objects.filter(
        sample_project_field__sample_projects_id=sample_project_obj
    ).exists():
        # Delete existing information
        s_p_option_objs = core.models.SamplesProjectsTableOptions.objects.filter(
            sample_project_field__sample_projects_id=sample_project_obj
        )
        for s_p_option_obj in s_p_option_objs:
            s_p_option_obj.delete()

    for row_line in excel_json_data:
        # ignore empty lines
        if row_line["Field name"] == "" and row_line["Change field name"] == "":
            continue

        sample_project_field_obj = None
        if row_line["Classification"] == "":
            row_line["SampleProjectFieldClassificationID"] = None
        else:
            # create classification object if not exists
            if not core.models.SampleProjectFieldClassification.objects.filter(
                classification_name__iexact=row_line["Classification"],
                sample_projects_id=sample_project_obj,
            ).exists():
                c_data = {}
                c_data["sample_project_id"] = sample_project_obj
                c_data["classification_name"] = row_line["Classification"]
                classification_obj = core.models.SampleProjectFieldClassification.objects.create_sample_project_field_classification(
                    c_data
                )
            else:
                classification_obj = (
                    core.models.SampleProjectFieldClassification.objects.filter(
                        classification_name__iexact=row_line["Classification"],
                        sample_projects_id=sample_project_obj,
                    ).last()
                )
            row_line["SampleProjectFieldClassificationID"] = classification_obj
        # Check if new parameter is added
        if row_line["Field name"] == "" and row_line["Change field name"] != "":
            # Add new field
            row_line["Field name"] = row_line["Change field name"]
            row_line["sample_project_id"] = sample_project_obj
            sample_project_field_obj = (
                core.models.SampleProjectsFields.objects.create_sample_project_fields(
                    row_line
                )
            )

        # check if field is a list to create the new opt fields in database
        elif row_line["Field name"] != "" and row_line["Change field name"] != "":
            # rename field name
            row_line["Field name"] = row_line["Change field name"]
        # update values for the existing parameters
        if sample_project_field_obj is None:
            sample_project_field_obj = get_sample_project_field_obj_from_id(
                row_line["s_p_field_id"]
            )
            if not sample_project_field_obj:
                # Unable to find the object class. Skipping this change
                continue
            sample_project_field_obj.update_sample_project_fields(row_line)
        if row_line["Field type"] == "Options List":
            option_list_values = row_line["Option Values"].split(",")
            for opt_value in option_list_values:
                value = opt_value.strip()
                if value == "":
                    continue
                data = {"s_proj_obj": sample_project_field_obj}
                data["opt_value"] = value
                core.models.SamplesProjectsTableOptions.objects.create_new_s_proj_table_opt(
                    data
                )
        saved_fields["fields"].append(
            sample_project_field_obj.get_sample_project_fields_name(include_search=True)
        )

    return saved_fields


def record_molecule_use(from_data, app_name):
    """
    Description:    The function collect the name for the molecule use field and the tag if massive
                and store it on database .
                Returns the molecule use object created.
    Input:
        form_data   # form from the user
        app_name    # application name to assign the right molecule use
    Return:
        molecule_use_information #
    """
    molecule_use_information = {}
    if core.models.MoleculeUsedFor.objects.filter(
        used_for__exact=from_data["moleculeUseName"]
    ).exists():
        molecule_use_information["ERROR"] = (
            core.core_config.ERROR_MOLECULE_USE_FOR_EXISTS
        )
        return molecule_use_information
    molecule_use_data = {}
    molecule_use_data["usedFor"] = from_data["moleculeUseName"]
    molecule_use_data["apps_name"] = app_name
    if "requiresMassive" in from_data:
        molecule_use_data["massiveUse"] = True
    else:
        molecule_use_data["massiveUse"] = False
    core.models.MoleculeUsedFor.objects.create_molecule_use_for(molecule_use_data)
    molecule_use_information["new_defined_molecule_use"] = from_data["moleculeUseName"]
    return molecule_use_information


def record_molecules(samples, excel_data, heading, user, app_name):
    """Recored the molecues defined in excel_data.  If information is missing
        returns the data to display again for correcting.

    Args:
        samples (list): _description_
        excel_data (list): _description_
        heading (list): _description_
        user (string): _description_
        app_name (string): _description_

    Returns:
        dict: return a dictionary with used protols and the list of molecule ids
    """

    incompleted = []
    for idx in range(len(samples)):
        if "" in list(excel_data[idx].values()):
            incompleted.append(samples[idx])
    if len(incompleted) > 0:
        # restore user imput data
        return_data = []
        for row in excel_data:
            r_data = []
            for item in heading:
                r_data.append(row[item])
            return_data.append(r_data)
        # collect data for dropdown selection
        protocol_filter_selection = []
        (protocols_dict, protocol_list) = get_molecule_protocols(app_name)
        for key, value in protocols_dict.items():
            protocol_filter_selection.append([key, value])
        return {
            "incomplete": True,
            "data": return_data,
            "protocol_filter_selection": protocol_filter_selection,
            "protocol_list": protocol_list,
            "type_of_molecules": get_modules_type(),
            "headings": heading[:-1],
        }

    prot_in_molecule = {}
    for row in excel_data:
        molecule_data = {}
        sample_obj = get_sample_obj_from_id(row["s_id"])
        molecule_data["app_name"] = app_name
        molecule_data["sample"] = sample_obj
        molecule_data["user"] = user
        molecule_data["molecule_type"] = row["Molecule type"]
        molecule_data["extraction_type"] = row["Type of Extraction"]
        molecule_data["molecule_extraction_date"] = row["Extraction date"]
        molecule_data["protocol_used"] = row["Protocol to be used"]
        if core.models.MoleculePreparation.objects.filter(sample=sample_obj).exists():
            last_molecule_code = (
                core.models.MoleculePreparation.objects.filter(sample=sample_obj)
                .last()
                .get_molecule_code_id()
            )
            code_split = re.search(r"(.*_E)(\d+)$", last_molecule_code)
            number_code = int(code_split.group(2))
            number_code += 1
            molecule_data["molecule_code_id"] = code_split.group(1) + str(number_code)
        else:
            molecule_data["molecule_code_id"] = sample_obj.get_sample_code() + "_E1"

        molecule_obj = core.models.MoleculePreparation.objects.create_molecule(
            molecule_data
        )
        mol_prot_obj = molecule_obj.get_protocol_obj()
        # Update Sample state to "Extracted molecule"
        sample_obj.set_state("Extract molecule")
        # create the list with protocols and molecule ids
        if mol_prot_obj not in prot_in_molecule:
            prot_in_molecule[mol_prot_obj] = []

        prot_in_molecule[mol_prot_obj].append(molecule_obj.get_molecule_id())

    return prot_in_molecule


def pending_sample_summary(req_user=None, friend_list=None):
    """_summary_

    Parameters
    ----------
    req_user : _type_, optional
        _description_, by default None
    friend_list : _type_, optional
        _description_, by default None
    """

    pending_data = {}
    pending_data["state"] = OrderedDict()
    pending_data["state_number"] = {}
    pending_state = [
        "Pre-defined",
        "Defined",
        "Molecule extraction",
        "Library preparation",
        "Pool preparation",
        "Sequencing",
    ]
    for p_state in pending_state:
        sample_objs = get_sample_objs_in_state(p_state, req_user, friend_list)
        if len(sample_objs) == 0:
            continue
        pending_data["state_number"][p_state] = len(sample_objs)
        pending_data["state"][p_state] = list(
            sample_objs.values_list(
                "pk",
                "sample_name",
                "lab_request__lab_name_coding",
                "sample_type__sample_type",
                "sample_project__sample_project_name",
                "sample_user__username",
            )
        )

    if req_user is None:
        pending_data["users"] = {}
        pend_users = list(
            core.models.Samples.objects.filter(
                sample_state__sample_state_name__in=pending_state
            )
            .values(User=F("sample_user__username"))
            .annotate(value=Count("sample_name"))
        )
        for pend_user in pend_users:
            pending_data["users"][pend_user["User"]] = pend_user["value"]
    return pending_data


def save_type_of_sample(form_data, app_name):
    """
    Description:
        The function store the new type of sample, together with the names of the mandatory fields

    Input:
        form_data # information collected from the form
        app_name # application name where are sample type are defined
    Return:
        save_s_type
    """
    save_s_type = {}
    mandatory_fields = []
    if core.models.SampleType.objects.filter(
        sample_type__exact=form_data["sampleTypeName"], apps_name__exact=app_name
    ).exists():
        save_s_type["ERROR"] = core.core_config.ERROR_TYPE_OF_SAMPLE_EXISTS
        return save_s_type
    # select the optional fields and get the indexes
    for field in form_data:
        if form_data[field] == "on":
            mandatory_fields.append(field)

    data = {}
    data["sampleType"] = form_data["sampleTypeName"]
    data["apps_name"] = app_name
    data["mandatory_fields"] = ",".join(mandatory_fields)

    # Store in database
    sample_type_obj = core.models.SampleType.objects.create_sample_type(data)
    save_s_type = {}
    save_s_type["new_defined_sample_type"] = form_data["sampleTypeName"]
    save_s_type["new_defined_id"] = sample_type_obj.get_sample_type_id()

    return save_s_type


def search_samples(sample_name, user_name, sample_state, start_date, end_date):
    sample_list = []

    if core.models.Samples.objects.all().exists():
        sample_founds = core.models.Samples.objects.all()
    else:
        return sample_list
    if user_name != "":
        if User.objects.filter(username__exact=user_name).exists():
            user_name_obj = core.models.User.objects.filter(
                username__exact=user_name
            ).last()
            user_friend_list = core.utils.common.get_friend_list(user_name_obj)
            if not sample_founds.filter(sample_user__in=user_friend_list).exists():
                return sample_list
            else:
                sample_founds = sample_founds.filter(sample_user__in=user_friend_list)
        else:
            return sample_list
    if sample_name != "":
        if sample_founds.filter(sample_name__exact=sample_name).exists():
            sample_founds = sample_founds.filter(sample_name__exact=sample_name)
            if len(sample_founds) == 1:
                sample_list.append(sample_founds[0].pk)
                return sample_list

        elif sample_founds.filter(sample_name__icontains=sample_name).exists():
            sample_founds = sample_founds.filter(sample_name__icontains=sample_name)
        else:
            return sample_list
    if sample_state != "":
        sample_founds = sample_founds.filter(
            sample_state__sample_state_name__exact=sample_state
        )

    if start_date != "" and end_date != "":
        sample_founds = sample_founds.filter(generated_at__range=(start_date, end_date))

    if start_date != "" and end_date == "":
        sample_founds = sample_founds.filter(generated_at__gte=start_date)

    if start_date == "" and end_date != "":
        sample_founds = sample_founds.filter(generated_at__lte=end_date)

    if len(sample_founds) == 1:
        sample_list.append(sample_founds[0].pk)
        return sample_list

    for sample in sample_founds:
        sample_list.append(sample.get_info_for_searching())
    return sample_list


def set_molecule_use(molecule_use_data, app_name):
    """_summary_

    Args:
        molecule_use_data (_type_): _description_
        app_name (_type_): _description_

    Returns:
        _type_: _description_
    """
    molecule_update = {
        "data": [],
        "heading": core.core_config.HEADING_FOR_SELECTING_MOLECULE_USE,
    }
    for molecule in molecule_use_data:
        molecule_obj = get_molecule_obj_from_id(molecule["m_id"])
        molecule_obj.set_molecule_use(molecule["Molecule use for"], app_name)
        sample_obj = molecule_obj.get_sample_obj()
        if molecule_obj.get_used_for_massive():
            sample_obj.set_state("Library preparation")
        else:
            sample_obj.set_state("Completed")
        molecule_update["data"].append(
            [
                molecule["Sample Name"],
                molecule["Molecule CodeID"],
                molecule["Molecule use for"],
            ]
        )
    return molecule_update


def set_sample_project_fields(data_form):
    sample_project_id = data_form["sample_project_id"]
    sample_project_obj = core.models.SampleProjects.objects.get(
        pk__exact=sample_project_id
    )
    json_data = json.loads(data_form["table_data1"])
    fields = core.core_config.HEADING_FOR_SAMPLE_PROJECT_FIELDS
    excel_json_data = core.utils.common.jspreadsheet_to_dict(fields, json_data)

    saved_fields = []
    stored_fields = {}
    for row_line in excel_json_data:
        if row_line["Field name"] == "":
            continue

        row_line["sample_project_id"] = sample_project_obj

        # check classification
        classification_name = row_line["Classification"]
        if classification_name != "":
            if core.models.SampleProjectFieldClassification.objects.filter(
                sample_projects_id=sample_project_obj,
                classification_name__iexact=classification_name,
            ).exists():
                classification_obj = (
                    core.models.SampleProjectFieldClassification.objects.filter(
                        sample_projects_id=sample_project_obj,
                        classification_name__iexact=classification_name,
                    ).last()
                )
            else:
                c_data = {
                    "sample_project_id": sample_project_obj,
                    "classification_name": classification_name,
                }
                classification_obj = core.models.SampleProjectFieldClassification.objects.create_sample_project_field_classification(
                    c_data
                )
            row_line["SampleProjectFieldClassificationID"] = classification_obj
        else:
            row_line["SampleProjectFieldClassificationID"] = None
        sample_project_field_obj = (
            core.models.SampleProjectsFields.objects.create_sample_project_fields(
                row_line
            )
        )
        if row_line["Field type"] == "Options List":
            option_list_values = row_line["Option Values"].split(",")
            if len(option_list_values) > 0:
                for opt_value in option_list_values:
                    value = opt_value.strip()
                    if value == "":
                        continue
                    data = {"s_proj_obj": sample_project_field_obj}
                    data["opt_value"] = value
                    core.models.SamplesProjectsTableOptions.objects.create_new_s_proj_table_opt(
                        data
                    )

        saved_fields.append(
            sample_project_field_obj.get_sample_project_fields_name(include_search=True)
        )

    stored_fields["fields"] = saved_fields
    stored_fields["heading"] = core.core_config.HEADING_FOR_SAMPLE_PROJECT_FIELDS
    stored_fields["sample_project_name"] = sample_project_obj.get_sample_project_name()
    return stored_fields


def update_molecule_reused(sample_id, molecule_code_id):
    """
    Description:    The function update the sample state to start reprocessing it and increments
            the number of reused.
    Input:
        sample_id     # sample id
        molecule_code_id # molecule id to be reprocessed
    Return:
        sample_obj #
    """

    sample_obj = get_sample_obj_from_id(sample_id)
    try:
        molecule_obj = core.models.MoleculePreparation.objects.get(
            sample=sample_obj, molecule_code_id__exact=molecule_code_id
        )
    except Exception:
        return None
    molecule_obj.set_increase_reuse()
    return molecule_obj


def update_sample_reused(reprocess_id):
    """
    Description:    The function update the sample state to start reprocessing it and increments
            the number of reused.
    Input:
        reprocess_id     # sample id to be reprocessed
    Return:
        sample_obj #
    """
    sample_obj = core.models.Samples.objects.get(pk__exact=reprocess_id)
    sample_obj.set_increase_reuse()

    return sample_obj


def defined_sample_type(sample_type, app_name):
    """_summary_

    Parameters
    ----------
    sample_type
        sample's sample type obtained from jspreadsheet
    app_name
        application name (wetlab, drylab, core, etc.)
    Returns
    -------
        error_cause
            String with the error explanation:
            "No Type of Samples are defined yet. Check documentation to define them"
    """
    if not core.models.SampleType.objects.filter(apps_name=app_name).exists():
        error_cause = core.core_config.ERROR_NO_DEFINED_TYPE_OF_SAMPLES
        return error_cause

    if not core.models.SampleType.objects.filter(
        sample_type__iexact=sample_type, apps_name__exact=app_name
    ).exists():
        error_cause = core.core_config.ERROR_NO_SAMPLE_TYPE.copy()
        error_cause.insert(1, sample_type)
        return " ".join(error_cause)


def defined_lab_request(lab_request, app_name):
    """_summary_

    Parameters
    ----------
    lab_request
        sample's laboratory requested obtained from jspreadsheet
    app_name
        application name (wetlab, drylab, core, etc.)
    Returns
    -------
        error_cause
            String with the error explanation:
            "No Laboratory is defined yet. Check documentation to define the Laboratory"
    """
    if not core.models.LabRequest.objects.filter(apps_name=app_name).exists():
        error_cause = core.core_config.ERROR_NO_DEFINED_LAB_REQUESTED
        return error_cause

    if not core.models.LabRequest.objects.filter(
        lab_name_coding__iexact=lab_request, apps_name__exact=app_name
    ).exists():
        error_cause = core.core_config.ERROR_NO_LAB_REQUESTED.copy()
        error_cause.insert(1, lab_request)
        return " ".join(error_cause)


def defined_species(sample_species, app_name):
    """_summary_

    Parameters
    ----------
    lab_request
        sample's laboratory requested obtained from jspreadsheet
    app_name
        application name (wetlab, drylab, core, etc.)
    Returns
    -------
        error_cause
            String with the error explanation:
            "No Laboratory is defined yet. Check documentation to define the Laboratory"
    """
    if not core.models.Species.objects.filter(apps_name=app_name).exists():
        error_cause = core.core_config.ERROR_NO_DEFINED_SPECIES
        return error_cause

    if not core.models.Species.objects.filter(
        species_name__iexact=sample_species, apps_name__exact=app_name
    ).exists():
        error_cause = core.core_config.ERROR_NO_SPECIES.copy()
        error_cause.insert(1, sample_species)
        return " ".join(error_cause)
