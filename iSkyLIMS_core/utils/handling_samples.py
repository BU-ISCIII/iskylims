import json
import re
import datetime
from iSkyLIMS_core.core_config import *
from iSkyLIMS_core.models import *
from iSkyLIMS_core.utils.generic_functions import get_friend_list
from iSkyLIMS_core.utils.handling_commercial_kits import get_lot_commercial_kits
from iSkyLIMS_core.utils.handling_protocols import *
from django.contrib.auth.models import User


def create_table_to_select_molecules(samples_list):
    """
    Description:
        The function return a dictionary with the information to display to user to
        select the molecules.
    Input:
        samples_list  : sample list object to get the information
    Variables:

    Return:
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
    sample_information["sample_heading"] = HEADING_FOR_DEFINED_SAMPLES_STATE
    sample_information["sample_code_ids"] = ",".join(sample_code_id)

    return sample_information


def display_molecule_protocol_parameters(molecule_ids, user_obj):
    """
    Description:
        The function return the quality parameters defined for the
        selected protocol.

    Input:
        molecule_ids
    Functions:
        get_protocol_parameters_and_type  # located at iSkyLIMS_core.handling_protocols.py
        get_lot_commercial_kits  # located at iSkyLIMS_core.handling_commercial_kits.py
    Return:
        laboratories.
    """

    molecule_recorded = {}
    showed_molecule = []
    molecule_recorded["data"] = []
    pending_molecule = []
    molecule_code_ids = []
    selected_protocol = ""
    parameter_list = []
    for molecule in molecule_ids:
        if not MoleculePreparation.objects.filter(pk=int(molecule)).exists():
            continue
        molecule_obj = MoleculePreparation.objects.get(pk=int(molecule))
        protocol_used = molecule_obj.get_protocol()
        protocol_used_obj = molecule_obj.get_protocol_obj()

        if selected_protocol == "":
            if ProtocolParameters.objects.filter(
                protocol_id__exact=protocol_used_obj
            ).exists():
                selected_protocol = protocol_used
                protocol_parameters = ProtocolParameters.objects.filter(
                    protocol_id__exact=protocol_used_obj, parameterUsed=True
                ).order_by("parameterOrder")

                for parameter in protocol_parameters:
                    parameter_list.append(parameter.get_parameter_name())
                length_heading = len(
                    HEADING_FOR_MOLECULE_ADDING_PARAMETERS + parameter_list
                )
                molecule_recorded[
                    "fix_heading"
                ] = HEADING_FOR_MOLECULE_ADDING_PARAMETERS
                molecule_recorded["param_heading"] = parameter_list
                molecule_recorded[
                    "protocol_parameters_heading_type"
                ] = get_protocol_parameters_and_type(protocol_used_obj)
                # if Protocols.objects.filter(name__exact = selected_protocol).exists():
                # protocol_obj = Protocols.objects.get(name__exact = selected_protocol)
                molecule_recorded["lot_kit"] = get_lot_commercial_kits(
                    protocol_used_obj
                )
                # else:
                #     molecule_recorded['lot_kit'] = ''
        if protocol_used == selected_protocol:
            showed_molecule.append(molecule)
            data = [""] * length_heading
            data[0] = molecule_obj.get_molecule_code_id()
            # data[1] = protocol_used
            molecule_recorded["data"].append(data)
            molecule_code_ids.append(molecule_obj.get_molecule_code_id())
        else:
            pending_molecule.append(molecule_obj.get_molecule_id())

    molecule_recorded["molecule_id"] = ",".join(showed_molecule)
    molecule_recorded["molecule_code_ids"] = ",".join(molecule_code_ids)
    molecule_recorded["pending_id"] = ",".join(pending_molecule)
    molecule_recorded["heading_in_excel"] = "::".join(parameter_list)

    return molecule_recorded


def add_molecule_protocol_parameters(form_data):
    """
    Description:
        The function will store in database the molecule parameters.
        Return the list of the molecules updated
    Input:
        request
    Variables:
        laboratories # list containing all laboratory names
    Return:
        molecule_updated_list.
    """

    molecule_parameter_value = {}
    molecule_updated_list = []
    sample_updated_list = []

    molecule_json_data = json.loads(form_data["parameters_data"])
    molecule_ids = form_data["molecules"].split(",")
    molecule_code_ids = form_data["molecule_code_ids"].split(",")
    parameter_heading = form_data["heading_in_excel"].split("::")
    parameters_length = len(molecule_json_data[0])
    fixed_heading_length = len(HEADING_FOR_MOLECULE_ADDING_PARAMETERS)

    # user_lot_commercial_kit_obj = UserLotCommercialKits.objects.get(chipLot__exact = molecule_json_data[0][1])
    # protocol_used_obj = user_lot_commercial_kit_obj.get_protocol_obj_for_kit()

    # protocol_used_obj = Protocols.objects.get(name__exact = molecule_json_data[0][1])
    for row_index in range(len(molecule_json_data)):
        user_lot_commercial_kit_obj = UserLotCommercialKits.objects.get(
            chipLot__exact=molecule_json_data[row_index][1]
        )
        # increase the number of use and updated the last use date
        user_lot_commercial_kit_obj.set_increase_use()
        user_lot_commercial_kit_obj.set_latest_use(datetime.datetime.now())
        right_id = molecule_ids[
            molecule_code_ids.index(molecule_json_data[row_index][0])
        ]

        molecule_obj = get_molecule_obj_from_id(right_id)
        molecule_obj.set_user_lot_kit_obj(user_lot_commercial_kit_obj)
        molecule_obj.set_state("Completed")
        molecule_updated_list.append(molecule_obj.get_molecule_code_id())
        protocol_used_obj = molecule_obj.get_protocol_obj()

        for p_index in range(fixed_heading_length, parameters_length):
            molecule_parameter_value[
                "moleculeParameter_id"
            ] = ProtocolParameters.objects.get(
                protocol_id=protocol_used_obj,
                parameterName__exact=parameter_heading[p_index - fixed_heading_length],
            )
            molecule_parameter_value["molecule_id"] = molecule_obj
            molecule_parameter_value["parameterValue"] = molecule_json_data[row_index][
                p_index
            ]
            new_parameters_data = (
                MoleculeParameterValue.objects.create_molecule_parameter_value(
                    molecule_parameter_value
                )
            )

        # Update sample state
        sample_obj = molecule_obj.get_sample_obj()
        sample_obj.set_state("Pending for use")
        # Update sample list
        # sample_updated_list.append(sample_obj.get_sample_id())

    return molecule_updated_list


def analyze_input_samples(request, app_name):
    """
    Description:
        The function will get the samples data that user filled in the form.
        defined_samples are the samples that either has no sample project or for
            the sample projects that no requires additional data
        pre_definde_samples are the ones that requires additional Information
        For already defined samples, no action are done on them and  they are included in not_valid_samples.
        it will return a dictionary which contains the processed samples.
    Input:
        request
        app_name
    Functions:
        check_if_sample_already_defined : located at this file
        check_empty_fields :            located at this file
        check_patient_code_exists :     located at this file
        create_empty_patient :          located at this file
        increase_unique_value :         located at this file
    Constants:
        HEADING_FOR_DISPLAY_RECORDED_SAMPLES
        HEADING_FOR_RECORD_SAMPLES
        OPTIONAL_SAMPLES_FIELDS
    Variables:
        defined_samples  # contains the list of sample in defined state
        samples_continue  # samples id's from the samples in defined state
        pre_defined_samples  # contains the list of sample in pre-defined state
        pre_defined_samples_id  # samples id's from the samples in pre-defined state
        invalid_samples     # samples that already exists on database
        invalid_samples_id  # sample id's from the samples that already exists on database
        incomplete_samples  # samples that contain missing information
    Return:
        sample_recorded # Dictionnary with all samples cases .
    """

    na_json_data = json.loads(request.POST["table_data"])
    heading_in_form = HEADING_FOR_RECORD_SAMPLES

    sample_recorded = {}

    defined_samples, samples_continue = [], []
    pre_defined_samples, pre_defined_samples_id = [], []
    invalid_samples, invalid_samples_id = [], []
    incomplete_samples = []
    sample_recorded["all_samples_defined"] = True

    reg_user = request.user.username

    for row in na_json_data:
        sample_data = {}
        sample_name = str(row[heading_in_form.index("Sample Name")])
        if sample_name == "":
            continue

        if not check_if_sample_already_defined(
            row[heading_in_form.index("Sample Name")], reg_user
        ):
            sample_type = str(row[heading_in_form.index("Type of Sample")])

            if sample_type == "":
                incomplete_samples.append(row)
                continue

            for i in range(len(heading_in_form)):
                sample_data[MAPPING_SAMPLE_FORM_TO_DDBB[i][1]] = row[i]
            # optional_fields = []

            # for opt_field in OPTIONAL_SAMPLES_FIELDS:
            #    optional_fields.append(HEADING_FOR_RECORD_SAMPLES.index(opt_field))
            optional_fields = SampleType.objects.get(
                sampleType__exact=sample_type, apps_name__exact=app_name
            ).get_optional_values()

            # check_empty_fields does not consider if the optional values are empty
            if check_empty_fields(row, optional_fields):
                incomplete_samples.append(row)
                sample_recorded["all_samples_defined"] = False
                continue
            # Check if patient code  already exists on database, If not if will be created giving a sequencial dummy value
            if sample_data["p_code_id"] != "":
                patient_obj = check_patient_code_exists(sample_data["p_code_id"])
                if patient_obj is False:
                    # Define the new patient only Patient code is defined
                    patient_obj = create_empty_patient(sample_data["p_code_id"])
            else:
                patient_obj = None
            sample_data["patient"] = patient_obj
            sample_data["user"] = reg_user
            sample_data["sample_id"] = str(reg_user + "_" + sample_name)
            if not Samples.objects.exclude(uniqueSampleID__isnull=True).exists():
                sample_data["new_unique_value"] = "AAA-0001"
            else:
                last_unique_value = (
                    Samples.objects.exclude(uniqueSampleID__isnull=True)
                    .last()
                    .uniqueSampleID
                )
                sample_data["new_unique_value"] = increase_unique_value(
                    last_unique_value
                )
            # set to Defined state the sample if not required to add more additional data

            if sample_data["project_service"] == "None":
                sample_data["sampleProject"] = None
                if sample_data["onlyRecorded"]:
                    sample_data["sampleState"] = "Completed"
                    sample_data["completedDate"] = datetime.datetime.now()
                else:
                    sample_data["sampleState"] = "Defined"
            else:
                sample_data["sampleProject"] = SampleProjects.objects.get(
                    sampleProjectName__exact=sample_data["project_service"]
                )
                if SampleProjectsFields.objects.filter(
                    sampleProjects_id=sample_data["sampleProject"]
                ).exists():
                    sample_recorded["all_samples_defined"] = False
                    sample_data["sampleState"] = "Pre-Defined"
                else:
                    sample_data["sampleState"] = "Defined"
            sample_data["app_name"] = app_name
            new_sample = Samples.objects.create_sample(sample_data)

            if (
                sample_data["sampleState"] == "Defined"
                or sample_data["sampleState"] == "Completed"
            ):
                defined_samples.append(new_sample.get_sample_definition_information())
                samples_continue.append(new_sample.get_sample_id())
            else:
                # select the samples that requires to add additional Information
                pre_defined_samples.append(new_sample.get_sample_name())
                pre_defined_samples_id.append(new_sample.get_sample_id())
        else:  # get the invalid sample to displays information to user
            sample_recorded["all_samples_defined"] = False
            sample_id = Samples.objects.get(
                sampleName__exact=sample_name
            ).get_sample_id()
            if not "sample_id_for_action" in sample_recorded:
                # get the first no valid sample to ask user for new action on the sample
                sample_recorded["sample_data_for_action"] = Samples.objects.get(
                    sampleName__exact=sample_name
                ).get_sample_definition_information()
                sample_recorded["sample_id_for_action"] = sample_id
                invalid_samples.append(
                    Samples.objects.get(
                        sampleName__exact=sample_name
                    ).get_sample_definition_information()
                )
            else:
                invalid_samples_id.append(sample_id)
                invalid_samples.append(
                    Samples.objects.get(
                        sampleName__exact=sample_name
                    ).get_sample_definition_information()
                )
    ## Add already recorded sample in Pre-defined that were not processed because incomplete informatio in samples
    if "pre_defined_id" in request.POST:
        old_pre_defined_list = request.POST["pre_defined_id"].split(",")
        for old_pre_defined in old_pre_defined_list:
            sample_obj = get_sample_obj_from_id(old_pre_defined)
            pre_defined_samples.append(sample_obj.get_sample_name())
            pre_defined_samples_id.append(old_pre_defined)
    if "pending_pre_defined" in request.POST:
        old_pending_pre_defined_list = request.POST["pending_pre_defined"].split(",")
        for old_pending_pre_defined in old_pending_pre_defined_list:
            sample_obj = get_sample_obj_from_id(old_pending_pre_defined)
            pre_defined_samples.append(sample_obj.get_sample_name())
            pre_defined_samples_id.append(old_pending_pre_defined)
    ##   collect data into sample_recorded
    if len(defined_samples) > 0:
        sample_recorded["defined_samples"] = defined_samples
    if len(invalid_samples) > 0:
        sample_recorded["invalid_samples"] = invalid_samples
        sample_recorded["invalid_samples_id"] = ",".join(invalid_samples_id)
        sample_recorded["invalid_heading"] = HEADING_FOR_DISPLAY_RECORDED_SAMPLES
    if len(incomplete_samples) > 0:
        sample_recorded["incomplete_samples"] = incomplete_samples
    if len(pre_defined_samples) > 0:
        sample_recorded["pre_defined_samples"] = pre_defined_samples
        sample_recorded["pre_defined_samples_id"] = pre_defined_samples_id

    if sample_recorded["all_samples_defined"]:
        sample_recorded["samples_to_continue"] = ",".join(samples_continue)
    sample_recorded["recorded_sample_heading"] = HEADING_FOR_DISPLAY_RECORDED_SAMPLES
    sample_recorded["valid_samples_ids"] = samples_continue
    return sample_recorded


def analyze_input_molecules(request):
    """
    Description:
        The function analyze the user data to assign samples to the molecule protocol.
        Molecule is created for the sample and sample state is updated to "Extracted molecule"

    Input:
        request
    Variables:

    Return:
        molecule_recorded.
    """
    # excel_data = request.POST['table_data']
    molecule_json_data = json.loads(request.POST["molecule_data"])
    samples = request.POST["samples"].split(",")
    heading_in_excel = [
        "sampleID",
        "molecule_type",
        "type_extraction",
        "extractionDate",
        "protocol_type",
        "protocol_used",
    ]

    molecule_recorded = {}
    showed_molecule = []
    pending_molecule = []
    prot_used_in_display = ""
    molecule_recorded["molecule_id"] = []
    molecule_recorded["data"] = []
    incomplete_molecules = []
    incomplete_molecules_ids = []
    for row_index in range(len(molecule_json_data)):

        molecule_data = {}
        if not Samples.objects.filter(pk=int(samples[row_index])).exists():
            continue
        sample_obj = Samples.objects.get(pk=int(samples[row_index]))
        # check_empty_fields does not consider if the optional values are empty
        if check_empty_fields(molecule_json_data[row_index], [""]):
            incomplete_samples.append(molecule_json_data[row_index])
            incomplete_molecules_ids.append(int(samples[row_index]))
            continue
        if MoleculePreparation.objects.filter(sample=sample_obj).exists():
            last_molecule_code = (
                MoleculePreparation.objects.filter(sample=sample_obj)
                .last()
                .get_molecule_code_id()
            )
            code_split = re.search(r"(.*_E)(\d+)$", last_molecule_code)
            number_code = int(code_split.group(2))
            number_code += 1
            molecule_code_id = code_split.group(1) + str(number_code)
        else:
            number_code = 1
            molecule_code_id = sample_obj.get_sample_code() + "_E1"
        protocol_used = molecule_json_data[row_index][
            heading_in_excel.index("protocol_used")
        ]
        protocol_used_obj = Protocols.objects.get(name__exact=protocol_used)
        molecule_type_obj = MoleculeType.objects.get(
            moleculeType__exact=molecule_json_data[row_index][
                heading_in_excel.index("molecule_type")
            ]
        )

        molecule_data["protocolUsed"] = protocol_used_obj
        molecule_data["sample"] = sample_obj
        molecule_data["moleculeType"] = molecule_type_obj
        molecule_data["moleculeCodeId"] = molecule_code_id
        molecule_data["extractionType"] = molecule_json_data[row_index][
            heading_in_excel.index("type_extraction")
        ]
        molecule_data["moleculeExtractionDate"] = molecule_json_data[row_index][
            heading_in_excel.index("extractionDate")
        ]
        molecule_data["numberOfReused"] = str(number_code - 1)

        new_molecule = MoleculePreparation.objects.create_molecule(molecule_data)
        # Update Sample state to "Extracted molecule"
        sample_obj.set_state("Extracted Molecule")
        if prot_used_in_display == "":
            if ProtocolParameters.objects.filter(
                protocol_id__exact=protocol_used_obj
            ).exists():
                prot_used_in_display = protocol_used
                protocol_parameters = ProtocolParameters.objects.filter(
                    protocol_id__exact=protocol_used_obj, parameterUsed=True
                ).order_by("parameterOrder")
                parameter_list = []
                for parameter in protocol_parameters:
                    parameter_list.append(parameter.get_parameter_name())
                length_heading = len(
                    HEADING_FOR_MOLECULE_ADDING_PARAMETERS + parameter_list
                )
                molecule_recorded[
                    "fix_heading"
                ] = HEADING_FOR_MOLECULE_ADDING_PARAMETERS
                molecule_recorded["param_heading"] = parameter_list

        if protocol_used == prot_used_in_display:
            showed_molecule.append(new_molecule.get_id())
            data = [""] * length_heading
            data[0] = molecule_code_id
            data[1] = protocol_used
            molecule_recorded["data"].append(data)
        else:
            pending_molecule.append(new_molecule.get_id())
    molecule_recorded["molecule_id"] = ",".join(showed_molecule)
    molecule_recorded["pending_id"] = ",".join(pending_molecule)
    molecule_recorded["heading_in_excel"] = ",".join(parameter_list)

    return molecule_recorded


def analyze_input_sample_project_fields(form_data):
    """
    Description:
        The function analyze the user data to assign values to the sample project field.


    Input:
        form_data
    Functions:
        get_sample_obj_from_id   : located at this file

    Return:
        sample_to_display.
    """
    sample_recorded = {}
    field_value_json_data = json.loads(form_data["table_data"])
    samples_name = form_data["pre_defined_samples"].split(",")
    samples_ids = form_data["pre_defined_id"].split(",")
    pending_ids = form_data["pending_pre_defined"].split(",")
    heading_list = form_data["pre_defined_heading"].split(",")
    sample_to_display = []
    for i in range(len(samples_ids)):
        right_id = samples_ids[samples_name.index(field_value_json_data[i][0])]
        sample_obj = get_sample_obj_from_id(right_id)
        sample_project_obj = sample_obj.get_sample_project_obj()
        for j in range(len(heading_list)):
            sample_project_field = SampleProjectsFields.objects.get(
                sampleProjects_id=sample_project_obj,
                sampleProjectFieldName__exact=heading_list[j],
            )
            field_value = {}
            field_value["sample_id"] = sample_obj
            field_value["sampleProjecttField_id"] = sample_project_field
            field_value["sampleProjectFieldValue"] = field_value_json_data[i][j + 1]
            new_sample_project_f_value = (
                SampleProjectsFieldsValue.objects.create_project_field_value(
                    field_value
                )
            )
        sample_to_display.append([field_value_json_data[i][0], right_id])
        if sample_obj.is_only_recorded():
            sample_obj.set_state("Completed")
        else:
            # Update Sample state to defined
            sample_obj.set_state("Defined")
    sample_recorded["display_samples"] = sample_to_display

    return sample_recorded


def build_record_sample_form(app_name):
    """
    Description:
        The function collect the stored information of  species, sample origin and sample type to use in the
        selected form.
    Input:

    Functions:
        get_species             located at this file
        get_lab_requested       located at this file
        get_sample_type         located at this file
    Variables:
        sample_information:     Dictionnary to collect the information
    Return:
        sample_information
    """

    sample_information = {}
    sample_information["species"] = get_species()
    sample_information["lab_requested"] = get_lab_requested()
    sample_information["sampleType"] = get_sample_type(app_name)
    sample_information["sample_project"] = get_defined_sample_projects(app_name)
    sample_information["sample_project"].insert(0, "None")
    return sample_information


def check_if_sample_already_defined(sample_name, reg_user):
    if Samples.objects.filter(
        sampleName__exact=sample_name, sampleUser__username__exact=reg_user
    ).exists():
        return True
    else:
        return False


def check_if_sample_project_id_exists(sample_project_id):
    if SampleProjects.objects.filter(pk__exact=sample_project_id).exists():
        return True
    return False


def check_if_sample_project_exists(sample_project, app_name):
    """
    Description:
        The function check if sample project name is defined in database.
    Input:
        sample_project:       sample project name
        app_name :           # application name
    Return:
        False is sample project does not exists. True if sample project exists
    """
    if SampleProjects.objects.filter(
        sampleProjectName__iexact=sample_project, apps_name__exact=app_name
    ).exists():
        return True
    return False


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
    if PatientCore.objects.filter(patientCode__iexact=p_code_id).exists():
        patient_obj = PatientCore.objects.get(patientCode__iexact=p_code_id)
    else:
        return False

    return patient_obj


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

    patient_obj = PatientCore.objects.create_patient(patient_data)

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
    new_sample_project = SampleProjects.objects.create_sample_project(s_project_data)

    new_sample_project_id = new_sample_project.get_id()

    return new_sample_project_id


def create_table_pending_use(sample_list, app_name):
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
    use_type["types"] = []
    use_type["data"] = []
    molecule_code_ids = []
    molecule_ids = []

    if MoleculeUsedFor.objects.filter(apps_name__exact=app_name).exists():
        m_used = MoleculeUsedFor.objects.filter(apps_name__exact=app_name)
        for used in m_used:
            use_type["types"].append(used.get_molecule_use_name())

    use_type["heading"] = HEADING_FOR_SELECTING_MOLECULE_USE
    length_heading = len(HEADING_FOR_SELECTING_MOLECULE_USE)
    if MoleculePreparation.objects.filter(
        moleculeUsedFor=None, state__moleculeStateName__exact="Completed"
    ).exists():
        molecules = MoleculePreparation.objects.filter(
            moleculeUsedFor=None, state__moleculeStateName__exact="Completed"
        )
        for molecule in molecules:
            molecule_code_id = molecule.get_molecule_code_id()
            data = [""] * length_heading
            data[0] = molecule.get_sample_name()
            data[1] = molecule_code_id
            molecule_ids.append(molecule.get_molecule_id())
            molecule_code_ids.append(molecule_code_id)
            use_type["data"].append(data)

    use_type["molecule_ids"] = ",".join(molecule_ids)
    use_type["molecule_code_ids"] = ",".join(molecule_code_ids)
    return use_type


def create_table_user_molecules(user_owner_molecules):
    """
    Description:
        The function prepare the user molecule information to display.
    Input:
        user_owner_molecules:  list of molecules belongs to user
    Return:
        molecule_data
    """
    molecule_data = {}
    molecule_data["data"] = []
    for molecule in user_owner_molecules:
        data = []
        data.append(molecule.get_sample_name())
        data += molecule.get_molecule_information()
        molecule_data["data"].append(data)

    molecule_data["molecule_heading"] = HEADING_FOR_USER_PENDING_MOLECULES

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
    sample_project_obj = SampleProjects.objects.get(pk__exact=sample_project_id)

    sample_project_data[
        "sample_project_name"
    ] = sample_project_obj.get_sample_project_name()
    sample_project_data["sample_project_id"] = sample_project_id
    sample_project_data["heading"] = HEADING_FOR_SAMPLE_PROJECT_FIELDS
    return sample_project_data


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

    if SampleType.objects.filter(apps_name__exact=app_name).exists():
        s_types = SampleType.objects.filter(apps_name__exact=app_name)
        for s_type in s_types:
            defined_sample_types.append([s_type.get_sample_type_id, s_type.get_name])
        sample_types["defined_sample_types"] = defined_sample_types
    sample_types["optional_values"] = HEADING_FOR_OPTIONAL_FIELD_SAMPLES
    return sample_types


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
    sample_type_data = {}
    sample_type_data["optional_data"] = []
    if SampleType.objects.filter(pk__exact=sample_type_id).exists():
        sample_type_obj = SampleType.objects.get(pk__exact=sample_type_id)
        opt_list = sample_type_obj.get_optional_values()
        sample_type_data["sample_type_name"] = sample_type_obj.get_name()
        for i in range(len(HEADING_FOR_RECORD_SAMPLES)):
            if i in opt_list:
                sample_type_data["optional_data"].append(
                    [HEADING_FOR_RECORD_SAMPLES[i], "Not Required"]
                )
            else:
                sample_type_data["optional_data"].append(
                    [HEADING_FOR_RECORD_SAMPLES[i], "Mandatory"]
                )
    else:
        sample_type_data["ERROR"] = ERROR_TYPE_OF_SAMPLE_ID_DOES_NOT_EXISTS
    return sample_type_data


def save_type_of_sample(form_data, app_name):
    """
    Description:
        The function store the new type of sample, together with the index of the optional fields
        that can be empty

    Input:
        form_data # information collected from the form
        app_name # application name where are sample type are defined
    Return:
        save_s_type
    """
    save_s_type = {}
    mandatory_index_field_list = []
    if SampleType.objects.filter(
        sampleType__exact=form_data["sampleTypeName"], apps_name__exact=app_name
    ).exists():
        save_s_type["ERROR"] = ERROR_TYPE_OF_SAMPLE_EXISTS
        return save_s_type
    # select the optional fields and get the indexes
    # get the main mandatory fields when recording a new sample
    for m_field in HEADING_FOR_RECORD_SAMPLES:
        if m_field in HEADING_FOR_OPTIONAL_FIELD_SAMPLES:
            continue
        mandatory_index_field_list.append(HEADING_FOR_RECORD_SAMPLES.index(m_field))
    for field in HEADING_FOR_RECORD_SAMPLES:
        if not field in form_data:
            continue
        mandatory_index_field_list.append(HEADING_FOR_RECORD_SAMPLES.index(field))
    # get the index that are optional
    index_opt_field_list = list(
        map(
            str,
            list(
                set(list(range(0, len(HEADING_FOR_RECORD_SAMPLES), 1)))
                - set(mandatory_index_field_list)
            ),
        )
    )

    data = {}
    data["sampleType"] = form_data["sampleTypeName"]
    data["apps_name"] = app_name
    data["optional_fields"] = ",".join(index_opt_field_list)

    sample_type = SampleType.objects.create_sample_type(data)

    save_s_type["new_defined_sample_type"] = form_data["sampleTypeName"]
    save_s_type["new_defined_id"] = sample_type.get_sample_type_id()
    return save_s_type


def get_sample_project_information(sample_project_obj, sample_obj):
    """Get the sample project, fields and value"""
    s_project_info = {}
    s_project_info[
        "sample_project_name"
    ] = sample_project_obj.get_sample_project_name()
    if SampleProjectsFields.objects.filter(
        sampleProjects_id=sample_project_obj
    ).exists():
        s_project_info["sample_project_field_heading"] = []
        s_project_info["sample_project_field_value"] = []
        sample_project_fields = SampleProjectsFields.objects.filter(
            sampleProjects_id=sample_project_obj
        )
        for s_p_field in sample_project_fields:
            s_project_info["sample_project_field_heading"].append(
                s_p_field.get_field_name()
            )
            if SampleProjectsFieldsValue.objects.filter(
                sample_id=sample_obj, sampleProjecttField_id=s_p_field
            ).exists():
                field_value = SampleProjectsFieldsValue.objects.get(
                    sample_id=sample_obj, sampleProjecttField_id=s_p_field
                ).get_field_value()
                if s_p_field.get_field_type() == "Date":
                    field_value = field_value.replace(" 00:00:00", "")
            else:
                field_value = VALUE_NOT_PROVIDED
            s_project_info["sample_project_field_value"].append(field_value)
    return s_project_info


def get_sample_definition_heading():
    """Function to return the sample fields if other apps need them"""
    return HEADING_FOR_SAMPLE_DEFINITION

def get_all_sample_information(sample_id, massive):
    sample_information = {}
    sample_information["sample_id"] = sample_id
    parameter_heading_values = []
    if not Samples.objects.filter(pk__exact=sample_id).exists():
        return "Error"
    sample_obj = Samples.objects.get(pk__exact=sample_id)
    sample_information["sample_definition"] = sample_obj.get_info_for_display()
    sample_information["sample_definition_heading"] = HEADING_FOR_SAMPLE_DEFINITION
    # get the sample project information fields
    sample_project_obj = sample_obj.get_sample_project_obj()
    if sample_project_obj is not None:
        sample_information.update(get_sample_project_information(sample_project_obj, sample_obj))

    # check if molecule information exists for the sample
    if MoleculePreparation.objects.filter(sample=sample_obj).exists():
        molecules = MoleculePreparation.objects.filter(sample=sample_obj)
        sample_information[
            "molecule_definition_heading"
        ] = HEADING_FOR_MOLECULE_DEFINITION
        sample_information["molecule_definition"] = []
        sample_information["molecule_parameter_values"] = []
        sample_information["molecule_definition_data"] = []

        for molecule in molecules:
            molecule_definition_data = []
            molecule_definition_data.append(molecule.get_info_for_display())
            protocol_used_obj = molecule.get_protocol_obj()
            if ProtocolParameters.objects.filter(
                protocol_id=protocol_used_obj
            ).exists():
                parameter_names = ProtocolParameters.objects.filter(
                    protocol_id=protocol_used_obj
                ).order_by("parameterOrder")
                molecule_param_heading = ["Molecule CodeID"]
                mol_param_value = [molecule.get_molecule_code_id()]
                for p_name in parameter_names:
                    molecule_param_heading.append(p_name.get_parameter_name())
                    if MoleculeParameterValue.objects.filter(
                        molecule_id=molecule
                    ).exists():
                        try:
                            mol_param_value.append(
                                MoleculeParameterValue.objects.get(
                                    molecule_id=molecule, moleculeParameter_id=p_name
                                ).get_param_value()
                            )
                        except MoleculeParameterValue.DoesNotExist:
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
    defined_samples["heading"] = HEADING_FOR_DEFINED_SAMPLES_STATE
    sample_information = []
    if Samples.objects.filter(
        sampleState__sampleStateName__exact="Defined",
        sampleUser__username__exact=register_user,
    ).exists():
        sample_list = Samples.objects.filter(
            sampleState__sampleStateName__exact="Defined",
            sampleUser__username__exact=register_user,
        ).order_by("generated_at")
        for sample in sample_list:
            sample_information.append(sample.get_info_in_defined_state())
        defined_samples["sample_information"] = sample_information
    return defined_samples


def get_extraction_kits(username):
    """
    Description:
        The function will return the kits that are associated to the user.
    Input:
        username
    Variables:
        laboratories # list containing all laboratory names
    Return:
        laboratories.
    """
    kits = 1
    return


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
    if LabRequest.objects.filter().exists():
        lab_requesteds = LabRequest.objects.all()

        for lab_requested in lab_requesteds:
            lab_requested_places.append(lab_requested.get_lab_request_code())
    return lab_requested_places


def get_info_to_display_sample_project(sample_project_id):
    """
    Description:
        The function return the information for the requested sample project
    Return:
        info_s_project.
    """
    info_s_project = {}
    if SampleProjects.objects.filter(pk__exact=sample_project_id).exists():
        sample_project_obj = SampleProjects.objects.get(pk__exact=sample_project_id)
        # collect data from project
        info_s_project["sample_project_id"] = sample_project_id
        info_s_project["main_data"] = list(
            zip(SAMPLE_PROJECT_MAIN_DATA, sample_project_obj.get_full_info_to_display())
        )

        if SampleProjectsFields.objects.filter(
            sampleProjects_id=sample_project_obj
        ).exists():
            sample_project_fields = SampleProjectsFields.objects.filter(
                sampleProjects_id=sample_project_obj
            ).order_by("sampleProjectFieldOrder")
            s_project_fields_list = []
            for sample_project_field in sample_project_fields:
                s_project_fields_list.append(
                    sample_project_field.get_sample_project_fields_name()
                )
            info_s_project["fields"] = s_project_fields_list
        info_s_project["heading"] = HEADING_FOR_SAMPLE_PROJECT_FIELDS
    else:
        return "ERROR"

    return info_s_project


def get_only_recorded_samples_and_dates():
    """
    Description:
        The function api return a list of list of samples which are defined as only recorded,
        project name, sample type, species, recorded date, and sample id
    Return:
        samples_data
    """
    samples_data = []
    if Samples.objects.filter(onlyRecorded=True).exists():
        sample_objs = Samples.objects.filter(onlyRecorded=True).order_by("generated_at")
        for sample_obj in sample_objs:
            data = [sample_obj.get_sample_name()]
            data.append(sample_obj.get_sample_project())
            data.append(sample_obj.get_sample_type())
            data.append(sample_obj.get_species())
            data.append(sample_obj.get_extraction_date())
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
    if SampleProjects.objects.filter(pk__exact=sample_project_id).exists():
        sample_project_obj = SampleProjects.objects.get(pk__exact=sample_project_id)
        if SampleProjectsFields.objects.filter(
            sampleProjects_id=sample_project_obj
        ).exists():
            sample_project_fields = SampleProjectsFields.objects.filter(
                sampleProjects_id=sample_project_obj
            ).order_by("sampleProjectFieldOrder")
            s_project_fields_list = []
            parameter_ids = []
            parameter_names = []
            for sample_project_field in sample_project_fields:
                parameter_data = (
                    sample_project_field.get_sample_project_fields_for_javascript()
                )
                parameter_names.append(parameter_data[0])
                parameter_ids.append(sample_project_field.get_field_id())
                parameter_data.insert(1, "")
                s_project_fields_list.append(parameter_data)
            parameters_s_project["fields"] = s_project_fields_list
        parameters_s_project["heading"] = HEADING_FOR_MODIFY_SAMPLE_PROJECT_FIELDS
        parameters_s_project["sample_project_id"] = sample_project_id
        parameters_s_project[
            "sample_project_name"
        ] = sample_project_obj.get_sample_project_name()
        parameters_s_project["parameter_names"] = ",".join(parameter_names)
        parameters_s_project["parameter_ids"] = ",".join(parameter_ids)
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
    if SampleProjects.objects.filter(apps_name__exact=app_name).exists():
        s_projects = SampleProjects.objects.filter(apps_name__exact=app_name)
        for s_project in s_projects:
            s_project_data = s_project.get_info_to_display()
            if SampleProjectsFields.objects.filter(
                sampleProjects_id=s_project
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
    if SampleProjects.objects.filter(apps_name__exact=app_name).exists():
        s_projects = SampleProjects.objects.filter(apps_name__exact=app_name)
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
    if MoleculePreparation.objects.filter(pk__exact=molecule_id).exists():
        molecules_obj = MoleculePreparation.objects.get(pk__exact=molecule_id)
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
    if MoleculePreparation.objects.filter(sample=sample_obj).exists():
        molecules_obj = MoleculePreparation.objects.filter(sample=sample_obj)
        return molecules_obj
    else:
        return ""


def get_molecule_in_state(state, user):
    """
    Description:
        The function will return a list with moelcules which are in the state defined
        in state variable
        If user is included the result are limited to this user and their friend list.
        If user is empty then no user restriction applies.
    Input:
        state       # state name to match
        user        # user object to limit the results
    Functions:
        get_friend_list
    Variables:
        molecule_state # Dictionnary with the heading and the molecule information
    Return:
        molecule_objs.
    """
    molecule_objs = ""
    if user != "":
        user_friend_list = get_friend_list(user)
        if MoleculePreparation.objects.filter(
            state__moleculeStateName__exact=state, moleculeUser__in=user_friend_list
        ).exists():
            molecule_objs = MoleculePreparation.objects.filter(
                state__moleculeStateName__exact=state, moleculeUser__in=user_friend_list
            ).order_by("generated_at")

    else:
        if MoleculePreparation.objects.filter(
            state__moleculeStateName__exact=state
        ).exists():
            molecule_objs = (
                MoleculePreparation.objects.filter(
                    state__moleculeStateName__exact=state
                )
                .order_by("moleculeUser")
                .order_by("generated_at")
            )

    return molecule_objs


##### For each state get samples per user
def get_samples_in_defined_state(user):
    """
    Description:
        The function will return a list with samples which are in defined state.
    Input:
        state  # string of the state to be matched. If empty then all samples in
                defined state is returned
    Variables:
        sample_type_names # list containing all sample types names
    Return:
        samples_in_state.
    """
    sample_information = []
    samples_in_state = {}
    if user != "":
        user_friend_list = get_friend_list(user)

        if Samples.objects.filter(
            sampleState__sampleStateName__exact="Defined",
            sampleUser__in=user_friend_list,
        ).exists():
            # if Samples.objects.filter(sampleState__sampleStateName__exact = 'Defined').exists():
            # samples_obj = Samples.objects.filter(sampleState__sampleStateName__exact = 'Defined')
            samples_obj = Samples.objects.filter(
                sampleState__sampleStateName__exact="Defined",
                sampleUser__in=user_friend_list,
            )
            for sample_obj in samples_obj:
                sample_information.append(sample_obj.get_info_in_defined_state())
            samples_in_state["sample_information"] = sample_information
            samples_in_state["sample_heading"] = HEADING_FOR_DEFINED_SAMPLES_STATE
            samples_in_state["length"] = len(sample_information)
            return samples_in_state

        else:
            samples_in_state["length"] = 0
            return samples_in_state
    else:
        if Samples.objects.filter(
            sampleState__sampleStateName__exact="Defined"
        ).exists():
            samples_obj = (
                Samples.objects.filter(sampleState__sampleStateName__exact="Defined")
                .order_by("sampleUser")
                .order_by("sampleEntryDate")
            )
            for sample_obj in samples_obj:
                sample_data = sample_obj.get_info_in_defined_state()
                sample_data.append(sample_obj.get_register_user())
                sample_information.append(sample_data)
            samples_in_state["sample_information"] = sample_information
            samples_in_state[
                "sample_heading"
            ] = HEADING_FOR_DEFINED_SAMPLES_STATE_WETLAB_MANAGER
            samples_in_state["length"] = len(sample_information)
            return samples_in_state

        else:
            samples_in_state["length"] = 0
            return samples_in_state


##### For each state get molecules per user
def get_samples_in_extracted_molecule_state(user):
    """
    Description:
        The function will return a list with samples which are in extracted_molecule state,
        excluding the ones that molecule state are Completed.
    Input:

    Variables:
        molecule_state # Dictionnary with the heading and the molecule information
    Return:
        molecule_state.
    """
    molecule_state = {}
    molecule_information = []
    if user != "":
        user_friend_list = get_friend_list(user)
        if Samples.objects.filter(
            sampleState__sampleStateName__exact="Extract molecule",
            sampleUser__in=user_friend_list,
        ).exists():
            # if Samples.objects.filter(sampleState__sampleStateName__exact = 'Extract molecule').exists():
            samples_obj = Samples.objects.filter(
                sampleState__sampleStateName__exact="Extract molecule",
                sampleUser__in=user_friend_list,
            )
            for sample_obj in samples_obj:
                molecules = MoleculePreparation.objects.filter(sample=sample_obj)

                for molecule in molecules:
                    if molecule.get_state() == "Completed":
                        continue
                    sample_information = []
                    sample_information.append(sample_obj.get_extraction_date())
                    sample_information.append(sample_obj.get_sample_name())
                    molecule_data = molecule.get_molecule_information()
                    molecule_information.append(sample_information + molecule_data)

            molecule_state["molecule_information"] = molecule_information
            molecule_state["molecule_heading"] = HEADING_FOR_EXTRACTED_MOLECULES_STATE
            molecule_state["length"] = len(molecule_information)
            return molecule_state

        else:
            molecule_state["length"] = 0
            return molecule_state
    else:
        if Samples.objects.filter(
            sampleState__sampleStateName__exact="Extract molecule"
        ).exists():
            samples_obj = (
                Samples.objects.filter(
                    sampleState__sampleStateName__exact="Extract molecule"
                )
                .order_by("sampleUser")
                .order_by("sampleEntryDate")
            )
            for sample_obj in samples_obj:
                molecules = MoleculePreparation.objects.filter(sample=sample_obj)
                for molecule in molecules:
                    if molecule.get_state() == "Completed":
                        continue
                    sample_information = sample_obj.get_info_in_defined_state()
                    sample_information.append(sample_obj.get_register_user())
                    molecule_data = molecule.get_molecule_information()
                    molecule_information.append(sample_information + molecule_data)

            molecule_state["molecule_information"] = molecule_information
            molecule_state[
                "molecule_heading"
            ] = HEADING_FOR_EXTRACTED_MOLECULES_STATE_WETLAB_MANAGER
            molecule_state["length"] = len(molecule_information)
            return molecule_state
        else:
            molecule_state["length"] = 0
            return molecule_state


##### End of getting samples by state


def get_sample_obj_from_sample_name(sample_name):
    if Samples.objects.filter(sampleName__exact=sample_name).exists():
        sample_obj = Samples.objects.get(sampleName__exact=sample_name)
        return sample_obj
    return


def get_sample_obj_from_id(sample_id):
    """
    Description:
        The function will return the class object from id number of the class.
    Input:
        sample_id
    Return:
        sample_obj.
    """

    if Samples.objects.filter(pk__exact=sample_id).exists():
        sample_obj = Samples.objects.get(pk__exact=sample_id)
    else:
        sample_obj = None
    return sample_obj


def get_sample_type(app_name):
    """
    Description:
        The function will return the type of samples defined in database.
    Input:
        none
    Variables:
        sample_type_names # list containing all sample types names
    Return:
        sample_type_names.
    """
    sample_type_names = []
    if SampleType.objects.filter(apps_name__exact=app_name).exists():
        sample_types = SampleType.objects.filter(apps_name__exact=app_name)

        for sample in sample_types:
            sample_type_names.append(sample.get_name())
    return sample_type_names


def get_sample_states():
    """
    Description:
        The function will return the sample states defined in database.
    Return:
        sample_states.
    """
    sample_states = []
    if StatesForSample.objects.all().exists():
        states = StatesForSample.objects.all()
        for state in states:
            sample_states.append(state.get_sample_state())
    return sample_states


def get_samples_in_state(state):
    """
    Description:
        The function returns a object list with the samples in the requested state.
    Return:
        sample_objs. False if no samples found in the requested state.
    """
    if Samples.objects.filter(sampleState__sampleStateName__exact=state).exists():
        sample_objs = Samples.objects.filter(sampleState__sampleStateName__exact=state)
        return sample_objs
    else:
        return False


def get_selected_recorded_samples(form_data):
    """
    Description:    The function get the selected recorded samples to add molecule information.
    Input:
        form_data     # form data from user
    Return:
        selected_samples
    """
    selected_samples = []
    samples_json_data = json.loads(form_data)
    for row_index in range(len(samples_json_data)):
        if samples_json_data[row_index][-1] == True:
            selected_samples.append(samples_json_data[row_index][-2])

    return selected_samples


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
    if Species.objects.filter().exists():
        all_species = Species.objects.all()

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
    if SampleProjects.objects.filter(pk__exact=s_project_id).exists():
        s_project_obj = SampleProjects.objects.get(pk__exact=s_project_id)
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
    if SampleProjectsFields.objects.filter(pk__exact=sample_project_field_id).exists():
        sample_project_field_obj = SampleProjectsFields.objects.get(
            pk__exact=sample_project_field_id
        )
    return sample_project_field_obj


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
    types = MoleculeType.objects.all()
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
    p_types = ProtocolType.objects.filter(
        molecule__isnull=False, apps_name__exact=apps_name
    )
    molecule_types = MoleculeType.objects.filter()
    for molecule in molecule_types:
        protocols[molecule.get_name()] = []
    for p_type in p_types:
        protocol_types.append(p_type.get_name())
    for protocol_type in p_types:
        # if protocol_type.get_name() not in protocols :
        #    protocols[protocol_type.get_molecule_type()] = []
        # protocols[protocol_type.get_name()] = []
        protocols_in_type = Protocols.objects.filter(type=protocol_type)
        # protocols_in_type = Protocols.objects.filter(type__protocol_type__exact = protocol_type)
        for p_in_type in protocols_in_type:
            protocol_name = p_in_type.get_name()
            protocols[protocol_type.get_molecule_type()].append(protocol_name)
            # protocols[protocol_type.get_name()].append(protocol_name)
            protocol_list.append(protocol_name)

    return protocols, protocol_list


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


def prepare_sample_input_table(app_name):
    """
    Description: The function collect the species, Lab request, type of
        samples, and heading used in the input table.
        Return a dictionary with collected information.

    Functions:
        build_record_sample_form  : located at this file
    Variables:
        s_information # dictionary which collects all info
    Return:
        s_information #
    """
    # get the choices to be included in the form
    s_information = build_record_sample_form(app_name)
    s_information["heading"] = HEADING_FOR_RECORD_SAMPLES
    s_information["table_size"] = len(HEADING_FOR_RECORD_SAMPLES)
    sample_objs = get_samples_in_state("Pre-defined")
    if sample_objs:
        s_information["pre_defined_samples"] = []
        s_information[
            "pre_defined_heading"
        ] = HEADING_FOR_COMPLETION_SAMPLES_PRE_DEFINED
        for sample_obj in sample_objs:
            s_information["pre_defined_samples"].append(
                sample_obj.get_info_in_defined_state()
            )

    return s_information


def check_if_molecule_use_defined(app_name):
    """
    Description:    The function check if there are defined the use for molecules

    Input:
        app_name    # application name to assign the right molecule use
    Return:
        True or False #
    """
    if MoleculeUsedFor.objects.filter(apps_name__exact=app_name).exists():
        return True
    return False


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
    if MoleculeUsedFor.objects.filter(apps_name__exact=app_name).exists():
        molecule_uses = MoleculeUsedFor.objects.filter(apps_name__exact=app_name)
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
    if MoleculeUsedFor.objects.filter(
        usedFor__exact=from_data["moleculeUseName"]
    ).exists():
        molecule_use_information["ERROR"] = ERROR_MOLECULE_USE_FOR_EXISTS
        return molecule_use_information
    molecule_use_data = {}
    molecule_use_data["usedFor"] = from_data["moleculeUseName"]
    molecule_use_data["apps_name"] = app_name
    if "requiresMassive" in from_data:
        molecule_use_data["massiveUse"] = True
    else:
        molecule_use_data["massiveUse"] = False
    new_molecule_use = MoleculeUsedFor.objects.create_molecule_use_for(
        molecule_use_data
    )
    molecule_use_information["new_defined_molecule_use"] = from_data["moleculeUseName"]
    return molecule_use_information


def record_molecules(form_data, user, app_name):
    """
    Description:    The function store in database the new molecule and molecule_updated_list
                    the sample state to Extracted molecule.
                    When user did not write any of the fields , the sample is added to incomplete_samples
    Input:
        form_data   # form from the user
        user        # logged user
        app_name    # application name to assign the right protocol
    Functions:
        check_empty_fields  : located at this file
    Constant:
        HEADING_FOR_MOLECULE_PROTOCOL_DEFINITION
    Variables:
        molecule_information # dictionary which collects all info
        molecules_code_ids
    Return:
        molecules_recorded with the list of the recorded molecules and the heading to
        display them
    """
    molecule_json_data = json.loads(form_data["molecule_data"])
    samples_ids = form_data["samples"].split(",")
    samples_code_ids = form_data["samples_code_ids"].split(",")
    molecules_recorded = {}
    molecules_ids, molecules_code_ids = [], []
    molecule_list = []
    incomplete_sample_data = []
    incomplete_sample_ids = []
    incomplete_sample_code_ids = []

    heading_in_excel = HEADING_FOR_MOLECULE_PROTOCOL_DEFINITION
    for row_index in range(len(molecule_json_data)):
        right_id = samples_ids[samples_code_ids.index(molecule_json_data[row_index][0])]
        if not Samples.objects.filter(pk__exact=right_id).exists():
            continue
        sample_obj = get_sample_obj_from_id(right_id)

        # check_empty_fields does not consider if the optional values are empty
        if check_empty_fields(molecule_json_data[row_index], [""]):
            incomplete_sample_data.append(molecule_json_data[row_index])
            incomplete_sample_ids.append(right_id)
            incomplete_sample_code_ids.append(molecule_json_data[row_index][0])
            continue

        molecule_data = {}
        protocol_used = molecule_json_data[row_index][
            heading_in_excel.index("Protocol to be used")
        ]
        # if MoleculePreparation.objects.filter(sample = sample_obj, moleculeCodeId__icontains = protocol_used).exists():
        if MoleculePreparation.objects.filter(sample=sample_obj).exists():
            last_molecule_code = (
                MoleculePreparation.objects.filter(sample=sample_obj)
                .last()
                .get_molecule_code_id()
            )
            code_split = re.search(r"(.*_E)(\d+)$", last_molecule_code)
            number_code = int(code_split.group(2))
            number_code += 1
            molecule_code_id = code_split.group(1) + str(number_code)
        else:
            protocol_code = protocol_used.replace(" ", "-")
            molecule_code_id = sample_obj.get_sample_code() + "_E1"

        # protocol_used_obj = Protocols.objects.get(name__exact = protocol_used)
        molecule_used = molecule_json_data[row_index][
            heading_in_excel.index("Molecule type")
        ]

        molecule_data["protocolUsed"] = protocol_used
        molecule_data["app_name"] = app_name
        molecule_data["sample"] = sample_obj
        molecule_data["moleculeType"] = molecule_used
        molecule_data["moleculeCodeId"] = molecule_code_id
        molecule_data["extractionType"] = molecule_json_data[row_index][
            heading_in_excel.index("Type of Extraction")
        ]
        molecule_data["moleculeExtractionDate"] = molecule_json_data[row_index][
            heading_in_excel.index("Extraction date")
        ]
        molecule_data["user"] = user
        # molecule_data['usedForMassiveSequencing'] = massive
        # molecule_data['numberOfReused'] = str(number_code - 1)

        new_molecule = MoleculePreparation.objects.create_molecule(molecule_data)

        molecule_list.append([molecule_code_id, protocol_used])
        # Update Sample state to "Extracted molecule"
        sample_obj.set_state("Extract molecule")
        # Include index key to allow adding quality parameter data
        molecules_ids.append(new_molecule.get_molecule_id())
        molecules_code_ids.append(molecule_code_id)
    if len(molecules_ids) > 0:
        molecules_recorded["heading"] = HEADING_CONFIRM_MOLECULE_RECORDED
        molecules_recorded["molecule_list"] = molecule_list
        molecules_recorded["molecule_ids"] = ",".join(molecules_ids)
        molecules_recorded["molecule_code_ids"] = ",".join(molecules_code_ids)
    if len(incomplete_sample_ids) > 0:
        molecules_recorded["incomplete_sample_data"] = incomplete_sample_data
        molecules_recorded["incomplete_sample_ids"] = ",".join(incomplete_sample_ids)
        molecules_recorded["incomplete_sample_code_ids"] = ",".join(
            incomplete_sample_code_ids
        )

    return molecules_recorded


def prepare_sample_project_input_table(pre_defined_samples_id):
    """
    Description:    The function collects the sample project fields for the samples in the input variable.
            It return the fields heading of the sample project
            In case that samples do not have the same project then it grouped store in database the new molecule and molecule_updated_list
                    the sample state to Extracted molecule.
    Input:
        pre_defined_samples_id  # sample_id list to be processed
    Variables:
        molecule_information # dictionary which collects all info
    Return:
        molecules_recorded with the list of the recorded molecules and the heading to
        display them
    """
    sample_projects = {}
    selected_sample_project = ""
    sample_project_field_heading = []
    updated_pre_defined_samples_id = []
    pending_pre_defined_samples_id = []
    only_field_heading_name = []
    sample_projects["pre_defined_sample_data"] = []
    pre_defined_samples_name = []
    for sample_id in pre_defined_samples_id:
        sample_obj = get_sample_obj_from_id(sample_id)
        s_project_obj = sample_obj.get_sample_project_obj()
        if selected_sample_project == "":
            selected_sample_project = s_project_obj.get_sample_project_name()
            s_project_fields = (
                SampleProjectsFields.objects.filter(sampleProjects_id=s_project_obj)
                .exclude(sampleProjectFieldUsed=None)
                .order_by("sampleProjectFieldOrder")
            )
            for s_project_field in s_project_fields:
                heading_item = []
                heading_item.append(s_project_field.get_field_name())
                heading_item.append(s_project_field.get_field_type())
                heading_item.append(s_project_field.get_field_options_list())
                sample_project_field_heading.append(heading_item)
                only_field_heading_name.append(s_project_field.get_field_name())
            heading_length = len(sample_project_field_heading)

        if selected_sample_project == s_project_obj.get_sample_project_name():
            updated_pre_defined_samples_id.append(sample_id)
            data = [""] * (heading_length + 1)
            sample_name = sample_obj.get_sample_name()
            data[0] = sample_name
            pre_defined_samples_name.append(sample_name)
            sample_projects["pre_defined_sample_data"].append(data)
        else:
            pending_pre_defined_samples_id.append(sample_id)

    sample_projects["updated_pre_defined_samples_id"] = ",".join(
        updated_pre_defined_samples_id
    )
    sample_projects["pending_pre_defined_samples_id"] = ",".join(
        pending_pre_defined_samples_id
    )
    sample_projects["pre_defined_fields_heading_type"] = sample_project_field_heading
    sample_projects["pre_defined_fields_heading_list"] = ",".join(
        only_field_heading_name
    )
    sample_projects["pre_defined_fields_length"] = heading_length + 1
    sample_projects["pre_defined_samples_length"] = len(updated_pre_defined_samples_id)
    sample_projects["pre_defined_samples_name"] = ",".join(pre_defined_samples_name)
    return sample_projects


def get_info_for_reprocess_samples(sample_ids, sample_in_action):
    sample_recorded = {}
    invalid_samples = []

    for sample_id in sample_ids:
        if sample_id == sample_in_action:
            sample_recorded["sample_data_for_action"] = Samples.objects.get(
                pk__exact=sample_id
            ).get_sample_definition_information()
        invalid_samples.append(
            Samples.objects.get(pk__exact=sample_id).get_sample_definition_information()
        )

    sample_recorded["invalid_samples"] = invalid_samples
    # sample_recorded['invalid_samples_id'] = ','.join(sample_ids)
    sample_recorded["invalid_heading"] = HEADING_FOR_DISPLAY_RECORDED_SAMPLES

    return sample_recorded


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
    molecule_information["headings"] = HEADING_FOR_MOLECULE_PROTOCOL_DEFINITION
    sample_code_ids = []
    valid_samples = []
    for sample in samples:
        try:
            if Samples.objects.filter(pk__exact=int(sample)).exists():
                valid_samples.append(int(sample))
        except:
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
        # sample_code_id.append(Samples.objects.get(pk__exact = sample).get_sample_code())
        data = [""] * len(HEADING_FOR_MOLECULE_PROTOCOL_DEFINITION)
        data[0] = sample_code_id
        data[1] = sample_obj.get_sample_type()
        molecule_information["data"].append(data)

    molecule_information["type_of_molecules"] = get_modules_type()
    (
        molecule_information["protocols_dict"],
        molecule_information["protocol_list"],
    ) = get_molecule_protocols(apps_name)
    molecule_information["number_of_samples"] = len(valid_samples)
    molecule_information["table_length"] = len(HEADING_FOR_MOLECULE_PROTOCOL_DEFINITION)
    molecule_information["protocol_type"] = list(
        molecule_information["protocols_dict"].keys()
    )
    molecule_information["protocol_filter_selection"] = []
    molecule_information["sample_code_ids"] = ",".join(sample_code_ids)
    for key, value in molecule_information["protocols_dict"].items():
        molecule_information["protocol_filter_selection"].append([key, value])

    return molecule_information


def search_samples(sample_name, user_name, sample_state, start_date, end_date):
    sample_list = []

    if Samples.objects.all().exists():
        sample_founds = Samples.objects.all()
    else:
        return sample_list
    if user_name != "":
        if User.objects.filter(username__exact=user_name).exists():
            user_name_obj = User.objects.filter(username__exact=user_name).last()
            user_friend_list = get_friend_list(user_name_obj)
            if not sample_founds.filter(sampleUser__in=user_friend_list).exists():
                return sample_list
            else:
                sample_founds = sample_founds.filter(sampleUser__in=user_friend_list)
        else:
            return sample_list
    if sample_name != "":
        if sample_founds.filter(sampleName__exact=sample_name).exists():
            sample_founds = sample_founds.filter(sampleName__exact=sample_name)
            if len(sample_founds) == 1:
                sample_list.append(sample_founds[0].pk)
                return sample_list

        elif sample_founds.filter(sampleName__icontains=sample_name).exists():
            sample_founds = sample_founds.filter(sampleName__icontains=sample_name)
        else:
            return sample_list
    if sample_state != "":
        sample_founds = sample_founds.filter(
            sampleState__sampleStateName__exact=sample_state
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


def set_molecule_use(form_data, app_name):
    """
    Description:    The function get the molecule use decided by user.
            Sample state is changed to Library preparation, and molecule is updated with the use value.
    Input:
        form_data     # form data from user
    Functions:
        get_modules_type         # located at this file
        get_molecule_protocols   # located at this file
    Variables:
        molecule_update # dictionary which collects all info
    Return:
        molecule_update #
    """
    molecule_json_data = json.loads(form_data["molecule_used_for"])
    molecule_ids = form_data["molecule_ids"].split(",")
    molecule_code_ids = form_data["molecule_code_ids"].split(",")
    molecule_update = {}
    molecule_update["data"] = []
    molecules_length = len(molecule_json_data[0])
    for row_index in range(len(molecule_json_data)):
        if molecule_json_data[row_index][2] != "":
            data = []

            right_id = molecule_ids[
                molecule_code_ids.index(molecule_json_data[row_index][1])
            ]
            molecule_obj = get_molecule_obj_from_id(right_id)
            molecule_obj.set_molecule_use(molecule_json_data[row_index][2], app_name)
            sample_obj = molecule_obj.get_sample_obj()
            if molecule_obj.get_used_for_massive():
                sample_obj.set_state("Library preparation")
            else:
                sample_obj.set_state("Completed")
            molecule_update["data"].append(molecule_json_data[row_index])

    if len(molecule_update["data"]) > 0:
        molecule_update["heading"] = HEADING_FOR_SELECTING_MOLECULE_USE
    return molecule_update


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
    parameter_ids = form_data["parameter_ids"].split(",")
    parameter_names = form_data["parameter_names"].split(",")
    json_data = json.loads(form_data["table_data1"])
    sample_project_obj = SampleProjects.objects.get(pk__exact=sample_project_id)
    fields = HEADING_FOR_MODIFY_SAMPLE_PROJECT_FIELDS
    saved_fields = {}
    saved_fields["fields"] = []
    saved_fields["heading"] = HEADING_FOR_SAMPLE_PROJECT_FIELDS
    saved_fields["sample_project_name"] = sample_project_obj.get_sample_project_name()
    # Delete existing optionns to add new values. Even they are the same,
    # is easier to remmmve all and created again instaed of infividual checking

    if SamplesProjectsTableOptions.objects.filter(sampleProjectField__sampleProjects_id=sample_project_obj).exists():
        # Delete existing information
        s_p_option_objs = SamplesProjectsTableOptions.objects.filter(sampleProjectField__sampleProjects_id=sample_project_obj)
        for s_p_option_obj in s_p_option_objs:
            s_p_option_obj.delete()
    for row_data in json_data:
        if row_data[0] == "" and row_data[1] == "":
            continue
        s_p_fields = {}
        for i in range(len(fields)):
            s_p_fields[fields[i]] = row_data[i]

        if row_data[0] == "" and row_data[1] != "":
            # Add new field
            s_p_fields["Field name"] = row_data[1]
            s_p_fields["sample_project_id"] = sample_project_obj
            sample_project_field_obj = SampleProjectsFields.objects.create_sample_project_fields(
                s_p_fields
            )
            saved_fields["fields"].append(sample_project_field_obj.get_sample_project_fields_name())
            # check if field is a list to create the new opt fields on database

        elif row_data[0] != "" and row_data[1] != "":
            # rename field name
            s_p_fields["Field name"] = row_data[1]
        else:
            s_p_fields["Field name"] = row_data[0]
        # Update  Field
            right_id = parameter_ids[parameter_names.index(row_data[0])]
            sample_project_field_obj = get_sample_project_field_obj_from_id(right_id)
            if not sample_project_field_obj:
                # Unable to find the object class. Skipping this change
                continue

            sample_project_field_obj.update_sample_project_fields(s_p_fields)
            saved_fields["fields"].append(
                sample_project_field_obj.get_sample_project_fields_name()
            )

        if row_data[fields.index("Field type")] == "Options List":
            option_list_values = row_data[fields.index("Option Values")].split(",")
            for opt_value in option_list_values:
                value = opt_value.strip()
                if value == "":
                    continue
                data = {"s_proj_obj": sample_project_field_obj}
                data["opt_value"] = value
                SamplesProjectsTableOptions.objects.create_new_s_proj_table_opt(data)
    return saved_fields


def set_sample_project_fields(data_form):
    sample_project_id = data_form["sample_project_id"]
    json_data = json.loads(data_form["table_data1"])
    fields = HEADING_FOR_SAMPLE_PROJECT_FIELDS

    sample_project_obj = SampleProjects.objects.get(pk__exact=sample_project_id)

    saved_fields = []
    stored_fields = {}
    for row_data in json_data:
        if row_data[0] == "":
            continue
        s_p_fields = {}

        s_p_fields["sample_project_id"] = sample_project_obj
        for i in range(len(fields)):
            s_p_fields[fields[i]] = row_data[i]

        if row_data[fields.index("Field type")] == "Option List":
            option_list_values = row_data[fields.index("Option Values")].split(",")
            clean_value_list = []
            for opt_value in option_list_values:
                value = opt_value.strip()
                if value != "":
                    clean_value_list.append(value)

            s_p_fields["Option Values"] = ",".join(clean_value_list)
        else:
            s_p_fields["Option Values"] = ""
        saved_fields.append(
            SampleProjectsFields.objects.create_sample_project_fields(
                s_p_fields
            ).get_sample_project_fields_name()
        )

    stored_fields["fields"] = saved_fields
    stored_fields["heading"] = HEADING_FOR_SAMPLE_PROJECT_FIELDS
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
        molecule_obj = MoleculePreparation.objects.get(
            sample=sample_obj, moleculeCodeId__exact=molecule_code_id
        )
    except:
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
    sample_obj = Samples.objects.get(pk__exact=reprocess_id)
    sample_obj.set_increase_reuse()

    return sample_obj
