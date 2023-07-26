# Generic imports
import json
import re
from datetime import datetime

import jsonschema
import pandas as pd
from jsonschema import Draft202012Validator

# Local imports
import core.core_config
import core.models

import core.utils.samples


def get_sample_projects_names(sample_df):
    """
    """
    sample_df["Project/Service"] = sample_df["Project/Service"].fillna("")
    return sample_df["Project/Service"].unique().tolist()


def check_samples_belongs_to_same_type_and_molecule_protocol(sample_batch_data):
    """
    Description:
        The Function check if only one type of sample is in data and if the mandatory parameters
        for the type of sample are included.
        In case that are more than 1 type of sample it checks one by one the mandatory parameters
    Input:
        sample_batch_data
    Return:
        True or False
    """
    sample_type = sample_batch_data["Type of Sample"].unique().tolist()
    sample_project = sample_batch_data["Project/Service"].unique().tolist()
    protocol = sample_batch_data["ProtocolName"].unique().tolist()
    if len(sample_type) != 1 or len(sample_project) != 1 or len(protocol) != 1:
        return False
    return True


def check_defined_option_values_in_samples(sample_batch_df, package):
    """
    Description:
        The Function checks if values (related to new sample creation) inside data frame are already defined
        in database.
        It also checks if the mandatory parameters in the type of sample are in dataframe
    Input:
        sample_batch_df     # sample data in dataframe
        package     # name of the apps that request the checking

    Return:
        OK or error code
    """
    # Check if option values are already defined
    unique_lab_request = sample_batch_df["Lab requested"].unique().tolist()
    if not core.models.LabRequest.objects.all().exists():
        return (
            core.core_config.ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_DEFINED_LAB_REQUESTED
        )
    lab_requested_values = list(
        core.models.LabRequest.objects.all().values_list("lab_name_coding", flat=True)
    )
    for lab_request in unique_lab_request:
        if lab_request not in lab_requested_values:
            error_cause = (
                core.core_config.ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_LAB_REQUESTED.copy()
            )
            error_cause.insert(1, lab_request)
            return error_cause
    if not core.models.SampleType.objects.filter(apps_name=package).exists():
        return (
            core.core_config.ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_DEFINED_TYPE_OF_SAMPLES
        )
    sample_type_values = list(
        core.models.SampleType.objects.filter(apps_name=package).values_list(
            "sample_type", flat=True
        )
    )
    unique_sample_types = sample_batch_df["Type of Sample"].unique().tolist()
    for sample_type in unique_sample_types:
        if sample_type not in sample_type_values:
            error_cause = (
                core.core_config.ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_SAMPLE_TYPE.copy()
            )
            error_cause.insert(1, sample_type)
            return error_cause
    if not core.models.Species.objects.filter(apps_name=package).exists():
        return core.core_config.ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_DEFINED_SPECIES
    species_values = list(
        core.models.Species.objects.filter(apps_name=package).values_list(
            "species_name", flat=True
        )
    )
    unique_species = sample_batch_df["Species"].unique().tolist()
    for specie in unique_species:
        if specie not in species_values:
            error_cause = (
                core.core_config.ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_SPECIES.copy()
            )
            error_cause.insert(1, specie)
            return error_cause
    # check if sample projects are defined
    if not core.models.SampleProjects.objects.filter(apps_name=package).exists():
        return (
            core.core_config.ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_DEFINED_SAMPLE_PROJECTS
        )
    sample_project_values = list(
        core.models.SampleProjects.objects.filter(apps_name=package).values_list(
            "sample_project_name", flat=True
        )
    )
    unique_sample_projects = sample_batch_df["Project/Service"].unique().tolist()
    for sample_project in unique_sample_projects:
        if sample_project not in sample_project_values:
            error_cause = (
                core.core_config.ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_SAMPLE_PROJECTS.copy()
            )
            error_cause.insert(1, sample_project)
            return error_cause

    # check if additional sample Project parameters are include in bathc file
    if core.models.SampleProjectsFields.objects.filter(
        sample_projects_id__sample_project_name__exact=sample_project,
        sample_projects_id__apps_name=package,
    ).exists():
        sample_project_fields_objs = core.models.SampleProjectsFields.objects.filter(
            sample_projects_id__sample_project_name__exact=sample_project,
            sample_projects_id__apps_name=package,
        )
        for sample_project_fields_obj in sample_project_fields_objs:
            if not sample_project_fields_obj.get_field_name() in sample_batch_df:
                error_cause = (
                    core.core_config.ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NOT_SAMPLE_FIELD_DEFINED.copy()
                )
                error_cause.insert(1, sample_project_fields_obj.get_field_name())
                return error_cause

    return "OK"


def check_molecule_has_same_data_type(sample_batch_df, package):
    """
    Description:
        The Function checks if values (related to molecule creation) inside data frame are already defined in database.
    Input:
        sample_batch_df     # sample data in dataframe
        package     # name of the apps that request the checking
    Constants:
        ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_DEFINED_MOLECULE_TYPES

        ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_MOLECULE_PROTOCOL_NAME
    """
    mandatory_fields = [
        "Molecule Type",
        "Type of Extraction",
        "Extraction date",
        "ProtocolName",
    ]
    for m_field in mandatory_fields:
        if m_field not in sample_batch_df:
            error_cause = (
                core.core_config.ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_MOLECULE_NOT_DEFINED.copy()
            )
            error_cause.insert(1, m_field)
            return error_cause
    if not core.models.MoleculeType.objects.filter(apps_name__exact=package).exists():
        return (
            core.core_config.ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_DEFINED_MOLECULE_TYPES
        )
    molecule_values = list(
        core.models.MoleculeType.objects.filter(apps_name__exact=package).values_list(
            "molecule_type", flat=True
        )
    )
    unique_molecule_types = sample_batch_df["Molecule Type"].unique().tolist()
    for molecule_type in unique_molecule_types:
        if molecule_type not in molecule_values:
            error_cause = (
                core.core_config.ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_MOLECULE_TYPE.copy()
            )
            error_cause.insert(1, molecule_type)
            return error_cause
    if not core.models.Protocols.objects.filter(
        type__molecule__molecule_type__exact=molecule_type,
        type__molecule__apps_name__exact=package,
    ).exists():
        return core.core_config.ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_DEFINED_PROTOCOL
    protocol_values = list(
        core.models.Protocols.objects.filter(
            type__molecule__molecule_type__exact=molecule_type,
            type__molecule__apps_name__exact=package,
        ).values_list("name", flat=True)
    )
    unique_protocol_names = sample_batch_df["ProtocolName"].unique().tolist()
    for protocol_name in unique_protocol_names:
        if protocol_name not in protocol_values:
            error_cause = (
                core.core_config.ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_MOLECULE_PROTOCOL_NAME.copy()
            )
            error_cause.insert(1, protocol_name)
            return error_cause

    return "OK"


def create_sample_from_batch_file(sample_data, reg_user, package):
    """
    Description:
        The function create new sample
    Inputs:
        sample_data     # panda dataframe
        reg_user        # user name
        package         # name of the apps that request the checking
    Functions:
        check_if_sample_already_defined     # located
    Return:
        new_sample
    """

    sample_new = {}
    # Check if sample alredy  defined for this user
    if not core.utils.samples.check_if_sample_already_defined(
        sample_data["Sample Name"], reg_user
    ):
        if sample_data["Type of Sample"] == "":
            print("Type of Sample is null for the Sample ", sample_data["Sample Name"])

    sample_new["sampleType"] = sample_data["Type of Sample"]
    sample_new["sampleName"] = sample_data["Sample Name"]

    #  Check if patient code  already exists on database, If not if will be created giving a sequencial dummy value
    if sample_data["Patient Code ID"] != "":
        patient_obj = core.utils.samples.check_patient_code_exists(
            sample_data["Patient Code ID"]
        )
        if patient_obj is False:
            # Define the new patient only Patient code is defined
            sample_new["patient"] = core.utils.samples.create_empty_patient(
                sample_data["Patient Code ID"]
            )
        else:
            sample_new["patient"] = patient_obj
    else:
        sample_new["patient"] = None

    sample_new["app_name"] = package
    sample_new["onlyRecorded"] = sample_data["Only Record"]
    sample_new["labRequest"] = sample_data["Lab requested"]
    sample_new["species"] = sample_data["Species"]
    sample_new["sampleLocation"] = sample_data["Sample storage"]
    sample_new["sampleEntryDate"] = str(sample_data["Date sample reception"])
    sample_new["sampleProject"] = core.models.SampleProjects.objects.filter(
        sample_project_name__exact=sample_data["Project/Service"], apps_name=package
    ).last()
    sample_new["user"] = reg_user
    sample_new["sample_id"] = str(reg_user + "_" + sample_data["Sample Name"])
    if not core.models.Samples.objects.exclude(unique_sample_id__isnull=True).exists():
        sample_new["new_unique_value"] = "AAA-0001"
    else:
        sample_new["new_unique_value"] = core.utils.samples.increase_unique_value(
            core.models.Samples.objects.exclude(unique_sample_id__isnull=True)
            .last()
            .get_unique_sample_id()
        )
    if sample_data["Only Record"]:
        sample_new["sampleState"] = "Completed"
    else:
        sample_new["sampleState"] = "Defined"
    sample_new["onlyRecorded"] = sample_data["Only Record"]
    new_sample = core.models.Samples.objects.create_sample(sample_new)

    return new_sample


def create_sample_project_fields_value(sample_obj, sample_data, package):
    """
    Description:
        The function set the additional parameters described in the sampleProject
    Inputs:
        sample_obj      # object of the sample
        sample_data     # panda dataframe
        package         # name of the apps that request the checking
    Functions:

    Return:
        None
    """
    s_project_field_data = {}
    s_project_field_data["sample_id"] = sample_obj
    sample_project_fields_objs = core.models.SampleProjectsFields.objects.filter(
        sample_projects_id__sample_project_name__exact=sample_data["Project/Service"],
        sample_projects_id__apps_name=package,
    )
    for sample_project_fields_obj in sample_project_fields_objs:
        s_project_field_data["sampleProjecttField_id"] = sample_project_fields_obj
        s_project_field_data["sampleProjectFieldValue"] = sample_data[
            sample_project_fields_obj.get_field_name()
        ]
        core.models.SampleProjectsFieldsValue.objects.create_project_field_value(
            s_project_field_data
        )
    return


def create_molecule_from_file(sample_obj, sample_data, reg_user, package):
    """
    Description:
        The function create new molecule from the data frame
    Inputs:
        sample_obj      # object of the sample
        sample_data     # panda dataframe
        reg_user        # user name
        package         # name of the apps that request the checking
    Return:
        new_molecule
    """
    molecule_data = {}
    molecule_data["protocolUsed"] = sample_data["ProtocolName"]
    molecule_data["moleculeType"] = sample_data["Molecule Type"]
    molecule_data["sample"] = sample_obj
    molecule_data["moleculeCodeId"] = sample_obj.get_sample_code() + "_E1"
    molecule_data["extractionType"] = sample_data["Type of Extraction"]
    molecule_data["moleculeExtractionDate"] = str(sample_data["Extraction date"])
    molecule_data["numberOfReused"] = str(0)
    molecule_data["app_name"] = package
    molecule_data["user"] = reg_user
    new_molecule = core.models.MoleculePreparation.objects.create_molecule(
        molecule_data
    )
    return new_molecule


def create_molecule_parameter_from_file(molecule_obj, sample_data):
    """
    Description:
        The function create new molecule from the data frame
    Inputs:
        molecule_obj    # object of the molecule
        sample_data     # panda dataframe
    Return:
        None
    """
    # add user commercial kit and increase usage
    user_lot_commercial_kit_obj = core.models.UserLotCommercialKits.objects.filter(
        chip_lot__exact=sample_data["Lot Commercial Kit"]
    ).last()
    molecule_obj.set_user_lot_kit_obj(user_lot_commercial_kit_obj)
    user_lot_commercial_kit_obj.set_increase_use()
    try:
        user_lot_commercial_kit_obj.set_latest_use(
            datetime.strptime(sample_data["Extraction data"], "%Y-%m-%d").date()
        )
    except Exception:
        user_lot_commercial_kit_obj.set_latest_use(datetime.datetime.now())

    protocol_used_obj = molecule_obj.get_protocol_obj()
    protocol_parameters_objs = core.models.ProtocolParameters.objects.filter(
        protocol_id__exact=protocol_used_obj, parameter_used=True
    )

    for p_parameter_obj in protocol_parameters_objs:
        molecule_parameter_value = {}
        molecule_parameter_value["moleculeParameter_id"] = p_parameter_obj
        molecule_parameter_value["molecule_id"] = molecule_obj
        molecule_parameter_value["parameterValue"] = sample_data[
            p_parameter_obj.get_parameter_name()
        ]
        core.models.MoleculeParameterValue.objects.create_molecule_parameter_value(
            molecule_parameter_value
        )
    return


def get_sample_project_fields(project_names):
    """
    """
    p_fields = []
    for p_name in project_names:
        if p_name == "" or p_name.lower() == "none":
            continue
        if core.models.SampleProjects.objects.filter(sample_project_name__iexact=p_name).exists():
            s_proj_obj = core.models.SampleProjects.objects.filter(sample_project_name__iexact=p_name).last()
            s_pro_field_objs = core.models.SampleProjectsFields.objects.filter(sample_projects_id=s_proj_obj)
            for s_pro_field_obj in s_pro_field_objs:
                p_fields.append(s_pro_field_obj.get_field_name())
        else:
            return {"ERROR":core.core_config.ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_SAMPLE_PROJECT_NO_DEFINED}
    return p_fields


def heading_validation(sample_batch_df):
    """
    """
    batch_head = list(sample_batch_df.head(0))
    for item in core.core_config.HEADING_FOR_RECORD_SAMPLES:
        if item not in batch_head:
            return {"ERROR": core.core_config.ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_INVALID_FORMAT}
    project_names = get_sample_projects_names(sample_batch_df)
    project_len = len(project_names)
    if "" in project_names or "None" in project_names:
        project_len -= 1
    if project_len > 1 :
        return {"ERROR": core.core_config.ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_TOO_MANY_PROJECTS}
    s_project_fields = get_sample_project_fields(project_names)
    if "ERROR" in s_project_fields:
        return s_project_fields
    for field in s_project_fields:
        if field not in batch_head:
            return {"ERROR": core.core_config.ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_INVALID_FORMAT}
    return "OK"

   

def read_batch_sample_file(batch_file):
    """
    Description:
        The Function read the batch file and return a panda data frame with the information
    Constants:
        ERROR_MESSAGE_FOR_EMPTY_SAMPLE_BATCH_FILE
    Return:
        sample_batch_data
    """

    sample_batch_df = pd.read_excel(batch_file, sheet_name=0)
    num_rows, _ = sample_batch_df.shape
    if num_rows == 0:
        sample_data = {}
        sample_data[
            "ERROR"
        ] = core.core_config.ERROR_MESSAGE_FOR_EMPTY_SAMPLE_BATCH_FILE
        return sample_data

    return sample_batch_df


def valid_sample_batch_file(sample_batch_df, package):
    """
    Description:
        The Function check if all parameters required in the samples are included in file
    Input:
        sample_batch_df     # sample data in dataframe
        package             # name of the apps that request the checking
    Return:
        sample_batch_data
    """

    # Check if dataframe has empty values
    validate_heading = heading_validation(sample_batch_df)
    if "ERROR" in validate_heading:
        return validate_heading
    # if sample_batch_df.isnull().values.any():
    #     return core.core_config.ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_EMPTY_VALUE
    import pdb; pdb.set_trace()
    if not check_samples_belongs_to_same_type_and_molecule_protocol(sample_batch_df):
        return (
            core.core_config.ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NOT_SAME_SAMPLE_PROTOCOL
        )
    check_opt_values = check_defined_option_values_in_samples(sample_batch_df, package)
    if check_opt_values != "OK":
        return (
            core.core_config.ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NOT_SAME_SAMPLE_PROTOCOL
        )
    # check molecule columns data
    check_molecule_par = check_molecule_has_same_data_type(sample_batch_df, package)
    if check_molecule_par != "OK":
        return check_molecule_par
    return "OK"


def save_samples_in_batch_file(sample_batch_df, reg_user, package):
    """
    Description:
        The Function save the sample and the molecule information in database
    Input:
        sample_batch_df     # sample data in dataframe
        reg_user            # user name
        package             # name of the apps that request the checking
    Functions:
        create_sample_from_batch_file   # located at this file
        create_molecule_from_file        # located at this file
        create_molecule_parameter_DNA_from_file     # located at this file
    """

    for index, row_data in sample_batch_df.iterrows():
        new_sample = create_sample_from_batch_file(row_data, reg_user, package)
        create_sample_project_fields_value(new_sample, row_data, package)
        new_molecule = create_molecule_from_file(
            new_sample, row_data, reg_user, package
        )
        create_molecule_parameter_from_file(new_molecule, row_data)
        new_sample.set_state("Library preparation")
        new_molecule.set_molecule_use("SNV/CNV", package)
        new_molecule.set_state("Completed")
    return "OK"


def read_json_schema(json_schema):
    """
    Description:
        Function reads the json and validate the schema.
    Input:
        json_schema     # File containing the schema
    Output:
        schema
    """
    schema = json.load(json_schema)
    try:
        Draft202012Validator.check_schema(schema)
    except jsonschema.ValidationError:
        return {"ERROR": core.core_config.ERROR_MESSAGE_INVALID_JSON_SCHEMA}
    return {"schema": schema}


def store_schema(schema, field, valid_fields, s_project_id):
    """ """
    if valid_fields == "":
        all_fields = True
    else:
        all_fields = False
        v_fields = valid_fields.split("\r\n")
    property_list = []
    # get the ontology value defined for sample. Do not include them in the
    # sample project list
    ont_list = []
    if core.models.OntologyMap.objects.all().exists():
        ont_list = list(
            core.models.OntologyMap.objects.values_list("ontology", flat=True)
        )
    else:
        ont_list = []
    for property in schema["properties"].keys():
        try:
            schema_dict = {"Field name": property}
            schema_dict["Description"] = schema["properties"][property]["label"]
            if field not in schema["properties"][property]:
                print(property)
                continue

            if not all_fields and schema["properties"][property][field] not in v_fields:
                continue
            # do not include the fields that are included in the main recorded
            # sample form
            if schema["properties"][property]["ontology"] in ont_list:
                continue
            if "format" in schema["properties"][property]:
                schema_dict["Field type"] = "Date"
            elif "enum" in schema["properties"][property]:
                schema_dict["Field type"] = "Options List"
                schema_dict["enum"] = []
                for item in schema["properties"][property]["enum"]:
                    enum = re.search(r"(.+) \[(.*)\]", item)
                    if enum:
                        schema_dict["enum"].append(enum.group(1))
                    else:
                        schema_dict["enum"].append(item)
            else:
                schema_dict["Field type"] = "String"
            if "classification" in schema["properties"][property]:
                schema_dict["classification"] = schema["properties"][property][
                    "classification"
                ]
            else:
                schema_dict["classification"] = None
            property_list.append(schema_dict)
        except KeyError as e:
            error = core.core_config.ERROR_FIELD_NOT_EXIST_IN_SCHEMA.copy()
            error.append(e)
            return {"ERROR": error}

    s_project_obj = core.utils.samples.get_sample_project_obj_from_id(s_project_id)
    # get classification

    counter = 0
    for property in property_list:
        # create classification if not exists
        if property["classification"] != "":
            if core.models.SampleProjectFieldClassification.objects.filter(
                sample_projects_id=s_project_obj,
                classification_name__iexact=property["classification"],
            ).exists():
                classification_obj = (
                    core.models.SampleProjectFieldClassification.objects.filter(
                        sample_projects_id=s_project_obj,
                        classification_name__iexact=property["classification"],
                    ).last()
                )
            else:
                c_data = {
                    "sample_project_id": s_project_obj,
                    "classification_name": property["classification"],
                }
                classification_obj = core.models.SampleProjectFieldClassification.objects.create_sample_project_field_classification(
                    c_data
                )
            property["SampleProjectFieldClassificationID"] = classification_obj
        else:
            property["SampleProjectFieldClassificationID"] = None
        counter += 1
        property["Used"] = True
        property["Searchable"] = False
        property["Option Values"] = ""
        property["sample_project_id"] = s_project_obj
        property["Order"] = counter
        new_s_project_prop = (
            core.models.SampleProjectsFields.objects.create_sample_project_fields(
                property
            )
        )
        if "enum" in property:
            for opt in property["enum"]:
                data = {"s_proj_obj": new_s_project_prop}
                data["opt_value"] = opt
                try:
                    core.models.SamplesProjectsTableOptions.objects.create_new_s_proj_table_opt(
                        data
                    )
                except Exception:
                    raise Exception
    return {"Success": core.core_config.SUCCESSFUL_JSON_SCHEMA}
