# Generic imports
import json
import re

import jsonschema
import pandas as pd
from jsonschema import Draft202012Validator

# Local imports
import core.core_config
import core.models

import core.utils.samples


def get_sample_project_fields(project_name):
    """_summary_

    Parameters
    ----------
    project_names
        Name of the project from which you want to obtain the fields

    Returns
    -------
    p_fields
        List of fields of the project
    """
    p_fields = []
    if project_name == "" or project_name.lower() == "none":
        return p_fields
    if core.models.SampleProjects.objects.filter(
        sample_project_name__iexact=project_name
    ).exists():
        s_proj_obj = core.models.SampleProjects.objects.filter(
            sample_project_name__iexact=project_name
        ).last()
        s_pro_field_objs = core.models.SampleProjectsFields.objects.filter(
            sample_projects_id=s_proj_obj
        )
        for s_pro_field_obj in s_pro_field_objs:
            p_fields.append(s_pro_field_obj.get_field_name())
    else:
        return {"ERROR": core.core_config.ERROR_SAMPLE_PROJECT_NO_DEFINED}
    return p_fields


def heading_refactor(sample_batch_data):
    """_summary_
    Description:
    ----------
        The Function read the excel file, checks if the header has proper format and changes the header to the no verbose field name

    Parameters
    ----------
    sample_batch_data
        Pandas dataframe with the batch data

    Returns
    -------
    sample_batch_data
        Pandas dataframe with the header in non verbose names, if applicable.
    """
    fields_info = {}
    for field in core.models.Samples._meta.get_fields():
        try:
            fields_info[field.name] = field.verbose_name
        except Exception:
            fields_info[field.name] = field.name

    for field in fields_info:
        if fields_info[field] in list(sample_batch_data.columns.values):
            sample_batch_data.rename(columns={fields_info[field]: field}, inplace=True)

    return sample_batch_data


def validate_header(columns, project_name):
    """_summary_

    Parameters
    ----------
    sample_batch_data
        Pandas dataframe with the batch data

    Returns
    -------
    error_cause
        If header was not validated, returns the error.
    """
    fields = {}
    heading_sample = []
    for field in core.models.Samples._meta.get_fields():
        try:
            fields[field.name] = field.verbose_name
        except Exception:
            fields[field.name] = field.name
    for item in core.core_config.HEADING_BATCH:
        heading_sample.append(fields[item])
    validate_result = {"result": "OK"}
    if project_name == "None":
        projects_fields = []
    else:
        projects_fields = get_sample_project_fields(project_name)
    right_heading = heading_sample + projects_fields
    if len(right_heading) != len(columns):
        validate_result = {"result": "NOK"}
    else:
        invalid_col_name = [i for i in right_heading if i not in columns]
        if invalid_col_name:
            validate_result = {"result": "NOK"}
    if validate_result["result"] == "NOK":
        error_message = core.core_config.ERROR_BATCH_INVALID_HEADER
        error_message += ", ".join(right_heading)
        validate_result["error_message"] = error_message
        return validate_result
    return validate_result


def check_and_format_date(sample_batch_data):
    """_summary_

    Description:
    ----------
        This function reads the batch dataframe and checks the date format is correct for date colummns

    Parameters
    ----------
    sample_batch_data
        Pandas dataframe with the batch excell information

    Returns
    -------
    error_cause
        If date was not validated, returns the error.
    """
    columns = list(sample_batch_data.columns.values)

    date_columns = []

    for value in columns:
        if "date" in value:
            date_columns.append(value)

    for column in date_columns:
        try:
            sample_batch_data[column] = pd.to_datetime(
                sample_batch_data[column]
            ).dt.strftime("%Y-%m-%d")
        except Exception:  # Unknown string format: string present at position x
            return core.core_config.ERROR_DATE_FORMAT_FIELD
    return sample_batch_data


def read_batch_sample_file(batch_file):
    """_summary_

    Description:
    ----------
        The Function reads the batch file and returns a json list of dictionaries with the information

    Parameters
    ----------
    batch_file
        Pandas dataframe with batch data

    Returns
    -------
    batch_data
        json object with the information of the batch dataframe
    """

    batch_file.fillna("", inplace=True)
    batch_data = batch_file.to_json(orient="records")
    batch_data = json.loads(batch_data)

    return batch_data


def project_validation(sample_batch_data, app_name):
    """_summary_

    Args:
        sample_batch_data (pandas dataframe): data from user input file
        app_name (string): name of the application

    Returns:
        dict: Contains the result and optional the error message or the project
        name
    """
    validation_result ={"result": "NOK"}
    if sample_batch_data["Sample Project"].isnull().values.any():
        validation_result["error_message"] = core.core_config.ERROR_EMPTY_PROJECT
        return validation_result, None

    project_list = sample_batch_data["Sample Project"].unique().tolist()

    if len(project_list) > 1:
        validation_result["error_message"] = core.core_config.ERROR_TOO_MANY_PROJECTS
        return validation_result, None

    # Check if project exist in the DB
    if "None" not in project_list:
        if not core.models.SampleProjects.objects.filter(apps_name=app_name, sample_project_name__iexact=project_list[0]).exists():
            validation_result["error_message"] = core.core_config.ERROR_NO_DEFINED_SAMPLE_PROJECTS
            return validation_result, None
    validation_result["result"] = "OK"

    return validation_result, project_list[0]

def read_json_schema(json_schema):
    """_summary_

    Description:
    ----------
        The function reads the json and validate the schema.

    Parameters
    ----------
    json_schema
        File containing the json schema

    Returns
    -------
    schema
        Json schema validation errors
    """

    schema = json.load(json_schema)
    try:
        Draft202012Validator.check_schema(schema)
    except jsonschema.ValidationError:
        return {"ERROR": core.core_config.ERROR_INVALID_JSON_SCHEMA}
    return {"schema": schema}


def store_schema(schema, field, valid_fields, s_project_id):
    """_summary_

    Parameters
    ----------
    schema
        Json schema
    field
        _description_
    valid_fields
        _description_
    s_project_id
        _description_

    Returns
    -------
        _description_

    Raises
    ------
    Exception
        _description_
    """
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
