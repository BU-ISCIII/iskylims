import pandas as pd

from iSkyLIMS_clinic.clinic_config import *
from iSkyLIMS_clinic.models import *
from iSkyLIMS_core.models import PatientCore, PatientProjects, PatientSex
from iSkyLIMS_core.utils.patient_core import *
from iSkyLIMS_core.utils.patient_projects import *
from iSkyLIMS_core.utils.samples import *


def add_additional_information(form_data):
    """
    Description:
        The function store in patient profile the additional information
    Input:
        form_data  #  information collected from form data
    Return:
        additional_data.
    """
    patient_core_obj = get_patient_core_obj_from_id(form_data["patient_id"])
    p_opt_data = {}
    p_opt_data["patienCore"] = patient_core_obj
    for item in FORM_OPT_DATA_PATIENT_DEFINITION:
        p_opt_data[item] = form_data[item]
    opt_data_obj = PatientData.objects.create_patient_opt_data(p_opt_data)

    return opt_data_obj


def create_new_patient(form_data, app_name):
    """
    Description:
        The function gets the patient data and stores on database.
        Returns info to display back to user and project fields data
        if project was selected
    Input:
        form_data  #  information collected from form data

    Return:
        p_main_data.
    """
    # p_recorded_info = {}
    if PatientCore.objects.filter(
        patient_code__iexact=form_data["patientCode"]
    ).exists():
        return "ERROR"
    p_main_data = {}
    p_main_data["p_name"] = (
        ConfigSetting.objects.filter(configuration_name__exact="PATIENT_NAME")
        .last()
        .get_configuration_value()
    )
    p_main_data["p_surname"] = (
        ConfigSetting.objects.filter(configuration_name__exact="PATIENT_SURNAME")
        .last()
        .get_configuration_value()
    )
    for item in FORM_MAIN_DATA_PATIENT_DEFINITION:
        p_main_data[item] = form_data[item]
    new_patient_core = PatientCore.objects.create_patient(p_main_data)
    p_opt_data = {}
    p_opt_data["patienCore"] = new_patient_core
    for item in FORM_OPT_DATA_PATIENT_DEFINITION:
        p_opt_data[item] = form_data[item]

    if form_data["patientProject"] != "None":
        project_obj = get_project_obj(form_data["patientProject"], app_name)
        new_patient_core.patientProjects.add(project_obj)
        p_main_data["patient_id"] = new_patient_core.get_patient_id()
        p_main_data["project_id"] = get_project_id(
            form_data["patientProject"], app_name
        )
        p_main_data["project_name"] = project_obj.get_project_name()
        fields = get_project_fields(p_main_data["project_id"])
        if fields:
            p_main_data["fields"] = fields

    return p_main_data


def display_one_patient_info(p_id, app_name):
    """
    Description:
        The function will return the patinent information.
    Input:
        p_id  # id of the patienteCore object
    Constants:
        HEADING_FOR_DISPLAY_PATIENT_BASIC_INFORMATION
        HEADING_FOR_DISPLAY_PATIENT_ADDITIONAL_INFORMATION
    Return:
        patient_info.
    """

    if not PatientCore.objects.filter(pk__exact=p_id).exists():
        return "ERROR"
    patient_info = {}
    fields_name_in_basic_info = HEADING_FOR_DISPLAY_PATIENT_BASIC_INFORMATION.copy()
    fields_name_in_basic_info[0] = (
        ConfigSetting.objects.filter(configuration_name__exact="PATIENT_NAME")
        .last()
        .get_configuration_value()
    )
    fields_name_in_basic_info[1] = (
        ConfigSetting.objects.filter(configuration_name__exact="PATIENT_SURNAME")
        .last()
        .get_configuration_value()
    )

    patient_core_obj = get_patient_core_obj_from_id(p_id)

    patient_info["patient_name"] = []
    patient_info["patient_name"].append(patient_core_obj.get_patient_name())
    patient_info["patient_name"].append(patient_core_obj.get_patient_surname())
    p_main_info = patient_info["patient_name"][:]
    p_main_info.append(patient_core_obj.get_patient_code())
    p_main_info.append(patient_core_obj.get_patient_sex())

    patient_info["patient_basic_info"] = list(
        zip(fields_name_in_basic_info, p_main_info)
    )
    patient_info["patient_id"] = p_id
    if PatientData.objects.filter(patien_core__exact=patient_core_obj).exists():
        p_data_obj = PatientData.objects.get(patien_core__exact=patient_core_obj)
        p_data_info = p_data_obj.get_patient_full_data()
        patient_info["patient_data"] = list(
            zip(HEADING_FOR_DISPLAY_PATIENT_ADDITIONAL_INFORMATION, p_data_info)
        )

    # get project information for the patient
    if patient_core_obj.patientProjects.all().exists():
        pat_projects = patient_core_obj.patientProjects.all()
        patient_info["project_information"] = {}

        for pat_project in pat_projects:
            project_name = pat_project.get_project_name()
            patient_info["project_information"][
                project_name
            ] = get_project_field_values(pat_project.get_project_id(), patient_core_obj)

    # get other projects that patient could belongs to
    new_projects = get_available_projects_for_patient(patient_core_obj, app_name)
    if new_projects:
        patient_info["available_projects"] = []
        for new_project in new_projects:
            patient_info["available_projects"].append(new_project.get_project_name())

    # get Samples belongs to Patient
    if ClinicSampleRequest.objects.filter(patient_core=patient_core_obj).exists():
        clinic_samples = ClinicSampleRequest.objects.filter(
            patient_core=patient_core_obj
        )
        patient_info[
            "samples_heading"
        ] = HEADING_FOR_DISPLAY_SAMPLE_DATA_IN_PATIENT_INFO
        patient_info["samples_data"] = []
        for clinic_sample in clinic_samples:
            sample_obj = clinic_sample.get_core_sample_obj()
            patient_info["samples_data"].append(sample_obj.get_info_for_patient())

    return patient_info


def display_patient_list(p_list):
    """
    Description:
        The function collect basic information about patient to display
    Input:
        p_list      # list of the patient_ids
    Return:
        p_list_info
    """
    p_list_info = {}
    p_list_data = []
    p_name = (
        ConfigSetting.objects.filter(configuration_name__exact="PATIENT_NAME")
        .last()
        .get_configuration_value()
    )
    p_surname = (
        ConfigSetting.objects.filter(configuration_name__exact="PATIENT_SURNAME")
        .last()
        .get_configuration_value()
    )
    p_list_info["heading"] = ["Patient Code", p_name, p_surname]
    for p_id in p_list:
        p_obj = get_patient_core_obj_from_id(p_id)
        if p_obj:
            p_list_data.append(
                [
                    p_obj.get_patient_code(),
                    p_obj.get_patient_name(),
                    p_obj.get_patient_surname(),
                    p_id,
                ]
            )
    p_list_info["patient_data"] = p_list_data
    return p_list_info


def fields_for_new_patient(app_name):
    """
    Description:
        The function will get the option values to display in the select menus.
    Return:
        patient_definition_data.
    """
    patient_definition_data = {}
    patient_definition_data["p_name"] = (
        ConfigSetting.objects.filter(configuration_name__exact="PATIENT_NAME")
        .last()
        .get_configuration_value()
    )
    patient_definition_data["p_surname"] = (
        ConfigSetting.objects.filter(configuration_name__exact="PATIENT_SURNAME")
        .last()
        .get_configuration_value()
    )
    if PatientSex.objects.filter().exists():
        sex_objs = PatientSex.objects.filter()
        patient_definition_data["sex_values"] = []
        for sex_obj in sex_objs:
            patient_definition_data["sex_values"].append(sex_obj.get_patient_sex())
    if PatientProjects.objects.filter(apps_name__exact=app_name).exists():
        projects = PatientProjects.objects.filter(apps_name__exact=app_name)
        patient_definition_data["project_names"] = []
        for project in projects:
            patient_definition_data["project_names"].append(project.get_project_name())

    return patient_definition_data


def get_patients_in_search(data_request):
    """
    Description:
        The function will return a patinent list which matches with the input conditions.
    Input:
        data_request  # input data from the user form
    Functions:
        get_friend_list # located at core.utils.common
    Return:
        patient_list.
    """
    patient_list = []
    if not PatientCore.objects.filter().exists():
        return patient_list
    patient_objs_found = PatientCore.objects.filter()
    if data_request["p_name"] != "":
        patient_objs_found = patient_objs_found.filter(
            patient_name__icontains=data_request["p_name"]
        )
    if data_request["p_surname"] != "":
        patient_objs_found = patient_objs_found.filter(
            patient_surname__icontains=data_request["p_surname"]
        )
    if data_request["p_code"] != "":
        patient_objs_found = patient_objs_found.filter(
            patient_code__icontains=data_request["p_code"]
        )

    for patient_found in patient_objs_found:
        patient_list.append(patient_found.pk)

    return patient_list


def from_data_for_search_patient():
    """
    Description:
        The function get the information to get information to display in search patient form
    Return:
        s_patient_data
    """
    s_patient_data = {}
    s_patient_data["p_name"] = (
        ConfigSetting.objects.filter(configuration_name__exact="PATIENT_NAME")
        .last()
        .get_configuration_value()
    )
    s_patient_data["p_surname"] = (
        ConfigSetting.objects.filter(configuration_name__exact="PATIENT_SURNAME")
        .last()
        .get_configuration_value()
    )

    return s_patient_data


def read_batch_patient_file(batch_file):
    """
    Description:
        The function read the batch file and return a list of dictionnary with
        the extracted information
    Constants:
        PATIENT_BATCH_FILE_HEADING
        PATIENT_NAME
        PATIENT_SURNAME
    Return:
        patient_batch_data
    """
    patient_batch_data = []
    error_data = {}
    patient_batch_heading = PATIENT_BATCH_FILE_HEADING.copy()

    patient_batch_heading[1] = (
        ConfigSetting.objects.filter(configuration_name__exact="PATIENT_NAME")
        .last()
        .get_configuration_value()
    )
    patient_batch_heading[2] = (
        ConfigSetting.objects.filter(configuration_name__exact="PATIENT_SURNAME")
        .last()
        .get_configuration_value()
    )
    patient_batch_df = pd.read_excel(batch_file, sheet_name=0)
    num_rows, num_cols = patient_batch_df.shape

    if num_cols != 5:
        error_data["ERROR"] = ERROR_MESSAGE_FOR_INVALID_PATIENT_BATCH_FILE
        return error_data
    for field in patient_batch_heading:
        if field not in patient_batch_df.columns:
            error_data["ERROR"] = ERROR_MESSAGE_FOR_INVALID_PATIENT_BATCH_FILE
            return error_data
    if num_rows == 0:
        error_data["ERROR"] = ERROR_MESSAGE_FOR_EMPTY_PATIENT_BATCH_FILE
        return error_data
    p_projects = set(patient_batch_df["patientProjects"])
    for p_project in p_projects:
        if p_project == "" or p_project == "None":
            continue
        if not PatientProjects.objects.filter(project_name__exact=p_project).exists():
            error_data["ERROR"] = ERROR_MESSAGE_NOT_VALID_PROJECT_IN_BATCH_FILE + [
                p_project
            ]
            return error_data

    for index, row_data in patient_batch_df.iterrows():
        patient_data = {}
        for index in range(len(PATIENT_BATCH_FILE_HEADING)):
            patient_data[PATIENT_BATCH_FILE_HEADING[index]] = row_data[
                patient_batch_heading[index]
            ]
        patient_batch_data.append(patient_data)

    return patient_batch_data


def store_batch_patient(batch_data):
    """
    Description:
        The function store the new patient, including the patient project
    Input:
        batch_data      # list of patient dictionnary data
    Constants:
        PATIENT_BATCH_FILE_HEADING
        PATIENT_NAME
        PATIENT_SURNAME
    Return:
        summary
    """
    summary = {}
    not_valid_projects = []

    for patient_data in batch_data:
        new_patient = PatientCore.objects.create_patient(patient_data)

        if (
            patient_data["patientProjects"] != ""
            or patient_data["patientProjects"] != "None"
        ):
            try:
                patient_project_obj = PatientProjects.objects.get(
                    projectName__exact=patient_data["patientProjects"]
                )
                new_patient.patientProjects.add(patient_project_obj)
            except Exception:
                not_valid_projects.append(patient_data["patientProjects"])
    summary["patients"] = len(batch_data)
    if len(not_valid_projects) > 0:
        summary["no_valid_projects"] = list(set(not_valid_projects))
    return summary
