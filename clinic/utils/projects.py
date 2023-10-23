# Local imports
import core.models
import core.utils.patient
import core.utils.patient_projects


def add_project_fields(form_data):
    """
    Description:
        The function store the patient information that it is related to the project
    Input:
        form_data  #  information collected from form data
    Return:
        p_fields.
    """
    p_fields = {}
    patient_obj = core.utils.patient.get_patient_core_obj_from_id(
        form_data["patient_id"]
    )
    project_id = form_data["project_id"]
    fields = core.utils.patient_projects.get_project_fields(form_data["project_id"])
    field_value = {}
    field_value["patientCore_id"] = patient_obj
    field_ids = core.utils.patient_projects.get_project_field_ids(
        form_data["project_id"]
    )

    for item in range(len(field_ids)):
        field_value["projectField_id"] = core.models.PatientProjectsFields.objects.get(
            pk__exact=field_ids[item]
        )
        field_value["projectFieldValue"] = form_data[fields[item]]
        core.models.PatientProjectFieldValue.objects.create_project_field_value(
            field_value
        )
    p_fields["project_name"] = core.utils.patient_projects.get_project_obj_from_id(
        project_id
    ).get_project_name()
    return p_fields


def assign_project_patient(form_data, app_name):
    """
    Description:
        The function add the new project to the patient information.
        Returns the project fields to get the information
    Input:
        form_data  #  information collected from form data
        app_name   # application name
    Return:
        project_fields.
    """
    patient_id = form_data["patient_id"]
    patient_obj = core.utils.patient.get_patient_core_obj_from_id(patient_id)
    project_id = core.utils.patient_projects.get_project_id(
        form_data["projectName"], app_name
    )
    project_obj = core.utils.patient_projects.get_project_obj_from_id(project_id)

    # Add new project to patient
    core.utils.patient_projects.add_project_to_patient(project_obj, patient_obj)

    # get the project fields to get their values
    project_fields = core.utils.patient_projects.get_project_fields(project_id)

    return project_fields, project_id
