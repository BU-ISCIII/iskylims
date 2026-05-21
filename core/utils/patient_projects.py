# Generic imports
import json

# Local imports
import core.models
import core.core_config


def add_project_to_patient(project_obj, patient_obj):
    """
    Description:
        The function add the project to the patient
    Input:
        project_obj  # New project object to be added
        patient_obj  #
    Return:
        None.
    """
    patient_obj.patientProjects.add(project_obj)
    return None


def check_if_project_exists(project_id, app_name):
    """
    Description:
        The function return True if protocol exists. False if not
    Input:
        protocol # "protocol name" or "protocol id" to be checked
    Return:
        True/False.
    """
    if core.models.PatientProjects.objects.filter(
        pk__exact=project_id, apps_name__exact=app_name
    ).exists():
        return True
    else:
        return False


def create_patient_project(form_data, app_name):
    """
    Description:
        The function gets the project definition  data and stores on database.
        Returns info to display back to user the project name and project id
    Input:
        form_data  #  information collected from form data

    Return:
        project_data.
    """
    if core.models.PatientProjects.objects.filter(
        project_name__iexact=form_data["projectName"], apps_name__exact=app_name
    ).exists():
        return "ERROR"
    project_data = {}
    for item in core.core_config.FORM_PROJECT_CREATION:
        project_data[item] = form_data[item]
    project_data["apps_name"] = app_name
    new_patient_core = core.models.PatientProjects.objects.create_project(project_data)
    project_data["project_id"] = new_patient_core.get_project_id()
    project_data["project_name"] = new_patient_core.get_project_name()
    return project_data


def get_all_project_info(proyect_id):
    """
    Description:
        The function return a dictionary with all definition fields fro the given project
    Return:
        project_data.
    """
    project_data = {}
    project_data["fields"] = []
    project_obj = core.models.PatientProjects.objects.get(pk__exact=proyect_id)

    if core.models.PatientProjectsFields.objects.filter(
        patient_projects_id=project_obj
    ).exists():
        project_data["field_heading"] = (
            core.core_config.HEADING_FOR_DEFINING_PROJECT_FIELDS
        )
        project_data["project_name"] = project_obj.get_project_name()
        project_fields = core.models.PatientProjectsFields.objects.filter(
            patient_projects_id=project_obj
        ).order_by("project_field_order")
        for field in project_fields:
            project_data["fields"].append(field.get_all_fields_info())

    return project_data


def get_available_projects_for_patient(patient_obj, app_name):
    """
    Description:
        The function gets the projects that are not defined for the patient but could be,
        because they are defined for the application
        Returns a list with the available project
    Input:
        patient_obj  # Patient class
        app_name     # application name to filter for the available projects
    Return:
        project_list.
    """
    if (
        core.models.PatientProjects.objects.filter(apps_name__exact=app_name)
        .exclude(id__in=patient_obj.patient_projects.values_list("id", flat=True))
        .exists()
    ):
        available_projects = core.models.PatientProjects.objects.filter(
            apps_name__exact=app_name
        ).exclude(id__in=patient_obj.patient_projects.values_list("id", flat=True))
        return available_projects
    else:
        return None


def get_defined_patient_projects(app_name):
    """
    Description:
        The function will return the patient projects that a patient could belongs to.
    Variables:
        patient_projects # list containing the patient projects defined in database
    Return:
        patient_projects.
    """
    patient_projects = []
    if core.models.PatientProjects.objects.filter(apps_name__exact=app_name).exists():
        pat_projects = core.models.PatientProjects.objects.filter(
            apps_name__exact=app_name
        )
        for pat_project in pat_projects:
            data = pat_project.get_patient_project_data()
            if core.models.PatientProjectsFields.objects.filter(
                patient_projects_id=pat_project
            ).exists():
                data.append(True)
            else:
                data.append(False)
            patient_projects.append(data)
    return patient_projects


def get_project_id(project_name, app_name):
    """
    Description:
        The function return the project id from the input project name
    Input:
        project_name    # project name
        app_name       # application that created the project
    Return:
        project_id
    """
    if core.models.PatientProjects.objects.filter(
        project_name__exact=project_name, apps_name__exact=app_name
    ):
        project_id = core.models.PatientProjects.objects.get(
            project_name__exact=project_name, apps_name__exact=app_name
        ).get_project_id()
        return project_id
    else:
        return None


def get_project_obj_from_id(project_id):
    """
    Description:
        The function return the project obj from the input project id
    Input:
        project_id    # project id

    Return:
        project_obj
    """
    if core.models.PatientProjects.objects.filter(pk__exact=project_id).exists():
        project_obj = core.models.PatientProjects.objects.get(pk__exact=project_id)
        return project_obj
    else:
        return None


def get_project_fields(project_id):
    """
    Description:
        The function gets the fields that user was defined in the project.
        Returns a list of the fields requested for the project
    Return:
        project_fields.
    """
    project_fields = []

    if core.models.PatientProjectsFields.objects.filter(
        patient_projects_id__pk=project_id
    ).exists():
        p_fields_obj = core.models.PatientProjectsFields.objects.filter(
            patient_projects_id__pk=project_id
        ).order_by("project_field_order")
        for p_field in p_fields_obj:
            project_fields.append(p_field.get_field_name())
    return project_fields


def get_project_field_ids(project_id):
    """
    Description:
        The function gets the fields id tha user was defined in the project.
        Returns a list of the field ids requested for the project
    Return:
        project_fields.
    """
    project_field_ids = []

    if core.models.PatientProjectsFields.objects.filter(
        patient_projects_id__pk=project_id
    ).exists():
        p_fields_obj = core.models.PatientProjectsFields.objects.filter(
            patient_projects_id__pk=project_id
        ).order_by("project_field_order")
        for p_field in p_fields_obj:
            project_field_ids.append(p_field.get_field_id())
    return project_field_ids


def get_project_field_values(project_id, patient_obj):
    """
    Description:
        The function gets the fields id tha user was defined in the project.
        Returns a list of the field ids requested for the project
    Return:
        project_field_values.
    """
    project_field_values = []

    p_fields_obj = core.models.PatientProjectsFields.objects.filter(
        patient_projects_id__pk=project_id
    ).order_by("project_field_order")
    for p_field in p_fields_obj:
        if core.models.PatientProjectFieldValue.objects.filter(
            patient_core_id=patient_obj, project_field_id=p_field
        ).exists():
            field_value = core.models.PatientProjectFieldValue.objects.get(
                patient_core_id=patient_obj, project_field_id=p_field
            ).get_field_value()
        else:
            field_value = "Not defined"
        project_field_values.append((p_field.get_field_name(), field_value))
    return project_field_values


def get_project_obj(project_name, app_name):
    """
    Description:
        The function return the project instance from the input project name
    Input:
        project_name    # project name
        app_name       # application that created the project
    Return:
        project_obj
    """
    if core.models.PatientProjects.objects.filter(
        project_name__exact=project_name, apps_name__exact=app_name
    ):
        project_obj = core.models.PatientProjects.objects.get(
            project_name__exact=project_name, apps_name__exact=app_name
        )
        return project_obj
    else:
        return None


def define_table_for_project_fields(project_id):
    """
    Description:
        The function return a dictionary with the information to create the table
        for defining the fields used in patient protocol
    Input:
        project_id # project id  to get project name
    Return:
        project_fields
    """
    project_fields = {}
    project = core.models.PatientProjects.objects.get(pk__exact=project_id)

    project_fields["project_name"] = project.get_project_name()
    project_fields["project_id"] = project_id
    project_fields["heading"] = core.core_config.HEADING_FOR_DEFINING_PROJECT_FIELDS
    return project_fields


def set_project_fields(data_form):
    project_id = data_form["project_id"]
    json_data = json.loads(data_form["table_data1"])
    fields = core.core_config.HEADING_FOR_DEFINING_PROJECT_FIELDS
    project_id_obj = core.models.PatientProjects.objects.get(pk__exact=project_id)
    saved_fields = []
    stored_fields = {}
    for row_data in json_data:
        if row_data[0] == "":
            continue
        project_fields = {}

        project_fields["project_id"] = project_id_obj
        for i in range(len(fields)):
            project_fields[fields[i]] = row_data[i]
        new_p_field_obj = (
            core.models.PatientProjectsFields.objects.create_project_fields(
                project_fields
            )
        )
        saved_fields.append(
            [new_p_field_obj.get_field_name(), new_p_field_obj.get_description()]
        )
    stored_fields["fields"] = saved_fields
    stored_fields["project_name"] = project_id_obj.get_project_name()

    return stored_fields
