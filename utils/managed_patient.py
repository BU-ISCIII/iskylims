from iSkyLIMS_clinic.models import *
from iSkyLIMS_clinic.clinic_config import *
from iSkyLIMS_core.models import PatientCore, PatientSex, PatientProjects


def create_new_patient(form_data):
    '''
    Description:
        The function gets the patient data and stores on database.
        Returns info to display back to user and project fields data
        if project was selected
    Input:
        request  #  information collected from form data

    Return:
        p_main_data.
    '''
    #p_recorded_info = {}
    if PatientCore.objects.filter(patientCode__iexact = form_data['patientCode']).exists():
        return 'ERROR'
    main_data = ['patientName', 'patientSurname', 'patientCode', 'patientSex']
    p_main_data = {}
    for item in main_data:
        p_main_data[item] = form_data[item]
    new_patient_core = PatientCore.objects.create_patient(p_main_data)
    opt_data = ['address', 'phone', 'birthday', 'smoker', 'notificationPreference' , 'comments']
    p_opt_data = {}
    p_opt_data['patienCore'] = new_patient_core
    for item in opt_data:
        p_opt_data[item] = form_data[item]
    if form_data['patientProject'] != None:
        p_fields = get_project_fields(form_data['patientProject'])
    #project_request =
    import pdb; pdb.set_trace()

    return p_main_data

def display_one_patient_info (p_id):
    '''
    Description:
        The function will return the patinent information.
    Input:
        p_id  # id of the patienteCore object

    Return:
        patient_info.
    '''

    if not PatientCore.objects.filter(pk__exact = p_id).exists():
        return 'ERROR'
    patient_info = {}

    patient_core_obj = get_patient_core_obj_from_id(p_id)

    patient_info['patient_name'] = []
    patient_info['patient_name'].append(patient_core_obj.get_patient_name())
    patient_info['patient_name'].append(patient_core_obj.get_patient_surname())
    p_main_info = patient_info['patient_name'][:]
    p_main_info.append(patient_core_obj.get_patient_code())

    patient_info['patient_basic_info'] = list(zip(HEADING_FOR_DISPLAY_PATIENT_BASIC_INFORMATION, p_main_info))

    if PatientData.objects.filter(patienCore__exact = patient_core_obj).exists():
        p_data_obj = PatientData.objects.get(patienCore__exact = patient_core_obj)
        p_data_info = p_data_obj.get_patient_full_data()
        patient_info['patient_data'] = list(zip(HEADING_FOR_DISPLAY_PATIENT_ADDITIONAL_INFORMATION, p_data_info))
    return patient_info

def fields_for_new_patient ():
    '''
    Description:
        The function will get the option values to display in the select menus.
    Return:
        patient_definition_data.
    '''
    patient_definition_data = {}
    if PatientSex.objects.filter().exists():
        sex_objs = PatientSex.objects.filter()
        patient_definition_data['sex_values'] = []
        for sex_obj in sex_objs:
            patient_definition_data['sex_values'].append(sex_obj.get_patient_sex())
    if PatientProjects.objects.filter().exists():
        projects =  PatientProjects.objects.filter()
        patient_definition_data['project_names'] = []
        for project in projects:
            patient_definition_data['project_names'].append(project.get_project_name())



    return patient_definition_data

def get_patients_in_search(data_request):
    '''
    Description:
        The function will return a patinent list which matches with the input conditions.
    Input:
        data_request  # input data from the user form
    Functions:
        get_friend_list # located at core.utils.generic_functions
    Return:
        patient_list.
    '''
    patient_list = []
    if not PatientCore.objects.filter().exists():
        return patient_list
    patient_objs_found = PatientCore.objects.filter()
    if data_request['p_name'] != '':
        patient_objs_found = patient_objs_found.filter(patientName__icontains = data_request['p_name'])
    if data_request['p_surname'] != '':
        patient_objs_found = patient_objs_found.filter(patientSurname__icontains = data_request['p_surname'])
    if data_request['p_code'] != '':
        patient_objs_found = patient_objs_found.filter(patientCode__icontains = data_request['p_code'])

    for patient_found in patient_objs_found :
        patient_list.append(patient_found.pk)

    return patient_list

def get_patient_core_obj_from_id (p_id):
    patient_obj = PatientCore.objects.get(pk__exact = p_id)
    return patient_obj

def get_project_fields(project_name):
    if PatientProjects.objects.filter(projectName__exact = project_name).exists():
        p_project = PatientProjects.objects.get(projectName__exact = project_name)


def set_additional_patient_data(patient_obj):

    return patient_data
