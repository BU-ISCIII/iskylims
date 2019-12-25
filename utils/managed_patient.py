from iSkyLIMS_clinic.models import *
from iSkyLIMS_clinic.clinic_config import *
from iSkyLIMS_core.models import PatientCore


def define_new_patient ():

    return

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

def set_additional_patient_data(patient_obj):

    return patient_data
