from iSkyLIMS_clinic.models import *
from iSkyLIMS_clinic.clinic_config import *
from iSkyLIMS_clinic.utils.generic_functions import *

import json

'''
Order', 'Confirmation Code', 'Priority', 'Requested Service Date',
                'History Number','Requested Service by', 'Doctor', 'Suspicious', 'Comments']
'''
def analyze_and_store_patient_data (user_post, user):

    stored_samples = []
    analyze_data = {}
    not_match = []
    #not_match_samples_ids = []
    incomplete_clinic_samples = []
    incomplete_clinic_samples_ids = []
    json_data = json.loads(user_post['patient_data'])
    clinic_samples = user_post['clinic_samples'].split(',')
    heading =  ADDITIONAL_HEADING_FOR_RECORDING_SAMPLES
    for c_samples_id in range(len(clinic_samples )):
        # check if patient history number matches

        if check_empty_fields(json_data[c_samples_id]):
            incomplete_clinic_samples.append(json_data[c_samples_id])
            incomplete_clinic_samples_ids.append(clinic_samples[c_samples_id])
            continue

        history_number = json_data[c_samples_id][heading.index('History Number')]

        patient_obj = get_patient_obj(history_number)

        if not patient_obj :
            not_match.append(json_data[c_samples_id][0:6])
            incomplete_clinic_samples.append(json_data[c_samples_id])
            incomplete_clinic_samples_ids.append(clinic_samples[c_samples_id])
            continue
        else:
            patient_data = {}
            patient_data['patient_id'] = patient_obj
            for map_column in  MAP_ADDITIONAL_HEADING_TO_DATABASE:

                patient_data[map_column[1]] = json_data[c_samples_id][heading.index(map_column[0])]
                #match_ids.append(clinic_samples[c_samples_id])
            patient_data['doctor_id'] = get_doctor_obj(json_data[c_samples_id][heading.index('Doctor')])
            patient_data['serviceUnit_id'] = get_service_unit_obj((json_data[c_samples_id][heading.index('Requested Service by')]))

            clinic_obj = get_clinic_sample_obj_from_id(clinic_samples[c_samples_id])
            clinic_obj.update(patient_data)
            stored_samples.append([json_data[c_samples_id][0], history_number , json_data[c_samples_id][heading.index('Requested Service by')]])
            # create suspicious data for this clinic sample
            suspicious_data = {}
            suspicious_data['patient_id'] = patient_obj
            suspicious_data['clinicSample_id'] = clinic_obj
            suspicious_data['description'] = json_data[c_samples_id][heading.index('Suspicious')]
            suspicious_obj = SuspiciousHistory.objects.create_suspicious_history(suspicious_data)
    if stored_samples :
        analyze_data['stored_samples_heading'] = HEADING_FOR_STORED_PATIENT_DATA
    if incomplete_clinic_samples_ids:
        analyze_data['heading'] = ADDITIONAL_HEADING_FOR_RECORDING_SAMPLES
        analyze_data['heading_length'] = len(ADDITIONAL_HEADING_FOR_RECORDING_SAMPLES)
        analyze_data['incomplete_clinic_samples'] = incomplete_clinic_samples
        analyze_data['data_length'] = len(incomplete_clinic_samples)
        analyze_data['incomplete_clinic_samples_ids'] = ','.join(incomplete_clinic_samples_ids)
        #analyze_data['data'] = not_match
        analyze_data['s_request_by'] =  get_service_units()
        analyze_data['doctor'] = get_available_doctor()


    if not_match :
        analyze_data['not_match'] = not_match
        analyze_data['heading_not_match'] = HEADING_FOR_NOT_MATCH
    #analyze_data['not_match_samples_ids'] =not_match_samples_ids
    if stored_samples :
        analyze_data['stored_samples'] = stored_samples

    return analyze_data

def check_if_sample_c_exists(sample_c_id):
    if  ClinicSampleRequest.objects.filter(pk__exact = sample_c_id).exists():
        return True
    else:
        return False

def display_one_sample_info(id):
    sample_info = {}
    sample_obj = get_clinic_sample_obj_from_id(id)
    sample_info['s_name'] = sample_obj.get_sample_name()
    p_info = sample_obj.get_patient_information()
    sample_info['patient_info'] = list(zip(HEADING_FOR_DISPLAY_PATIENT_INFORMATION, p_info))
    r_by_info = sample_obj.get_requested_by_information()
    sample_info['requested_by'] = list(zip(HEADING_FOR_DISPLAY_REQUESTED_BY_INFORMATION, r_by_info))
    sample_info['sample_core_info'] = sample_obj.get_sample_core_info()
    sample_info['sample_core_heading'] = HEADING_FOR_DISPLAY_SAMPLE_CORE_INFORMATION
    if SuspiciousHistory.objects.filter(clinicSample_id__exact = sample_obj).exists():
        sample_info['suspicious'] = SuspiciousHistory.objects.get(clinicSample_id__exact = sample_obj).get_suspicious_text()
    else:
        sample_info['suspicious'] = 'Information not available'
    sample_info['comments'] = sample_obj.get_comments()
    return sample_info



def display_sample_list(sample_c_list):
    display_sample_list_info = {}
    display_sample_list_info['heading'] = HEADING_SEARCH_LIST_SAMPLES_CLINIC
    sample_c_data = []
    for sample_c in sample_c_list:
        sample_c_obj = ClinicSampleRequest.objects.get(pk__exact = sample_c)
        sample_c_data.append(sample_c_obj.get_sample_info_for_list())
        #import pdb; pdb.set_trace()
    display_sample_list_info['sample_c_data'] = sample_c_data
    return display_sample_list_info

def get_clinic_sample_obj_from_id(id):
    clinic_sample_obj = ClinicSampleRequest.objects.get(pk__exact = id)
    return clinic_sample_obj

def get_clinic_samples_by_state(state):

    if ClinicSampleRequest.objects.filter(clinicSampleState__clinicState__exact = state).exists():
        c_samples  =  ClinicSampleRequest.objects.filter(clinicSampleState__clinicState__exact = state)
        c_samples_ids = []
        for c_sample in c_samples:
            c_samples_ids.append(c_sample.get_id())
        return c_samples_ids
    else :
        return None

def get_available_doctor ():
    doctor_list = []
    if Doctor.objects.all().exists():
        doctors = Doctor.objects.all().order_by('doctorName')
        for doctor in doctors:
            doctor_list.append(doctor.get_name())
    return doctor_list

def get_doctor_obj (doctor_name):
    if Doctor.objects.filter(doctorName__exact = doctor_name).exists():
        return Doctor.objects.get(doctorName__exact = doctor_name)
    else:
        return None

def get_patient_obj (history_number):
    if Patient.objects.filter(numberOfHistory__exact = history_number).exists():
        return Patient.objects.get(numberOfHistory__exact = history_number)
    else:
        return None
def get_service_unit_obj(unit_name):
    if ServiceUnits.objects.filter(serviceUnitName__exact = unit_name).exists():
        return ServiceUnits.objects.get(serviceUnitName__exact = unit_name)
    else:
        return None

def get_service_units():
    service_unit_list = []
    if ServiceUnits.objects.all().exists():
        services = ServiceUnits.objects.all().order_by('serviceUnitName')
        for service in services:
            service_unit_list.append(service.get_name())
    return service_unit_list

def get_samples_clinic_in_search (data_request):
    clinic_s_list = []
    if ClinicSampleRequest.objects.all().exists():
        clinic_s_found = ClinicSampleRequest.objects.all()
    else:
        return clinic_s_list
    if data_request['history_number'] != '':
        clinic_s_found = clinic_s_found.filter(patient_id__numberOfHistory__icontains = data_request['history_number'])

    if data_request['patient_name'] != '':
        clinic_s_found = clinic_s_found.filter(patient_id__patientName__icontains = data_request['patient_name'])
    if data_request['sample_name'] != '':
        clinic_s_found = clinic_s_found.filter(sampleCore__sampleName__icontains = data_request['sample_name'])
        if len(clinic_s_found) == 1:
                clinic_s_list.append(clinic_s_found[0].pk)
                return clinic_s_list
    if data_request['doctor_name'] != '':
        clinic_s_found = clinic_s_found.filter(doctor_id__doctorName__exact = data_request['doctor_name'])
    if data_request['requested_service_by'] != '':
        clinic_s_found = clinic_s_found.filter(serviceUnit_id__serviceUnitName__exact = data_request['requested_service_by'])

    if data_request['start_date'] !='' and data_request['end_date'] != '':
        clinic_s_found = clinic_s_found.filter(generated_at___range=(data_request['start_date'], data_request['end_date'] ))

    if data_request['start_date'] !='' and data_request['end_date']  == '':
        clinic_s_found = clinic_s_found.filter(generated_at__gte = data_request['start_date'])

    if data_request['start_date'] =='' and data_request['end_date']  != '':
            clinic_s_found = clinic_s_found.filter(run_date__lte = data_request['end_date'] )

    for clinic_s in clinic_s_found :
        clinic_s_list.append(clinic_s.pk)

    return clinic_s_list

def prepare_patient_form (clinic_samples_ids):

    patient_info = {}
    patient_info['data'] = []
    patient_info['heading'] = ADDITIONAL_HEADING_FOR_RECORDING_SAMPLES
    heading_length = len(ADDITIONAL_HEADING_FOR_RECORDING_SAMPLES)
    for clinic_s_id in clinic_samples_ids:

        clinic_sample_obj = get_clinic_sample_obj_from_id(clinic_s_id)
        data_sample = ['']*heading_length
        data_sample[0] = clinic_sample_obj.get_sample_name()
        patient_info['data'].append(data_sample)
    patient_info['clinic_samples_ids'] = ','.join(clinic_samples_ids)
    patient_info['heading_length'] = heading_length
    patient_info['data_length'] = len(patient_info['data'])
    patient_info['doctor'] = get_available_doctor()
    patient_info['s_request_by'] = get_service_units()
    patient_info['clinic_samples'] = ','.join(clinic_samples_ids)
    return patient_info
