from iSkyLIMS_clinic.models import *
from iSkyLIMS_clinic.clinic_config import *
import json

'''
Order', 'Confirmation Code', 'Priority', 'Requested Service Date',
                'History Number','Requested Service by', 'Doctor', 'Suspicious', 'Comments']
'''
def analyze_and_store_patient_data (user_post, user):

    stored_samples = []
    analyze_data = {}
    not_match = []
    not_match_samples_ids = []

    json_data = json.loads(user_post['patient_data'])
    clinic_samples = user_post['clinic_samples'].split(',')
    heading =  ADDITIONAL_HEADING_FOR_RECORDING_SAMPLES
    for c_samples_id in range(len(clinic_samples )):
        # check if patient history number matches

        history_number = json_data[c_samples_id][heading.index('History Number')]
        if history_number == '':
            continue
        patient_obj = get_patient_obj(history_number)

        if not patient_obj :
            not_match.append(json_data[c_samples_id])
            not_match_samples_ids.append(clinic_samples[c_samples_id])
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
    if not_match:
        analyze_data['heading'] = ADDITIONAL_HEADING_FOR_RECORDING_SAMPLES
        analyze_data['data'] = not_match
        analyze_data['s_request_by'] =  get_service_units()
        analyze_data['doctor'] = get_available_doctor()
    analyze_data['not_match'] = not_match
    analyze_data['not_match_samples_ids'] =not_match_samples_ids
    analyze_data['stored_samples'] = stored_samples

    return analyze_data



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
