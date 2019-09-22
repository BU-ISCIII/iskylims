from iSkyLIMS_clinic.models import *
from iSkyLIMS_clinic.clinic_config import *

def define_patient_information (clinic_samples_ids):

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
    return patient_info


def get_clinic_sample_obj_from_id(id):
    clinic_sample_obj = ClinicSampleRequest.objects.get(pk__exact = id)
    return clinic_sample_obj
