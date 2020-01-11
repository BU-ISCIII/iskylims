from iSkyLIMS_core.models import *
from iSkyLIMS_core.core_config import *

def get_patient_core_obj_from_id (p_id):
    '''
    Description:
        The function return the patient core object from the id value
    Input:
        p_id # "    patient id
    Return:
    patient_obj
    '''
    patient_obj = PatientCore.objects.get(pk__exact = p_id)
    return patient_obj
