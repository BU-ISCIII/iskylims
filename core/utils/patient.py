# Local imports
import core.models

def get_patient_core_obj_from_id (p_id):
    '''
    Description:
        The function return the patient core object from the id value
    Input:
        p_id #    patient id
    Return:
        patient_obj
    '''
    try:
        patient_obj = core.models.PatientCore.objects.get(pk__exact = p_id)
    except:
        patient_obj = None
    return patient_obj
