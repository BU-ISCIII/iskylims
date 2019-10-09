from iSkyLIMS_clinic.models import *
from iSkyLIMS_clinic.clinic_config import *
import json




def define_table_for_result_parameters():
    '''
    Description:
        The function return a dictionary with the information to create the table
        for defining the result parameters used
    Input:
        protocol_id # protocol id  to get protocol information
    Return:
        prot_parameters
    '''
    result_parameters = {}
    #protocol_obj = Protocols.objects.get(pk__exact = protocol_id)

    #prot_parameters['protocol_name'] = protocol_obj.get_name()
    #prot_parameters['protocol_id'] = protocol_id
    result_parameters['heading'] = HEADING_FOR_DEFINING_RESULT_PARAMETERS
    return result_parameters



def set_result_parameters(request):
    #protocol_id = request.POST['protocol_id']
    json_data = json.loads(request.POST['table_data1'])
    parameters = HEADING_FOR_DEFINING_RESULT_PARAMETERS
    #protocol_id_obj = Protocols.objects.get(pk__exact = protocol_id)

    saved_parameters = []
    stored_parameters = {}
    for row_data in json_data:

        if row_data[0] == '':
            continue
        result_parameters = {}

        #prot_parameters['protocol_id'] = protocol_id_obj
        for i in range(len(parameters)):
            result_parameters[parameters[i]] = row_data[i]

        saved_parameters.append(ResultParameters.objects.create_result_parameter(result_parameters).get_parameter_name())
    stored_parameters['parameters'] = saved_parameters
    #stored_parameters['protocol_name'] = protocol_id_obj.get_name()

    return stored_parameters
