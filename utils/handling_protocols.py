import json, re
from iSkyLIMS_core.core_config import *
from iSkyLIMS_core.models import *
from django.contrib.auth.models import User



def create_new_protocol (new_protocol, protocol_type, description):
    '''
    Description:
        The function create a new item in database.
    Input:
        new_protocol # new protocol name
        protocol_type # protocol type for the new protocol
    Return:
        ID of the new created object.
    '''
    protocol_type_obj = ProtocolType.objects.get(protocol_type__exact = protocol_type)
    new_protocol_object = Protocols(type = protocol_type_obj, name = new_protocol,
                        description = description)
    new_protocol_object.save()
    return new_protocol_object.pk


def check_if_protocol_exists (protocol, app_name):
    '''
    Description:
        The function return True if protocol exists. False if not
    Input:
        protocol # "protocol name" or "protocol id" to be checked
    Return:
        True/False.
    '''
    if type(protocol) is int:
        if Protocols.objects.filter(pk__exact = protocol, type__apps_name__exact = app_name).exists():
            return True
        else:
            return False
    else:
        if Protocols.objects.filter(name__exact = protocol,  type__apps_name__exact = app_name).exists():
            return True
        else:
            return False

def define_table_for_prot_parameters(protocol_id):
    '''
    Description:
        The function return a dictionary with the information to create the table
        for defining the parameters used in the protocol
    Input:
        protocol_id # protocol id  to get protocol information
    Return:
        prot_parameters
    '''
    prot_parameters = {}
    protocol_obj = Protocols.objects.get(pk__exact = protocol_id)

    prot_parameters['protocol_name'] = protocol_obj.get_name()
    prot_parameters['protocol_id'] = protocol_id
    prot_parameters['heading'] = HEADING_FOR_DEFINING_PROTOCOL_PARAMETERS
    return prot_parameters

def display_available_protocols (app_name):
    '''
    Description:
        The function return a list with all defined protocols that contains
        molecule definition. This means to exclude any other protocol that their
        parameters are not stored in iSkyLIMS_core.
    Return:
        protocol_list.
    '''

    molecule_protocol_list = []
    if ProtocolType.objects.filter(apps_name__exact = app_name).exclude(molecule = None).exists():
        protocol_types = ProtocolType.objects.filter(apps_name__exact = app_name).exclude(molecule = None).order_by('molecule')
        for protocol_type in protocol_types :
            molecule_type = protocol_type.get_molecule_type()
            prot_type_str =  protocol_type.get_name()
            if Protocols.objects.filter(type = protocol_type).exists():
                protocols = Protocols.objects.filter(type = protocol_type).order_by('type')
                for protocol in protocols:
                    data_prot = []
                    data_prot.append(molecule_type)
                    data_prot.append(prot_type_str)
                    data_prot.append(protocol.get_name())
                    data_prot.append(protocol.pk)
                    if ProtocolParameters.objects.filter(protocol_id = protocol).exists():
                        data_prot.append(True)
                    else:
                        data_prot.append(False)
                    molecule_protocol_list.append(data_prot)
    other_protocol_list = []
    if ProtocolType.objects.filter(molecule = None, apps_name__exact = app_name).exists():
        protocol_types = ProtocolType.objects.filter(molecule = None, apps_name__exact = app_name)
        for protocol_type in protocol_types:
            prot_type_str =  protocol_type.get_name()
            if Protocols.objects.filter(type = protocol_type).exists():
                protocols = Protocols.objects.filter(type = protocol_type).order_by('type')
                for protocol in protocols:
                    data_prot = []
                    data_prot.append(prot_type_str)
                    data_prot.append(protocol.get_name())
                    data_prot.append(protocol.pk)
                    if ProtocolParameters.objects.filter(protocol_id = protocol).exists():
                        data_prot.append(True)
                    else:
                        data_prot.append(False)
                    other_protocol_list.append(data_prot)

    return molecule_protocol_list , other_protocol_list

def get_defined_protocols():
    defined_protocols = []
    if Protocols.objects.all().exists():
        protocols_obj = Protocols.objects.all().order_by('type')
        for protocol_obj in protocols_obj:
            defined_protocols.append(protocol_obj.get_name())
    return defined_protocols

def display_protocol_types (app_name):
    '''
    Description:
        The function return a list with all protocol types defined.
    Return:
        protocol_types_list.
    '''
    protocol_types_list = []
    if ProtocolType.objects.filter(apps_name__exact = app_name).exists():
        protocol_types = ProtocolType.objects.filter(apps_name__exact = app_name).order_by('molecule')
        for protocol_type in protocol_types:
            protocol_types_list.append(protocol_type.get_name())
    return protocol_types_list

def get_all_protocol_info(protocol_id):
    '''
    Description:
        The function return a dictionary with all definition parameters fro the given protocol.
    Return:
        protocol_data.
    '''
    protocol_data = {}
    protocol_data['parameters'] = []
    protocol_obj = Protocols.objects.get(pk__exact = protocol_id)

    if ProtocolParameters.objects.filter(protocol_id = protocol_obj).exists():
        protocol_data['parameter_heading'] = HEADING_FOR_DEFINING_PROTOCOL_PARAMETERS
        protocol_parameters = ProtocolParameters.objects.filter(protocol_id = protocol_obj).order_by('parameterOrder')
        for parameter in protocol_parameters:
            protocol_data['parameters'].append(parameter.get_all_parameter_info())

    return protocol_data

def get_protocol_obj_from_name(protocol_name):
    if Protocols.objects.filter(name__exact = protocol_name).exists():
        return Protocols.objects.get(name__exact = protocol_name)
    else:
        return None

def get_protocol_parameters(protocol_obj):
    '''
    Description:
        The function return a list of the used parameters .
    Return:
        protocol_parameter_list.
    '''
    protocol_parameter_list = []
    if ProtocolParameters.objects.filter(protocol_id = protocol_obj).exists():
        protocol_parameters = ProtocolParameters.objects.filter(protocol_id = protocol_obj, parameterUsed = True). order_by('parameterOrder')
        for protocol_parameter in protocol_parameters :
            protocol_parameter_list.append(protocol_parameter.get_parameter_name())
    return protocol_parameter_list

def set_protocol_parameters(request):
    protocol_id = request.POST['protocol_id']
    json_data = json.loads(request.POST['table_data1'])
    parameters = HEADING_FOR_DEFINING_PROTOCOL_PARAMETERS
    protocol_id_obj = Protocols.objects.get(pk__exact = protocol_id)

    saved_parameters = []
    stored_parameters = {}
    for row_data in json_data:

        if row_data[0] == '':
            continue
        prot_parameters = {}

        prot_parameters['protocol_id'] = protocol_id_obj
        for i in range(len(parameters)):
            prot_parameters[parameters[i]] = row_data[i]

        saved_parameters.append(ProtocolParameters.objects.create_protocol_parameter(prot_parameters).get_parameter_name())
    stored_parameters['parameters'] = saved_parameters
    stored_parameters['protocol_name'] = protocol_id_obj.get_name()

    return stored_parameters
