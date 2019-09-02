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


def check_if_protocol_exists (protocol_name):
    '''
    Description:
        The function return True if protocol exists. False if not
    Input:
        protocol_name # protocol name to be checked
    Return:
        True/False.
    '''
    if Protocols.objects.filter(name__exact = protocol_name).exists():
        return True
    else:
        return False

def display_availble_protocols ():
    '''
    Description:
        The function return a list with all defined protocols that contains
        molecule definition. This means to exclude any other protocol that their
        parameters are not stored in iSkyLIMS_core.
    Return:
        protocol_list.
    '''

    protocol_list = []
    if ProtocolType.objects.all().exclude(molecule = None).exists():
        protocol_types = ProtocolType.objects.all().exclude(molecule = None).order_by('molecule')
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
                    protocol_list.append(data_prot)

    return protocol_list


def display_protocol_types ():
    '''
    Description:
        The function return a list with all protocol types defined.
    Return:
        protocol_types_list.
    '''
    protocol_types_list = []
    if ProtocolType.objects.all().exists():
        protocol_types = ProtocolType.objects.all().order_by('molecule')
        for protocol_type in protocol_types:
            protocol_types_list.append(protocol_type.get_name())
    return protocol_types_list
