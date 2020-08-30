import json
from iSkyLIMS_wetlab.models import *
from iSkyLIMS_core.models import Protocols, CommercialKits
from django.contrib.auth.models import User
from iSkyLIMS_core.utils.handling_protocols import get_protocol_obj_from_id
from iSkyLIMS_wetlab.wetlab_config import *



def define_table_for_additional_kits(protocol_id):
    '''
    Description:
        The function get the commercial protocols defined for the library preparation
        protocol.
    Input:
        protocol_id__exact  # Id of the protocol
    Constant:
        HEADING_ADDING_COMMERCIAL_KITS_TO_PROTOCOL
    Return:
        kit_data
    '''
    kit_data = {}
    protocol_obj = get_protocol_obj_from_id(protocol_id)
    if CommercialKits.objects.filter(protocolKits = protocol_obj).exists():
        c_kits = CommercialKits.objects.filter(protocolKits = protocol_obj)
        c_kit_names = []
        for c_kit in c_kits:
            c_kit_names.append(c_kit.get_name())
        kit_data['kits'] = c_kit_names
        kit_data['protocol_name'] = protocol_obj.get_name()
        kit_data['protocol_id'] = protocol_id
        kit_data['heading'] = HEADING_ADDING_COMMERCIAL_KITS_TO_PROTOCOL
    return kit_data

def get_additional_kits_list (app_name):
    '''
    Description:
        The function get protocol list with information if additional kits are:
        Defined; Pending; or Not requested
    Return:
        data_prot with "protocol_name" "protocol_id" and status of the
        additional index
    '''
    additional_kits = []
    if ProtocolType.objects.filter(molecule = None, apps_name__exact = app_name).exists():
        protocol_types = ProtocolType.objects.filter(molecule = None, apps_name__exact = app_name)
        for protocol_type in protocol_types:
            prot_type_str =  protocol_type.get_name()
            if Protocols.objects.filter(type = protocol_type).exists():
                protocols = Protocols.objects.filter(type = protocol_type).order_by('type')
                for protocol in protocols:
                    data_prot = []
                    data_prot.append(protocol.get_name())
                    data_prot.append(protocol.pk)
                    if AdditionaKitsLibraryPreparation.objects.filter( protocol_id = protocol).exists():
                        data_prot.append('Defined')
                    else :
                        data_prot.append('')
                    additional_kits.append(data_prot)
    return additional_kits

def set_additional_kits (form_data, user):
    '''
    Description:
        The function get additional kits used in the library preparation and save
        them in database
    Input:
        form_data   # user form data
        user        # requested user
    Return:
        kit_names
    '''
    kit_names = []
    protocol_id = form_data['protocol_id']
    kit_heading_names = ['kitName','kitOrder', 'kitUsed', 'commercial_kit', 'description']
    json_data = json.loads(form_data['kits_data'])
    for row_index in range(len(json_data)):
        if json_data[row_index][0] == '':
            continue
        kit_data = {}
        for i in range(len(kit_heading_names)):
            kit_data[kit_heading_names[i]] = json_data[row_index][i]
        kit_data['protocol_id'] = get_protocol_obj_from_id(protocol_id)
        kit_data['commercialKit_id'] = CommercialKits.objects.filter(name__exact = kit_data['commercial_kit']).last()
        kit_data['user'] = user
        new_additional_kit = AdditionaKitsLibraryPreparation.objects.create_additional_kit(kit_data)
        kit_names.append(kit_data['kitName'])
    return kit_names
