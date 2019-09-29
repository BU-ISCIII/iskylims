
from iSkyLIMS_core.models import *
from iSkyLIMS_core.utils.handling_protocols import get_protocol_obj_from_name
from iSkyLIMS_core.core_config import * 
from django.contrib.auth.models import User


def get_commercial_kit_id(kit_name):
    if ComercialKits.objects.filter(name__exact = kit_name).exists():
        return ComercialKits.objects.get(name__exact = kit_name)
    else:
        return None

def get_comercial_kit_basic_data(kit_obj):
    kit_data = {}
    kit_data['data'] = kit_obj.get_basic_data()
    kit_data['heading'] = HEADING_FOR_COMMERCIAL_KIT_BASIC_DATA
    return  kit_data

def get_lot_comercial_kits(register_user_obj, protocol_obj):
    user_kit_list = []

    if UserLotComercialKits.objects.filter(user = register_user_obj, basedComercial__protocol_id = protocol_obj.type).exists():
        user_kits =UserLotComercialKits.objects.filter(user = register_user_obj, basedComercial__protocol_id = protocol_obj.type)
        for user_kit in user_kits:
            user_kit_list.append(user_kit.get_nick_name)

    return user_kit_list

def store_commercial_kit (kit_data):
    commercial_kit_values = {}
    commercial_kit_values['protocol_id']= get_protocol_obj_from_name(kit_data['protocol'])
    commercial_kit_values['name'] = kit_data['kitName']
    commercial_kit_values['provider'] = kit_data['provider']
    commercial_kit_values['cat_number'] = kit_data ['catNo']
    commercial_kit_values['description'] = kit_data['description']
    #commercial_kit_values['maximumUses'] = kit_data['usesNumber']
    commercial_kit_values['maximumUses'] =  0

    new_kit = ComercialKits.objects.create_commercial_kit(commercial_kit_values )
    return new_kit
