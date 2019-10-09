
from iSkyLIMS_core.models import *
from iSkyLIMS_core.utils.handling_protocols import get_protocol_obj_from_name
from iSkyLIMS_core.core_config import *
from django.contrib.auth.models import User


def get_commercial_kit_id(kit_name):
    if CommercialKits.objects.filter(name__exact = kit_name).exists():
        return CommercialKits.objects.get(name__exact = kit_name)
    else:
        return None

def get_defined_commercial_kits():
    commercial_kit_list = []
    if CommercialKits.objects.exists():
        kits = CommercialKits.objects.all().order_by('protocol_id')
        for kit in kits:
            commercial_kit_list.append(kit.get_name())
    return commercial_kit_list

def get_lot_user_commercial_kit_id(nick_name):
    if UserLotCommercialKits.objects.filter(nickName__exact = nick_name).exists():
        return UserLotCommercialKits.objects.get(nickName__exact = nick_name)
    else:
        return None

def get_commercial_kit_obj_from_name(kit_name):
    if CommercialKits.objects.filter(name__exact = kit_name).exists():
        return CommercialKits.objects.get(name__exact = kit_name)
    return None

def get_data_for_commercial_kits():
    data_kits = []
    data_commercial_kits = {}
    if CommercialKits.objects.exists():
        kits = CommercialKits.objects.all().order_by('protocol_id')
        for kit in kits:
            data_kits.append(kit.get_name())
        data_commercial_kits['heading'] = HEADING_FOR_COMMERCIAL_KIT_BASIC_DATA
        data_commercial_kits['c_kits_data'] = data_kits
    return data_commercial_kits

def get_commercial_kit_basic_data(kit_obj):
    kit_data = {}
    kit_data['data'] = kit_obj.get_basic_data()
    kit_data['heading'] = HEADING_FOR_COMMERCIAL_KIT_BASIC_DATA
    return  kit_data

def get_lot_user_commercial_kit_basic_data(kit_obj):
    lot_kit_data = {}
    lot_kit_data['data'] = kit_obj.get_basic_data()
    lot_kit_data['heading'] = HEADING_FOR_LOT_USER_COMMERCIAL_KIT_BASIC_DATA
    return  lot_kit_data

def get_lot_commercial_kits(register_user_obj, protocol_obj):
    user_kit_list = []

    if UserLotCommercialKits.objects.filter(user = register_user_obj, basedCommercial__protocol_id = protocol_obj).exists():
        user_kits =UserLotCommercialKits.objects.filter(user = register_user_obj, basedCommercial__protocol_id = protocol_obj)
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

    new_kit = CommercialKits.objects.create_commercial_kit(commercial_kit_values )
    return new_kit

def store_lot_user_commercial_kit (kit_data, user_name):
    commercial_kit_obj = get_commercial_kit_obj_from_name(kit_data['commercialKit'])
    lot_kit_values = {}
    lot_kit_values['user'] = user_name
    lot_kit_values['basedCommercial']= commercial_kit_obj
    lot_kit_values['nickName'] = kit_data['nickName']
    lot_kit_values['chipLot'] = kit_data['txtCode']
    lot_kit_values['expirationDate'] = kit_data ['expirationDate']
    lot_kit_values['maximumUses'] = commercial_kit_obj.get_maximum_uses()


    new_kit = UserLotCommercialKits.objects.create_user_lot_commercial_kit(lot_kit_values )
    return new_kit
