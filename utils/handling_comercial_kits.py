
from iSkyLIMS_core.models import *
from django.contrib.auth.models import User

def get_lot_comercial_kits(register_user_obj, protocol_obj):
    user_kit_list = []

    if UserComercialKits.objects.filter(user = register_user_obj, basedComercial__protocol_id = protocol_obj.type).exists():
        user_kits =UserComercialKits.objects.filter(user = register_user_obj, basedComercial__protocol_id = protocol_obj.type)
        for user_kit in user_kits:
            user_kit_list.append(user_kit.get_nick_name)

    return user_kit_list
