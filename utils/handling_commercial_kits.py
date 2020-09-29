
from iSkyLIMS_core.models import *
from iSkyLIMS_core.utils.handling_protocols import get_protocol_obj_from_name
from iSkyLIMS_core.core_config import *
from django.contrib.auth.models import User
from datetime import date

def get_commercial_kit_id(kit_name):
    if CommercialKits.objects.filter(name__iexact = kit_name).exists():
        return CommercialKits.objects.get(name__iexact = kit_name)
    else:
        return None

def get_defined_commercial_kits():
    commercial_kit_list = []
    if CommercialKits.objects.exists():
        kits = CommercialKits.objects.all().order_by('name')
        for kit in kits:
            commercial_kit_list.append(kit.get_name())
    return commercial_kit_list

def get_lot_user_commercial_kit_id(lot_number):
    if UserLotCommercialKits.objects.filter(chipLot__iexact = lot_number).exists():
        return UserLotCommercialKits.objects.get(chipLot__iexact = lot_number)
    else:
        return None

def get_commercial_kit_obj_from_name(kit_name):
    if CommercialKits.objects.filter(name__exact = kit_name).exists():
        return CommercialKits.objects.get(name__exact = kit_name)
    return None

def get_data_for_commercial_kits(platform):
    '''
    Description:
        The function get the commercial kit data. If platform is included it returns
        the commercial kits defined for platform
    Input:
        platform    # platform name
    Return:
        data_commercial_kits
    '''

    data_commercial_kits = {}
    data_commercial_kits['protocol'] = {}
    data_commercial_kits['protocol']['data'] = {}
    if CommercialKits.objects.all().exclude(protocolKits = None).exists():
        kits = CommercialKits.objects.all().exclude(protocolKits = None).order_by('name')
        for kit in kits:
            data_kits = []
            commercial_kit_name = kit.get_name()

            protocol_objs = kit.protocolKits.all()
            protocols = []
            for protocol_obj in protocol_objs:
                protocols.append(protocol_obj.get_name())

            data_kits.append(protocols)
            #data_kits.append(kit.get_name())
            data_kits.append(kit.get_provider_kit_name())
            data_kits.append(kit.get_cat_number())

            #if not protocol in data_commercial_kits['data']:
            #   data_commercial_kits['data'][protocols] = []
            data_commercial_kits['protocol']['data'][commercial_kit_name] = [data_kits]
        data_commercial_kits['protocol']['headings'] = HEADING_FOR_COMMERCIAL_PROTOCOL_KIT_BASIC_DATA
    # get the platform CommercialKits if platorm is not empty
    if platform != '':
        if CommercialKits.objects.all().exclude(platformKits = None).exists():
            kits = CommercialKits.objects.all().exclude(platformKits = None).order_by('name')
            data_commercial_kits['platform'] = {}
            data_commercial_kits['platform']['data'] = {}
            for kit in kits:
                data_kits = []
                commercial_kit_name = kit.get_name()
                data_kits.append(kit.get_platform_name())
                data_kits.append(kit.get_provider_kit_name())
                data_kits.append(kit.get_cat_number())
                data_commercial_kits['platform']['data'][commercial_kit_name] = [data_kits]
            data_commercial_kits['platform']['headings'] = HEADING_FOR_COMMERCIAL_PLATFORM_KIT_BASIC_DATA

    return data_commercial_kits

def get_commercial_kit_basic_data(kit_obj):
    kit_data = {}
    if kit_obj.platform_kit_obj():

        kit_data['data'] = kit_obj.get_commercial_platform_basic_data()
        kit_data['heading'] = HEADING_FOR_NEW_SAVED_COMMERCIAL_PLATFORM_KIT
    else:
        kit_data['data'] = kit_obj.get_commercial_protocol_basic_data()
        kit_data['heading'] = HEADING_FOR_NEW_SAVED_COMMERCIAL_PROTOCOL_KIT
        kit_data['protocol_kit'] = True
    return  kit_data


def display_user_lot_kit_information_from_query_list(user_kits_objs):
    '''
    Description:
        The function gets the user kits objects list and return a list with basic
        information
    Input:
        user_kits_objs  # list of the user kit objects
    Constant:
        HEADING_FOR_USER_LOT_SEARCH_RESULTS
    Return:
        user_lot
    '''
    user_lot = {}
    user_lot['kit_data'] = []
    user_lot['heading'] = HEADING_FOR_USER_LOT_SEARCH_RESULTS

    for user_kit_obj in user_kits_objs:
        data = user_kit_obj.get_basic_data()

        commercial_kit_obj = user_kit_obj.get_commercial_obj()
        protocols = commercial_kit_obj. get_protocol_objs()
        # protocols = commercial_kit_obj.protocolKits.all()
        protocols_name = []
        for protocol in protocols:
            protocols_name.append(protocol.get_name())

        data.append(protocols_name)
        data.append(user_kit_obj.get_user_lot_kit_id())
        user_lot['kit_data'].append(data)
    return user_lot


def get_expired_lot_user_kit (register_user_obj):

    user_expired_kits = {}
    user_expired_kits['data'] = {}
    if UserLotCommercialKits.objects.filter(user = register_user_obj, expirationDate__lt = date.today()).exists():

        user_kits = UserLotCommercialKits.objects.filter(user = register_user_obj, expirationDate__lt = date.today()).order_by('basedCommercial')
        for user_kit in user_kits:
            data_kit = []
            c_kit = user_kit.get_commercial_kit()
            data_kit.append(user_kit.get_lot_number())
            data_kit.append(user_kit.get_expiration_date())
            data_kit.append(user_kit.get_number_of_uses())
            if not c_kit in user_expired_kits['data']:
                user_expired_kits['data'][c_kit] = []
            user_expired_kits['data'][c_kit].append(data_kit)
    user_expired_kits['headings'] = HEADING_FOR_USER_LOT_INVENTORY
    return user_expired_kits

def get_valid_lot_user_kit (register_user_obj):

    valid_kits = {}
    valid_kits['data'] = {}
    if UserLotCommercialKits.objects.filter(user = register_user_obj, runOut = True).exists():

        user_kits = UserLotCommercialKits.objects.filter(user = register_user_obj, runOut = True).order_by('basedCommercial')
        for user_kit in user_kits:
            data_kit = []
            c_kit = user_kit.get_commercial_kit()
            data_kit.append(user_kit.get_lot_number())
            data_kit.append(user_kit.get_expiration_date())
            data_kit.append(user_kit.get_number_of_uses())
            if not c_kit in valid_kits['data']:
                valid_kits['data'][c_kit] = []
            valid_kits['data'][c_kit].append(data_kit)
    valid_kits['headings'] = HEADING_FOR_USER_LOT_INVENTORY
    return valid_kits

def get_lot_user_commercial_kit_basic_data(kit_obj):
    lot_kit_data = {}
    lot_kit_data['data'] = kit_obj.get_basic_data()
    lot_kit_data['heading'] = HEADING_FOR_LOT_USER_COMMERCIAL_KIT_BASIC_DATA
    return  lot_kit_data


def get_lot_commercial_kits(protocol_obj):
    '''
    Description:
        The function get the user commercial kits that are defined for using
        for the protocol.
        Because of the sharing lot commercial kits between the investigators
        the result is not longer filtered by the user whom record the kit.
    Input:
        protocol_obj  # protocol object
    Return
        user_kit_list
    '''
    user_kit_list = []

    if CommercialKits.objects.filter(protocolKits = protocol_obj).exists():
        commercial_kits = CommercialKits.objects.filter(protocolKits = protocol_obj)
        if UserLotCommercialKits.objects.filter(basedCommercial__in = commercial_kits, runOut = True).exists():
            user_kits = UserLotCommercialKits.objects.filter(basedCommercial__in = commercial_kits, runOut = True).order_by('expirationDate')
            for user_kit in user_kits:
                user_kit_list.append(user_kit.get_lot_number())
    return user_kit_list

def get_lot_reagent_commercial_kits(platform):
    '''
    Description:
        The function get the user commercial kits that are defined for using
        platform.
    Input:
        platform  # platform name
    Return
        user_platform_kit_list and commercial_list
    '''
    user_platform_kit_dict = {}
    user_platform_kit_list = []
    commercial_kit_names = []
    if CommercialKits.objects.filter(platformKits__platformName__exact = platform).exists():
        commercial_objs = CommercialKits.objects.filter(platformKits__platformName__exact = platform).order_by('name')
        for commercial_obj in commercial_objs:
            commercial_name = commercial_obj.get_name()
            commercial_kit_names.append(commercial_name)
            user_platform_kit_dict[commercial_name] = []
            if UserLotCommercialKits.objects.filter(basedCommercial = commercial_obj, runOut = True).exists():
                user_kits = UserLotCommercialKits.objects.filter(basedCommercial = commercial_obj, runOut = True).order_by('expirationDate')
                for user_kit in user_kits:
                    user_platform_kit_dict[commercial_name].append(user_kit.get_lot_number())
        user_platform_kit_list =list([(k,v) for k, v in  user_platform_kit_dict.items()])
        commercial_list = ','.join(commercial_kit_names)
    return user_platform_kit_list, commercial_list


def get_molecule_lot_kit_in_sample(sample_id):
    '''
    Description:
        The function get the user Lot commercial kits that used during molecule
        extraction.
    Input:
        sample_id  # sample id
    Return
        user_kit_list
    '''
    extraction_kits = {}
    extraction_kits['protocols_user_kits'] ={}
    if MoleculePreparation.objects.filter(sample__pk__exact = sample_id).exists():
        extraction_kits['heading_lot_kits'] = HEADING_FOR_DISPLAY_IN_SAMPLE_INFO_USER_KIT_DATA
        molecule_objs = MoleculePreparation.objects.filter(sample__pk__exact = sample_id).order_by('protocolUsed')
        for molecule_obj in molecule_objs:
            protocol_name = molecule_obj.get_protocol()
            if protocol_name not in extraction_kits['protocols_user_kits']:
                extraction_kits['protocols_user_kits'][protocol_name] = []
            kit_used_obj = molecule_obj.get_user_lot_kit_obj()
            if kit_used_obj:
                data = kit_used_obj.get_basic_data()
                data.append(molecule_obj.get_molecule_code_id())
                extraction_kits['protocols_user_kits'][protocol_name].append(data)

    return extraction_kits


def store_commercial_kit (kit_data):
    commercial_kit_values = {}
    #commercial_kit_values['protocol_id']= get_protocol_obj_from_name(kit_data['protocol'])
    commercial_kit_values['name'] = kit_data['kitName']
    commercial_kit_values['provider'] = kit_data['provider']
    commercial_kit_values['cat_number'] = kit_data ['catNo']
    commercial_kit_values['description'] = kit_data['description']
    if 'platform' in kit_data:
        commercial_kit_values['platform'] = kit_data['platform']

    new_kit = CommercialKits.objects.create_commercial_kit(commercial_kit_values )
    if not 'platform' in kit_data:
        for protocol in kit_data.getlist('protocol'):
            new_kit.protocolKits.add(get_protocol_obj_from_name(protocol))
    return new_kit

def store_lot_user_commercial_kit (kit_data, user_name):
    commercial_kit_obj = get_commercial_kit_obj_from_name(kit_data['commercialKit'])
    lot_kit_values = {}
    lot_kit_values['user'] = user_name
    lot_kit_values['basedCommercial']= commercial_kit_obj
    lot_kit_values['chipLot'] = kit_data['barCode']
    lot_kit_values['expirationDate'] = kit_data ['expirationDate']

    new_kit = UserLotCommercialKits.objects.create_user_lot_commercial_kit(lot_kit_values )
    return new_kit


def update_usage_user_lot_kit (lot_number, commercial_kit):
    '''
    Description:
        The function fetch the user lot kit filtering the lot number and the commercial
        kit name (to avoid that 2 user lot kit could have the same lot number)
        It steps in one the number of use.
        Return the user lot kit object.
    Input:
        lot_number    # number of the user lot
        commercial_kit # name of the commercial kit
    Return:
        user_lot_obj
    '''
    user_lot_obj = ''
    if UserLotCommercialKits.objects.filter(chipLot__exact = lot_number, basedCommercial__name__exact = commercial_kit).exists():
        user_lot_obj= UserLotCommercialKits.objects.filter(chipLot__exact = lot_number, basedCommercial__name__exact = commercial_kit).last()
        user_lot_obj.set_increase_use()
    return user_lot_obj
