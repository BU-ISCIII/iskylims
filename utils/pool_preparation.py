import json
from datetime import date
from iSkyLIMS_wetlab.models import *
from iSkyLIMS_wetlab.wetlab_config import *

def check_lib_prep_adapter(lib_prep_ids):
    adapters = []
    for lib_prep_id in lib_prep_ids :
        if not LibraryPreparation.objects.filter(pk__exact = lib_prep_id).exists():
            continue
        lib_prep_obj =  LibraryPreparation.objects.get(pk__exact = lib_prep_id)
        adapters.append(lib_prep_obj.get_adapter())
    adap_list = list(set(adapters))
    if len(adap) == 1:
        return True, ''.join(adap_list)
    else:
        return False, adap_list



def check_index_compatible(lib_prep_ids):
    index_values ={}
    s_paired_values = {}
    incompatible_samples = {}
    lib_prep_count = 0
    for lib_prep_id in lib_prep_ids :
        if not LibraryPreparation.objects.filter(pk__exact = lib_prep_id).exists():
            continue
        lib_prep_obj =  LibraryPreparation.objects.get(pk__exact = lib_prep_id)
        lib_prep_count += 1
        lib_index = lib_prep_obj.get_indexes()
        combined_index = str(lib_index['i7_indexID'] + '_' + lib_index['i5_indexID'])
        lib_prep_paired_end = lib_prep_obj.get_single_paired()
        if not lib_prep_paired_end in s_paired_values:
            s_paired_values[lib_prep_paired_end] = []
        s_paired_values[lib_prep_paired_end].append(lib_prep_id)
        if combined_index not in index_values:
            index_values[combined_index] = []
        index_values[combined_index].append(lib_prep_id)
    if (len(index_values) == lib_prep_count and len(s_paired_values) == 1):
        return True
    else:
        incompatible_index = []

        for key, values in index_values.items():
            if len(values) > 1:
                s_name  = []
                for value in values:
                    s_name.append(LibraryPreparation.objects.get(pk__exact = value).get_sample_name())
                incompatible_index.append([' and  '.join(s_name), key])
        incompatible_samples['incompatible_index'] =incompatible_index
    if len(s_paired_values) > 1:
        incompatible_s_p_end= []
        key_values = []
        for key in s_paired_values.keys():
            key_values.append([len(s_paired_values[key]), key])
        if key_values[0][0] > key_values[1][0]:
            incompatible_s_p_end = s_paired_values[key_values[1][1]]
        else:
            incompatible_s_p_end = s_paired_values[key_values[0][1]]
        incompatible_samples['incompatible_s_p_end'] = incompatible_s_p_end
    return incompatible_samples

def check_single_paired(lib_prep_ids):

    for lib_prep_id in lib_prep_ids :
        if not LibraryPreparation.objects.filter(pk__exact = lib_prep_id).exists():
            continue
        lib_prep_obj =  LibraryPreparation.objects.get(pk__exact = lib_prep_id)
        if lib_prep_obj.get_single_paired() == 'Paired End' :
            return 'Paired End'
    return 'Single Read'


def generate_pool_code_id():

    today_date = date.today().strftime("%Y_%m_%d")
    if LibraryPool.objects.filter(poolCodeID__icontains = today_date).exists():
        lib_pool = LibraryPool.objects.filter(poolCodeID__icontains = today_date).last()
        latest_code_id = lib_pool.get_pool_code_id()
        last_seq_number =  int(latest_code_id.split('_')[-1])
        return str(today_date + '_' + str(last_seq_number + 1))
    else:
        return str(today_date + '_1')


def get_info_to_display_created_pool(pool_obj):
    information_for_created_pool = {}
    information_for_created_pool['data'] = pool_obj.get_info()
    information_for_created_pool['pool_name'] = pool_obj.get_pool_name()
    information_for_created_pool['heading_pool'] = HEADING_FOR_DISPLAY_CREATED_POOL
    lib_prep_data = []
    if LibraryPreparation.objects.filter(pools = pool_obj).exists():
        lib_prep_ids = LibraryPreparation.objects.filter(pools = pool_obj)
        for lib_prep_obj in lib_prep_ids :
            lib_prep_data.append(lib_prep_obj.get_info_for_display_pool())
    information_for_created_pool['lib_prep_data'] = lib_prep_data
    information_for_created_pool['heading_library_pool'] = HEADING_FOR_DISPLAY_LIB_PREP_IN_POOL

    return information_for_created_pool

def get_lib_prep_collection_index(lib_prep_ids):
    collection_dict = {}
    for lib_prep_id in lib_prep_ids :
        if not LibraryPreparation.objects.filter(pk__exact = lib_prep_id).exists():
            continue
        lib_prep_obj =  LibraryPreparation.objects.get(pk__exact = lib_prep_id)
        collection_name = lib_prep_obj.get_collection_index_name()
        if collection_name not in collection_dict:
            collection_dict[collection_name] = 0
        collection_dict[collection_name] += 1
    if len(collection_dict) == 1:
        return collection_name
    else:
        max_value = 0
        for key, value in collection_dict.items():
            if value > max_value :
                max_value = value
                collection_name = key
        return collection_name

def get_lib_prep_assay(lib_prep_ids):
    assay_dict = {}
    for lib_prep_id in lib_prep_ids :
        if not LibraryPreparation.objects.filter(pk__exact = lib_prep_id).exists():
            continue
        lib_prep_obj =  LibraryPreparation.objects.get(pk__exact = lib_prep_id)
        assay_value = lib_prep_obj.get_assay()
        if assay_value not in assay_dict :
            assay_dict[assay_value] = 0
        assay_dict[assay_value] += 1
    if len(assay_dict) == 1:
        return assay_value
    else:
        max_value = 0
        for key, value in assay_dict.items():
            if value > max_value :
                max_value = value
                assay = key
        return assay
