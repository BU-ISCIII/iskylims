import json
from iSkyLIMS_wetlab.models import *
from iSkyLIMS_wetlab.wetlab_config import *



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




def get_info_to_display_created_pool(lib_prep_ids, pool_name, pool_codeID):
    information_for_created_pool = {}
    information_for_created_pool['pool_name'] = pool_name
    information_for_created_pool['pool_codeID'] = pool_codeID
    information_for_created_pool['heading'] = HEADING_FOR_DISPLAY_CREATED_POOL
    lib_prep_data = []
    for lib_prep_id in lib_prep_ids :
        lib_prep_obj = LibraryPreparation.objects.get(pk__exact = lib_prep_id)
        lib_prep_data.append(lib_prep_obj.get_info_for_display_pool())
    information_for_created_pool['lib_prep_data'] = lib_prep_data
    information_for_created_pool['samples_in_pool'] = len(lib_prep_ids)
    return information_for_created_pool
