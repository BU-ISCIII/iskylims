from iSkyLIMS_wetlab.models import *
from iSkyLIMS_wetlab.wetlab_config import *

def check_index_compatible(lib_prep_ids):
    index_values ={}
    lib_prep_count = 0
    for lib_prep_id in lib_prep_ids :
        if not libraryPreparation.objects.filter(pk__exact = lib_prep_id).exists():
            continue
        lib_prep_obj =  libraryPreparation.objects.get(pk__exact = lib_prep_id)
        lib_prep_count += 1
        lib_index = lib_prep_obj.get_indexes()
        if lib_index not in index_values:
            index_values[lib_index] = []
        index_values[lib_index].append(lib_prep_id)
    if len(index_values) == lib_prep_count :
        return True
    else:
        incompatible_index = []
        
        for key, values in index_values.items():
            if len(values) > 1:
                s_name  = []
                for value in values:
                    s_name.append(libraryPreparation.objects.get(pk__exact = value).get_sample_name())
                incompatible_index.append([' and  '.join(s_name), key])

        return incompatible_index

def get_info_for_create_pool(lib_prep_ids):
    information_for_selected_run = {}
    information_for_selected_run['data'] = []
    information_for_selected_run['heading'] = HEADING_FOR_POOL_SELECTED_IN_RUN
    length_heading = len(HEADING_FOR_POOL_SELECTED_IN_RUN)
    for lib_prep_id in lib_prep_ids :
        lib_prep_obj = libraryPreparation.objects.get(pk__exact = lib_prep_id)
        info_lib_data = lib_prep_obj.get_info_for_pool()
        data = ['']*length_heading
        for i in range(len(info_lib_data)):
            data[i] = info_lib_data[i]
        information_for_selected_run['data'].append(data)
    information_for_selected_run['lib_prep_id'] = ','.join(lib_prep_ids)
    return information_for_selected_run
