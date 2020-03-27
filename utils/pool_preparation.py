import json
from datetime import date
from iSkyLIMS_wetlab.models import *
from iSkyLIMS_wetlab.wetlab_config import *

def get_lib_prep_adapter(lib_prep_ids):
    '''
    Description:
        The function gets the different adapters used for the library preparations instance.
        Returns the adapter value
    Input:
        lib_prep_ids  # library preparation id list
    Return:
        returns a list of adapters in adapter_list
    '''
    adapters = []
    for lib_prep_id in lib_prep_ids :
        if not LibraryPreparation.objects.filter(pk__exact = lib_prep_id).exists():
            continue
        lib_prep_obj =  LibraryPreparation.objects.get(pk__exact = lib_prep_id)
        adapters.append(lib_prep_obj.get_adapters()[0])
    adapter_list = list(set(adapters))
    return adapter_list



def check_single_paired_compatible(lib_prep_ids):
    '''
    Description:
        The function check if the library preparations instance has the same index, (Single Read or Paired End).
        If there is some library ids that are not in database the function return an error message.
        If the library preparation contains Single reads and Paired End , return False.
    Input:
        lib_prep_ids  # library preparation id list
    Return:
        error message if library id does not exists.
        False when incompatible
        True all library preparations belong to the same index.
    '''
    single_paired = ''
    not_defined_library_preparation_ids = []

    for lib_prep_id in lib_prep_ids :
        if not LibraryPreparation.objects.filter(pk__exact = lib_prep_id).exists():
            not_defined_library_preparation_ids.append(lib_prep_id)
            return 'ERROR'
        lib_prep_obj =  LibraryPreparation.objects.get(pk__exact = lib_prep_id)

        lib_prep_paired_end = lib_prep_obj.get_single_paired()
        if single_paired == '':
            single_paired = lib_prep_paired_end
            continue
        if single_paired != lib_prep_paired_end:
            return 'False'
    return 'True'


def check_if_duplicated_index (lib_prep_ids):
    '''
    Description:
        The function check if some of the library preparations instances has the same indexes.
        If there is not duplicated returns True else returns the list of the library preparation ids which
        have duplicated index .
    Input:
        lib_prep_ids  # library preparation id list
    Return:
        True if unique indexes are used.
        a dictionary with 'incompatible_index' as key and the list of samples names as value.
    '''
    index_values ={}
    for lib_prep_id in lib_prep_ids :
        lib_prep_obj =  LibraryPreparation.objects.get(pk__exact = lib_prep_id)
        #lib_prep_count += 1
        lib_index = lib_prep_obj.get_indexes()
        combined_index = str(lib_index['i7_indexID'] + '_' + lib_index['i5_indexID'])
        '''
        lib_prep_paired_end = lib_prep_obj.get_single_paired()
        if not lib_prep_paired_end in s_paired_values:
            s_paired_values[lib_prep_paired_end] = []
        s_paired_values[lib_prep_paired_end].append(lib_prep_id)
        '''
        if combined_index not in index_values:
            index_values[combined_index] = []
        index_values[combined_index].append(lib_prep_id)
    if len(index_values) == len(lib_prep_ids) :
        return 'True'
    else:
        incompatible_samples ={}
        incompatible_index = []
        # check which index
        for key, values in index_values.items():
            if len(values) > 1:
                s_name  = []
                for value in values:
                    s_name.append(LibraryPreparation.objects.get(pk__exact = value).get_sample_name())
                incompatible_index.append([' and  '.join(s_name), key])
        incompatible_samples['incompatible_index'] =incompatible_index
        return incompatible_samples

def get_single_paired(lib_prep_id):

    lib_prep_obj =  LibraryPreparation.objects.get(pk__exact = lib_prep_id)
    if lib_prep_obj.get_single_paired() == 'Paired End' :
        return 'Paired End'
    return 'Single Read'

def define_new_pool(form_data, user_obj):
    '''
    Description:
        The function performs some checks to verify the new pool can be created.

    Input:
        lib_prep_ids  # library preparation id list
    Constants:
        ERROR_LIBRARY_PREPARATION_NOT_EXISTS
    Functions:
        check_single_paired_compatible      # located at this file
        check_if_duplicated_index      # located at this file
        get_lib_prep_adapter      # located at this file
        generate_pool_code_id      # located at this file

    Return:
        error if some of the checks are not ok.
            ERROR if some of library preparation id does not run_exists
            incompatible_s_p_end if single reads and paired end are both in the selected library preparation ids
            incompatible_index with tuple of sample name and index values
            not_defined_lib_prep with error message in case that when reaching this steps they were not available

        protocol_parameter_list.

    '''

    if  'lib_prep_in_list' in form_data:
        lib_prep_ids = form_data.getlist('lib_prep_id')
        if len('lib_prep_in_list') == 0:
            lib_prep_ids = list(form_data['lib_prep_id'])
    else:
        lib_prep_ids = form_data['lib_prep_id'].split(',')

    error ={}
    # check if index are not duplicate in the library preparation
    single_paired_compatible = check_single_paired_compatible(lib_prep_ids)
    if 'ERROR' in single_paired_compatible:
        error['ERROR'] = ERROR_LIBRARY_PREPARATION_NOT_EXISTS
        return error
    if 'False' == single_paired_compatible :
        error['incompatible_s_p_end'] = True
        return error
    duplicated_index = check_if_duplicated_index(lib_prep_ids)
    if 'incompatible_index' in duplicated_index :
        error['duplicated_index'] = duplicated_index
        return error
    # check that all library preparation use the same adapter
    adapters = get_lib_prep_adapter(lib_prep_ids)
    if len(adapters) > 1 :
        error['multiple_adapters'] = adapters
        return error

    pool_data = {}
    pool_data['poolName'] = form_data['poolName']
    pool_data['poolCodeID'] = generate_pool_code_id()
    pool_data['registerUser'] = user_obj
    pool_data['adapter'] = adapters[0]
    pool_data['pairedEnd'] = get_single_paired(lib_prep_ids[0])
    pool_data['n_samples'] = len(lib_prep_ids)

    new_pool = LibraryPool.objects.create_lib_pool(pool_data)
    # update pool_id in each library_preparation belongs the new pool
    for lib_prep in lib_prep_ids:
        lib_prep_obj = LibraryPreparation.objects.get(pk__exact = lib_prep)
        lib_prep_obj.set_pool(new_pool)
        # return back to state completed if a library prepatarion was used for reused pool
        if lib_prep_obj.get_state() == 'Reused pool':
            lib_prep_obj.set_state('Completed')

    # update the number of samples
    new_pool.set_pool_state('Selected')

    return new_pool

def generate_pool_code_id():
    '''
    Description:
        The function get the date of today ( in format yyyy_mm_dd) and check if exists a pool_code_id with the same value.
        If yes the function increment the last digits else is set to "_1"
    Return:
        today_date + subindex
    '''
    today_date = date.today().strftime("%Y_%m_%d")
    if LibraryPool.objects.filter(poolCodeID__icontains = today_date).exists():
        lib_pool = LibraryPool.objects.filter(poolCodeID__icontains = today_date).last()
        latest_code_id = lib_pool.get_pool_code_id()
        last_seq_number =  int(latest_code_id.split('_')[-1])
        return str(today_date + '_' + str(last_seq_number + 1))
    else:
        return str(today_date + '_1')


def get_lib_prep_to_select_in_pool():
    '''
    Description:
        The function get the library preparation instances that do not have linked to a pool and
        collect the information to display.
    Constants:
        HEADING_FOR_DISPLAY_SAMPLES_IN_POOL
    Return:
        display_list
    '''
    display_list = {}
    display_list['data'] = []
    #import pdb; pdb.set_trace()
    from django.db.models import Q
    if LibraryPreparation.objects.filter(Q (libPrepState__libPrepState__exact = 'Completed', pools = None )| Q(libPrepState__libPrepState__exact = 'Reused pool' )).exists():
        lib_preparations =  LibraryPreparation.objects.filter(Q (libPrepState__libPrepState__exact = 'Completed', pools = None )| Q(libPrepState__libPrepState__exact = 'Reused pool' )).order_by('registerUser')
        #lib_preparations =  LibraryPreparation.objects.filter(libPrepState__libPrepState = 'Completed', pools = None).order_by('registerUser')
        for lib_prep in lib_preparations :
            display_list['data'].append( lib_prep.get_info_for_selection_in_pool())

        display_list['heading'] = wetlab_config.HEADING_FOR_DISPLAY_SAMPLES_IN_POOL

    return display_list

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
    '''
    Description:
        The function return a list of the used parameters .
    Return:
        protocol_parameter_list.

    '''
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
'''
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
'''
