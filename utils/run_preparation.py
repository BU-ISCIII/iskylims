import json, os
import random
from datetime import datetime
import string

from Bio.Seq import Seq
from iSkyLIMS_wetlab.models import *
from iSkyLIMS_wetlab.wetlab_config import *

from iSkyLIMS_wetlab.utils.pool_preparation import  check_if_duplicated_index
from iSkyLIMS_core.utils.handling_protocols import *
from iSkyLIMS_core.utils.handling_commercial_kits import *
from django.conf import settings

def handle_input_samples_for_run (user_post, user):
    '''
    Description:
        Function will check the information to create the new run.
        Based on the selected pool
        Return True if all are in , False if not
    Input:
        project_list    #list of the project to check
    variables:
        logger # logging object to write in the log file
    Return:
        Error if experiment name already exists
        dictionary wit the information to display as response
    '''
    lib_prep_ids = user_post['lib_prep_ids'].split(',')
    run_id = user_post['run_process_id']

    plate_name = user_post['plateName']
    container_id = user_post['containerID']
    paired = user_post['pairedEnd']
    json_data = json.loads(user_post['s_sheet_data'])
    # check if run exists
    run_obj = RunProcess.objects.get(pk__exact = run_id)
    exp_name = run_obj.get_run_name()
    record_data = {}
    if RunProcess.objects.filter(runName__exact = exp_name).exclude(state__runStateName__exact = 'Pre-Recorded').exists():
        record_data['Error'] = 'Run name already exists'
        return record_data
    run_data= {}
    run_data['experiment_name'] = exp_name
    run_data['plateName'] = user_post['plateName']
    run_data['containerID'] = user_post['containerID']

    if paired == 'True':
        heading = HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_PAIREDEND
        mapping = MAPPING_BASESPACE_SAMPLE_SHEET_TWO_INDEX
        heading_base_space = BASESPACE_FILE_TWO_INDEX
        heading_sample_sheet = HEADING_FOR_SAMPLE_SHEET_TWO_INDEX
    else:
        heading = HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_SINGLEREAD
        mapping = MAPPING_BASESPACE_SAMPLE_SHEET_ONE_INDEX
        heading_base_space = BASESPACE_FILE_ONE_INDEX
        heading_sample_sheet = HEADING_FOR_SAMPLE_SHEET_ONE_INDEX
    sample_sheet_data = {}
    # collect the data to prepare the sample Sheet
    for lib_id  in range(len(lib_prep_ids)) :
        sample_sheet_data[lib_prep_ids[lib_id]] = {}
        for column_index in range(len(heading)):
            sample_sheet_data[lib_prep_ids[lib_id]][heading[column_index]] = json_data[lib_id][column_index]

    #parsing data for Base Space
    base_space_lib = parsing_data_for_bs_file(sample_sheet_data, mapping, paired, heading_base_space )
    # Collect the information to prepare and save the file for Base Space
    project_bs_files = {}
    for bs_lib, projects in base_space_lib.items():
        project_bs_files.update(create_base_space_file(projects, bs_lib, plate_name, container_id , exp_name, paired))
    #import pdb; pdb.set_trace()
    # Handle the input information to create the sample sheet file

    #sample_sheet_file = handle_sample_sheet(sample_sheet_data)
    new_sample_sheet_data = update_sample_sheet(sample_sheet_data, lib_prep_ids)

    data_for_sample_sheet_file = parsing_data_for_sample_sheet_file(new_sample_sheet_data, mapping, heading_sample_sheet)
    # sample_sheet_file_name = create_sample_sheet_file(data_for_sample_sheet_file, user,reads, adapter, exp_name,
    import pdb; pdb.set_trace()
    sample_sheet_file_name = create_sample_sheet_file(data_for_sample_sheet_file, user, user_post.reads, user_post.adapter,  exp_name, user_post.assay , user_post.colection_index, paired)



    record_data['sample_sheet'] = sample_sheet_file_name
    record_data['exp_name'] = exp_name
    record_data['lib_prep_ids'] = lib_prep_ids
    record_data['projects_in_lib']= project_bs_files
    return record_data
    '''

        pool_name = json_data[row_index][HEADING_FOR_CREATING_POOL.index('Pool Name')]
        bs_lib_name = json_data[row_index][HEADING_FOR_CREATING_POOL.index('BaseSpace Library')]
        if bs_lib_name == '':
            bs_lib_name = 'default'
        if not bs_lib_name in sample_sheet_data:
            sample_sheet_data[bs_lib_name] = {}
        if not pool_name in sample_sheet_data[bs_lib_name]:
            sample_sheet_data[bs_lib_name][pool_name] = {}

        sample_sheet_data[bs_lib_name][pool_name][lib_prep_ids[row_index]] = {}

        # collect the index sequence information to allow to set values different that created by IEM
        i7_index_id = json_data[row_index][HEADING_FOR_CREATING_POOL.index('I7 Index')]
        i7_seq = json_data[row_index][HEADING_FOR_CREATING_POOL.index('I7 Sequence')]

        i5_index = json_data[row_index][HEADING_FOR_CREATING_POOL.index('I5 Index')]
        i5_seq = json_data[row_index][HEADING_FOR_CREATING_POOL.index('I5 Sequence')]

        sample_sheet_data[bs_lib_name][pool_name][lib_prep_ids[row_index]] ['i7_seq'] = i7_seq
        sample_sheet_data[bs_lib_name][pool_name][lib_prep_ids[row_index]] ['i5_seq'] = i5_seq
        sample_sheet_data[bs_lib_name][pool_name][lib_prep_ids[row_index]] ['project_name']= json_data[row_index][HEADING_FOR_CREATING_POOL.index('Project Name')]
        import pdb; pdb.set_trace()
    sample_sheet_file = ''
    '''
    return project_bs_files

def get_pool_objs_from_ids (pool_ids):
    '''
    Description:
        The function get the pool instance from the list of pool ids, and return a lisf of pool objs
    Return:
        pool_obj
    '''
    pool_objs = []
    for pool in pool_ids :
        pool_objs.append(LibraryPool.objects.get(pk__exact = pool))
    return pool_objs

def get_pool_adapters(pool_objs):
    '''
    Description:
        The function get the adapters used for each pool.
        return a dictionary with adapter as a key and the pool name list as value
    Return:
        adapters
    '''
    adapters = {}
    for pool in pool_obj :
        adapter = pool.get_adapter()
        if not adapter in adapters:
            adapters[adapter] = []
        adapters[adapter].append(pool.get_pool_name())

    return adapters

def get_pool_single_paired(pool_objs):
    '''
    Description:
        The function get the single read and paired ened  used for each pool.
        return a dictionary with single_paired as a key and the pool name list as value
    Return:
        adapters
    '''
    single_paired = {}
    for pool in pool_obj :
        s_p = pool.get_pool_single_paired()
        if not s_p in single_paired:
            single_paired[s_p] = []
        single_paired[s_p].append(pool.get_pool_name())

    return single_paired

def get_pool_duplicated_index(pool_objs):
    '''
    Description:
        The function get the single read and paired ened  used for each pool.
        return a dictionary with single_paired as a key and the pool name list as value
    Functions:
        check_if_duplicated_index   # located at iSkyLIMS_wetlab/utils/pool_preparation.py
    Return:
        False if no duplicated found or the result_index cotaining the duplicated samples index
    '''
    library_preparation_ids = []
    for pool in pool_objs :
        lib_prep_objs = LibraryPreparation.objects.filter(pools = pool)
        for lib_prep_obj in lib_prep_objs :
            library_preparation_ids.append(lib_prep_obj.get_id)

        result_index = check_if_duplicated_index(library_preparation_ids)
    if 'True' in result_index:
        return 'False'
    else:
        return result_index

def check_pools_compatible(data_form):
    '''
    Description:
        The function in case that more than one pools are used for a new run it will check if, single_paired, adapters are the same.
        It checks that there are not duplicate_index.
    Functions:
        get_pool_objs_from_ids       # located at this file
        get_pool_adapters            # located at this file
        get_single_paired            # located at this file
        get_pool_duplicated_index   # located at this file
    Return:
        True if all cheks are ok, or error message to display to user
    '''
    pool_ids = data_form.getlist('poolID')
    if len(pool_ids) == 1 :
        return 'True'
    pool_objs = get_pool_objs_from_ids (pool_ids)
    # get adapters used in the pools
    adapters = get_pool_adapters(pool_objs)
    if len(adapters) > 1:
        error = {}
        error['ERROR'] = adapters
        return error
    single_paired = get_single_paired (pool_objs)
    if len(single_paired) > 1:
        error['single_paired'] = single_paired
        return error
    duplicated_index = get_pool_duplicated_index(pool_objs)
    if not 'False' in duplicate_index:
        error['duplicated_index'] = duplicate_index
        return error
    return 'True'

def create_base_space_file(project_data, bs_lib, plate, container_id, experiment_name, paired):
    data = []
    project_dict ={}

    today_date = datetime.datetime.today().strftime("%Y%m%d_%H%M%S")
    file_name = str(experiment_name + today_date + '_for_basespace_'+ bs_lib + '.csv')
    bs_file_relative_path = os.path.join( MIGRATION_DIRECTORY_FILES, file_name)
    bs_file_full_path = os.path.join(settings.MEDIA_ROOT, bs_file_relative_path)
    for project, values in project_data.items():
        for value in values :
            data.append(value)
        project_dict[project] = (os.path.join(settings.MEDIA_URL, bs_file_relative_path)).replace('/', '', 1)

    sample_data = '\n'.join(data)
    if container_id == '' :
        today_date = datetime.datetime.today().strftime("%Y%m%d")
        container_id = str('B' + today_date + id_generator())
    if paired:
        template_file = os.path.join(settings.MEDIA_ROOT,TEMPLATE_FILES_DIRECTORY,BASE_SPACE_TWO_INDEX_TEMPLATE_NAME)
    else:
        template_file = os.path.join(settings.MEDIA_ROOT,TEMPLATE_FILES_DIRECTORY,BASE_SPACE_ONE_INDEX_TEMPLATE_NAME)

    with open (template_file, 'r') as filein:
        #filein = open(template_file, 'r')
        bs_template = string.Template (filein.read())
    d = {'bs_lib': bs_lib , 'plate': plate, 'container_id': container_id,  'sample_data': sample_data}
    updated_info = bs_template.substitute(d)

    fh = open(bs_file_full_path, 'w')
    fh.write(updated_info)
    fh.close()

    return project_dict


def create_sample_sheet_file(data, user, reads, adapter, exp_name, colection_index, assay, paired):
    today_date = datetime.datetime.today().strftime("%Y%m%d_%H%M%S")
    file_name = str(exp_name + today_date + SAMPLE_SHEET)
    ss_file_relative_path = os.path.join( RUN_SAMPLE_SHEET_DIRECTORY, file_name)
    ss_file_full_path = os.path.join(settings.MEDIA_ROOT, ss_file_relative_path)

    today_date = datetime.datetime.today().strftime("%d/%m/%Y")

    if paired:
        template_file = os.path.join(settings.MEDIA_ROOT,TEMPLATE_FILES_DIRECTORY,SAMPLE_SHEET_TWO_INDEX_TEMPLATE_NAME)
    else:
        template_file = os.path.join(settings.MEDIA_ROOT,TEMPLATE_FILES_DIRECTORY,SAMPLE_SHEET_ONE_INDEX_TEMPLATE_NAME)
    sample_data = '\n'.join(data)

    with open (template_file, 'r') as filein:
        #filein = open(template_file, 'r')
        ss_template = string.Template (filein.read())

    d = {'exp_name': exp_name ,'user':user, 'reads': reads, 'date': today_date, 'colection_index': colection_index, 'assay':assay ,'adapter':adapter, 'sample_data': sample_data}
    updated_info = ss_template.substitute(d)

    fh = open(ss_file_full_path, 'w')
    fh.write(updated_info)
    fh.close()

    return ss_file_relative_path


def get_library_preparation_data_in_run (lib_prep_ids, pool_ids):
    '''
    Description:
        The function returns the value of th pool singlePaired
    Input:
        lib_prep_ids    # list having the library preparation id
        pool_ids        # pool id list
    Constant:
        HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_PAIREDEND
        HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_SINGLEREAD
    Function:
        get_type_read_sequencing             # located at this file
        collect_lib_prep_data_for_new_run    # located at this file
    Return:
        display_sample_information
    '''
    display_sample_information ={}
    single_paired = get_type_read_sequencing(pool_ids[0])
    if single_paired == 'Paired End':
        paired = True
        display_sample_information['heading'] = HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_PAIREDEND
    else:
        paired = False
        display_sample_information['heading'] = HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_SINGLEREAD
    import pdb; pdb.set_trace()
    display_sample_information['data'] = collect_lib_prep_data_for_new_run(lib_prep_ids, paired)
    display_sample_information['lib_prep_ids'] = ','.join(lib_prep_ids)
    display_sample_information['paired_end'] = paired
    display_sample_information['reads'] = ''
    display_sample_information['assay'] = ''

    return display_sample_information

def get_type_read_sequencing(pool_id):
    '''
    Description:
        The function returns the value of th pool singlePaired
    Input:
        pool_id         # pool id
    Return:
        single_paired
    '''
    single_paired = LibraryPool.objects.get(pk__exact = pool_id).get_pool_single_paired()
    return single_paired

def collect_lib_prep_data_for_new_run(lib_prep_ids, paired):
    '''
    Description:
        The function returns the library preparation data for each one in the lib_prep_ids
    Input:
        lib_prep_ids        # list of library preparations
        pairedEnd           # Boolean variable of True is "paired end" False is "single read"
    Return:
        data
    '''
    data = []
    for lib_prep_id in lib_prep_ids:
        if paired:
            data.append( LibraryPreparation.objects.get(pk__exact = lib_prep_id).get_info_for_run_paired_end())
        else:
            data.append(LibraryPreparation.objects.get(pk__exact = lib_prep_id).get_info_for_run_single_read())
    return data



def update_index_in_sample_sheet(sample_sheet_data, lib_prep_ids) :
    # check in index hve been changed to adapt them to base space.
    # if changes were made the reverte the changes in the sample sheet but
    # keep the changes in the sequencing in case the number of the adapter
    # is intencianality shorter to get increase te length reads
    for  index_lib in range(len(lib_prep_ids)):
        import pdb; pdb.set_trace()
        lib_prep_obj = LibraryPreparation.objects.get(pk__exact = lib_prep_ids[index_lib])
        sample_sheet_data[lib_prep_ids[index_lib]]['I7_Index_ID'] = lib_prep_obj.get_i7_index()
        lib_prep_obj.update_i7_index(sample_sheet_data[lib_prep_ids[index_lib]]['index'])
        if 'I5_Index_ID' in sample_sheet_data[lib_prep_ids[index_lib]]:
            sample_sheet_data[lib_prep_ids[index_lib]]['I5_Index_ID'] = lib_prep_obj.get_i5_index()
            lib_prep_obj.update_i5_index(sample_sheet_data[lib_prep_ids[index_lib]]['index2'])

    return sample_sheet_data


def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))


def get_library_prep_in_pools (pool_ids):
    '''
    Description:
        The function get the pool id list and returns the the ids for the LibraryPreparation
    Return:
        lib_prep_ids
    '''
    lib_prep_ids = []

    for pool_id in pool_ids:
        if LibraryPreparation.objects.filter(pools__exact = pool_id).exists():
            lib_prep_objs = LibraryPreparation.objects.filter(pools__exact = pool_id).order_by('registerUser')
            for lib_prep_obj in lib_prep_objs:
                lib_prep_ids.append(lib_prep_obj.get_id())
    return lib_prep_ids





def get_available_pools_for_run():
    '''
    Description:
        The function get the pool which are not used in a run. It split in 2 keys. Pool which are not assigned
        yet to a run and the pools that are associated but there is information missing to be filled.
    Return:
        pools_to_update
    '''
    pools_to_update = {}
    if LibraryPool.objects.filter(poolState__poolState__exact = 'Selected', runProcess_id = None).exists():
        pools_to_update['pools_available'] = LibraryPool.objects.filter(poolState__poolState__exact = 'Selected', runProcess_id = None)
    # get the pools that are associated to a run but not yet completed
    if LibraryPool.objects.filter(poolState__poolState__exact = 'Selected').exclude(runProcess_id = None).exists():
        pools_to_update['defined_runs'] = LibraryPool.objects.filter(poolState__poolState__exact = 'Selected').exclude(runProcess_id = None).order_by('runProcess_id')

    return pools_to_update

def get_pool_info (pools_to_update):
    '''
    Description:
        The function get the information for pool which are not used in a run.
        And information for pools that are only defined the run name
    Constant:
        HEADING_FOR_SELECTING_POOLS
        HEADING_FOR_INCOMPLETED_SELECTION_POOLS
    Return:
        pools_to_update
    '''
    pool_info = {}
    if 'pools_available' in pools_to_update:
        pool_data = {}
        pool_data['heading'] = wetlab_config.HEADING_FOR_SELECTING_POOLS
        pool_data['data'] = []
        pool_ids = []
        for pool in pools_to_update['pools_available']:
            data = pool.get_info()
            if int(pool.get_number_of_samples()) == len(LibraryPreparation.objects.filter(pools = pool)) :
                data.append(pool.get_id())
                pool_data['data'].append(data)
                pool_ids.append(pool.get_id())
            else:
                if not 'invalid_run_data' in pool_info :
                    pool_info ['invalid_run_data'] = {}
                    pool_info['invalid_run_data']['data']= []
                pool_info['invalid_run_data']['data'].append(data)
        pool_data['pool_ids'] = ','.join(pool_ids)
        pool_info['pool_data'] = pool_data
    if 'defined_runs' in pools_to_update:
        run_data = {}
        tmp_data = {}
        #run_data['r_name'] = {}
        for pool in pools_to_update['defined_runs']:
            if int(pool.get_number_of_samples()) == len(LibraryPreparation.objects.filter(pools = pool)) :
                run_name = pool.get_run_name()
                if not run_name in tmp_data:
                    tmp_data[run_name] = {}
                    tmp_data[run_name]['data'] = []
                tmp_data[run_name]['run_id'] = pool.get_run_id()
                pool_name = pool.get_pool_name()
                pool_code = pool.get_pool_code_id()
                pool_numbers = pool.get_number_of_samples()
                tmp_data[run_name]['data'].append([pool_name, pool_code, pool_numbers])
            else:
                if not 'invalid_run_data' in pool_info :
                    pool_info['invalid_run_data'] ={}
                    pool_info['invalid_run_data']['data']= []
                pool_info['invalid_run_data']['data'].append(data)
        run_info_data = []
        for r_name, values in tmp_data.items():
            run_info_data.append([r_name,values['data'],values['run_id']])
        if len(run_info_data) > 0:
            run_data['heading'] = wetlab_config.HEADING_FOR_INCOMPLETED_SELECTION_POOLS
            run_data['data'] = run_info_data
            pool_info['run_data'] = run_data
        if 'invalid_run_data' in pool_info :
            pool_info['invalid_run_data']['heading'] = wetlab_config.HEADING_FOR_INCOMPLETED_SELECTION_POOLS

    return pool_info

def display_available_pools():
    '''
    Description:
        The function call 2 functions:
        - get_available_pools_for_run to get pool instances that are available
        - get_pool_info for adding the information to display
    Functions:
        get_available_pools_for_run     # located at this file
        get_pool_info                    # located at this file
    Return:
        display_pools_for_run
    '''
    display_pools_for_run = {}
    pools_to_update = get_available_pools_for_run()
    if pools_to_update:
        display_pools_for_run = get_pool_info(pools_to_update)
    return display_pools_for_run

def get_pool_instance_from_id (pool_id):
    pool_obj = LibraryPool.objects.get(pk__exact = pool_id)
    return pool_obj

def parsing_data_for_bs_file(sample_sheet_data, mapping, paired, heading_base_space):
    base_space_lib ={}
    for sample_id in sample_sheet_data:
        lib_name = sample_sheet_data[sample_id]['Base Space Library']
        if not lib_name in base_space_lib :
            base_space_lib[lib_name] = {}
        project_name = sample_sheet_data[sample_id]['Sample_Project']

        if not project_name in base_space_lib[lib_name]:
            base_space_lib[lib_name][project_name] = []
        bs_sample_data = {}
        for item in range(len(mapping)):
            bs_sample_data[mapping[item][0]] = sample_sheet_data[sample_id][mapping[item][1]]
        # modify the index 5 if paired End

        if paired :
            if bs_sample_data['Index2Sequence'] != '':
              seq=Seq(bs_sample_data['Index2Sequence'])
              bs_sample_data['Index2Sequence']=str(seq.reverse_complement())
        # add the default values
        bs_sample_data['Species'] = ''
        bs_sample_data['NucleicAcid'] = 'DNA'
        bs_sample_data['Project']= project_name
        # build the data row for base space file
        bs_row_data = []
        for item in heading_base_space:
            bs_row_data.append(bs_sample_data[item])
        string_row_data = ','.join(bs_row_data)
        base_space_lib[lib_name][project_name].append(string_row_data)
    return base_space_lib

def parsing_data_for_sample_sheet_file(new_sample_sheet_data, mapping, heading_sample_sheet):
    data = []

    for key, values in new_sample_sheet_data.items():
        row_data = []
        #for item in range(len(mapping)):
        #    ss_sample_data[mapping[item][0]] = new_sample_sheet_data[key][mapping[item][1]]

        for item in heading_sample_sheet:
            row_data.append(values[item])
        data.append(','.join(row_data))
    return data


def prepare_lib_prep_table_new_run (index_adapters, request, extracted_data_list, file_name, assay, adapter1, adapter2):

    '''

    BORRAR
    BORRAR esta funcion
    '''
    protocol = request.POST['lib_protocols']
    single_paired = request.POST['singlePairedEnd']
    read_length = request.POST['readlength']
    user_sample_sheet_data = {}
    stored_lib_prep = {}
    stored_lib_prep['data'] = []

    if CollectionIndexKit.objects.filter(collectionIndexName__exact = index_adapters).exists():
        collection_index_kit_id = CollectionIndexKit.objects.get(collectionIndexName__exact = index_adapters)
    else:
        collection_index_kit_id = None
    register_user_obj = User.objects.get(username__exact = request.user.username)
    user_sample_sheet_data['registerUser'] = register_user_obj
    protocol_obj = Protocols.objects.get(name__exact = protocol)
    user_sample_sheet_data['collectionIndexKit_id'] = collection_index_kit_id

    user_sample_sheet_data['sampleSheet'] = file_name
    user_sample_sheet_data['assay'] = assay
    user_sample_sheet_data['adapter1'] = adapter1
    user_sample_sheet_data['adapter2'] = adapter2
    new_user_s_sheet_obj = libPreparationUserSampleSheet.objects.create_lib_prep_user_sample_sheet(user_sample_sheet_data)


    parameter_heading = get_protocol_parameters(protocol_obj)
    length_heading = len(HEADING_FIX_FOR_ADDING_LIB_PARAMETERS) + len (parameter_heading)
    stored_lib_prep['heading'] = HEADING_FIX_FOR_ADDING_LIB_PARAMETERS
    stored_lib_prep['par_heading'] = parameter_heading
    stored_lib_prep['heading_in_excel'] = ','.join(HEADING_FIX_FOR_ADDING_LIB_PARAMETERS + parameter_heading)
    lib_prep_id = []
    samples_not_available = []
    stored_lib_prep['reagents_kits'] = get_lot_commercial_kits(register_user_obj, protocol_obj)
    for extracted_data in extracted_data_list :

        if Samples.objects.filter(sampleName__exact = extracted_data['sample_id'], sampleUser = register_user_obj,
                        sampleState__sampleStateName = 'Library preparation').exists():

            sample_obj = Samples.objects.get(sampleName__exact = extracted_data['sample_id'])
            extracted_data['sample_id'] = sample_obj
            #samples_id.append(sample_obj.get_sample_id())

            extracted_data['protocol_obj'] = protocol_obj
            extracted_data['collection_index_kit_id'] = collection_index_kit_id

            molecule_obj = MoleculePreparation.objects.filter(sample = sample_obj).last()
            if LibraryPreparation.objects.filter(sample_id = sample_obj, libPrepState__libPrepState__exact = 'Created for Reuse').exists():
                lib_prep_obj = LibraryPreparation.objects.get(sample_id = sample_obj, libPrepState__libPrepState__exact = 'Created for Reuse')
                molecule_obj = lib_prep_obj.get_molecule_obj()

                last_lib_prep_for_molecule = LibraryPreparation.objects.filter(sample_id = sample_obj, molecule_id = molecule_obj).exclude(libPrepState__libPrepState__exact = 'Created for Reuse').last()
                if last_lib_prep_for_molecule :
                    last_lib_prep_code_id = last_lib_prep_for_molecule.get_lib_prep_code()
                    split_code = re.search('(.*_)(\d+)$',last_lib_prep_code_id)
                    index_val = int(split_code.group(2))
                    new_index = str(index_val +1).zfill(2)
                    lib_prep_code_id = split_code.group(1) + new_index
                    unique_s_id_split = lib_prep_obj.get_unique_sample_id().split('-')
                    inc_value = int(unique_s_id_split[-1]) + 1
                    unique_s_id_split[-1] = str(inc_value)
                    extracted_data['uniqueID'] = '-'.join(unique_s_id_split)

                else:
                    lib_prep_code_id = molecule_obj.get_molecule_code_id() + '_LIB_01'
                    split_code = lib_prep_code_id.split('_')
                    extracted_data['uniqueID'] = sample_obj.get_unique_sample_id() +'-' + split_code[-3][1:] + '-' + split_code[-1]

                extracted_data['lib_code_id'] = lib_prep_code_id
                #lib_prep_obj.update_lib_preparation_info_in_reuse_state(extracted_data)

                new_library_preparation = lib_prep_obj.update_lib_preparation_info_in_reuse_state(extracted_data, new_user_s_sheet_obj, single_paired , read_length)
                #new_library_preparation = LibraryPreparation.objects.update_library_preparation(extracted_data)
            else:
                lib_prep_code_id = molecule_obj.get_molecule_code_id() + '_LIB_01'
                extracted_data['lib_code_id'] = lib_prep_code_id
                split_code = lib_prep_code_id.split('_')
                extracted_data['uniqueID'] = sample_obj.get_unique_sample_id() +'-' + split_code[-3][1:] + '-' + split_code[-1]
                extracted_data['assay'] = assay
                # Create the new library preparation object
                new_library_preparation = LibraryPreparation.objects.create_lib_preparation(extracted_data, new_user_s_sheet_obj, register_user_obj,
                                        molecule_obj,  single_paired , read_length)
            lib_prep_id.append(new_library_preparation.get_id())
            data = ['']*length_heading
            data[0] = extracted_data['sample_id']
            data[1] = lib_prep_code_id

            if not collection_index_kit_id :
                data[2] = 'collection Index not defined'
            else:
                data[2] = collection_index_kit_id.get_collection_index_name()
            stored_lib_prep['data'].append(data)

        else:
            samples_not_available.append(extracted_data['sample_id'])

    stored_lib_prep['lib_prep_id'] = ','.join(lib_prep_id)
    stored_lib_prep['samples_not_available'] = samples_not_available
    return stored_lib_prep
