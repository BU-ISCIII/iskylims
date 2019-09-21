import json, os
import random
from datetime import datetime
import string

from Bio.Seq import Seq
from iSkyLIMS_wetlab.models import *
from iSkyLIMS_wetlab.wetlab_config import *
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
    exp_name = user_post['experimentName']
    plate_name = user_post['plateName']
    container_id = user_post['containerID']
    paired = user_post['pairedEnd']
    json_data = json.loads(user_post['s_sheet_data'])
    # check if run exists
    record_data = {}
    if RunProcess.objects.filter(runName__exact = exp_name).exists():
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
    sample_sheet_file_name = create_sample_sheet_file(data_for_sample_sheet_file, user, '151', 'adapter-ccc',  exp_name, 'colection_index-22', 'Nextera', paired)



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
        project_dict[project] = (bs_file_relative_path)

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

def check_type_read_sequencing(lib_prep_ids):
    single_read = 0
    paired_end = 0
    for lib_prep_id in lib_prep_ids :
        single_paired = LibraryPreparation.objects.get(pk__exact = lib_prep_id).get_single_paired()
        if single_paired == 'Paired End':
            paired_end +=1
        else:
            single_read += 1
    if paired_end == 0:
        return 'Single Read'
    else:
        return 'Paired End'

def collect_lib_prep_data_for_new_run(lib_prep_ids, paired):
    data = []
    for lib_prep_id in lib_prep_ids:
        if paired:
            data.append( LibraryPreparation.objects.get(pk__exact = lib_prep_id).get_info_for_run_paired_end())
        else:
            data.append(LibraryPreparation.objects.get(pk__exact = lib_prep_id).get_info_for_run_single_read())
    return data



def update_sample_sheet(sample_sheet_data, lib_prep_ids) :
    # check in index hve been changed to adapt them to base space.
    # if changes were made the reverte the changes in the sample sheet but
    # keep the changes in the sequencing in case the number of the adapter
    # is intencianality shorter to get increase te length reads
    for  index_lib in range(len(lib_prep_ids)):
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
    lib_prep_ids = []

    for pool_id in pool_ids:
        if LibraryPreparation.objects.filter(pool_id__exact = pool_id).exists():
            lib_prep_objs = LibraryPreparation.objects.filter(pool_id__exact = pool_id).order_by('registerUser')
            for lib_prep_obj in lib_prep_objs:
                lib_prep_ids.append(lib_prep_obj.get_id())
    return lib_prep_ids



def get_available_pools_for_run():
    if not LibraryPool.objects.filter(poolState__poolState__exact = 'Defined').exists():
        return
    pools_available = LibraryPool.objects.filter(poolState__poolState__exact = 'Defined')
    for pool in pools_available :
        pass
    return pools_available

def get_pool_info (pools_available):
    pool_data = {}
    pool_data['heading'] = wetlab_config.HEADING_FOR_SELECTING_POOLS
    pool_data['data'] = []
    pool_ids = []
    for pool in pools_available:

        data = pool.get_info()
        data.append(pool.get_id())
        pool_data['data'].append(data)
        pool_ids.append(pool.get_id())
    pool_data['pool_ids'] = ','.join(pool_ids)
    return pool_data


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
