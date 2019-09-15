import json
from iSkyLIMS_wetlab.models import *
from iSkyLIMS_wetlab.wetlab_config import *

def analyze_input_pool (user_post, user):
    '''
    Description:
        Function will check the information to create the pool.
        Ch  if all projects given in the project list
        are defined on database
        Return True if all are in , False if not
    Input:
        project_list    #list of the project to check
    variables:
        logger # logging object to write in the log file
    Return:
        Error if experiment name already exists
        dictionary wit the information to display as response
    '''

    if  'lib_prep_in_list' in user_post:
        lib_prep_ids = user_post.getlist('lib_prep_id')
        if len('lib_prep_in_list') == 0:
            lib_prep_ids = list(user_post['lib_prep_id'])
    else:
        lib_prep_ids = user_post['lib_prep_id'].split(',')
    exp_name = user_post['experimentName']
    plate_name = user_post['plateName']
    container_id = user_post['containerID']
    json_data = json.loads(user_post['pool_data'])
    # check if run exists
    if RunProcess.objects.filter(runName__exact = exp_name).exists():
        return 'Error'
    pool_run_data= {}
    pool_run_data['experiment_name'] = exp_name
    pool_run_data['plateName'] = user_post['plateName']
    pool_run_data['containerID'] = user_post['containerID']

    sample_sheet_data = {}

    for row_index in range(len(json_data)) :

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
    return

def create_sample_sheet(data):
    filein = open('template_sample_sheet.csv')
    src = Template (filein.read())
    d = {'exp_name': 'CNM_0023' , 'read':151, 'project': 'Project mine'}
    result = src.substitute(d)
    fh = open('recorded_sample.csv', 'w')
    fh.write(result)
    fh.close()
    return


def get_available_pools_for_run():
    if not LibraryPool.objects.filter(poolState__poolState__exact = 'Defined').exists():
        return
    pools_available = LibraryPool.objects.filter(poolState__poolState__exact = 'Defined')

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
