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
    import pdb; pdb.set_trace()
    if  'lib_prep_in_list' in user_post:
        lib_prep_ids = user_post.getlist('lib_prep_id')
        if len('lib_prep_in_list') == 0:
            lib_prep_ids = list(user_post['lib_prep_id'])
    else:
        lib_prep_ids = user_post['lib_prep_id'].split(',')
    exp_name = user_post['experimentName']
    pool_name = user_post['poolName']
    plate_name = user_post['plateName']
    container_id = user_post['containerID']
    if RunProcess.objects.filter(runName__exact = exp_name).exists():
        return 'Error'
    pool_run_data= {}
    pool_run_data['poolName'] = user_post['poolName']
    pool_run_data['plateName'] = user_post['plateName']
    pool_run_data['containerID'] = user_post['containerID']
    pool_run_data['libUsedInBaseSpace'] = user_post['basespace']

    sample_sheet_file = ''
    return

def check_index_compatible(lib_prep_ids):
    index_values ={}
    s_paired_values = {}
    incompatible_samples = {}
    lib_prep_count = 0
    for lib_prep_id in lib_prep_ids :
        if not libraryPreparation.objects.filter(pk__exact = lib_prep_id).exists():
            continue
        lib_prep_obj =  libraryPreparation.objects.get(pk__exact = lib_prep_id)
        lib_prep_count += 1
        lib_index = lib_prep_obj.get_indexes()
        lib_prep_paired_end = lib_prep_obj.get_single_paired()
        if not lib_prep_paired_end in s_paired_values:
            s_paired_values[lib_prep_paired_end] = []
        s_paired_values[lib_prep_paired_end].append(lib_prep_id)
        if lib_index not in index_values:
            index_values[lib_index] = []
        index_values[lib_index].append(lib_prep_id)
    if (len(index_values) == lib_prep_count and len(s_paired_values) == 1):
        return True
    else:
        incompatible_index = []

        for key, values in index_values.items():
            if len(values) > 1:
                s_name  = []
                for value in values:
                    s_name.append(libraryPreparation.objects.get(pk__exact = value).get_sample_name())
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

def create_sample_sheet(data):
    filein = open('template_sample_sheet.csv')
    src = Template (filein.read())
    d = {'exp_name': 'CNM_0023' , 'read':151, 'project': 'Project mine'}
    result = src.substitute(d)
    fh = open('recorded_sample.csv', 'w')
    fh.write(result)
    fh.close()
    return



def get_info_to_create_pool(lib_prep_ids):
    information_for_selected_run = {}
    base_space_library = []
    information_for_selected_run['data'] = []
    information_for_selected_run['heading'] = HEADING_FOR_CREATING_POOL
    length_heading = len(HEADING_FOR_CREATING_POOL)
    for lib_prep_id in lib_prep_ids :
        lib_prep_obj = libraryPreparation.objects.get(pk__exact = lib_prep_id)
        info_lib_data = lib_prep_obj.get_info_for_pool()
        data = ['']*length_heading
        for i in range(len(info_lib_data)):
            data[i] = info_lib_data[i]
        information_for_selected_run['data'].append(data)
    information_for_selected_run['lib_prep_id'] = ','.join(lib_prep_ids)
    if BaseSpaceLibraryName.objects.all().exists():
        bs_names = BaseSpaceLibraryName.objects.all()
        for bs_name in bs_names:
            base_space_library.append(bs_name.get_bs_lib_name())
    information_for_selected_run['base_space_library'] = base_space_library
    return information_for_selected_run
