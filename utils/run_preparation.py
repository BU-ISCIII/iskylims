import json, os, shutil
import random
from datetime import datetime
import string

from Bio.Seq import Seq
from iSkyLIMS_wetlab.models import *
from iSkyLIMS_wetlab.wetlab_config import *

from iSkyLIMS_wetlab.utils.pool_preparation import  check_if_duplicated_index
from iSkyLIMS_wetlab.utils.library_preparation import  get_lib_prep_obj_from_id, update_library_preparation_for_reuse, get_iem_version_for_library_prep_ids
from iSkyLIMS_core.utils.handling_samples import update_sample_reused, get_sample_obj_from_sample_name, get_molecule_objs_from_sample, update_molecule_reused
from iSkyLIMS_core.utils.handling_protocols import *
from iSkyLIMS_core.utils.handling_commercial_kits import *
from django.conf import settings
from django_utils.models import Profile, Center


def check_valid_data_for_creation_run(form_data,user_obj):
    '''
    Description:
        Function checks if polls are compatible and if run name is not defined
    Input:
        form_data    # user form
        user_obj    # user who is requesting the creation run
    Functions:
        check_pools_compatible   # located at this file
        get_library_prep_in_pools    # located at this file
    Return:
        True if all ckecks are ok, ERROR and error message if not
    '''
    error = {}
    pool_ids = form_data.getlist('poolID')
    experiment_name = form_data['experimentName']
    if RunProcess.objects.filter(runName__exact = experiment_name).exists():
        error_message = ERROR_RUN_NAME_ALREADY_DEFINED.copy()
        error_message[0]  = error_message[0] + experiment_name
        error['ERROR'] = error_message
        return error
        #return render (request,'iSkyLIMS_wetlab/CreateNewRun.html',{'invalid_exp_name':experiment_name, 'display_pools_for_run':display_pools_for_run})
    result_compatibility = check_pools_compatible (form_data)

    if 'ERROR' in result_compatibility :
        return result_compatibility
        #display_pools_for_run.update(result_compatibility)
        #return  render(request, 'iSkyLIMS_wetlab/CreateNewRun.html',{'display_pools_for_run': display_pools_for_run})

    lib_prep_ids = get_library_prep_in_pools (pool_ids)
    if len(lib_prep_ids) == 0:
        error['ERROR'] = wetlab_config.ERROR_POOLS_WITH_NO_LIBRARY
        #return  render(request, 'iSkyLIMS_wetlab/CreateNewRun.html',{'display_pools_for_run': display_pools_for_run})
    # return an error message if logged user does not have assigned either Profile or Center
    try:
        center_requested_id = Profile.objects.filter(profileUserID = user_obj).last().profileCenter.id
        center_requested_by = Center.objects.get(pk__exact = center_requested_id)
    except:
        error['ERROR'] = wetlab_config.ERROR_NO_PROFILE_OR_CENTER_FOR_USER
        return  error
    return 'OK'

def create_run_in_pre_recorded_and_get_data_for_confirmation (form_data, user_obj):
    '''
    Description:
        Function get the user data and create a new run instance. Function also gets
        the pool information to confirm the data
    Input:
        form_data    # user form
        user_obj    # user who is requesting the creation run
    Functions:
        get_library_preparation_data_in_run    # located at this file
        get_stored_user_sample_sheet   # located at this file
        fetch_reagent_kits_used_in_run    # located at this file
    Return:
        display_sample_information
    '''
    display_sample_information ={}
    pool_ids = form_data.getlist('poolID')
    lib_prep_ids = get_library_prep_in_pools (pool_ids)

    center_requested_id = Profile.objects.filter(profileUserID = user_obj).last().profileCenter.id
    center_requested_by = Center.objects.get(pk__exact = center_requested_id)
    reagent_kit_objs = fetch_reagent_kits_used_in_run(form_data)

    display_sample_information = get_library_preparation_data_in_run(lib_prep_ids, pool_ids)
    display_sample_information.update(get_stored_user_sample_sheet(lib_prep_ids))

    # update Reagents kits
    new_run_obj =  RunProcess(runName=form_data['experimentName'], sampleSheet= '',
                            state = RunStates.objects.get(runStateName__exact = 'Pre-Recorded'),
                            centerRequestedBy = center_requested_by)
    new_run_obj.save()
    for reagent_kit_obj in reagent_kit_objs :
        new_run_obj.reagent_kit.add(reagent_kit_obj)

    for pool in pool_ids:
        pool_obj = get_pool_instance_from_id(pool)
        pool_obj.update_run_name(new_run_obj)

    display_sample_information['experiment_name'] = form_data['experimentName']
    display_sample_information['run_process_id'] = new_run_obj.get_run_id()
    return display_sample_information


def create_or_add_project_to_run(run_data, user_obj, project_name, proj_data):
    '''
    Description:
        Function will check if project is already defined. If proyect is new, Then
        it is created but if already exists then run is added to the project
    Input:
        run_data    # data related to the run
        user_obj    # owner of the project
        project_name    # name of the project
        proj_data   # data related to project
    Functions:
        get_project_obj_from_project_name   # located at this file
    Return:
        None
    '''
    if not Projects.objects.filter(projectName__exact = project_name).exists():
        project_data = {}
        project_data['user_id'] = user_obj
        project_data['projectName'] = proj_name
        project_data['libraryKit'] = run_data['collection_index']
        project_data['baseSpaceFile'] = proj_data[0]
        project_data['BaseSpaceLibrary'] = proj_data[1]
        project_obj = Projects.objects.create_new_project(project_data)
    else:
        project_obj = get_project_obj_from_project_name(project_name)
    project_obj.runProcess.add(run_data['run_obj'])
    return

def get_project_obj_from_project_name(project_name):
    '''
    Description:
        The function get the project obj for the given project name
    Input:
        project_name    # project name
    Return:
        project_obj
    '''
    project_obj = None
    if Projects.objects.filter(projectName__exact = project_name).exists():
        project_obj = Projects.objects.filter(projectName__exact = project_name).last()
    return project_obj

def collect_data_and_update_library_preparation_samples_for_run (data_form, user):
    '''
    Description:
        Function collect the information in the form and update the library preparation
        data with the value confirmed bu user.
        Return a dictionary with the user form
    Input:
        data_form    # user data form
    Functions:
        prepare_fields_to_create_sample_sheet_from_template

    Constants:
        HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_MISEQ_SINGLE_READ
        HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_MISEQ_PAIRED_END
        HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_NEXTSEQ_SINGLE_READ
        HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_NEXTSEQ_PAIRED_END
        MAP_USER_SAMPLE_SHEET_TO_DATABASE_MISEQ_SINGLE_READ
        MAP_USER_SAMPLE_SHEET_TO_DATABASE_MISEQ_PAIRED_END
        MAP_USER_SAMPLE_SHEET_TO_DATABASE_NEXTSEQ_SINGLE_READ
        MAP_USER_SAMPLE_SHEET_TO_DATABASE_NEXTSEQ_PAIRED_END
    Return:
        record_data with the information collected from the user form
    '''
    record_data = {}
    lib_prep_ids = data_form['lib_prep_ids'].split(',')
    lib_prep_unique_ids = data_form['lib_prep_unique_ids'].split(',')
    run_id = data_form['run_process_id']
    json_data = json.loads(data_form['s_sheet_data'])
    single_read = data_form['single_read']
    projects = []
    iem_version = get_iem_version_from_user_sample_sheet(lib_prep_ids)

    if data_form['platform_type'] == 'MiSeq':
        if single_read == 'TRUE':
            if iem_version == '4':
                heading = HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_MISEQ_SINGLE_READ_VERSION_4
                mapping_fields = MAP_USER_SAMPLE_SHEET_TO_DATABASE_MISEQ_SINGLE_READ_VERSION_4
            else:
                heading = HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_MISEQ_SINGLE_READ_VERSION_5
                mapping_fields = MAP_USER_SAMPLE_SHEET_TO_DATABASE_MISEQ_SINGLE_READ_VERSION_5
        else:
            if iem_version == '4':
                heading = HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_MISEQ_PAIRED_END_VERSION_4
                mapping_fields = MAP_USER_SAMPLE_SHEET_TO_DATABASE_MISEQ_PAiRED_END_VERSION_4
            else:
                heading = HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_MISEQ_PAIRED_END_VERSION_5
                mapping_fields = MAP_USER_SAMPLE_SHEET_TO_DATABASE_MISEQ_PAiRED_END_VERSION_5
    else:
        if single_read == 'TRUE':
            heading = HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_NEXTSEQ_SINGLE_READ
            mapping_fields = MAP_USER_SAMPLE_SHEET_TO_DATABASE_NEXTSEQ_SINGLE_READ
        else:
            heading = HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_NEXTSEQ_PAIRED_END
            mapping_fields = MAP_USER_SAMPLE_SHEET_TO_DATABASE_NEXTSEQ_PAIRED_END

    #sample_sheet_data = {}
    sample_sheet_data_field = []
    # Store the confirmation information form user form
    for row_index in range(len(json_data)) :
        confirmation_data = {}
        right_id = lib_prep_ids[lib_prep_unique_ids.index(json_data[row_index][0])]
        lib_prep_obj = get_lib_prep_obj_from_id(right_id)

        for mapping_index in range(len(mapping_fields)):
            confirmation_data[mapping_fields[mapping_index][1]] = json_data[row_index][heading.index(mapping_fields[mapping_index][0])]
        projects.append(json_data[row_index][heading.index('Sample_Project')])
        # keep the old value of the user sample sheet and update the library prepation
        confirmation_data['user_sample_sheet'] = lib_prep_obj.get_user_sample_sheet_obj()
        lib_prep_obj.update_library_preparation_with_indexes(confirmation_data)
        sample_sheet_data_field.append(json_data[row_index])
    record_data['sample_sheet_data_field'] = sample_sheet_data_field
    record_data['sample_sheet_data_heading'] = heading
    record_data['run_obj'] = get_run_obj_from_id(data_form['run_process_id'])
    # get data for sample Sheet

    record_data['investigator'] = user.username
    record_data['application'] = data_form['application']
    if 'instrument' in data_form:
        record_data['instrument'] = data_form['instrument']
    else:
        record_data['instrument'] = ''
    record_data['assay'] = data_form['assay']
    if 'collection_index' in data_form:
        record_data['collection_index'] = data_form['collection_index']
    else:
        record_data['collection_index'] = ''
    record_data['reads'] = data_form['reads']
    record_data['adapter'] = data_form['adapter']
    if 'adapter2' in data_form:
            record_data['adapter2'] = data_form['adapter2']

    record_data['exp_name'] = record_data['run_obj'].get_run_name()
    record_data['lib_prep_ids'] = lib_prep_ids
    record_data['projects'] = list(set(projects))

    return record_data

def create_new_projects_added_to_run(project_list, run_obj, user_obj):
    '''
    Description:
        The function create the project instance and join them to the run obj
    Input:
        project_list    # List of project to create
        run_obj     # run object
        user_obj    # user object
    Ouput:
        project_objs
    '''
    projects_objs = []
    duplicated_projects = []
    if ConfigSetting.objects.filter(configurationName__exact = 'PROJECTS_ALLOWED_IN_MULTIPLE_RUNS').exists():
        allow_projects_in_multi_run = ConfigSetting.objects.filter(configurationName__exact = 'PROJECTS_ALLOWED_IN_MULTIPLE_RUNS').last().get_configuration_value()
    else:
        allow_projects_in_multi_run = 'TRUE'
    for project in project_list :
        if Projects.objects.filter(projectName__iexact = project).exists():
            if allow_projects_in_multi_run == 'FALSE':
                duplicated_projects.append(project)
                continue
            project_obj = Projects.objects.filter(projectName__iexact = project).last()
            projects_objs.append(project_obj)
        else:
            project_data = {}
            project_data['user_id'] = user_obj
            project_data['projectName'] = project
            project_obj = Projects.objects.create_new_empty_project(project_data)
            projects_objs.append(project_obj)
    if len(duplicated_projects) > 0 :
        # delete defined projects
        for project_obj in projects_objs:
            project_obj.delete()
        error_message = ERROR_NOT_ALLOWED_REPEATED_PROJECTS.copy()
        error_message.append(', '.join(duplicated_projects))
        error = {'ERROR':error_message}
        return error
    # join project wirh run
    for project_obj in projects_objs:
        project_obj.runProcess.add(run_obj)

    return projects_objs

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
    for pool in pool_objs :
        adapter = pool.get_adapter()
        if not adapter in adapters:
            adapters[adapter] = []
        adapters[adapter].append(pool.get_pool_name())

    return adapters

#  def get_pool_single_paired(pool_objs):
    '''
    Description:
        The function get the single read and paired end  used for each pool.
        return a dictionary with single_paired as a key and the pool name list as value
    Return:
        adapters

    single_paired = {}
    for pool in pool_obj :
        s_p = pool.get_pool_single_paired()
        if not s_p in single_paired:
            single_paired[s_p] = []
        single_paired[s_p].append(pool.get_pool_name())

    return single_paired
    '''

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
            library_preparation_ids.append(lib_prep_obj.get_id())

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

        get_pool_duplicated_index   # located at this file
    Return:
        True if all cheks are ok, or error message
    '''
    error = {}
    pool_ids = data_form.getlist('poolID')
    if len(pool_ids) == 1 :
        return 'True'
    pool_objs = get_pool_objs_from_ids (pool_ids)
    # get adapters used in the pools
    adapters = get_pool_adapters(pool_objs)
    if len(adapters) > 1:
        error_message = ERROR_DIFFERENT_ADAPTERS_USED_IN_POOL.copy()
        error['ERROR'] = error_message.insert(1, ','.join(adapters))
        return error

    duplicated_index = get_pool_duplicated_index(pool_objs)
    if not 'False' in duplicated_index:
        error_message = ERROR_DUPLICATED_INDEXES_FOUND_IN_DIFFERENT_POOLS.copy()
        for duplicated in duplicated_index['incompatible_index']:
            error_message.insert(1, ','.join(duplicated))
        error['ERROR'] = error_message
        return error
    return 'True'
'''
def create_base_space_file(projects, bs_lib, plate, container_id, experiment_name, paired, ):
    data = []
    project_dict ={}

    today_date = datetime.datetime.today().strftime("%Y%m%d_%H%M%S")
    file_name = str(experiment_name + today_date + '_for_basespace_'+ bs_lib + '.csv')
    bs_file_relative_path = os.path.join( MIGRATION_DIRECTORY_FILES, file_name)
    bs_file_full_path = os.path.join(settings.MEDIA_ROOT, bs_file_relative_path)
    for project in projects.keys():
        for value in projects[project]['data'] :
            data.append(value)
        project_dict[project] = [(os.path.join(settings.MEDIA_URL, bs_file_relative_path)).replace('/', '', 1), bs_lib, projects[project]['project_user']]

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
'''

def store_confirmation_sample_sheet(fields):
    '''
    Description:
        The function store the sample sheet for the run, using the template and returning the file name including the relative path

    Functions:
        get_pool_objs_from_ids       # located at this file
        get_pool_adapters            # located at this file
        get_single_paired            # located at this file
        get_pool_duplicated_index   # located at this file
    Constants:
        SAMPLE_SHEET
        RUN_TEMP_DIRECTORY
        RUN_SAMPLE_SHEET_DIRECTORY
        MEDIA_ROOT
        TEMPLATE_FILES_DIRECTORY
        SAMPLE_SHEET_TWO_ADAPTERS_TEMPLATE_NAME
        SAMPLE_SHEET_ONE_ADAPTER_TEMPLATE_NAME
    Return:
        sample sheet name including the relative path
    '''
    exp_name_in_file = fields['exp_name'].replace(' ', '̣_')
    today_date = datetime.datetime.today().strftime("%Y%m%d_%H%M%S")
    file_name = str(exp_name_in_file + '_'+ today_date + '_' + SAMPLE_SHEET)
    tmp_file_relative_path = os.path.join( RUN_TEMP_DIRECTORY, file_name)
    ss_file_relative_path = os.path.join(RUN_SAMPLE_SHEET_DIRECTORY, file_name)
    ss_file_full_path = os.path.join(settings.MEDIA_ROOT, tmp_file_relative_path)

    today_date = datetime.datetime.today().strftime("%d/%m/%Y")

    d = {'investigator':fields['investigator'],'exp_name': fields['exp_name'] , 'date': today_date, 'application': fields['application'],
        'instrument':fields['instrument'], 'assay':fields['assay'] , 'collection_index': fields['collection_index'], 'reads': fields['reads'],
        'adapter':fields['adapter']}

    if 'adapter2' in fields:
        d['adapter2'] = fields['adapter2']
        template_file = os.path.join(settings.MEDIA_ROOT,TEMPLATE_FILES_DIRECTORY,SAMPLE_SHEET_TWO_ADAPTERS_TEMPLATE_NAME)
    else:

        template_file = os.path.join(settings.MEDIA_ROOT,TEMPLATE_FILES_DIRECTORY,SAMPLE_SHEET_ONE_ADAPTER_TEMPLATE_NAME)

    #sample_data = '\n'.join(data)

    with open (template_file, 'r') as filein:
        #filein = open(template_file, 'r')
        ss_template = string.Template (filein.read())


    updated_info = ss_template.substitute(d)
    fh = open(ss_file_full_path, 'w')

    fh.write(updated_info)
    fh.write(','.join(fields['sample_sheet_data_heading']) + '\n')
    for sample in fields['sample_sheet_data_field']:
        fh.write(','.join(sample) + '\n')
    fh.close()
    # store sample sheet in database
    fields['run_obj'].update_sample_sheet(ss_file_full_path , ss_file_relative_path, file_name)
    #os.remove(ss_file_full_path)
    import pdb; pdb.set_trace()
    return ss_file_relative_path

def get_experiment_name (run_id):
    '''
    Description:
        The function returns the experiment name
    Input:
        run_id # run process id
    Return:
        experiment_name
    '''
    run_obj= RunProcess.objects.get(pk__exact = run_id)
    experiment_name = run_obj.get_run_name()
    return experiment_name

def get_library_preparation_unique_id (lib_prep_ids):
    '''
    Description:
        The function returns the library preparation unique id for the input list
    Input:
        lib_prep_ids # list having the library preparation id
    Functions:
        get_lib_prep_obj_from_id   # located at iSkyLIMS_wetlab/utils/library_preparation.py
    Return:
        unique_id_list
    '''
    unique_id_list =[]
    for lib_prep_id in lib_prep_ids :
        unique_id_list.append(get_lib_prep_obj_from_id(lib_prep_id).get_unique_id())
    return unique_id_list

def get_library_preparation_data_in_run (lib_prep_ids, pool_ids):
    '''
    Description:
        The function returns the information for the library preparation
    Input:
        lib_prep_ids    # list having the library preparation id
        pool_ids        # pool id list
    Constant:
        HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_PAIREDEND
        HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_SINGLEREAD
    Function:
        get_type_read_sequencing             # located at this file
        collect_lib_prep_data_for_new_run    # located at this file
        get_library_preparation_unique_id    # located at this file
    Return:
        display_sample_information
    '''
    display_sample_information ={}

    pool_objs = get_pool_objs_from_ids(pool_ids)
    sequencers_platform = []
    for pool_obj in pool_objs:
        sequencers_platform.append(pool_obj.get_platform_name())
    if 'MiSeq' in sequencers_platform:
        platform_in_pool = 'MiSeq'
    else:
        platform_in_pool = ''
    #single_paired = get_type_read_sequencing(pool_ids)
    lib_prep_data  = collect_lib_prep_data_for_new_run(lib_prep_ids, platform_in_pool)
    display_sample_information['data'] = lib_prep_data['data']
    display_sample_information['heading'] = lib_prep_data['heading']
    display_sample_information['lib_prep_ids'] = ','.join(lib_prep_ids)
    display_sample_information['lib_prep_unique_ids'] = ','.join(lib_prep_data['uniqueID_list'])

    #display_sample_information['lib_prep_unique_ids'] = ','.join(get_library_preparation_unique_id(lib_prep_ids))

    display_sample_information['date'] = today_date = datetime.datetime.today().strftime("%Y%m%d")
    display_sample_information['single_read'] = lib_prep_data['single_read']
    display_sample_information['platform_type'] = platform_in_pool

    return display_sample_information


def get_iem_version_from_user_sample_sheet(lib_prep_id):
    '''
    Description:
        The function returns the latest IEM version used for the requested library
        preparation
    Input:
        lib_prep_id   # id of the library preparation
    Return:
        iem_version
    '''
    iem_versions =  get_iem_version_for_library_prep_ids(lib_prep_id)
    if len(iem_versions) == 1:
        iem_version = iem_versions[0]
    else:
        # multiple version used in IEM, select the latest version
        latest_version = 0
        for version in iem_versions:
            version = int(version)
            if version > latest_version:
                latest_version = version
        iem_version = str(latest_version)
    return iem_version


def get_stored_user_sample_sheet(lib_prep_ids):
    '''
    Description:
        The function returns the stored values of the user sample sheet
    Input:
        lib_prep_ids   # id of the library preparation
    Return:
        sample_sheet_data
    '''
    sample_sheet_data = {}
    iem_version = get_iem_version_from_user_sample_sheet(lib_prep_ids)
    if iem_version == '0':
        lib_prep_obj = get_lib_prep_obj_from_id(lib_prep_ids[0])
    else:
        for lib_prep_id in lib_prep_ids:
            lib_prep_obj = get_lib_prep_obj_from_id(lib_prep_id)
            if lib_prep_obj.get_iem_version() == iem_version:
                break
    u_sample_sheet_obj = lib_prep_obj.get_user_sample_sheet_obj()
    fields = ['collection_index','application', 'instrument', 'adapter1', 'adapter2','assay','reads']
    not_index_fields_in_version_4 = [0,2]
    data = u_sample_sheet_obj.get_all_data()

    for i in range(len(fields)) :
        if iem_version == '4' and i in not_index_fields_in_version_4:
            continue
        else:
            sample_sheet_data[fields[i]] = data[i]
    return sample_sheet_data

# def get_type_read_sequencing(pool_ids):
    '''
    Description:
        The function returns the value of th pool singlePaired
    Input:
        pool_ids        # pool id list
    Return:
        PairedEnd or SingleRead
    '''
    '''
    for pool_id in pool_ids :
        single_paired = LibraryPool.objects.get(pk__exact = pool_id).get_pool_single_paired()
        if single_paired == 'PairedEnd':
            return 'PairedEnd'
    return 'SingleRead'
    '''

def collect_lib_prep_data_for_new_run(lib_prep_ids, platform_in_pool):
    '''
    Description:
        The function returns the library preparation data for each one in the lib_prep_ids.
        Sample_uniqueID is modified by adding '-' and the number of reused for the library preparation
    Input:
        lib_prep_ids        # list of library preparations
        platform_in_pool    # platfrom used to add new fields if needed
    Return:
        data
    '''
    lib_data = {}
    data = []
    single_read = True
    uniqueID_list = []
    iem_version = get_iem_version_from_user_sample_sheet(lib_prep_ids)
    for lib_prep_id in lib_prep_ids:
        lib_prep_obj =LibraryPreparation.objects.get(pk__exact = lib_prep_id)
        # get index 7 and index 5. then index 5 is deleted if all samples are single end
        row_data = lib_prep_obj.get_info_for_run_paired_end()
        if single_read:
            if lib_prep_obj.get_i5_index() != '':
                single_read = False
        if platform_in_pool == 'MiSeq':
            if iem_version == '5':
                row_data.insert(8,lib_prep_obj.get_manifest())
                row_data.insert(9,lib_prep_obj.get_genome_folder())
            else:
                row_data.insert(8,lib_prep_obj.get_genome_folder())
        else:
            row_data.insert(4,lib_prep_obj.get_index_plate_well())
        row_data[0] = row_data[0] + '-' + lib_prep_obj.get_reused_value()
        uniqueID_list.append(row_data[0])
        data.append(row_data)
    lib_data['data'] = data
    lib_data['uniqueID_list'] = uniqueID_list
    lib_data['single_read'] = 'TRUE' if single_read else 'FALSE'

    # if there is not index 5 in any of the samples. Then delete the I5 index colunm for all samples
    if single_read:
        for item in data:
            if platform_in_pool == 'MiSeq':
                del item[7]
                del item[6]
            else:
                # sample well was inserted in position 3, so delete index for i5 are stepped in one
                del item[8]
                del item[7]
    if platform_in_pool == 'MiSeq':
        if single_read:
            if iem_version == '4':
                lib_data['heading'] = HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_MISEQ_SINGLE_READ_VERSION_4
            else:
                lib_data['heading'] = HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_MISEQ_SINGLE_READ_VERSION_5
        else:
            if iem_version == '4':
                lib_data['heading'] = HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_MISEQ_PAIRED_END_VERSION_4
            else:
                lib_data['heading'] = HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_MISEQ_PAIRED_END_VERSION_5
    else:
        if single_read:

            lib_data['heading'] = HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_NEXTSEQ_SINGLE_READ
        else:
            lib_data['heading'] = HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_NEXTSEQ_PAIRED_END
    return lib_data

def increase_reuse_if_samples_exists(sample_list):
    '''
    Description:
        The function check if sample was exists (only this function is requested when
        run is reused, by uploading the sample sheet
    Input:
        sample_list        # list of samples fetched from sample sheet
    Return:
        samples_reused
    '''
    samples_reused = []
    for sample in sample_list:
        sample_obj = get_sample_obj_from_sample_name(sample)
        if sample_obj:
            update_sample_reused(sample_obj.get_sample_id())
            molecules = get_molecule_objs_from_sample(sample_obj)
            last_molecule_obj = molecules.reverse()[0]
            update_molecule_reused(sample_obj.get_sample_id(), last_molecule_obj.get_molecule_code_id())
            update_library_preparation_for_reuse(samples_reused)
            samples_reused.append(sample)

    return samples_reused

#def update_index_in_sample_sheet(sample_sheet_data, lib_prep_ids) :
    '''
    Description:
        The function check if in index have been changed.
        The function does not change the original values recorded at the time that
        sample sheet was uploaded.
    Input:
        sample_sheet_data        # Data inside sample sheet
        lib_prep_ids        #
    Return:
        samples_reused
    '''
    '''
    Borrar
    # check in index hve been changed
    #
    for  index_lib in range(len(lib_prep_ids)):

        lib_prep_obj = LibraryPreparation.objects.get(pk__exact = lib_prep_ids[index_lib])
        sample_sheet_data[lib_prep_ids[index_lib]]['I7_Index_ID'] = lib_prep_obj.get_i7_index()
        lib_prep_obj.update_i7_index(sample_sheet_data[lib_prep_ids[index_lib]]['index'])
        if 'I5_Index_ID' in sample_sheet_data[lib_prep_ids[index_lib]]:
            sample_sheet_data[lib_prep_ids[index_lib]]['I5_Index_ID'] = lib_prep_obj.get_i5_index()
            lib_prep_obj.update_i5_index(sample_sheet_data[lib_prep_ids[index_lib]]['index2'])

    return sample_sheet_data
    '''

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
        The function get the pool which are not used in a run. It splits in 2 keys. Pool which are not assigned
        yet to a run and the pools that are associated but there is information missing to be filled.
    Return:
        pools_to_update
    '''
    pools_to_update = {}
    if LibraryPool.objects.filter(poolState__poolState__exact = 'Selected', runProcess_id = None).exists():
        pool_objs = LibraryPool.objects.filter(poolState__poolState__exact = 'Selected', runProcess_id = None).order_by('platform')
        pools_to_update['pools_available'] = {}
        for pool_obj in pool_objs :
            platform = pool_obj.get_platform_name()
            if not platform in pools_to_update['pools_available']:
                pools_to_update['pools_available'][platform] = []
            pools_to_update['pools_available'][platform].append(pool_obj)
    # get the pools that are associated to a run but not yet completed
    if LibraryPool.objects.filter(poolState__poolState__exact = 'Selected').exclude(runProcess_id = None).exists():
        pools_to_update['defined_runs'] = LibraryPool.objects.filter(poolState__poolState__exact = 'Selected').exclude(runProcess_id = None).order_by('runProcess_id')

    return pools_to_update

def get_pool_info (pools_to_update):
    '''
    Description:
        The function get the information for pool which are not used in a run.
        And information for pools that are only defined the run name
    Input:
        pools_to_update     # contains the pool objects split in pools_available and defined_runs
                        pools_available is a dictionary which contains as key the
                        platform and value a list of the pool objects
    Constant:
        HEADING_FOR_SELECTING_POOLS
        HEADING_FOR_INCOMPLETED_SELECTION_POOLS
    Functions:
        get_lot_reagent_commercial_kits # located at iSkyLIMS_core/utils/handling_commercial_kits
    Return:
        pool_info
    '''
    pool_info = {}
    if 'pools_available' in pools_to_update:
        pool_data = {}
        pool_data['heading'] = wetlab_config.HEADING_FOR_SELECTING_POOLS
        pool_data['platform'] = {}
        reagents_kits = {}
        commercial_kits = {}
        for platform, protocol_objs in  pools_to_update['pools_available'].items():
            pool_data['platform'][platform] = []
            for pool in protocol_objs:
                data = pool.get_info()
                # compare the number of samples to check that no samples are deleted
                if int(pool.get_number_of_samples()) == len(LibraryPreparation.objects.filter(pools = pool)) :
                    data.append(pool.get_id())
                    pool_data['platform'][platform].append(data)

                    # get the reagents kits used for the platform
                    reagents_kits[platform], commercial_kits[platform] = get_lot_reagent_commercial_kits(platform)
                else:
                    if not 'invalid_run_data' in pool_info :
                        pool_info ['invalid_run_data'] = {}
                        pool_info['invalid_run_data']['data']= []
                    pool_info['invalid_run_data']['data'].append(data)

        #pool_data['pool_ids'] = ','.join(pool_ids)
        pool_info['pool_data'] = pool_data
        pool_info['reagents_kits'] = reagents_kits
        pool_info['commercial_kits'] = commercial_kits
    if 'defined_runs' in pools_to_update:
        run_data = {}
        tmp_data = {}
        #run_data['r_name'] = {}
        for pool in pools_to_update['defined_runs']:
            run_name = pool.get_run_name()
            if int(pool.get_number_of_samples()) == len(LibraryPreparation.objects.filter(pools = pool)) :

                if not run_name in tmp_data:
                    tmp_data[run_name] = {}
                    tmp_data[run_name]['data'] = []
                tmp_data[run_name]['run_id'] = pool.get_run_id()
                pool_name = pool.get_pool_name()
                pool_code = pool.get_pool_code_id()
                pool_numbers = pool.get_number_of_samples()
                tmp_data[run_name]['data'].append([pool_name, pool_code, pool_numbers])
            else:
                data = pool.get_info()
                data.append(pool.get_id())

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

def get_run_obj_from_id(run_id):
    '''
    Description:
        The function return the runprocess object for the run_id
    Input:
        run_id     # runProcess id
    Return:
        run_obj
    '''
    run_obj = RunProcess.objects.get(pk__exact = run_id)
    return run_obj

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
            base_space_lib[lib_name][project_name] = {}
            base_space_lib[lib_name][project_name]['data'] = []
            base_space_lib[lib_name][project_name]['project_user'] = sample_sheet_data[sample_id]['Description']
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
        base_space_lib[lib_name][project_name]['data'].append(string_row_data)
    return base_space_lib
'''
def parsing_data_for_sample_sheet_file(new_sample_sheet_data, heading_sample_sheet):
    data = []

    for key, values in new_sample_sheet_data.items():
        row_data = []
        #for item in range(len(mapping)):
        #    ss_sample_data[mapping[item][0]] = new_sample_sheet_data[key][mapping[item][1]]

        for item in heading_sample_sheet:
            row_data.append(values[item])
        data.append(','.join(row_data))
    return data
'''

def prepare_fields_to_create_sample_sheet_from_template(data_form, user):
    '''
    Description:
        The function collect data from user form to build the dictionary with the right values to replace on the tmeplate
    Input:
        data_for    # information from the user form
        user        # userid
    Functions:
        get_run_obj_from_id                    # located at this file
    Return:
        fields  # dictionary having the key and values used to replace in the template
    '''
    fields={}
    data = []
    temp_data =  json.loads(data_form['s_sheet_data'])
    for line in temp_data:
        # remove the last item in data (base space field)
        line.pop()
        data.append(','.join(line))

    fields['data'] = '\n'.join(data)
    fields['exp_name']= get_run_obj_from_id(data_form['run_process_id']).get_run_name()
    fields['user'] = user
    fields['application'] = data_form['application']
    fields['instrument'] = data_form['instrument']
    fields['assay'] = data_form['assay']
    fields['collection_index'] = data_form['collection_index']
    fields['pairedEnd'] = data_form['pairedEnd']
    fields['reads'] = data_form['reads']
    fields['adapter'] = data_form['adapter']
    if 'adapter2' in data_form:
        fields['adapter2'] = data_form['adapter2']

    return fields


#def store_sample_sheet_in_tmp_folder(run_process_id):
    '''
    Description:
        The function get the orignal sample sheet and copy it into the tmp/run_process_id folder
    Input:
        run_process_id    # run process id
    Constants:
        MEDIA_ROOT
        RUN_TEMP_DIRECTORY_RECORDED
    Return:
        sample_sheet_relative
    '''
    '''
    run_obj = get_run_obj_from_id(run_process_id)
    sample_sheet_relative = run_obj.get_sample_file()
    sample_sheet_original = os.path.join(settings.MEDIA_ROOT ,sample_sheet_relative )
    temp_directory = os.path.join(settings.MEDIA_ROOT , wetlab_config.RUN_TEMP_DIRECTORY_RECORDED, run_process_id)
    if not os.path.exists(temp_directory) :
        os.mkdir(temp_directory)
    # set group writing permission to the temporary directory
    os.chmod(temp_directory, 0o774)

    sample_sheet_copy= os.path.join(temp_directory, 'samplesheet.csv' )
    shutil.copy(sample_sheet_original,sample_sheet_copy)
    # set the group write permission to the Sample Sheet File
    os.chmod(sample_sheet_copy, 0o664)

    return sample_sheet_relative
    '''

def fetch_reagent_kits_used_in_run (form_data):
    '''
    Description:
        The function fetch the reagent kits in the form
        Return an object list with the reagents user kits.
    Input:
        form_data    # data from the user form
    Return:
        user_reagents_kit_objs
    '''
    user_reagents_kit_objs = []
    commercial_kit_names = form_data['commercialKits'].split(',')
    user_reagentKit_id_kit_dict = {}
    for kit_name in commercial_kit_names:
        user_reagents_kit_objs.append(update_usage_user_lot_kit(form_data[kit_name]))

    return user_reagents_kit_objs


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
