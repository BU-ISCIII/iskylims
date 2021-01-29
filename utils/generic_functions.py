import logging
import os, re
from datetime import datetime
from logging.config import fileConfig
from logging.handlers import RotatingFileHandler
from smb.SMBConnection import SMBConnection


from django.core.mail import send_mail
from django.contrib.auth.models import Group

from .sample_sheet_utils import get_projects_in_run
from django.conf import settings
from iSkyLIMS_wetlab import wetlab_config
from iSkyLIMS_wetlab.models import RunProcess, RunStates, Projects, RunningParameters, SambaConnectionData, EmailData, ConfigSetting
from iSkyLIMS_core.models import SequencerInLab, SequencingPlatform


def check_all_projects_exists (project_list):
    '''
    Description:
        Function will check if all projects given in the project list
        are defined on database
        Return True if all are in , False if not
    Input:
        project_list    #list of the project to check
    variables:
        logger # logging object to write in the log file
    Return:
        True /False
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function for check_all_projects_exists')
    for project in project_list :
        if not Projects.objects.filter(projectName__exact = project).exists():
            logger.debug ('End function for check_all_projects_exists: Found some projects not defined')
            return False
    logger.debug ('End function for check_all_projects_exists')
    return True

def get_configuration_value(parameter_name):
    '''
    Description:
        Function will get the parameter value defined in the configutration table
        if not exists return 'False'

    Input:
        parameter_name    #parameter name
    Return:
        parameter_value
    '''
    parameter_value = 'False'
    if ConfigSetting.objects.filter(configurationName__exact = parameter_name).exists():
        parameter_obj = ConfigSetting.objects.filter(configurationName__exact = parameter_name).last()
        parameter_value = parameter_obj.get_configuration_value()
    return parameter_value

def get_run_in_same_year_to_compare (run_object):
    '''
    Description:
        The function return a list with all runs that have been created on the same year and they
        are based on the same chemistry.
    Input:
        run_object    # runProcess object
    Return:
        same_run_in_year run_object list
    '''
    # get the chemistry type for the run, that will be used to compare runs with the same chemistry value
    chem_high_mid = RunningParameters.objects.get(runName_id__exact = run_object).get_run_chemistry()
    run_different_chemistry = RunningParameters.objects.all(). exclude(Chemistry__exact = chem_high_mid)
    run_year = run_object.get_run_year()

    start_date = str(run_year) + '-1-1'
    end_date = str(run_year) +'-12-31'
    same_run_in_year = RunProcess.objects.filter(run_date__range=(start_date, end_date)).exclude(runName__in = run_different_chemistry)
    return same_run_in_year

def check_valid_date_format (date):
    try:
        datetime.strptime(date, '%Y-%m-%d')
        return True
    except:
        return False


    '''
    def get_sequencer_lanes_number_from_file (input_file, experiment_name):
    '''
    '''
    Description:
        Function read the RunInfo.xml file to get the number of lanes
    Input:
        input_file        # input file to get the number of lanes
        experiment_name     # experiment name used for logging
    Return:
        number_of_lanes
    '''
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting function get_sequencer_lanes_number_from_file', experiment_name)
    number_of_lanes = ''
    fh = open (input_file, 'r')
    search_line = '.*LaneCount="(\d).*'
    for line in fh:
        lane_found = re.search('^\s+ %s' % search_line, line)
        if lane_found:
            number_of_lanes = lane_found.group(1)
            break
    fh.close()
    logger.info('%s : number of lanes %s' , experiment_name , number_of_lanes)
    logger.debug ('%s : End function get_sequencer_lanes_number_from_file', experiment_name)
    return number_of_lanes
    '''


def save_email_data(email_fields):
    '''
    Description:
        create the email configuration file . If exists the old information is deleted
    Input:
        email_fields    # Email fields settings
    Return:
        email_data_obj
    '''
    if EmailData.objects.all().exists():
        email_data_obj = EmailData.objects.last().update_data(email_fields)
    else:
        email_data_obj = EmailData.objects.create_email_data(email_fields)
    return email_data_obj


def get_email_data():
    '''
    Description:
        Fetch the email configuration file
    Return:
        email_data
    '''
    email_data = {}
    if EmailData.objects.all().exists():
        email_data_obj = EmailData.objects.last()
        email_data = email_data_obj.get_email_data()
    return email_data


def get_samba_connection_data():
    '''
    Description:
        Fetch the samba configuration from database
    Return:
        samba_data
    '''
    samba_data = {}
    if SambaConnectionData.objects.all().exists():
        samba_connection_obj = SambaConnectionData.objects.all().last()
        samba_data = samba_connection_obj.get_samba_data()

    return samba_data

def save_samba_connection_data(data):
    '''
    Description:
        store the information in database. If it is the first attempt to store data
        the instance is create if not the table is updated with the new data
    Input:
        data    # Dictionary having samba connection data
    Return:
        samba_connection_obj
    '''
    if not SambaConnectionData.objects.all().exists():
        samba_connection_obj = SambaConnectionData.objects.create()
    else:
        samba_connection_obj = SambaConnectionData.objects.all().last()
    samba_connection_obj.update_data(data)
    return samba_connection_obj


def open_samba_connection():
    '''
    Description:
        The function open a samba connection with the parameter settings
        defined in wetlab configuration file
    Functions:
        get_samba_connection_data   # located at this file
    Return:
        conn object for the samba connection
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function open_samba_connection')
    samba_data = get_samba_connection_data()
    if not samba_data :
        string_message = 'Samba connection data on database is empty'
        logging_errors (string_message, True, False)
    if samba_data['SAMBA_APPLICATION_FOLDER_NAME'] != '':
        samba_data['SAMBA_SHARED_FOLDER_NAME'] = os.path.join(samba_data['SAMBA_SHARED_FOLDER_NAME'] , samba_data['SAMBA_APPLICATION_FOLDER_NAME'])
    conn = SMBConnection(samba_data['SAMBA_USER_ID'], samba_data['SAMBA_USER_PASSWORD'],
        samba_data['SAMBA_SHARED_FOLDER_NAME'],samba_data['SAMBA_REMOTE_SERVER_NAME'],
        use_ntlm_v2=samba_data['SAMBA_NTLM_USED'], domain=samba_data['SAMBA_DOMAIN'],
        is_direct_tcp=samba_data['IS_DIRECT_TCP'] )
    #try:
    if samba_data['SAMBA_HOST_NAME'] :
        conn.connect(socket.gethostbyname(samba_data['SAMBA_HOST_NAME']), int(samba_data['SAMBA_PORT_SERVER']))
    else:
        conn.connect(samba_data['SAMBA_IP_SERVER'], int(samba_data['SAMBA_PORT_SERVER']))
    #except:
        #string_message = 'Unable to connect to remote server'
        #logging_errors (string_message, True, True)
    #    raise IOError ('Samba connection error')


    logger.debug ('End function open_samba_connection')
    return conn

def get_type_of_data (data):
    '''
    Description:
        The function get always as input a string class.
        By trying to convert the input data to int or bolean it will decide the type of data
        If not possible to conver it returns string
    Return:
        type_of_data
    '''
    boolean_values = ['True', 'False', 'None']
    if data in boolean_values :
        return 'boolean'
    try:
        integer = int(data)
        return 'integer'
    except:
        return 'string'


def find_xml_tag_text (input_file, search_tag):
    '''
    Description:
        The function will look for the xml element tag in the file and it will return the text value
    Input:
        input_file  # file to find the tag
        search_tag  # xml tag to be found in the input_file
    Return:
        tag value
    '''
    fh = open (input_file, 'r')
    search_line = '<' + search_tag+ '>(.*)</' + search_tag+'>'
    for line in fh:
        found_tag = re.search('^\s+ %s' % search_line, line)
        if found_tag:
            fh.close()
            return found_tag.group(1)
    fh.close()
    return 'NOT FOUND'


def get_attributes_remote_file (conn, run_dir, remote_file) :
    '''
    Description:
        Function will fetch the file from remote server and copy on local
        directory
    Input:
        conn    # Samba connection object
        run_dir # run folder to fetch the file
        remote_file # file name to fetch on remote server
    Constants:
        SAMBA_SHARED_FOLDER_NAME
    variables:
        file_attributes # contains the file attributes
        logger # logging object to write in the log file
    Return:
        file_attributes if sucessful .
        Exception if attributes could not be fetched
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function for getting remote attributes')
    try:
        file_attributes = conn.getAttributes(wetlab_config.SAMBA_SHARED_FOLDER_NAME , remote_file)
        logger.info('Got attributes from %s', remote_file)
    except Exception as e:
        string_message = 'Unable to get attributes for ' + remote_file
        logging_errors (string_message, True, False)
        raise Exception('Not get attributes')
    logger.debug ('End function for  getting remote attributes')
    return file_attributes

def get_available_platform ():
    '''
    Description:
        The function fetch the available run states by quering
        the run_state table
    Variable:
        available_states   # list containing the run state names
        platforms   # all platform objects
    Return:
        available_platforms
    '''
    available_platforms =[]

    platforms = SequencingPlatform.objects.all()
    for platform in platforms :
        available_platforms.append(platform.get_platform_name())
    return available_platforms

def get_available_run_state ():
    '''
    Description:
        The function fetch the available run states by quering
        the run_state table
    Variable:
        available_states   # list containing the run state names
    Return:
        available_states
    '''
    available_states = []
    run_states = RunStates.objects.all()
    for state in run_states :
        available_states.append(state.runStateName)
    return available_states




def get_experiment_name_from_file (l_run_parameter) :
    '''
    Description:
        The function will get the experiment name  for the xml element tag in the
        file and it will return the experiment name value
    Input:
        l_run_parameter  # file to find the tag
    Functions:
        find_xml_tag_text # located at this file
    Variables:
        experiment_name # name of the experiment found in runParameter file
    Return:
        experiment_name
    '''
    experiment_name = find_xml_tag_text (l_run_parameter, wetlab_config.EXPERIMENT_NAME_TAG)

    return experiment_name


def get_log_file_name(config_log_file) :
    '''
    Description:
        The function will get the log file name from the configuration
        file and it will return the fullpath log file name
    Input:
        config_log_file  # configuration log file
    Variables:
        log_file_name # name of found log file name
    Return:
        log_file_name
    '''
    log_file_name = ''
    with open(config_log_file) as fh :
        for line in fh :
            if '.log' in line :
                 log_file_name = line.split('\'')[1]
    return log_file_name







def is_wetlab_manager (request):
    '''
    Description:
        The function will check if the logged user belongs to wetlab
        manager group
    Input:
        request # contains the session information
    Return:
        Return True if the user belongs to Wetlab Manager, False if not
    '''
    try:
        groups = Group.objects.get(name = wetlab_config.WETLAB_MANAGER)
        if groups not in request.user.groups.all():
            return False
    except:
        return False

    return True



def get_run_disk_utilization (conn, run_folder):
    '''
    Description:
        Function to get the size of the run directory on the
        remote server.
        It will use the function get_size_dir to get the size utilization
        for each subfolder
    Input:
        conn # Connection samba object
        run_folder   # root folder to start the checking file size
        run_processing_id #
    Functions:
        get_size_dir    # Located on this file
    Variables:
        file_list # contains the list of file and subfolders
        count_file_size # partial size for the subfolder
        disk_utilization # dictionnary used to keep the disk space used
                        in the run
    Return:
        disk_utilization # in the last iteraction will return the total
                    size of the folder
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function get_run_disk_utilization')
    full_path_run = os.path.join(wetlab_config.SAMBA_APPLICATION_FOLDER_NAME, run_folder)
    try:
        get_full_list = conn.listPath(wetlab_config.SAMBA_SHARED_FOLDER_NAME ,full_path_run)
    except:
        string_message = 'Unable to get the folder ' + run_folder
        logging_errors (string_message, True, False)
        logger.debug ('End function get_run_disk_utilization with error')
        raise
    rest_of_dir_size = 0
    data_dir_size = 0
    images_dir_size = 0
    MEGA_BYTES = 1024*1024
    disk_utilization = {}
    logger.info('Start folder iteration ')
    for item_list in get_full_list:
        if item_list.filename == '.' or item_list.filename == '..':
            continue
        if item_list.filename == 'Data':
            logger.info('Starting getting disk space utilization for Data Folder')
            dir_data = os.path.join(full_path_run,'Data')
            data_dir_size = get_size_dir(dir_data , conn)

        elif item_list.filename == 'Images':
            logger.info('Starting getting disk space utilization for Images Folder')
            dir_images = os.path.join(full_path_run, 'Images')
            images_dir_size = get_size_dir(dir_images , conn)

        if item_list.isDirectory:
            item_dir = os.path.join(full_path_run, item_list.filename)
            rest_of_dir_size += get_size_dir(item_dir, conn)
        else:
            rest_of_dir_size += item_list.file_size

    disk_utilization ['useSpaceFastaMb'] = '{0:,}'.format(round(data_dir_size/MEGA_BYTES))
    disk_utilization ['useSpaceImgMb'] = '{0:,}'.format(round(images_dir_size/MEGA_BYTES))
    disk_utilization ['useSpaceOtherMb'] = '{0:,}'.format(round(rest_of_dir_size/MEGA_BYTES))
    '''
        run_be_updated.useSpaceImgMb= images_dir_size_formated
        run_be_updated.useSpaceFastaMb= data_dir_size_formated
        run_be_updated.useSpaceOtherMb= rest_of_dir_size_formated
        run_be_updated.save()
    '''
    logger.info('End  disk space utilization for runID  %s', run_folder)
    logger.debug ('End function get_run_disk_utilization')
    return disk_utilization


def get_size_dir (directory, conn):
    '''
    Description:
        Recursive function to get the size of the run directory on the
        remote server.
        Optional can send an email to inform about the issue
    Input:
        conn # Connection samba object
        directory   # root folder to start the checking file size
    Variables:
        file_list # contains the list of file and subfolders
        count_file_size # partial size for the subfolder
    Return:
        count_file_size # in the last iteraction will return the total
                    size of the folder
    '''
    count_file_size = 0
    file_list = conn.listPath(wetlab_config.SAMBA_SHARED_FOLDER_NAME, directory)
    for sh_file in file_list:
        if sh_file.isDirectory:
            if (sh_file.filename == '.' or sh_file.filename == '..'):
                continue
            sub_directory = os.path.join (directory,sh_file.filename)
            count_file_size += get_size_dir (sub_directory, conn)
        else:
            count_file_size += sh_file.file_size

    return count_file_size


def normalized_data (set_data, all_data) :
    '''
    Description:
        The function is used to normalized data from diferent range of values
    Input:
        set_data    # contains a gruop of data
        all_data    # contains all data to be used for the normalization
    Variables:
        normalized_set_data # to keep the normalized set data
        normalized_all_data # to keep the normalized value of all data
    Return:
        normalized_set_data
        normalized_all_data.
    '''
    normalized_set_data, normalized_all_data = [] , []

    min_value = min(min(set_data),min(all_data))
    max_value = max(max(set_data), max(all_data))
    for value in set_data :
        normalized_set_data.append(format((value - min_value)/max_value,'.2f'))
    for value in all_data :
        normalized_all_data.append(format((value - min_value)/max_value,'.2f'))

    return normalized_set_data, normalized_all_data


def send_error_email_to_user ( subject, body_message, from_user, to_user):
    '''
    Description:
        Send an email to users  defined in the user list "to_user"
    Input:
    '''
    send_mail (subject, body_message, from_user, to_user)



def get_project_search_fields_form():
    project_form_data = {}

    project_form_data['run_states'] = []
    project_form_data['available_platforms'] = []
    project_form_data['available_sequencers'] = []
    run_states = RunStates.objects.all()
    for r_state in run_states :
        project_form_data['run_states'].append(r_state.get_run_state_name())

    platforms = SequencingPlatform.objects.all()
    for platform in platforms :
        project_form_data['available_platforms'].append(platform.get_platform_name())
    sequencers = SequencerInLab.objects.all()
    for sequencer in sequencers :
        project_form_data['available_sequencers'].append(sequencer.get_sequencer_name())

    return project_form_data


def get_run_search_fields_form():
    run_form_data = {}

    run_form_data['run_states'] = []
    run_form_data['available_platforms'] = []
    run_form_data['available_sequencers'] = []
    run_states = RunStates.objects.all()
    for r_state in run_states :
        run_form_data['run_states'].append(r_state.get_run_state_name())

    platforms = SequencingPlatform.objects.all()
    for platform in platforms :
        run_form_data['available_platforms'].append(platform.get_platform_name())
    machines = SequencerInLab.objects.all()
    for machine in machines :
        run_form_data['available_sequencers'].append(machine.get_sequencer_name())

    return run_form_data
