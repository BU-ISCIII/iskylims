import logging
import os, re
from datetime import datetime
from logging.config import fileConfig
from logging.handlers import RotatingFileHandler
from smb.SMBConnection import SMBConnection
import socket

from django.core.mail import send_mail
from django.contrib.auth.models import Group

from .sample_sheet_utils import get_projects_in_run
from django.conf import settings
from iSkyLIMS_wetlab import wetlab_config
from iSkyLIMS_wetlab.models import RunProcess, RunStates, Projects


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

def check_valid_date_format (date):
    try:
        datetime.strptime(date, '%Y-%m-%d')
        return True
    except:
        return False

def copy_to_remote_file (conn, run_dir, remote_file, local_file) :
    '''
    Description:
        Function will fetch the file from remote server and copy on local
        directory
    Input:
        conn    # Samba connection object
        run_dir # run folder to fetch the file
        remote_file # file name to fetch on remote server
        local_file # local copy of the file fetched
    Constants:
        SAMBA_SHARED_FOLDER_NAME
    variables:
        logger # logging object to write in the log file
    Return:
        True if file was successfuly copy.
        Exception if file could not be fetched
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function for copy file to remote')
    with open(local_file ,'rb') as r_par_fp :
        try:
            conn.storeFile(wetlab_config.SAMBA_SHARED_FOLDER_NAME, remote_file, r_par_fp)
            logger.info('Saving the file %s to remote server', local_file)
        except Exception as e:
            string_message = 'Unable to copy the ' + local_file + 'file on folder ' + run_dir
            logging_errors (string_message, True, True)
            raise Exception('File not copied')
    logger.debug ('End function for copy file to remote')
    return True


def fetch_remote_file (conn, run_dir, remote_file, local_file) :
    '''
    Description:
        Function will fetch the file from remote server and copy on local
        directory
    Input:
        conn    # Samba connection object
        run_dir # run folder to fetch the file
        remote_file # file name to fetch on remote server
        local_file # local copy of the file fetched
    Constants:
        SAMBA_SHARED_FOLDER_NAME
    Import:
        datetime
    variables:
        logger # logging object to write in the log file
    Return:
        local_file if the file was successfuly copy.
        Exception if file could not be fetched
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function for fetching remote file')
    with open(local_file ,'wb') as r_par_fp :
        try:
            conn.retrieveFile(wetlab_config.SAMBA_SHARED_FOLDER_NAME, remote_file, r_par_fp)
            logger.info('Retrieving the remote %s file for %s', local_file, run_dir)
        except Exception as e:
            string_message = 'Unable to fetch the ' + local_file + 'file on folder ' + run_dir
            logging_errors (string_message, False, True)
            os.remove(local_file)
            raise Exception('File not found')
    logger.debug ('End function for fetching remote file')
    return local_file


def find_xml_tag_text (input_file, search_tag):
    '''
    Description:
        The function will look for the xml element tag in the
        file and it will return the text value
    Input:
        input_file  # file to find the tag
        search_tag  # xml tag to be found in the input_file
    Variables:
        found_tag   # line containing the tag
    '''
    fh = open (input_file, 'r')
    search_line = '<' + search_tag+ '>(.*)</' + search_tag+'>'
    for line in fh:
        found_tag = re.search('^\s+ %s' % search_line, line)
        if found_tag:
            fh.close()
            return found_tag.group(1)
    fh.close()
    return ''


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
    Import:
        platform object from iSkyLIMS_drylab
    Variable:
        available_states   # list containing the run state names
        platforms   # all platform objects
    Return:
        available_states
    '''
    from iSkyLIMS_drylab.models import Platform
    available_platforms =[]

    platforms = Platform.objects.all()
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


def get_new_runs_from_remote_server (processed_runs, conn, shared_folder):
    '''
    Description:
        The function fetch the folder names from the remote server and
        returns a list containing the folder names that have not been
        processed yet.
    Input:
        processed_runs  # full path and name of the file
        conn # samba connection object
        shared_folder   # shared folder in the remote server
        logger          # log object
    Variable:
        new_runs        # list containing the new folder run names
        run_data_root_folder # main folder where all directory runs are
        run_folder_list  # list of the folder names on remote server
    Return:
        new runs
    '''
    logger = logging.getLogger(__name__)
    logger.debug('Starting function get_new_runs_on_remote_server' )
    new_runs = []
    #run_data_root_folder = os.path.join(wetlab_config.SAMBA_APPLICATION_FOLDER_NAME, '/')
    run_data_root_folder = os.path.join('/', wetlab_config.SAMBA_APPLICATION_FOLDER_NAME )
    run_folder_list = conn.listPath( shared_folder, run_data_root_folder)
    for sfh in run_folder_list:
        if sfh.isDirectory:
            folder_run = sfh.filename
            if (folder_run == '.' or folder_run == '..'):
                continue
            # if the run folder has been already process continue searching
            if folder_run in processed_runs:
                logger.debug('folder run  %s already processed', folder_run)
                continue
            else:
                logger.info ('Found a new run  %s ',folder_run)
                new_runs.append(folder_run)
    logger.debug('End function get_new_runs_on_remote_server' )
    return new_runs


def get_experiment_name_from_file (l_run_parameter) :
    '''
    Description:
        The function will get the experiment name  for the xml element tag in the
        file and it will return the experiment name value
    Input:
        l_run_parameter  # file to find the tag
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


def get_run_platform_from_file (l_run_parameter) :
    '''
    Description:
        The function will get the run platform for the xml element tag in the
        file and it will return the platform used
    Input:
        l_run_parameter  # file to find the tag
    Variables:
        platform # name of the experiment found in runParameter file
    Return:
        platform
    '''
    platform = find_xml_tag_text (l_run_parameter, wetlab_config.APPLICATION_NAME_TAG)

    return platform


def handling_errors_in_run (experiment_name, error_code):
    '''
    Description:
        Function will manage the error situation where the run must be
        set to run state ERROR
    Input:
        experiment_name # name of the run to be updated
        error_code      # Error code
        state_when_error # Run state when error was found
    Constants:
        SAMBA_SHARED_FOLDER_NAME
    Import:
        datetime
        Projects
        RunProcess
    variables:
        logger  # logging object to write in the log file
        run_process     # RunProcess object to be updated
        project_name_list # list of all project related to the run
    Return:
        True
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function handling_errors_in_run')
    logger.info('Set run to ERROR state')
    run_process = RunProcess.objects.get(runName__exact = experiment_name)
    run_process.set_run_error_code(error_code)
    project_name_list = Projects.objects.filter(runprocess_id__exact = run_process)
    logger.info('Set projects to ERROR state')
    for project in project_name_list:
        project.set_project_state('Error')
    logger.debug ('End function handling_errors_in_run')
    return True


def is_wetlab_manager (request):
    '''
    Description:
        The function will check if the logged user belongs to wetlab
        manager group
    Input:
        request # contains the session information
    Variables:
        groups # wetlab manager object group
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

def logging_errors(string_text, showing_traceback , print_on_screen ):
    '''
    Description:
        The function will log the error information to file.
        Optional can send an email to inform about the issue
    Input:
        logger # contains the logger object
        string_text # information text to include in the log
    Functions:
        send_error_email_to_user # located on utils.wetlab_misc_utilities
    Constant:
        SENT_EMAIL_ON_ERROR
    Variables:
        subject # text to include in the subject email
    '''
    logger = logging.getLogger(__name__)
    logger.error('-----------------    ERROR   ------------------')
    logger.error(string_text )
    if showing_traceback :
        logger.error('Showing traceback: ',  exc_info=True)
    logger.error('-----------------    END ERROR   --------------')
    if wetlab_config.SENT_EMAIL_ON_ERROR :
        subject = 'Error found on wetlab'
        send_error_email_to_user (subject, string_text, FROM_EMAIL_ADDRESS,
                                TO_EMAIL_ADDRESS)
    if print_on_screen :
        from datetime import datetime
        print('********* ERROR **********')
        print(string_text)
        print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        print('When processing run . Check log for detail information')
        print('******* END ERROR ********')
    return ''

def logging_warnings(string_text, print_on_screen ):
    '''
    Description:
        The function will log the error information to file.
        Optional can send an email to inform about the issue
    Input:
        logger # contains the logger object
        string_text # information text to include in the log
    '''
    logger = logging.getLogger(__name__)
    logger.warning('-----------------    WARNING   ------------------')
    logger.warning(string_text )
    logger.warning('-----------------    END WARNING   --------------')
    if print_on_screen :
        from datetime import datetime
        print('******* WARNING ********')
        print(string_text)
        print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        print('When processing run . Check log for detail information')
        print('**** END WARNING *******')
    return ''


def need_to_wait_more (experiment_name, waiting_time):
    '''
    Description:
        The function get the time run was recorded to compare
        with the present time. If the value is less that the allowed time
        to wait  will return True.
        False is returned if the time is bigger
    Input:
        experiment_name  # experiment name to be checked
    Import:
        RunProccess     # from iSkyLIMS_wetlab.models
        datetime
    Variable:
        run_date  # date when the run was create on the sequencer
        number_of_days # difference of days between run_date and present
        today  # present day
    Return:
        True if the number of days is less that the maximum number of days
        to wait
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function need_to_wait_sample_sheet')
    run_date = RunProcess.objects.get(runName__exact = experiment_name).get_run_date()
    run_date =  datetime.strptime(run_date,"%B %d, %Y").date()
    today = datetime.now().date()
    number_of_days = abs((today - run_date).days)
    if number_of_days > int (waiting_time):
        logger.info('Waiting time already exceeded')
        logger.debug ('End function need_to_wait_sample_sheet')
        return False
    else:
        logger.info('It is allowed to waiting more time')
        logger.debug ('End function need_to_wait_sample_sheet')
        return True

def open_log(config_file):
    '''
    Description:
        The function will create the log object to write all logging information
    Input:
        logger_name    # contains the logger name that will be included
                        in the log file
    Constant:
        LOGGING_CONFIG_FILE
    Return:
        logger object
    '''
    fileConfig(config_file)
    logger = logging.getLogger(__name__)
    return logger


def open_samba_connection():
    '''
    Description:
        The function open a samba connection with the parameter settings
        defined in wetlab configuration file
    Return:
        conn object for the samba connection
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function open_samba_connection')
    conn=SMBConnection(wetlab_config.SAMBA_USER_ID, wetlab_config.SAMBA_USER_PASSWORD,
        wetlab_config.SAMBA_SHARED_FOLDER_NAME,wetlab_config.SAMBA_REMOTE_SERVER_NAME,
        use_ntlm_v2=wetlab_config.SAMBA_NTLM_USED,domain=wetlab_config.SAMBA_DOMAIN,
        is_direct_tcp=wetlab_config.IS_DIRECT_TCP )
    try:
        if wetlab_config.SAMBA_HOST_NAME :
            conn.connect(socket.gethostbyname(wetlab_config.SAMBA_HOST_NAME), int(wetlab_config.SAMBA_PORT_SERVER))
        else:
            conn.connect(wetlab_config.SAMBA_IP_SERVER, int(wetlab_config.SAMBA_PORT_SERVER))
    except:
        string_message = 'Unable to connect to remote server'
        logging_errors (string_message, True, True)
        raise IOError ('Samba connection error')


    logger.debug ('End function open_samba_connection')
    return conn


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
