import logging
from logging.config import fileConfig
from logging.handlers import RotatingFileHandler

import os, re
from smb.SMBConnection import SMBConnection
import socket
import traceback
import xml.etree.ElementTree as ET
from datetime import datetime

from django.contrib.auth.models import User

from iSkyLIMS_wetlab.models import RunProcess, RunStates, Projects, RunningParameters, SambaConnectionData, EmailData, ConfigSetting
from .generic_functions import get_samba_connection_data, get_email_data, send_error_email_to_user, find_xml_tag_text
from iSkyLIMS_wetlab.wetlab_config import *
from .sample_sheet_utils import get_projects_in_run, get_index_library_name
from iSkyLIMS_core.models import SequencerInLab


def assign_used_library_in_run (run_process_obj, l_sample_sheet_path, experiment_name, ):
    '''
    Description:
        The function first check if the run has already assinged to projects.
        If not , it reads the sample sheet, fetch the projects and asign them to the run.
    Input:
        run_process_obj     # runProcess object
        l_sample_sheet_path   # path the sample sheet
        experiment_name     # experiment name
    Functions:
        get_projects_in_run # located at utils/sample_sheet_utils.py
    Return:
        None
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting function assign_used_library_in_run', experiment_name)
    if run_process_obj.get_index_library() == 'None':
        index_library_name = get_index_library_name(l_sample_sheet_path)
        run_process_obj.update_index_library(index_library_name)
        logger.info('%s : Defined index library name', experiment_name)
    else:
        logger.info('%s : Already defined index library name', experiment_name)
    logger.debug ('%s : End function assign_used_library_in_run', experiment_name)
    return

def assign_projects_to_run(run_process_obj, sample_sheet_file , experiment_name):
    '''
    Description:
        The function first check if the run has already assinged to projects.
        If not , it reads the sample sheet, fetch the projects and asign them to the run.
    Input:
        run_process_obj     # runProcess object
        sample_sheet_file   # path the sample sheet
        experiment_name     # experiment name
    Functions:
        get_projects_in_run # located at utils/sample_sheet_utils.py
    Return:
        None
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting function assign_projects_to_run', experiment_name)
    if run_process_obj.projects_set.all().exists():
        logger.info('%s  : Projects already defined in the run', experiment_name)
        logger.debug ('%s : End function assign_projects_to_run', experiment_name)
        return
    # fetch the project from sample sheet
    projects_objs = []
    project_with_users = get_projects_in_run(sample_sheet_file)
    for project, user in project_with_users.items():
        if not Projects.objects.filter(projectName__exact = project).exists():
            project_data = {}
            project_data['projectName'] = project
            if User.objects.filter(username__iexact = user).exists():
                project_data['user_id'] = User.objects.filter(username__iexact = user).last()
            else:
                project_data['user_id'] = None
            project_obj = Projects.objects.create_new_empty_project(project_data)
            projects_objs.append(project_obj)
            # assign project to runº
            project_obj.runProcess.add(run_process_obj)
            logger.info('%s : Project name  %s added ', experiment_name, project)
        else:
            project_obj = Projects.objects.filter(projectName__exact = project).last()
            if not project_obj in projects_objs:
                projects_objs.append(project_obj)
                project_obj.runProcess.add(run_process_obj)
                logger.info('%s : Project name  %s added ', experiment_name, project)
    logger.debug ('%s : End function assign_projects_to_run', experiment_name)
    return

def check_sequencer_status_from_log_file(log_file_content, log_cycles, number_of_cycles, experiment_name):
    '''
    Description:
        The function checks in the logs if run was canceled. If not canceled, then check
        the number of created log files are the same as they expected to have when
        run processing is completed
    Input:
        log_file_content    # content of the log file
        log_cycles          # number of cycles executed on sequencer
        number_of_cycles    # number of cycle to check in logs
        experiment_name     # experiment name
    Return:
        completed and run_completion_date in case the sequencer ends sucessfully.
        If not return status of the run (cancelled/still_running)
        ERROR returns if not able to fecth log files
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting function check_sequencer_status_from_log_file', experiment_name)
    run_completion_date = ''
    if 'Cancel' in log_file_content :
        logger.warning ('%s : Sequencer execution was canceled', experiment_name)
        status = 'Cancelled'
    elif log_cycles != number_of_cycles :
        logger.info('%s : run sequencer is still running', experiment_name)
        status = 'still_running'
    else:
        # Fetch the run completcion time
        status = 'completed'
        logger.info('%s : run in sequencer is completed', experiment_name)
        last_line_in_file = log_file_content.split('\n')[-2]
        last_log_time = last_line_in_file.split(' ')[0:2]
        last_log_time[0] = str('20'+ last_log_time[0])
        string_completion_date = ' '.join(last_log_time)
        run_completion_date = datetime.strptime(string_completion_date, '%Y-%m-%d %H:%M:%S.%f')
    logger.debug ('%s : End function check_sequencer_status_from_log_file', experiment_name)
    return status, run_completion_date

def check_sequencer_status_from_completion_file(l_run_completion, experiment_name):
    '''
    Description:
        The function will check if the run in sequencer was successful
    Input:
        l_run_completion  # local path for the run completion file
        experiment_name   # Contains the experiment name
    Functions:
        find_xml_tag_text # located at utils.generic_functions
    Constant:
        COMPLETION_TAG
    Return
        True if successfuly completed
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting function for check_sequencer_status_from_completion_file', experiment_name)
    # check if NextSEq run have been successful completed
    status_run = find_xml_tag_text (l_run_completion, COMPLETION_TAG )
    if  status_run != COMPLETION_SUCCESS:
        logger.info('%s : Run in sequencer was not completed but %s', experiment_name, stats_run)
        string_message = experiment_name + ' : Sequencer Run was not completed. Reason was ' + status_run
        logging_warnings (string_message, False)
        logger.debug ('%s : End function for check_sequencer_status_from_completion_file', experiment_name)
        return False
    else:
        logger.info ('%s : Sequencer Run successfuly completed ', experiment_name)
        logger.debug ('%s : End function for check_sequencer_status_from_completion_file', experiment_name)
        return True


def check_sequencer_run_is_completed(conn, run_folder , platform ,number_of_cycles, experiment_name):
    '''
    Description:

    Input:
        conn                # Connection Samba object
        run_folder          # run folder on the remote server
        platform            # platform name
        number_of_cycles    # number of cycle to check in logs
        experiment_name     # experiment name
    Functions:
        get_samba_application_shared_folder     # located at this file
        get_latest_run_procesing_log            # located at this file
        check_sequencer_status_from_log_file    # located at this file
    Constants:
        PLATFORM_WAY_TO_CHECK_RUN_COMPLETION
        RUN_LOG_FOLDER
        RUN_COMPLETION_FILE
    Return:
        status and run_completion_date
        ERROR returns if not able to fecth log files or not defined the way to
        check the termination of the run in sequencer
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting function check_sequencer_run_is_completed', experiment_name)
    way_to_check =''
    for method in PLATFORM_WAY_TO_CHECK_RUN_COMPLETION:
        if method[0] not in platform:
            continue
        way_to_check = method[1]

    if way_to_check == 'logs':

        log_folder = os.path.join(get_samba_application_shared_folder(), run_folder , RUN_LOG_FOLDER)
        try:
            log_cycles, log_file_content = get_latest_run_procesing_log(conn, log_folder ,  experiment_name)
        except :
            string_message = experiment_name + ' : Unable to fetch the log files on run folder ' +  run_folder + '/' + log_folder
            logging_errors( string_message, True, True)
            logger.debug ('%s : End function check_sequencer_run_is_completed with exeception', experiment_name)
            return {'ERROR':18}, ''
        status , run_completion_date = check_sequencer_status_from_log_file(log_file_content ,log_cycles,number_of_cycles, experiment_name)
        logger.debug ('%s : End function check_sequencer_run_is_completed', experiment_name)
        return status , run_completion_date

    elif way_to_check == 'xml_file':
        l_run_completion = os.path.join(RUN_TEMP_DIRECTORY, RUN_COMPLETION_FILE)
        s_run_completion = os.path.join(get_samba_application_shared_folder() , run_folder, RUN_COMPLETION_FILE)

        try:
            l_run_completion = fetch_remote_file (conn, run_folder, s_run_completion, l_run_completion)
            logger.info('%s : Sucessfully fetch of Completion status file',experiment_name)
        except Exception as e:
            logger.warning ('%s : Completion status file is not present on the run folder %s', experiment_name, run_folder)
            logger.debug ('%s : End function check_sequencer_run_is_completed', experiment_name)
            return 'still_running', ''

        if check_sequencer_status_from_completion_file (l_run_completion, experiment_name):
            os.remove(l_run_completion)
            logger.info('%s : Deleted Run Completion file', experiment_name)
            shared_folder = get_samba_shared_folder()
            conversion_attributes = conn.getAttributes( shared_folder, s_run_completion)
            run_completion_date = datetime.fromtimestamp(int(conversion_attributes.create_time)).strftime('%Y-%m-%d %H:%M:%S')
            logger.debug ('%s : End function for handling NextSeq run with exception', experiment_name)
            return 'completed' , run_completion_date
        os.remove(l_run_completion)
        logger.info('%s : Deleted Run Completion file', experiment_name)
        logger.debug ('%s : End function check_sequencer_run_is_completed with exception', experiment_name)
        return 'cancelled', ''
    else:
        string_message = experiment_name + ' : way to check the completion run is not defined for  ' +  platform
        logging_errors( string_message, False, True)
        logger.debug ('%s : End function check_sequencer_run_is_completed with exeception', experiment_name)
        return {'ERROR':25}, ''


def copy_sample_sheet_to_remote_folder(conn, sample_sheet_path, run_folder ,experiment_name):
    '''
    Description:
        The functions copy the sample sheet to the remote run folder
    Input:
        conn                # Connection Samba object
        sample_sheet_path   # path of the sample sheet
        run_folder          # run folder on the remote server
        experiment_name     # experiment name
    Functions:
        get_samba_application_shared_folder     # located at this file
        copy_to_remote_file                     # located at this file
    Constants:
        EMPTY_FIELDS_IN_SEQUENCER
    Return:
        new_sequencer_obj
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting function copy_sample_sheet_to_remote_folder', experiment_name)

    logger.info('%s : Copy sample sheet to remote folder %s', experiment_name, run_folder)
    s_sample= os.path.join(get_samba_application_shared_folder(), run_folder, SAMPLE_SHEET)

    try:
        sample_sheet_copied = copy_to_remote_file  (conn, run_folder, s_sample, sample_sheet_path)
        logger.info('%s : Sucessfully copy Sample sheet to remote folder', experiment_name)
    except Exception as e:
        string_message = experiment_name + ': Unable to copy the Sample Sheet to remote folder ' + run_folder
        logging_errors (string_message, True, False)
        handling_errors_in_run (experiment_name, '23')
        logger.debug ('%s : End function for copy_sample_sheet_to_remote_folder with exception', experiment_name)
        raise Exception
    logger.debug ('%s : End function copy_sample_sheet_to_remote_folder', experiment_name)
    return

def copy_to_remote_file (conn, run_dir, remote_file, local_file) :
    '''
    Description:
        Function will fetch the file from remote server and copy on local  directory
    Input:
        conn    # Samba connection object
        run_dir # run folder to fetch the file
        remote_file # file name to fetch on remote server
        local_file # local copy of the file fetched
    Constants:
        SAMBA_SHARED_FOLDER_NAME
    Return:
        True if file was successfuly copy.
        Exception if file could not be fetched
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function for copy file to remote')
    with open(local_file ,'rb') as r_par_fp :
        try:
            samba_folder = get_samba_shared_folder()
            conn.storeFile(samba_folder, remote_file, r_par_fp)
            logger.info('Saving the file %s to remote server', local_file)
        except Exception as e:
            string_message = 'Unable to copy the ' + local_file + 'file on folder ' + run_dir
            logging_errors (string_message, True, True)
            raise Exception('File not copied')
    logger.debug ('End function for copy file to remote')
    return True


def create_new_sequencer_lab_not_defined (sequencer_name,num_of_lanes, experiment_name):
    '''
    Description:

        creates a new entry in database wit only the sequencer name and the lane numbers
    Input:
        sequencer_name    # sequencer name
        num_of_lanes        # number of lanes
        experiment_name     # experiment name
    Functions:
        get_sequencer_lanes_number_from_file # located at this file
    Constants:
        EMPTY_FIELDS_IN_SEQUENCER
    Return:
        new_sequencer_obj
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting function create_new_sequencer_lab_not_defined', experiment_name)
    seq_data = {}
    for item in EMPTY_FIELDS_IN_SEQUENCER :
        seq_data[item] = None
    seq_data['sequencerNumberLanes'] = num_of_lanes
    seq_data['sequencerName'] = sequencer_name
    new_sequencer_obj = SequencerInLab.objects.create_sequencer_in_lab(seq_data)
    logger.info('%s : Created the new sequencer in database' , experiment_name )
    logger.debug ('%s : End function create_new_sequencer_lab_not_defined', experiment_name)
    return new_sequencer_obj




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
    Return:
        local_file if the file was successfuly copy.
        Exception if file could not be fetched
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting function for fetching remote file', run_dir)
    with open(local_file ,'wb') as r_par_fp :
        try:
            samba_folder = get_samba_shared_folder()
            conn.retrieveFile(samba_folder, remote_file, r_par_fp)
            logger.info('%s : Retrieving the remote %s file for %s ', run_dir, remote_file, local_file)
        except Exception as e:
            string_message = 'Unable to fetch the ' + local_file + ' file on folder : ' + run_dir
            logging_errors (string_message, False, True)
            os.remove(local_file)
            logger.debug ('%s : End function for fetching remote file', run_dir)
            raise Exception('File not found')
    logger.debug ('%s : End function for fetching remote file', run_dir)
    return local_file


def get_latest_run_procesing_log(conn, log_folder, experiment_name) :
    '''
    Description:
        The function will find the latest log file for the input folder
    Input:
        conn        # samba connection object
        log_folder  # folder name
        experiment_name     # experiment name
    Constant:
        RUN_TEMP_DIRECTORY
    Return:
        number_of_cycles
        file_content
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting function get_latest_run_procesing_log',  experiment_name)
    shared_folder = get_samba_shared_folder()
    folder_logs = os.path.join('/', get_samba_application_shared_folder(),log_folder )

    remote_file_list = conn.listPath( shared_folder, folder_logs)
    max_cycle = -1
    logger.info('%s : Succesful connection to fetch logs files', experiment_name)
    for sfh in remote_file_list:
        if sfh.isDirectory:
            continue
        file_remote = sfh.filename
        if file_remote.endswith('.log') :
            log_file = re.search('.*_Cycle(\d+)_.*',file_remote)

            cycle_number = int(log_file.group(1))
            if cycle_number > max_cycle :
                max_cycle = cycle_number
                latest_log = file_remote
    logger.info('%s : Fetching the latest log file  %s  ', experiment_name, latest_log)
    temporary_log = os.path.join(RUN_TEMP_DIRECTORY,'miseq_cycle.log')
    s_latest_log = os.path.join(log_folder,latest_log)

    temporary_log = fetch_remote_file (conn, log_folder, s_latest_log, temporary_log)
    logger.info('%s : copied to tmp folder the log is : %s',  experiment_name, s_latest_log)
    with open (temporary_log, 'r', encoding='utf8') as fh :
        log_file_content = fh.read()

    os.remove(temporary_log)
    logger.debug ('%s : End function get_latest_run_procesing_log',  experiment_name)
    return max_cycle, log_file_content

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
    Functions:
        get_samba_application_shared_folder     # located at this file
    Return:
        new runs
    '''
    logger = logging.getLogger(__name__)
    logger.debug('Starting function get_new_runs_on_remote_server' )
    new_runs = []
    run_data_root_folder = os.path.join('/', get_samba_application_shared_folder() )
    logger.debug('Shared folder  on remote server is : %s', run_data_root_folder)
    run_folder_list = conn.listPath( shared_folder, run_data_root_folder)
    for sfh in run_folder_list:
        if sfh.isDirectory:
            folder_run = sfh.filename
            if (folder_run == '.' or folder_run == '..'):
                continue
            # if the run folder has been already process continue searching
            if folder_run in processed_runs:
                continue
            else:
                logger.info (' %s  : Found new folder run ',folder_run)
                new_runs.append(folder_run)
    logger.debug('End function get_new_runs_on_remote_server' )
    return new_runs


def get_remote_sample_sheet(conn, new_run, experiment_name):
    '''
    Description:
        The function will get the sample sheet from remote server, and it return the
        path of the sample sheet
    Input:
        conn                # samba connection instance
        new_run             # run folder on the remote server
        experiment_name     # experiment name
    Constants:
        RUN_TEMP_DIRECTORY
        SAMPLE_SHEET
        SAMBA_APPLICATION_FOLDER_NAME
        SAMPLE_SHEET
    Functions:
        get_samba_application_shared_folder     # located at this file
    Return:
        l_sample_sheet_path
    '''
    logger = logging.getLogger(__name__)
    logger.debug('%s  : Starting function get_remote_sample_sheet',experiment_name )

    l_sample_sheet_path = os.path.join(RUN_TEMP_DIRECTORY, SAMPLE_SHEET)
    s_sample_sheet_path = os.path.join(get_samba_application_shared_folder(), new_run, SAMPLE_SHEET)
    try:
        l_sample_sheet = fetch_remote_file (conn, new_run, s_sample_sheet_path, l_sample_sheet_path)
        logger.info('%s : Sucessfully fetch of Sample Sheet file', experiment_name)
    except Exception as e:
        error_message = 'Unable to fetch Sample Sheet file for folder :' +  new_run
        logging_errors(error_message, True, True)
        logger.debug('%s  : End function get_remote_sample_sheet',experiment_name )
        return None

    logger.debug('%s  : End function get_remote_sample_sheet',experiment_name )
    return  l_sample_sheet_path


def get_run_platform_from_file (l_run_parameter) :
    '''
    Description:
        The function will get the run platform for the xml element tag in the
        file and it will return the platform used
    Input:
        l_run_parameter  # file to find the tag
    Return:
        platform
    '''
    platform = find_xml_tag_text (l_run_parameter, APPLICATION_NAME_TAG)

    return platform

def get_run_process_obj_or_create_if_not_exists(running_parameters, experiment_name):
    '''
    Description:
        The function get the run_proces obj or it is created if does not exists
    Input:
        running_parameters  # dictionary to collect information to create the new run process
        experiment_name     # experiment name
    Return:
        run_process_obj
    '''
    logger = logging.getLogger(__name__)
    logger.debug('Starting function get_run_process_obj_or_create_if_not_exists' )
    if RunProcess.objects.filter(runName__exact = experiment_name).exists():
        run_process_obj = RunProcess.objects.filter(runName__exact = experiment_name).last()
    else:
        run_data = {}
        run_data['experiment_name'] = experiment_name
        run_data['run_date'] = running_parameters['run_date']
        run_process_obj = RunProcess.objects.create_new_run_from_crontab(run_data)
        logger.info('%s  : New RunProcess instance created')
    logger.debug('End function get_run_process_obj_or_create_if_not_exists' )
    return run_process_obj

def get_sequencer_obj_or_create_if_no_exists(running_parameters, experiment_name):
    '''
    Description:
        The function get the sequencer obj or it is created if does not exists
    Input:
        running_parameters      # dictionnary from parsing to get the information in case a new sequencer must be defined
        experiment_name         # experiment name
    Functions:
        create_new_sequencer_lab_not_defined    # located at this file
    Return:
        sequencer_obj
    '''
    logger = logging.getLogger(__name__)
    logger.debug('Starting function get_sequencer_obj_or_create_if_no_exists' )
    if SequencerInLab.objects.filter(sequencerName__exact = running_parameters['instrument']).exists():
        sequencer_obj = SequencerInLab.objects.filter(sequencerName__exact = running_parameters['instrument']).last()

    else:
        string_message = experiment_name + ' : ' + running_parameters['instrument'] + ' no sequencer defined '
        logging_errors(string_message, False, True)
        sequencer_obj = create_new_sequencer_lab_not_defined ( running_parameters['instrument'],  running_parameters['running_data']['NumLanes'], experiment_name)
        logger.info('%s : Continue the proccess after creating the new sequencer' ,experiment_name)

    logger.debug('End function get_sequencer_obj_or_create_if_no_exists' )
    return sequencer_obj

def get_samba_application_shared_folder():
    '''
    Description:
        The function get in database the shared application folder used for samba connections
        Application folder is a sub_directory of the shared folder. In many cases this value is empty
    Return:
        samba_folder_name
    '''
    return SambaConnectionData.objects.last().get_samba_application_folder_name()


def get_samba_shared_folder():
    '''
    Description:
        The function get in database the shared folder used for samba connections
    Return:
        samba_folder_name
    '''
    return SambaConnectionData.objects.last().get_samba_folder_name()


def handling_errors_in_run (experiment_name, error_code):
    '''
    Description:
        Function will manage the error situation where the run must be
        set to run state ERROR
    Input:
        experiment_name # name of the run to be updated
        error_code      # Error code
    Return:
        True
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting function handling_errors_in_run' , experiment_name)
    logger.info('%s : Set run to ERROR state',  experiment_name)
    if  RunProcess.objects.filter(runName__exact = experiment_name).exists():
        run_process_obj = RunProcess.objects.filter(runName__exact = experiment_name).last()
        run_process_obj.set_run_error_code(error_code)
        logger.info ('%s : is now on ERROR state' , experiment_name)
    else:
        logger.info ('%s : experiment name is not defined yet in database' , experiment_name)
    logger.debug ('%s : End function handling_errors_in_run' , experiment_name)
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
        send_error_email_to_user # located on utils.generic_functions
    Constant:
        SENT_EMAIL_ON_ERROR
        EMAIL_USER_CONFIGURED
    Variables:
        subject # text to include in the subject email
    '''
    logger = logging.getLogger(__name__)
    logger.error('-----------------    ERROR   ------------------')
    logger.error(string_text )
    if showing_traceback :
        logger.error('################################')
        logger.error(traceback.format_exc())
        logger.error('################################')
    logger.error('-----------------    END ERROR   --------------')
    email_data = get_email_data()
    if email_data['SENT_EMAIL_ON_ERROR'] :
        subject = 'Error found on wetlab when running crontab'
        send_error_email_to_user (subject, string_text, email_data['USER_EMAIL'],
                            [email_data['USER_EMAIL']])
    if print_on_screen :
        from datetime import datetime
        print('********* ERROR **********')
        print(string_text)
        print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        print('Check log for detail information')
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
        print('Check log for detail information')
        print('**** END WARNING *******')
    return ''




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


def parsing_run_info_and_parameter_information(l_run_info, l_run_parameter, experiment_name) :
    '''
    Description:
        The function is called for parsing the RunInfo and RunParameter  files.
        After parsing the RunningParameters database table will be
        updated with a new row having the parsed data
        Empty values will be set for MiSeq runs that exist on NextSeq
        but not in MiSeq runs
    Input:
        run_info    # contains the path for RunInfo.xml file
        run_parameter # contains the path for RunParameter.xml file
        experiment_name      # contains the experiment name
    CONSTANTS:
        FIELDS_TO_COLLECT_FROM_RUN_INFO_FILE
        SETUP_TAG
        NUMBER_OF_LANES_TAG
        APPLICATION_NAME_TAG
     Return:
        parsing_data
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting function parsing_run_information', experiment_name)
    running_data={}
    parsing_data = {}
    image_channel=[]

    #################################################
    ## parsing RunInfo.xml file
    #################################################
    run_data=ET.parse(l_run_info)
    run_root=run_data.getroot()
    logger.info('%s : parsing the runInfo.xml file ', experiment_name)
    p_run=run_root[0]
    ## getting the common values NextSeq and MiSeq
    logger.info('%s  : Fetching Flowcell and FlowcellLayout data ', experiment_name)
    running_data['Flowcell']=p_run.find('Flowcell').text
    try:
        running_data['FlowcellLayout']=p_run.find('FlowcellLayout').attrib
    except:
        running_data['FlowcellLayout'] = ''
        string_message = experiment_name + ' : Parameter  FlowcellLayout  not found in RunParameter.xml'
        logging_warnings(string_message, False)

    for i in run_root.iter('Name'):
        image_channel.append(i.text)

    running_data['ImageChannel'] = image_channel
    try:
        running_data['ImageDimensions'] = p_run.find('ImageDimensions').attrib
    except:
        running_data['ImageDimensions'] = ''
        logger.debug('%s : There is no image dimesions on runInfo file', experiment_name)
    ## get the instrument for NextSeq run
    parsing_data['instrument'] = p_run.find('Instrument').text

    #################################################
    ## parsing RunParameter.xml file
    #################################################
    logger.info('%s : Parsing the runParameter.xml file  ',experiment_name )
    parameter_data=ET.parse(l_run_parameter)
    parameter_data_root=parameter_data.getroot()
    p_parameter=parameter_data_root[1]
    ## getting the common values NextSeq and MiSeq
    for field in FIELDS_TO_COLLECT_FROM_RUN_INFO_FILE:
        try:
            running_data[field] = parameter_data_root.find(field).text
        except:
            running_data[field] = ''
            string_message = experiment_name + ' : Parameter ' + field  + ' unable to fetch in RunParameter.xml'
            logging_warnings(string_message, False)

    ## get the nuber of lanes in case sequencer lab is not defined

    param_in_setup = ['ApplicationVersion', 'NumTilesPerSwath']
    for i in range(len(param_in_setup)):
        try:
            running_data[param_in_setup[i]] = parameter_data_root.find('Setup').find(param_in_setup[i]).text
        except:
            string_message = experiment_name + ' : Parameter ' + param_in_setup[i] + ' not found in RunParameter.xml'
            logging_warnings(string_message, False)
            continue
    for setup_field in FIELDS_TO_FETCH_FROM_SETUP_TAG:

        try:
            running_data[setup_field] = parameter_data_root.find(SETUP_TAG).find(setup_field).text

        except:
            running_data[setup_field] = ''
            string_message = experiment_name + ' : Parameter in Setup -- ' + setup_field + ' unable to fetch in RunParameter.xml'
            logging_errors(string_message, False, True)

    if 'MiSeq' in running_data[APPLICATION_NAME_TAG] :
        # initialize paramters in case there are not exists on runParameter file
        for i in range(len(READ_NUMBER_OF_CYCLES)):
            running_data[READ_NUMBER_OF_CYCLES[i]] = ''
        # get the length index number for reads and indexes for MiSeq Runs
        for run_info_read in parameter_data_root.iter(RUN_INFO_READ_TAG):
            try:
                index_number = int(run_info_read.attrib[NUMBER_TAG]) -1
                running_data[READ_NUMBER_OF_CYCLES[index_number]] = run_info_read.attrib[NUMBER_CYCLES_TAG]
            except:
                string_message = experiment_name + ' : Parameter RunInfoRead: Read Number not found in RunParameter.xml'
                logging_warnings(string_message, False)
                continue
    ##############################################
    ## updating the date fetched from the Date tag for run and project
    ##############################################
    date = p_run.find('Date').text
    logger.debug('%s : Found date that was recorded the Run %s', experiment_name , date)
    run_date = datetime.strptime(date, '%y%m%d')
    parsing_data['running_data'] = running_data
    parsing_data['run_date'] = run_date
    logger.debug('%s : End function nextseq_parsing_run_information', experiment_name)
    return parsing_data


def save_run_parameters_data_to_database(run_parameters, run_process_obj, experiment_name):
    '''
    Description:
        The function save the run parameters if they are not store yet
    Input:
        run_parameters      # dictionnary with the running parameters
        run_process_obj      # run process object
        experiment_name     # name of the run
    Return:
        run_parameter_obj

    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting function save_run_parameters_data_to_database', experiment_name)

    if RunningParameters.objects.filter(runName_id = run_process_obj).exists():
        run_parameter_objs = RunningParameters.objects.filter(runName_id = run_process_obj)
        for run_parameter_obj in run_parameter_objs:
            logger.info('%s  : Deleting RunParameters object from database', experiment_name)
            run_parameter_obj.delete()
    run_parameter_obj = RunningParameters.objects.create_running_parameters(run_parameters, run_process_obj)
    logger.info( '%s  : Created RunParameters object on database', experiment_name)
    logger.debug ('%s : End function save_run_parameters_data_to_database', experiment_name)
    return run_parameter_obj



def store_sample_sheet_if_not_defined_in_run (run_process_obj, l_sample_sheet_path, experiment_name ) :
    '''
    Description:
        The function will move the sample sheet from the local temporary
        folder to the folder destination defined in RUN_SAMPLE_SHEET_DIRECTORY
        It will update the field sampleSheet in database
    Input:
        l_sample_sheet_path  # local copy of sample sheet file
        experiment_name  # name of the run
    Import:
        datetime
        os
    Variable:
        run_updated     # RunProcess object
        new_sample_sheet_name  # renamed sample sheet including time to have
                                a unique file name
        now     # Present time of the server
        run_updated # RunProcess object for miseq run
        sample_sheet_on_database # sample sheet path used in database
        timestr     # Present time including 3 digits for miliseconds
    Return:
        sample_sheet_on_database
    '''
    logger = logging.getLogger(__name__)
    logger.debug('%s : Starting the function store_sample_sheet_in_run', experiment_name)
    # Get the present time in miliseconds to add to have an unique file name
    now = datetime.now()
    timestr = now.strftime("%Y%m%d-%H%M%S.%f")[:-3]
    new_sample_sheet_name = 'SampleSheet' + timestr + '.csv'

    new_sample_sheet_file = os.path.join (settings.MEDIA_ROOT, RUN_SAMPLE_SHEET_DIRECTORY, new_sample_sheet_name)
    logger.info('%s : new sample sheet name %s', experiment_name, new_sample_sheet_file)
    # Path to be included in database
    sample_sheet_on_database = os.path.join(RUN_SAMPLE_SHEET_DIRECTORY, new_sample_sheet_name)
    ## Move sample sheet to final folder
    os.rename(l_sample_sheet_path, new_sample_sheet_file)
    # Update the run with the sample sheet information  (full_path, relative_path, file_name)
    run_process_obj.update_sample_sheet(new_sample_sheet_file, sample_sheet_on_database, new_sample_sheet_name)

    logger.info('%s : Updated runProccess table with the sample sheet', experiment_name)
    logger.debug('%s : End function store_sample_sheet_in_run', experiment_name)
    return sample_sheet_on_database


def waiting_time_expired(run_process_obj,time_to_check, maximun_time , experiment_name):
    '''
    Description:
        The function get the time run was recorded to compare  with the present time.
        If the value is less that the allowed time to wait  will return False.
        True is returned if the time is bigger
    Input:
        run_process_obj     # run process object
        time_to_check       # reference time to be checked
        maximun_time        # maximm number of days to wait
        experiment_name     # experiment name to be checked
    Return:
        True if the number of days is bigger that the maximum number of days to wait
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting function waiting_time_expired', experiment_name)
    today = datetime.now().date()
    number_of_days = abs((today - time_to_check).days)
    if number_of_days > int (maximun_time):
        logger.info('%s  : Waiting time already exceeded', experiment_name)
        logger.debug ('%s  : End function waiting_time_expired', experiment_name)
        return True
    else:
        logger.info('%s  : It is allowed to waiting more time', experiment_name)
        logger.debug ('%s  : End function waiting_time_expired', experiment_name)
        return False
