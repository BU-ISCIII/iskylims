import logging
from logging.config import fileConfig
from logging.handlers import RotatingFileHandler
from smb.SMBConnection import SMBConnection
import os, re
from .sample_sheet_utils import get_projects_in_run
from django.conf import settings
from iSkyLIMS_wetlab import wetlab_config


def get_new_runs_from_remote_server (processed_runs, conn, shared_folder):
    '''
    Description:
        The function fetch the folder names from the remote server and 
        returns a list containing the folder names that have not been
        porcessed yet.
    Input:
        processed_runs  # full path and name of the file
        conn # samba connection object
        shared_folder   # shared folder in the remote server
        logger          # log object 
    Variable:
        new_runs        # list containing the new folder run names
        run_folder_list  # list of the folder names on remote server
    Return:
        new runs 
    '''
    logger = logging.getLogger(__name__)
    logger.debug('Starting function get_new_runs_on_remote_server' ) 
    new_runs = []
    run_folder_list = conn.listPath( shared_folder, '/')
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
            logging_errors (logger, string_message, True, True)
            os.remove(l_run_parameter)
            raise Exception('File not found')
    logger.debug ('End function for fetching remote file')
    return local_file

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
        logging_errors (logger, string_message, True, False)
        raise Exception('Not get attributes')
    logger.debug ('End function for  getting remote attributes')
    return file_attributes

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
    with open(local_file ,'wb') as r_par_fp :
        try:
            conn.storeFile(wetlab_config.SAMBA_SHARED_FOLDER_NAME, remote_file, r_par_fp)
            logger.info('Saving the file %s to remote server', local_file)
        except Exception as e:
            string_message = 'Unable to copy the ' + local_file + 'file on folder ' + run_dir
            logging_errors (logger, string_message, True, True)
            os.remove(l_run_parameter)
            raise Exception('File not copied')
    logger.debug ('End function for copy file to remote')
    return True

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

def handling_errors_in_run (experiment_name):
    '''
    Description:
        Function will manage the error situation where the run must be
        set to run state ERROR
    Input:
        experiment_name # name of the run to be updated
    Constants:
        SAMBA_SHARED_FOLDER_NAME
    Import:
        datetime
        Projects
    variables:
        logger  # logging object to write in the log file
        run     # RunProcess object to be updated
    Return:
        run
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function handling_errors_in_run')
    logger.info('Set run to ERROR state')
    run = RunProcess.objects.get(runName__exact = experiment_name).set_run_state('ERROR')
    
    project_name_list = Projects.objects.filter(runprocess_id__exact = run_found_id)
    logger.info('Set projects to ERROR state')
    for project in project_name_list:
        project.procState= 'ERROR'
        project.save()
    logger.debug ('End function handling_errors_in_run')
    return run


def open_log(logger_name):
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
    config_file = os.path.join(settings.BASE_DIR,'iSkyLIMS_wetlab',  wetlab_config.LOGGING_CONFIG_FILE )
    fileConfig(config_file)
    logger = logging.getLogger(logger_name)
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
        use_ntlm_v2=wetlab_config.SAMBA_NTLM_USED,domain=wetlab_config.SAMBA_DOMAIN)
    try:
        conn.connect(wetlab_config.SAMBA_IP_SERVER, int(wetlab_config.SAMBA_PORT_SERVER))
    except:
        #import pdb; pdb.set_trace()
        string_message = 'Unable to connect to remote server'
        logging_errors (logger, string_message, True, True)
        raise IOError ('Samba connection error')

    
    logger.debug ('End function open_samba_connection')
    return conn


def logging_errors(logger, string_text, showing_traceback , print_on_screen ):
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
        print('******* ERROR ********')
        print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        print('When processing run in recorded state. Check log for detail information')
    return ''


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



