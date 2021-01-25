import logging
from logging.config import fileConfig
from logging.handlers import RotatingFileHandler

import os
from smb.SMBConnection import SMBConnection
import socket
import traceback

from iSkyLIMS_wetlab.models import RunProcess, RunStates, Projects, RunningParameters, SambaConnectionData, EmailData, ConfigSetting
from .generic_functions import get_samba_connection_data, get_email_data, send_error_email_to_user


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
    Return:
        new runs
    '''
    logger = logging.getLogger(__name__)
    logger.debug('Starting function get_new_runs_on_remote_server' )
    new_runs = []
    run_data_root_folder = os.path.join('/', get_samba_application_shared_folder() )
    logger.debug('Shared folder  on remote server is : %s', run_data_root_folder)
    import pdb; pdb.set_trace()
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
                logger.info ('Found a new run  %s ',folder_run)
                new_runs.append(folder_run)
    logger.debug('End function get_new_runs_on_remote_server' )
    return new_runs

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
        import pdb; pdb.set_trace()
        send_error_email_to_user (subject, string_text, email_data['USER_EMAIL'],
                            [email_data['USER_EMAIL']])
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
