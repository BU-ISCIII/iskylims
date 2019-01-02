##General purpose utility functions used in different iSkyLIMS_wetlab modules

import datetime
import logging
from logging.config import fileConfig
from logging.handlers import RotatingFileHandler
import os, errno, re
from smb.SMBConnection import SMBConnection

from django.conf import settings
from iSkyLIMS_wetlab import wetlab_config
from django.core.mail import send_mail

def send_error_email_to_user ( subject, body_message, from_user, to_user):
    '''
    Description:
        Send an email to users  defined in the user list "to_user"
    Input:
    '''
    send_mail (subject, body_message, from_user, to_user)

def open_samba_connection():
    '''
    Description:
        The function open a samba connection with the parameter settings
        defined in wetlab configuration file
    Return:
        conn object for the samba connection
    '''
    conn=SMBConnection(wetlab_config.SAMBA_USER_ID, wetlab_config.SAMBA_USER_PASSWORD,
        wetlab_config.SAMBA_SHARED_FOLDER_NAME,wetlab_config.SAMBA_REMOTE_SERVER_NAME,
        use_ntlm_v2=wetlab_config.SAMBA_NTLM_USED,domain=wetlab_config.SAMBA_DOMAIN)

    conn.connect(wetlab_config.SAMBA_IP_SERVER, int(wetlab_config.SAMBA_PORT_SERVER))

    return conn




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
    print ('config = ', config_file )
    print ('directorio ', os.getcwd())
    
    fileConfig(config_file)
    logger = logging.getLogger(logger_name)
    return logger

'''
def fetch_samba_dir_filelist(logger,conn, smb_root_path='/'):

    ## If no exceptions the function will leave a SMB connection opened so that the user can
    ##interact with SMB server
    #timestamp_print('Starting process to fetch the list of elements of directory via SAMBA')
    logger.info('Starting process to to fetch the list of elements of directory via SAMBA')
    try:
        file_list= conn.listPath(wetlab_config.SAMBA_SHARED_FOLDER_NAME,smb_root_path)

        logger.debug('number of existing folder elements excluding . and ..= '+str(len(file_list)-2))

    except:
        logger.error('==>Exception when accessing SMB with listPath')
        #timestamp_print('==>Exception when accessing SMB with listPath')
        raise

    #timestamp_print('Leaving the process to fetch the list of elements of directory via SAMBA')
    logger.info('Leaving the process to fetch the list of elements of directory via SAMBA')

    return file_list
'''

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

def get_machine_lanes(run_id):
    from iSkyLIMS_wetlab.models import RunProcess
    number_of_lanes = RunProcess.objects.get(pk=run_id).sequencerModel.get_number_of_lanes()
    return int(number_of_lanes)


def logging_errors(logger, string_text, showing_traceback):
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
    return ''



def get_miseq_run_cycles (run_parameter_file):
    '''
    Description:
        The function will find the number of cycles for a miseq run
    Input:
        logger_name    # contains the logger name that will be included 
                        in the log file
    Variables:
        found_cycles # re object to get the cycle value
        number_of_cycles # will add the partial number of cycles
    Return:
        number_of_cycles 
    '''
    number_of_cycles = 0
    with open (run_parameter_file, 'r') as fh:
        for line in lines :
            if '<RunInfoRead' in line :
                found_cycles = re.search('.*NumCycles="(\d+)".*',line)
                number_of_cycles += int(found_cycles.group(1))
    return number_of_cycles

def get_latest_miseq_log(conn, log_folder) :
    '''
    Description:
        The function will find the latest log file for the input folder
    Input:
        conn        # samba connection object
        log_folder  # folder name
    Constant:
        RUN_TEMP_DIRECTORY
    Variables:
        file_content # text information in the latest log
        file_remote # file name on the remote server
        remote_file_list # samba object with the list of files 
        max_cycle   # Counter with the maximum cycle found 
        number_of_cycles # will add the partial number of cycles
    Return:
        number_of_cycles 
        file_content
    '''
    remote_file_list = conn.listPath( wetlab_config.SAMBA_SHARED_FOLDER_NAME, log_folder)
    #run_folder = log_folder.split('/')[0]
    max_cycle = -1
    for sfh in remote_file_list:
        if sfh.isDirectory:
            continue
        file_remote = sfh.filename
        # if the run folder has been already process continue searching
        if file_remote.endswith('*.log') :
            log_file = re.search('.*_Cycle(\d+)_.*',file_remote)
            cycle_number = int(log_file.groups(1))
            if cycle_number > max_cycle :
                max_cycle = cycle_number
                latest_log = file_remote
    temporary_log = os.path.join(wetlab_config.RUN_TEMP_DIRECTORY,'miseq_cycle.log')
    with open(temporary_log ,'wb') as log_fp :
        try: # get RunParameter.xml if NextSeq
            s_latest_log = os.path.join(log_folder,latest_log)
            conn.retrieveFile(wetlab_config.SAMBA_SHARED_FOLDER_NAME, s_latest_log, log_fp)
        except Exeption as e :  
            return 'Error', 'Error' 
    with open (temporary_log, 'r') as fh :
        file_content = fh.read()
    os.remove(temporary_log)
    
    return max_cycle, file_content
