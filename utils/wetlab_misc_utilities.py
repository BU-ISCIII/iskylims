##General purpose utility functions used in different iSkyLIMS_wetlab modules

import datetime
import logging
from logging.config import fileConfig
from logging.handlers import RotatingFileHandler
import os, errno
from smb.SMBConnection import SMBConnection
from django.conf import settings
from iSkyLIMS_wetlab import wetlab_config


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
    number_of_lanes = RunProcess.objects.get(pk=run_id).sequencerModel.get_number_of_lanes()
    return int(number_of_lanes)

def logging_exception_errors(logger, error, string_text):
    logger.error('-----------------    ERROR   ------------------')
    logger.error(string_text , error )
    logger.error('Showing traceback: ',  exc_info=True)
    logger.error('-----------------    END ERROR   --------------')
    return ''


