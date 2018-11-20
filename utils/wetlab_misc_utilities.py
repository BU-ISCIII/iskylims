##General purpose utility functions used in different iSkyLIMS_wetlab modules
## BE CAREFUL if modifying these resources as they are used by other functions in wetlab


import datetime
import logging
import os, errno
from smb.SMBConnection import SMBConnection
from iSkyLIMS_wetlab import wetlab_config

def timestamp_print(message):
    starting_time= datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(starting_time+' '+message)




def open_samba_connection():

# Lines for 'cuadrix' testing case:
#    timestamp_print('Starting the process for open_samba_connection() (cron.py- cuadrix testing)')
#    logger.info('user ID= '+wetlab_config.SAMBA_USER_ID+'. domain= '+wetlab_config.SAMBA_DOMAIN)
#     conn=SMBConnection(wetlab_config.SAMBA_USER_ID, wetlab_config.SAMBA_USER_PASSWORD,
#         wetlab_config.SAMBA_SHARED_FOLDER_NAME,wetlab_config.SAMBA_REMOTE_SERVER_NAME,
#         use_ntlm_v2=wetlab_config.SAMBA_NTLM_USED,domain=wetlab_config.SAMBA_DOMAIN)

    conn=SMBConnection(wetlab_config.SAMBA_USER_ID, wetlab_config.SAMBA_USER_PASSWORD,
        wetlab_config.SAMBA_SHARED_FOLDER_NAME,wetlab_config.SAMBA_REMOTE_SERVER_NAME,
        use_ntlm_v2=wetlab_config.SAMBA_NTLM_USED)
    if True != conn.connect(wetlab_config.SAMBA_IP_SERVER, int(wetlab_config.SAMBA_PORT_SERVER)):
        logger=open_log('open_samba_connection_testing.log')
        logger.error('Cannot set up SMB connection with '+wetlab_config.SAMBA_REMOTE_SERVER_NAME)
        timestamp_print('Cannot set up SMB connection with '+wetlab_config.SAMBA_REMOTE_SERVER_NAME)

    timestamp_print('Leaving open_samba_connection() (cron.py- cuadrix testing)')
    return conn



def fetch_samba_dir_filelist(logger,conn, smb_root_path='/'):
    ## By default, it returns the contents of the root ==> a list of the run directories (
    ## + '.' and '..')
    ## If no exceptions the function will leave a SMB connection opened so that the user can
    ##interact with SMB server
    timestamp_print('Starting process to fetch the list of elements of directory via SAMBA')
    logger.info('Starting process to to fetch the list of elements of directory via SAMBA')
    file_list=[]
    try:
        file_list= conn.listPath(wetlab_config.SAMBA_SHARED_FOLDER_NAME,smb_root_path)
        #file_list=(file_list[2:4]).append(file_list[9:10]) ##Debug
        #file_list=file_list[2:4] ##Debug
        file_list_filenames_debug=[x.filename for x in file_list] ##debug
        logger.debug(
            'number of existing folder elements excluding . and ..= '+str(len(file_list)-2))
        #logger.debug('run dir list=\n'+'\n'.join(file_list_filenames_debug))

    except: ##
        logger.error('==>Exception when accessing SMB with listPath')
        timestamp_print('==>Exception when accessing SMB with listPath')
        raise

    timestamp_print('Leaving the process to fetch the list of elements of directory via SAMBA')
    logger.info('Leaving the process to fetch the list of elements of directory via SAMBA')

    return file_list






