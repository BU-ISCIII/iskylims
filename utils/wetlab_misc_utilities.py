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




#TBD
def open_samba_connection():
    ## needed for testing in cuadrix
    ## to be commented-out when delivery


    timestamp_print('Starting the process for open_samba_connection() (cron.py- cuadrix testing)')

    ###logger.info('user ID= '+wetlab_config.SAMBA_USER_ID+'. domain= '+wetlab_config.SAMBA_DOMAIN)
    conn=SMBConnection(wetlab_config.SAMBA_USER_ID, wetlab_config.SAMBA_USER_PASSWORD,
        wetlab_config.SAMBA_SHARED_FOLDER_NAME,wetlab_config.SAMBA_REMOTE_SERVER_NAME,
        use_ntlm_v2=wetlab_config.SAMBA_NTLM_USED,domain=wetlab_config.SAMBA_DOMAIN)
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
    timestamp_print('Starting process to fetch the directory list via SAMBA')
    logger.info('Starting process to fetch the directory list via SAMBA')
    file_list=[]
    try:
        file_list= conn.listPath(wetlab_config.SAMBA_SHARED_FOLDER_NAME,smb_root_path)
        #file_list=file_list[2:3] ##TBDDebugEndDebug
        file_list_filenames_debug=[x.filename for x in file_list] ##debug
        logger.debug(
            'number of existing directory runs of any kind= '+str(len(file_list)-2))## -2 ->"." and ".."
        #logger.debug('run dir list=\n'+'\n'.join(file_list_filenames_debug))

    except: ##
        logger.error('==>Exception when accessing SMB with listPath')
        timestamp_print('==>Exception when accessing SMB with listPath')
        conn.close()
        raise

    timestamp_print('Leaving the process to fetch the run-directory list via SAMBA')
    logger.info('Leaving the process to fetch the run-directory list via SAMBA')

    return file_list

'''
def open_samba_connection():
    ## open samba connection
    # There will be some mechanism to capture userID, password, client_machine_name, server_name and server_ip
    # client_machine_name can be an arbitary ASCII string
    # server_name should match the remote machine name, or else the connection will be rejected
    conn=SMBConnection(wetlab_config.SAMBA_USER_ID, wetlab_config.SAMBA_USER_PASSWORD, wetlab_config.SAMBA_SHARED_FOLDER_NAME,wetlab_config.SAMBA_REMOTE_SERVER_NAME, use_ntlm_v2=wetlab_config.SAMBA_NTLM_USED)
    conn.connect(wetlab_config.SAMBA_IP_SERVER, int(wetlab_config.SAMBA_PORT_SERVER))
    #conn=SMBConnection('bioinfocifs', 'fCdEg979I-W.gUx-teDr', 'NGS_Data', 'quibitka', use_ntlm_v2=True)
    #conn.connect('172.21.7.11', 445)

    #conn=SMBConnection('Luigi', 'Apple123', 'NGS_Data_test', 'LUIGI-PC', use_ntlm_v2=True)
    #conn.connect('192.168.1.3', 139)
    #conn=SMBConnection('bioinfocifs', 'bioinfocifs', 'NGS_Data_test', 'barbarroja', use_ntlm_v2=True)
    #conn.connect('10.15.60.54', 139)


    ###conn = SMBConnection(userid, password, client_machine_name, remote_machine_name, use_ntlm_v2 = True)
    ###conn.connect(server_ip, 139)
    return conn
'''
### End TBD





