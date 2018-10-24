##General purpose utility functions used in different iSkyLIMS_wetlab modules
import datetime
import logging
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





