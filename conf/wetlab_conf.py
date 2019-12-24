################# USER SETTINGS ##############################
## Wetlab manager
WETLAB_MANAGER = 'WetlabManager'
##############################################################
################### SAMBA SETTINGS  ##########################
## SAMBA settings for connecting quibitka server to fetch the run files
SAMBA_USER_ID = 'luigi'
SAMBA_USER_PASSWORD = 'Apple123'
SAMBA_SHARED_FOLDER_NAME = 'NGS_Data'
#    Write the subfolder name in case that run folder are not under the
#    shared folder directory. Leave empty in other case
SAMBA_APPLICATION_FOLDER_NAME = ''
SAMBA_REMOTE_SERVER_NAME = 'Luigi-PC'
SAMBA_NTLM_USED = True
SAMBA_IP_SERVER = '192.168.1.3'
SAMBA_PORT_SERVER = '445'
SAMBA_HOST_NAME = ''
IS_DIRECT_TCP = True
## SAMBA_DOMAIN MUST be empty if domain value is not used for samba connection
SAMBA_DOMAIN=''
################### EMAIL SETTINGS  ##########################
SENT_EMAIL_ON_ERROR = False
TO_EMAIL_ADDRESS = ['bioinformatica@isciii.es']
FROM_EMAIL_ADDRESS = 'iSkyLIMS@isciii.es'
##############################################################
