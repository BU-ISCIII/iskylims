'''

'''
import sys
sys.path.append
try:
    from .drylab_samba_conf import *
    SAMBA_USER_CONFIGURED = True
except:
    SAMBA_USER_CONFIGURED = False
try:
    from .drylab_email_conf import *
    EMAIL_USER_CONFIGURED = True
except:
    EMAIL_USER_CONFIGURED = False


DRYLAB_MANAGER = 'ServiceManager'

INTERNAL_SEQUENCING_UNIT = "GENOMIC_SEQ_UNIT"

## CSS file to be used for creating the PDF files
CSS_FOR_PDF = '/documents/drylab/services_templates/css/print_services.css'

## template files for generating the PDF files
REQUESTED_CONFIRMATION_SERVICE = 'request_service_template.html'
RESOLUTION_TEMPLATE = 'resolution_template.html'
OUTPUT_DIR_TEMPLATE ='documents/drylab/' # Directory to store the pdf templates before moving to service folder


SAMBA_SERVICE_FOLDER = 'services'
## Folders to be created when service is accepted
FOLDERS_FOR_SERVICES = ['request', 'resolution', 'result'] # 0= request, 1= resolution, 2 = result (keep order as suggested)
#RESOLUTION_PREFIX = 'Resolution_'

CONFIRMATION_TEXT_MESSAGE = ['Your service request has been successfully recorded.',
                    'The sequence number assigned for your request is: SERVICE_NUMBER', 'Keep this number safe for refering your request',
                    'You will be contacted shortly.']

ERROR_USER_NOT_ALLOWED = ['You do have the enough privileges to see this page ','Contact with your administrator .']
ERROR_UNABLE_TO_RECORD_YOUR_SERVICE = ['Your service request cannot be recorded.', 'Check that all information is provided correctly.',
                                'If problem persist, contact your administrator']

ERROR_WRONG_SAMBA_CONFIGURATION_SETTINGS = ['Unsuccessful configuration settings for Samba connection']
ERROR_UNABLE_TO_SAVE_SAMBA_CONFIGURATION_SETTINGS = ['Unable to save the Samba configuration file ', 'check if folder iSkyLIMS_wetlab has write permision for apache user']
ERROR_UNABLE_TO_SAVE_EMAIL_CONFIGURATION_SETTINGS = ['Unable to save the email configuration file ', 'check if folder iSkyLIMS_wetlab has write permision for apache user']

ERROR_PIPELINE_ALREADY_EXISTS = ['Pipeline name and version is already defined']

AVAILABLE_ACTIONS_IN_PIPELINE = ['Copy','Symbolic link',]
HEADING_ACTIONS_PIPELINES = ['Given name for action', 'Order', 'Action', 'Fake Action']
HEADING_ACTIONS_PARAMETERS = ['Parameter1', 'Parameter2', 'Parameter3']


HEADING_MANAGE_PIPELINES = ['Service', 'User' ,'Pipeline Name', 'Pipeline Version', 'Date', 'Default', 'In use', 'id']

DISPLAY_NEW_DEFINED_PIPELINE = ['Service', 'Pipeline Name' , 'Pipeline Version']
DISPLAY_MULTYPLE_DEFINED_PIPELINE = ['Service', 'User', 'Pipeline Name' , 'Pipeline Version', 'Date', 'Default', 'In Use']

DISPLAY_DETAIL_PIPELINE_BASIC_INFO = ['Pipeline Name', 'Pipeline Version', 'Service']
DISPLAY_DETAIL_PIPELINE_ADDITIONAL_INFO = ['User', 'Creation Date', 'String Folder', 'Default', 'In Use', 'Automatic']

################ EMAIL TEXT   ##################################
SUBJECT_SERVICE_RECORDED = ['Service ', ' has been recorded']
BODY_SERVICE_RECORDED = ['Dear ','Your service  SERVICE_NUMBER has been recorded.','You will received the resolution of the request as soon as possible.',
                    'Kind regards', 'BU-ISCIII']

################ EMAIL CONFIGURATION FIELDS ###############################
EMAIL_CONFIGURATION_FIELDS = ['USER_NAME', 'USER_EMAIL','SENT_EMAIL_ON_ERROR']
EMAIL_CONFIGURATION_FILE_HEADING = '############# EMAIL CONFIGURATION FILE ########\n#DO NOT MODIFY MANUALLY THIS FILE\n#VALUES WILL BE MODIFIED WHEN USING THE CONFIGURATION FORM\n'
EMAIL_CONFIGURATION_FILE_END = '########## END EMAIL CONFIGURATION FILE'


################ SAMBA CONFIGURATION FIELDS ###############################
SAMBA_CONFIGURATION_FIELDS = ['SAMBA_USER_ID', 'SAMBA_USER_PASSWORD', 'SAMBA_SHARED_FOLDER_NAME', 'SAMBA_APPLICATION_FOLDER_NAME', 'SAMBA_REMOTE_SERVER_NAME',
                'SAMBA_NTLM_USED', 'SAMBA_IP_SERVER', 'SAMBA_HOST_NAME', 'SAMBA_PORT_SERVER', 'IS_DIRECT_TCP', 'SAMBA_DOMAIN']

SAMBA_CONFIGURATION_FILE_HEADING = '############# SAMBA CONFIGURATION FILE ########\n#DO NOT MODIFY MANUALLY THIS FILE\n#VALUES WILL BE MODIFIED WHEN USING THE CONFIGURATION FORM\n'
SAMBA_CONFIGURATION_FILE_END = '########## END SAMBA CONFIGURATION FILE'
