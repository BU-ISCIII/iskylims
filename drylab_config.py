'''

'''
import sys
sys.path.append
try:
    from .drylab_email_conf import *
    EMAIL_USER_CONFIGURED = True
except:
    EMAIL_USER_CONFIGURED = False


SERVICE_MANAGER = 'ServiceManager'

INTERNAL_SEQUENCING_UNIT = "GENOMIC_SEQ_UNIT"

USER_REQUESTED_SERVICE_FILE_DIRECTORY ='drylab/servicesRequest'

## CSS file to be used for creating the PDF files
CSS_FOR_PDF = '/documents/drylab/services_templates/css/print_services.css'

## template files for generating the PDF files
REQUESTED_CONFIRMATION_SERVICE = 'request_service_template.html'
RESOLUTION_TEMPLATE = 'resolution_template.html'
OUTPUT_DIR_SERVICE_REQUEST_PDF ='documents/drylab/service_request' # Directory to store service request pdf
OUTPUT_DIR_RESOLUTION_PDF ='documents/drylab/resolution' # Directory to store resolution pdf

RESOLUTION_FILES_DIRECTORY = 'documents/drylab/resolutions'

ABBREVIATION_USED_FOR_SERVICE_REQUEST = 'SRV'
USER_CENTER_USED_WHEN_NOT_PROVIDED = 'NO_CENTER'

#SAMBA_SERVICE_FOLDER = 'services'
## Folders to be created when service is accepted
#FOLDERS_FOR_SERVICES = ['request', 'resolution', 'result'] # 0= request, 1= resolution, 2 = result (keep order as suggested)
#RESOLUTION_PREFIX = 'Resolution_'

############## FOLDER SETTINGS ###############################
## Directory settings for processing the run data files ######
## Relative path from settings.BASE_DIR
LOG_DIRECTORY = 'logs/'


################# CONFIG FILE LOG NAME ###############################
LOGGING_CONFIG_FILE = 'logging_config.ini'

CONFIRMATION_TEXT_MESSAGE = ['Your service request has been successfully recorded.',
                    'The sequence number assigned for your request is: SERVICE_NUMBER', 'Keep this number safe for refering your request',
                    'You will be contacted shortly.']


################# ERROR  #########################
ERROR_USER_NOT_ALLOWED = ['You do not have the enough privileges to see this page ','Contact with your administrator .']
ERROR_UNABLE_TO_RECORD_YOUR_SERVICE = ['Your service request cannot be recorded.', 'Check that all information is provided correctly.',
                                'If problem persist, contact your administrator']

ERROR_SERVICE_ID_NOT_FOUND = ['Service ID not found']
'''
ERROR_WRONG_SAMBA_CONFIGURATION_SETTINGS = ['Unsuccessful configuration settings for Samba connection']
ERROR_UNABLE_TO_SAVE_SAMBA_CONFIGURATION_SETTINGS = ['Unable to save the Samba configuration file ', 'check if folder iSkyLIMS_wetlab has write permision for apache user']
ERROR_UNABLE_TO_SAVE_EMAIL_CONFIGURATION_SETTINGS = ['Unable to save the email configuration file ', 'check if folder iSkyLIMS_wetlab has write permision for apache user']
'''
ERROR_PIPELINE_ALREADY_EXISTS = ['Pipeline name and version is already defined']

ERROR_NO_SERVICES_ARE_SELECTED = ['Unable to process your request', 'No services have been selected']

HEADING_SELECT_SAMPLE_IN_SERVICE = ['Run Name', 'Run ID', 'Project Name',  'Project ID','Sample Name', 'sample id', 'Run finish date']

HEADING_ADDITIONAL_PIPELINE_PARAMETERS = ['Additional parameter name', 'Additional Parameter Value']

HEADING_SERVICE_DATES = ['Service Date Creation', 'Approval Service Date', 'Rejected Service Date']

HEADING_MANAGE_PIPELINES = ['Service', 'User' ,'Pipeline Name', 'Pipeline Version', 'Date', 'Default', 'In use', 'id']

HEADING_ADDITIONAL_RESOLUTION_PARAMETERS = ['Parameter name', 'Parameter value', 'Notes']

MAPPING_ADDITIONAL_RESOLUTION_PARAMETERS = [('resolutionParameter', 'Parameter name'),('resolutionParamValue', 'Parameter value'),('resolutionParamNotes', 'Notes')]

DATE_NOT_YET_DEFINED = 'Not Yet Defined'

DISPLAY_NEW_DEFINED_PIPELINE = ['Service', 'Pipeline Name' , 'Pipeline Version']
DISPLAY_MULTYPLE_DEFINED_PIPELINE = ['Service', 'User', 'Pipeline Name' , 'Pipeline Version', 'Date', 'Default', 'In Use']

DISPLAY_DETAIL_PIPELINE_BASIC_INFO = ['Pipeline Name', 'Pipeline Version', 'Service']
DISPLAY_DETAIL_PIPELINE_ADDITIONAL_INFO = ['User', 'Creation Date', 'String Folder', 'Default', 'In Use', 'Automatic']

################ EMAIL TEXT   ##################################
SUBJECT_SERVICE_RECORDED = ['Service ', ' has been recorded']
BODY_SERVICE_RECORDED = ['Dear USER_NAME','Your service  SERVICE_NUMBER has been recorded.','You will received the resolution of the request as soon as possible.',
                    'Kind regards', 'BU-ISCIII']

SUBJECT_RESOLUTION_RECORDED = ['Resolution ', ' has been updated']
BODY_RESOLUTION_ACCEPTED = ['Dear  USER_NAME','A new resolution has been added for your service:  SERVICE_NUMBER.', 'Your service has been STATUS',
                    'and your delivery estimated date is DATE','Your service is now queued and you will be notified when it is updated',
                     'Kind regards', 'BU-ISCIII']
BODY_RESOLUTION_REJECTED = ['Dear  USER_NAME','A new resolution has been added for your service:  SERVICE_NUMBER.', 'Your service has been STATUS',
                    'because it does not fullfil our requirements or is not in our services portfolio. If you have any question please contact us.',
                     'Kind regards', 'BU-ISCIII']


SUBJECT_SERVICE_ON_QUEUED = ['Service ', 'sent to preparation pipelines Jobs']
BODY_SERVICE_ON_QUEUED = ['Service  SERVICE_NUMBER is on queued ']

################ EMAIL CONFIGURATION FIELDS ###############################
EMAIL_CONFIGURATION_FIELDS = ['USER_NAME', 'USER_EMAIL','SENT_EMAIL_ON_ERROR']
EMAIL_CONFIGURATION_FILE_HEADING = '############# EMAIL CONFIGURATION FILE ########\n#DO NOT MODIFY MANUALLY THIS FILE\n#VALUES WILL BE MODIFIED WHEN USING THE CONFIGURATION FORM\n'
EMAIL_CONFIGURATION_FILE_END = '########## END EMAIL CONFIGURATION FILE'
