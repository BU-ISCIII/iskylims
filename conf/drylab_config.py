'''

'''
## CSS file to be used for creating the PDF files
CSS_FOR_PDF = '/documents/drylab/services_templates/css/print_services.css'

## template files for generating the PDF files
REQUESTED_CONFIRMATION_SERVICE = 'request_service_template.html'
RESOLUTION_TEMPLATE = 'resolution_template.html'
OUTPUT_DIR_TEMPLATE ='documents/drylab/' # Directory to store the pdf templates before moving to service folder

## SAMBA settings for connect to bioinfodoc server to create the folder request services
SAMBA_USER_ID = ''
SAMBA_USER_PASSWORD = ''
SAMBA_SHARED_FOLDER_NAME = ''
SAMBA_REMOTE_SERVER_NAME = ''
SAMBA_NTLM_USED = True
SAMBA_DOMAIN = ''
SAMBA_IP_SERVER = ''
SAMBA_PORT_SERVER = '445'

SAMBA_SERVICE_FOLDER = 'services'
## Folders to be created when service is accepted
FOLDERS_FOR_SERVICES = ['request', 'resolution', 'result'] # 0= request, 1= resolution, 2 = result (keep order as suggested)
#RESOLUTION_PREFIX = 'Resolution_'


