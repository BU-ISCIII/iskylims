'''

'''
## CSS file to be used for creating the PDF files
CSS_FOR_PDF = '/css/print_services.css'

## template files for generating the PDF files
REQUESTED_CONFIRMATION_SERVICE = 'request_service_template'
RESOLUTION_TEMPLATE = 'resolution_template.html'
OUTPUT_DIR_TEMPLATE ='documents/drylab/' # Directory to store the templates before moving to service folder

## SAMBA settings for connect to bioinfodoc server to create the folder request services
SAMBA_USER_ID = 'smonzon'
SAMBA_USER_PASSWORD = 's2r2m0nz0n'
SAMBA_SHARED_FOLDER_NAME = 'bioinfo_doc'
SAMBA_REMOTE_SERVER_NAME = 'panoramix'
SAMBA_NTLM_USED = True
SAMBA_IP_SERVER = '173.23.2.11'
SAMBA_PORT_SERVER = '139'

SAMBA_SERVICE_FOLDER = 'services'
## Folders to be created when service is accepted
FOLDERS_FOR_SERVICES = ['solicitud', 'resolucion', 'resultados']


