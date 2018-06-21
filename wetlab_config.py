

## SAMBA settings for connect to quibitka server to connect to the run folders
SAMBA_USER_ID = 'bioinfocifs'
SAMBA_USER_PASSWORD = 'fCdEg979I-W.gUx-teDr'
SAMBA_SHARED_FOLDER_NAME = 'NGS_Data'
SAMBA_REMOTE_SERVER_NAME = 'quibitka'
SAMBA_NTLM_USED = True
SAMBA_IP_SERVER = '172.21.7.11'
SAMBA_PORT_SERVER = '445'

## Directory settings for processing the run execution process
## Relative path from settings.BASE_DIR
LOG_DIRECTORY = 'logs/'
## Relative path from settings.MEDIA_ROOT
RUN_TEMP_DIRECTORY_RECORDED = 'wetlab/tmp/recorded/'
RUN_TEMP_DIRECTORY = 'wetlab/tmp'
RUN_TEMP_DIRECTORY_PROCESSING = 'wetlab/tmp/processing'
RUN_IMAGES_DIRECTORY = 'wetlab/images_plot'
RUN_SAMPLE_SHEET_DIRECTORY = 'wetlab/SampleSheets/'


MIGRATION_DIRECTORY_FILES = 'wetlab/BaseSpaceMigrationFiles/'
