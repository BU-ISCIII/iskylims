
## (PARAMETERS FOR TESTING) SAMBA settings for connecting quibitka server to fetch the run files
SAMBA_USER_ID = 'bioinfocifs'
SAMBA_USER_PASSWORD = 'fCdEg979I-W.gUx-teDr'
SAMBA_SHARED_FOLDER_NAME = 'NGS_Data_test'
#SAMBA_SHARED_FOLDER_NAME = 'NGS_Data' ##commented out for testing
SAMBA_REMOTE_SERVER_NAME = 'panoramix'
#SAMBA_REMOTE_SERVER_NAME = 'quibitka'  ##commented out for testing
SAMBA_NTLM_USED = True
SAMBA_IP_SERVER = '172.21.7.11'
SAMBA_PORT_SERVER = '445'

'''
## SAMBA settings for connecting quibitka server to fetch the run files
SAMBA_USER_ID = 'bioinfocifs'
SAMBA_USER_PASSWORD = 'fCdEg979I-W.gUx-teDr'
SAMBA_SHARED_FOLDER_NAME = 'NGS_Data'
SAMBA_REMOTE_SERVER_NAME = 'quibitka'
SAMBA_NTLM_USED = True
SAMBA_IP_SERVER = '172.21.7.11'
SAMBA_PORT_SERVER = '445'
'''

## Directory settings for processing the run execution process
## Relative path from settings.BASE_DIR
LOG_DIRECTORY = 'logs/'
## Relative path from settings.MEDIA_ROOT
RUN_TEMP_DIRECTORY_RECORDED = 'wetlab/tmp/recorded/'
RUN_TEMP_DIRECTORY = 'wetlab/tmp'
RUN_TEMP_DIRECTORY_PROCESSING = 'wetlab/tmp/processing'
RUN_IMAGES_DIRECTORY = 'wetlab/images_plot'
RUN_SAMPLE_SHEET_DIRECTORY = 'wetlab/SampleSheets/'
PROCESSED_RUN_FILE='processed_run_file'
## file with MiSeq runs whose samplesheets fail sanity checks
UNEXPECTED_SAMPLESHEET_MISEQRUNS_FILE='unexpected_samplesheet_miseq_runs'


MIGRATION_DIRECTORY_FILES = 'wetlab/BaseSpaceMigrationFiles/'

## Configuration for the sample sheet conversion file to BaseSpace format
# column names when sample sheet has only one index
BASESPACE_FILE_ONE_INDEX = ['SampleID','Name','Species','Project','NucleicAcid',
               'Well','Index1Name','Index1Sequence']
# colum names when sample sheet has two index
BASESPACE_FILE_TWO_INDEX = ['SampleID','Name','Species','Project','NucleicAcid',
               'Well','Index1Name','Index1Sequence','Index2Name','Index2Sequence']
# mapping structure when sample sheet contains only one index
MAP_BASESPACE_SAMPLE_SHEET_ONE_INDEX = [('SampleID','Sample_ID'),('Name','Sample_Name'), ('Project','Sample_Project'),
                ('Index1Name','I7_Index_ID'), ('Index1Sequence','index' ) ]
# mapping structure when sample sheet contains two index
MAP_BASESPACE_SAMPLE_SHEET_TWO_INDEX = [('SampleID','Sample_ID'),('Name','Sample_Name'), ('Project','Sample_Project'),
                ('Index1Name','I7_Index_ID'), ('Index1Sequence','index' ),('Index2Name','I5_Index_ID'),('Index2Sequence','index2') ]
