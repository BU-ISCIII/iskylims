import os
from django.conf import settings
################# USER SETTINGS ##############################
## Wetlab manager 
WETLAB_MANAGER = 'WetlabManager'
##############################################################

################### SAMBA SETTINGS  ##########################
## SAMBA settings for connecting quibitka server to fetch the run files
SAMBA_USER_ID = 'Luigi'
SAMBA_USER_PASSWORD = 'Apple123'
SAMBA_SHARED_FOLDER_NAME = 'NGS_Data'
SAMBA_REMOTE_SERVER_NAME = 'LUIGI-PC'
SAMBA_NTLM_USED = True
SAMBA_IP_SERVER = '192.168.1.3'
SAMBA_PORT_SERVER = '139'

################### SAMBA SETTINGS  ##########################
## SAMBA settings for connecting quibitka server to fetch the run files
#SAMBA_USER_ID = 'bioinfocifs'
#SAMBA_USER_PASSWORD = 'fCdEg979I-W.gUx-teDr'
#SAMBA_SHARED_FOLDER_NAME = 'NGS_Data'
#SAMBA_REMOTE_SERVER_NAME = 'quibitka'
#SAMBA_NTLM_USED = True
#SAMBA_IP_SERVER = '172.21.7.11'
#SAMBA_PORT_SERVER = '445'
## SAMBA_DOMAIN MUST be empty if domain value is not used for samba connection
SAMBA_DOMAIN=''
##############################################################

################### EMAIL SETTINGS  ##########################
SENT_EMAIL_ON_ERROR = False
TO_EMAIL_ADDRESS = ['bioinfo@isciii.es']
FROM_EMAIL_ADDRESS = 'iSkyLIMS@isciii.es'
##############################################################


############## FOLDER SETTINGS ###############################
## Directory settings for processing the run data files ######
## Relative path from settings.BASE_DIR
LOG_DIRECTORY = 'logs/'
## Relative path from settings.MEDIA_ROOT
RUN_TEMP_DIRECTORY_RECORDED = 'wetlab/tmp/recorded/'
RUN_TEMP_DIRECTORY = 'wetlab/tmp'
RUN_TEMP_DIRECTORY_PROCESSING = 'wetlab/tmp/processing'
RUN_IMAGES_DIRECTORY = 'wetlab/images_plot'
RUN_SAMPLE_SHEET_DIRECTORY = 'wetlab/SampleSheets/'

##############################################################

################# LOG NAMES #############################
LOG_NAME_RUN_IN_RECORDED_STATE = 'run_in_recorded'
LOG_NAME_MISEQ_FETCH_SAMPLE_SHEET = 'miseq_sample_sheet'

################# CONFIG FILE LOG NAME ###############################
LOGGING_CONFIG_FILE = 'logging_config.ini'


##############################################################

## FILE NAME CONTAINING THE PROCESSED RUNS ###################
PROCESSED_RUN_FILE='processed_run_file'
##############################################################

############# ILLUMINA OUTPUT FILES ##########################
RUN_PARAMETER_NEXTSEQ = 'RunParameters.xml'
RUN_PARAMETER_MISEQ = 'runParameters.xml'
RUN_INFO = 'RunInfo.xml'
RUN_COMPLETION = 'RunCompletionStatus.xml'
SAMPLE_SHEET = 'SampleSheet.csv'
## sample sheet to be copied on the remote folder 
COPY_SAMPLE_SHEET_TO_REMOTE = False # boolean constant True if NestSeq 
                                # sample sheet needs to be copied to remote server

##############################################################

############ VALUE TAG FOR XML FILES #########################
COMPLETION_TAG = 'CompletionStatus'
COMPLETION_SUCCESS = 'CompletedAsPlanned'
EXPERIMENT_NAME_TAG = 'ExperimentName'
APPLICATION_NAME_TAG = 'ApplicationName'
##############################################################

############ DEFAULT VALUES FOR MISEQ SAMPLE SHEET  ##########
DEFAULT_LIBRARY_KIT = 'Unknown'
DEFAULT_CENTER = 'CNM'

##############################################################
## 


'''
MISEQ_PROCESSED_RUN_FILE='miseq_processed_run_file'
MISEQ_PROCESSED_RUN_FILEPATH=os.path.join(settings.MEDIA_ROOT,
    RUN_TEMP_DIRECTORY,MISEQ_PROCESSED_RUN_FILE)




PROCESSED_RUN_FILEPATH=os.path.join(settings.MEDIA_ROOT,
    RUN_TEMP_DIRECTORY,PROCESSED_RUN_FILE)

## file with MiSeq runs whose samplesheets fail sanity checks
FAULTY_SAMPLESHEET_MISEQRUNS_FILE='faulty_samplesheet_miseq_runs'
FAULTY_SAMPLESHEET_MISEQRUNS_FILEPATH=os.path.join(settings.MEDIA_ROOT,
    RUN_TEMP_DIRECTORY,FAULTY_SAMPLESHEET_MISEQRUNS_FILE)

## file containing MiSeq runs in RECORDED state to check whether they have
## already reached SAMPLE SENT
RECORDED_MISEQRUNS_FILE='recorded_miseq_runs'
RECORDED_MISEQRUNS_FILEPATH=os.path.join(settings.MEDIA_ROOT,
    RUN_TEMP_DIRECTORY,RECORDED_MISEQRUNS_FILE)

'''

## Directory settings for processing the library kits
## Relative path from settings.BASE_DIR
LIBRARY_KITS_DIRECTORY = 'wetlab/library_kits/'
## Maximum file size allowed for the index library kits (in bytes)
LIBRARY_MAXIMUM_SIZE = '3145728'
## Configuration for index library file 
INDEX_LIBRARY_HEADING = ['[Version]','[Name]', '[PlateExtension]','[Settings]', '[I7]','[I5]']


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

