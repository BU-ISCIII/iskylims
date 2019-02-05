import os
from django.conf import settings
################# USER SETTINGS ##############################
## Wetlab manager
WETLAB_MANAGER = 'WetlabManager'
##############################################################
################### SAMBA SETTINGS  ##########################
## SAMBA settings for connecting quibitka server to fetch the run files
SAMBA_USER_ID = ''
SAMBA_USER_PASSWORD = ''
SAMBA_SHARED_FOLDER_NAME = ''
#    Write the subfolder name in case that run folder are not under the
#    shared folder directory. Leave empty in other case
SAMBA_APPLICATION_FOLDER_NAME = ''
SAMBA_REMOTE_SERVER_NAME = ''
SAMBA_NTLM_USED = True
SAMBA_IP_SERVER = ''
SAMBA_PORT_SERVER = '445'
## SAMBA_DOMAIN MUST be empty if domain value is not used for samba connection
SAMBA_DOMAIN=''
##############################################################

################### EMAIL SETTINGS  ##########################
SENT_EMAIL_ON_ERROR = False
TO_EMAIL_ADDRESS = ['bioinformatica@isciii.es']
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
LOG_NAME_FIND_UPDATE_RUNS = 'update_run'

################# CONFIG FILE LOG NAME ###############################
LOGGING_CONFIG_FILE = 'logging_config.ini'


##############################################################

## FILE NAME CONTAINING THE PROCESSED RUNS ###################
PROCESSED_RUN_FILE='processed_run_file'
##############################################################

############# ILLUMINA OUTPUT FILES ##########################
RUN_PARAMETER_NEXTSEQ = 'RunParameters.xml'
#RUN_PARAMETER_MISEQ = 'runParameters.xml'
RUN_INFO = 'RunInfo.xml'
RUN_COMPLETION = 'RunCompletionStatus.xml'
SAMPLE_SHEET = 'samplesheet.csv'
## sample sheet to be copied on the remote folder
COPY_SAMPLE_SHEET_TO_REMOTE = False # boolean constant True if NestSeq
                                # sample sheet needs to be copied to remote server
RUN_LOG_FOLDER = 'Logs'

DEMULTIPLEXION_BCL2FASTQ_FOLDER = 'Data/Intensities/BaseCalls'
REPORT_FOLDER = 'Reports'
STATS_FOLDER = 'Stats'


CONVERSION_STATS_FILE = 'ConversionStats.xml'
DEMULTIPLEXION_STATS_FILE = 'DemultiplexingStats.xml'
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

############ MAXIMUM TIME TO WAIT BEFORE MOVING TO ERROR #####
MAXIMUM_TIME_WAIT_SAMPLE_SHEET = '2' # in days
MAXIMUM_TIME_WAIT_RUN_COMPLETION = '2' # in days
##############################################################

############ RUN METRIC FOLDERS AND FILES ####################
INTEROP_PATH = '/opt/interop/bin/'
RUN_METRIC_FOLDER = 'InterOp'
PLOT_EXTENSION = '.png'

##############################################################

#########  FOLDER SETTINGS FOR PROCESSING LIBRARY KITS #######
## Relative path from settings.BASE_DIR
LIBRARY_KITS_DIRECTORY = 'wetlab/library_kits/'
## Maximum file size allowed for the index library kits (in bytes)
LIBRARY_MAXIMUM_SIZE = '3145728'
## Configuration for index library file
INDEX_LIBRARY_HEADING = ['[Version]','[Name]', '[PlateExtension]','[Settings]', '[I7]','[I5]']
##############################################################

MIGRATION_DIRECTORY_FILES = 'wetlab/BaseSpaceMigrationFiles/'


##### CONFIGURATION FOR SAMPLE SHEET CONVERSION ##############
######### FILE TO BASESPACE FORMAT     #######################
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

##############################################################

##### SETTINGS FOR CLEANUP RUNS AND NOT VALID FILES ##########

RETENTION_TIME = '7' # in days

##############################################################

