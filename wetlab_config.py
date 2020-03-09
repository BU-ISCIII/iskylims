import os
from django.conf import settings
##### Allow to import the configuration samba files from configuration folder
import sys
sys.path.append('../')
from conf.wetlab_conf import *
'''
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
'''
'''
## SAMBA settings for connecting quibitka server to fetch the run files
SAMBA_USER_ID = 'bioinfocifs'
SAMBA_USER_PASSWORD = 'fCdEg979I-W.gUx-teDr'
SAMBA_SHARED_FOLDER_NAME = 'NGS_Data'
#    Write the subfolder name in case that run folder are not under the
#    shared folder directory. Leave empty in other case
SAMBA_APPLICATION_FOLDER_NAME = ''
SAMBA_REMOTE_SERVER_NAME = 'galera'
SAMBA_NTLM_USED = True
SAMBA_PORT_SERVER = '445'
SAMBA_HOST_NAME = 'galera.isciii.es'
IS_DIRECT_TCP=True
## SAMBA_DOMAIN MUST be empty if domain value is not used for samba connection
SAMBA_DOMAIN='ISCIII
'''


##############################################################
'''
################### EMAIL SETTINGS  ##########################
SENT_EMAIL_ON_ERROR = False
TO_EMAIL_ADDRESS = ['bioinformatica@isciii.es']
FROM_EMAIL_ADDRESS = 'iSkyLIMS@isciii.es'
##############################################################
'''

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
#SAMPLE_SHEET_CREATED_ON_LAB = 'wetlab/SampleSheetsFromLab'
TEMPLATE_FILES_DIRECTORY ='wetlab/templates'

## Directory to store the imported user sampleSheets
LIBRARY_PREPARATION_SAMPLE_SHEET_DIRECTORY = 'wetlab/SampleSheets4LibPrep/'

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
COPY_SAMPLE_SHEET_TO_REMOTE = False # boolean constant True if NextSeq
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
LIBRARY_KITS_DIRECTORY = 'wetlab/collection_index_kits/'
COLLECTION_INDEX_KITS_DIRECTORY ='wetlab/collection_index_kits/'
## Maximum file size allowed for the index library kits (in bytes)
LIBRARY_MAXIMUM_SIZE = '3145728'
## Configuration for index library file
COLLECTION_INDEX_HEADING = ['[Version]','[Name]', '[PlateExtension]','[Settings]', '[I7]']
##############################################################

MIGRATION_DIRECTORY_FILES = 'wetlab/BaseSpaceMigrationFiles/'

##############################################################
BASE_SPACE_TWO_INDEX_TEMPLATE_NAME = 'base_space_two_index_template.csv'
BASE_SPACE_ONE_INDEX_TEMPLATE_NAME = 'base_space_one_index_template.csv'
SAMPLE_SHEET_TWO_INDEX_TWO_ADAPTERS_TEMPLATE_NAME = 'sample_sheet_two_index_two_adapters_template.csv'
SAMPLE_SHEET_TWO_INDEX_ONE_ADAPTER_TEMPLATE_NAME = 'sample_sheet_two_index_one_adapter_template.csv'
SAMPLE_SHEET_ONE_INDEX_TWO_ADAPTERS_TEMPLATE_NAME = 'sample_sheet_one_index_two_adapters_template.csv'
SAMPLE_SHEET_ONE_INDEX_ONE_ADAPTER_TEMPLATE_NAME = 'sample_sheet_one_index_one_adapter_template.csv'


##### CONFIGURATION FOR SAMPLE SHEET CONVERSION ##############
######### FILE TO BASESPACE FORMAT     #######################
# column names when sample sheet has only one index
BASESPACE_FILE_ONE_INDEX = ['SampleID','Name','Species','Project','NucleicAcid',
               'Well','Index1Name','Index1Sequence']
# colum names when sample sheet has two index
BASESPACE_FILE_TWO_INDEX = ['SampleID','Name','Species','Project','NucleicAcid',
               'Well','Index1Name','Index1Sequence','Index2Name','Index2Sequence']

HEADING_FOR_SAMPLE_SHEET_ONE_INDEX = ['Unique_Sample_ID','Sample_Name','Sample_Plate','Sample_Well','Index_Plate_Well','I7_Index_ID','index','Sample_Project','Description']

HEADING_FOR_SAMPLE_SHEET_TWO_INDEX = ['Unique_Sample_ID','Sample_Name','Sample_Plate','Sample_Well','Index_Plate_Well','I7_Index_ID','index','I5_Index_ID','index2','Sample_Project','Description']


# mapping structure when sample sheet contains only one index
MAP_BASESPACE_SAMPLE_SHEET_ONE_INDEX = [('SampleID','Sample_ID'),('Name','Sample_Name'), ('Project','Sample_Project'),
                ('Index1Name','I7_Index_ID'), ('Index1Sequence','index' ) ]
# mapping structure when sample sheet contains two index
MAP_BASESPACE_SAMPLE_SHEET_TWO_INDEX = [('SampleID','Sample_ID'),('Name','Sample_Name'), ('Project','Sample_Project'),
                ('Index1Name','I7_Index_ID'), ('Index1Sequence','index' ),('Index2Name','I5_Index_ID'),('Index2Sequence','index2') ]

######### MAPPIING COLUMNS IN SAMPLE SHEET FROM USER TO DATABASE   #############
MAP_USER_SAMPLE_SHEET_TO_DATABASE_TWO_INDEX = [('Sample_ID','userSampleID'), ('Sample_Name','sample_name'), ('Sample_Plate','samplePlate'),
            ('Sample_Well','sampleWell'),('Index_Plate_Well','indexPlateWell'),('I7_Index_ID','i7IndexID'),
            ('index','i7Index'),('I5_Index_ID','i5IndexID'),('index2','i5Index'),('Sample_Project','projectInSampleSheet'),('Description', 'registerUser')]

MAP_USER_SAMPLE_SHEET_TO_DATABASE_ONE_INDEX = [('Sample_ID','userSampleID'), ('Sample_Name','sample_name'), ('Sample_Plate','samplePlate'),
            ('Sample_Well','sampleWell'),('Index_Plate_Well','indexPlateWell'),('I7_Index_ID','i7IndexID'),
            ('index','i7Index'),('Sample_Project','projectInSampleSheet'),('Description', 'registerUser')]
###### mapping including WELL



# mapping structure when sample sheet contains only one index
MAPPING_BASESPACE_SAMPLE_SHEET_ONE_INDEX = [('SampleID','Unique_Sample_ID'),('Name','Sample_Name'), ('Project','Sample_Project'),('Well', 'Sample_Well'),
                ('Index1Name','I7_Index_ID'), ('Index1Sequence','index' ) ]
# mapping structure when sample sheet contains two index
MAPPING_BASESPACE_SAMPLE_SHEET_TWO_INDEX = [('SampleID','Unique_Sample_ID'),('Name','Sample_Name'), ('Project','Sample_Project'),('Well', 'Sample_Well'),
                ('Index1Name','I7_Index_ID'), ('Index1Sequence','index' ),('Index2Name','I5_Index_ID'),('Index2Sequence','index2') ]

# Sections to check in the IEM file created by user
SECTIONS_IN_IEM_SAMPLE_SHEET = ['[Header]', '[Reads]', '[Settings]', '[Data]']

##############################################################

##### SETTINGS FOR CLEANUP RUNS AND NOT VALID FILES ##########

RETENTION_TIME = '7' # in days

##############################################################

##### HEADINGS VALUES

## Heading for pending Library Preparation state

HEADING_FOR_LIBRARY_PREPARATION_STATE = ['Sample extraction date', 'Sample', 'Molecule Code ID', 'Molecule Extraction Date',
                                    'Used Protocol', 'UserID']

#######HEADING_FOR_ADD_LIBRARY_PREPARATION = ['Molecule Code ID', 'Protocol', 'Extraction Date', 'To be included']
HEADING_FOR_ADD_LIBRARY_PREPARATION_PARAMETERS = ['Library Preparation Code ID', 'Sample Name', 'Protocol used', 'Add Library']
HEADING_FIX_FOR_ADDING_LIB_PARAMETERS = ['Sample Name', 'Library Preparation Code ID',  'Lot Regents Kit used']
HEADING_FOR_CREATION_LIBRARY_PREPARATION = ['Molecule Code ID', 'Protocol used', 'Single/Paired end', 'Length read']

### Heading for display information on library Preparation definition
HEADING_FOR_LIBRARY_PREPARATION_DEFINITION = ['Library CodeID','Molecule CodeID ','Lib Preparation State', 'Protocol name',
                    'Project Name', 'I7 Index', 'I5 Index', 'Single/PairedEnd', 'Read Length', 'Number of reused']

### Heading for display pool with the samples belongs to
HEADING_FOR_DISPLAY_SAMPLES_IN_POOL = ['Library CodeID', 'Sample Name', 'User Name','Collection Index', 'I7 Index', 'I5 Index', 'Include in Pool']

### Heading for display pool with the samples belongs to
HEADING_FOR_DISPLAY_CREATED_POOL = ['Pool Name' , 'Pool Code', 'Number of Samples in Pool']

HEADING_FOR_DISPLAY_LIB_PREP_IN_POOL = ['Library Code', 'Sample Name', 'User Name']

### Heading for display samples to be selected for Run
HEADING_FOR_SELECTING_POOLS = ['Pool Name', 'Pool Code', 'Number of samples', 'Include in the Run']

HEADING_FOR_INCOMPLETED_SELECTION_POOLS = ['Pool Name', 'Pool Code', 'Number of Samples']

### Heading for display pool information when showing sample information
HEADING_FOR_DISPLAY_POOL_INFORMATION_IN_SAMPLE_INFO = ['Library Code', 'Pool Name', 'Pool Code', 'Run Name']

### Heading for creating the pool for selected samples
HEADING_FOR_CREATING_RUN = ['Library CodeID', 'Sample Name', 'Pool Name', 'Sample Well','I7 Index', 'I7 Sequence', 'I5 Index', 'I5 Sequence', 'BaseSpace Library','Project Name', 'User Name']


### Heading for getting information when creating a new Run
HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_PAIREDEND = ['Unique_Sample_ID','Sample_Name','Sample_Plate','Sample_Well','Index_Plate_Well','I7_Index_ID','index','I5_Index_ID','index2','Sample_Project','Description', 'Base Space Library']

HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_SINGLEREAD = ['Unique_Sample_ID','Sample_Name','Sample_Plate','Sample_Well','Index_Plate_Well','I7_Index_ID','index','Sample_Project','Description','Base Space Library']


########## ERROR MESSAGES  #########################
ERROR_INVALID_FILE_FORMAT = ['Invalid file format for the selected file', 'Select the valid file and submit it again']
ERROR_UNABLE_TO_DELETE_USER_FILE = 'Unable to delete user file form iSkyLIMS'
ERROR_SAMPLE_SHEET_CONTAINS_NOT_DEFINED_SAMPLES = ['Sample sheet cannot be uploaded because there are samples', 'which are not defined yet.']
ERROR_SAMPLES_INVALID_STATE_FOR_LIBRARY_PREPARATION = ['Sample sheet cannot be uploaded because there are samples', ' in a state from which cannot accept data from library preparation.']
ERROR_SAMPLES_INVALID_DUPLICATED_INDEXES = ['Sample sheet cannot be uploaded because there are samples', 'which have duplicated index']
ERROR_COLLECTION_INDEX_KIT_NOT_DEFINED = ['Sample sheet cannot be uploaded because collection Index Kit', 'is not defined']

ERROR_EMPTY_VALUES = ['Your request cannot be recorded because ', 'it contains empty values']

ERROR_SAMPLE_PROJECT_ALREADY_EXISTS = ['Sample Project is already defined']
ERROR_SAMPLE_PROJECT_DOES_NOT_EXISTS = ['The Sample Project that you requested ', 'Does not exist']

ERROR_LIBRARY_PREPARATION_NOT_EXISTS = ['Library preparation are not defined', '']

ERROR_POOLS_WITH_NO_LIBRARY  = ['The selected Pools have not assigned to any Library Preparation.' ]

ERROR_NO_PROFILE_OR_CENTER_FOR_USER = ['Unable to save your request.', 'Update your Profile/ Center first']
