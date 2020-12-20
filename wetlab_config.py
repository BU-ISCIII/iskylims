import os
from django.conf import settings
##### Allow to import the configuration samba files from configuration folder
import sys
sys.path.append
try:
    from .wetlab_samba_conf import *
    SAMBA_USER_CONFIGURED = True
except:
    SAMBA_USER_CONFIGURED = False
try:
    from .wetlab_email_conf import *
    EMAIL_USER_CONFIGURED = True
except:
    EMAIL_USER_CONFIGURED = False

##############################################################
####
WETLAB_MANAGER = 'WetlabManager'
####
##############################################################

##### SETTINGS FOR CLEANUP RUNS AND NOT VALID FILES ##########
####
RETENTION_TIME = '7' # in days
####
##############################################################

##############  Define in proyect names can be the same in different Runs ######
PROJECTS_ALLOWED_IN_MULTIPLE_RUNS = 'False' ## ('True'/'False')
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


##### CONFIGURATION FOR ADDING KITS FOR LIBRARY PREPARATION ########
HEADING_ADDING_COMMERCIAL_KITS_TO_PROTOCOL = ['Given name ', 'Order', 'Used', 'Commercial Kit Name', 'Description']


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

######### MAPPING COLUMNS IN SAMPLE SHEET FROM USER TO DATABASE   #############
MAP_USER_SAMPLE_SHEET_TO_DATABASE_TWO_INDEX = [('Sample_ID','userSampleID'), ('Sample_Name','sample_name'), ('Sample_Plate','samplePlate'),
            ('Sample_Well','sampleWell'),('I7_Index_ID','i7IndexID'),
            ('index','i7Index'),('I5_Index_ID','i5IndexID'),('index2','i5Index'),('Sample_Project','projectInSampleSheet'),('Description', 'userInSampleSheet')]

MAP_USER_SAMPLE_SHEET_TO_DATABASE_ONE_INDEX = [('Sample_ID','userSampleID'), ('Sample_Name','sample_name'), ('Sample_Plate','samplePlate'),
            ('Sample_Well','sampleWell'),('I7_Index_ID','i7IndexID'),
            ('index','i7Index'),('Sample_Project','projectInSampleSheet'),('Description', 'userInSampleSheet')]

######### MAPPING OPTIONAL COLUMNS THAT COULD BE IN SAMPLE SHEET FROM USER TO DATABASE   #############
MAP_USER_SAMPLE_SHEET_ADDITIONAL_FIELDS_FROM_TYPE_OF_SECUENCER = [('Index_Plate_Well','indexPlateWell'), ('Manifest', 'manifest'), ('GenomeFolder', 'genomeFolder')]


MAP_USER_SAMPLE_SHEET_TO_DATABASE_TWO_INDEX_WITH_WELL = [('Sample_ID','userSampleID'), ('Sample_Name','sample_name'), ('Sample_Plate','samplePlate'),
            ('Sample_Well','sampleWell'),('Index_Plate_Well','indexPlateWell'),('I7_Index_ID','i7IndexID'),
            ('index','i7Index'),('I5_Index_ID','i5IndexID'),('index2','i5Index'),('Sample_Project','projectInSampleSheet'),('Description', 'registerUser')]

MAP_USER_NEXTSEQ_SAMPLE_SHEET_TO_DATABASE_TWO_INDEX_WITH_WELL = [('Sample_ID','userSampleID'), ('Sample_Name','sample_name'), ('Sample_Plate','samplePlate'),
            ('Sample_Well','sampleWell'),('Index_Plate_Well','indexPlateWell'),('I7_Index_ID','i7IndexID'),
            ('index','i7Index'),('Sample_Project','projectInSampleSheet'),('Description', 'registerUser')]


# mapping structure when sample sheet contains only one index
MAPPING_BASESPACE_SAMPLE_SHEET_ONE_INDEX = [('SampleID','Unique_Sample_ID'),('Name','Sample_Name'), ('Project','Sample_Project'),('Well', 'Sample_Well'),
                ('Index1Name','I7_Index_ID'), ('Index1Sequence','index' ) ]
# mapping structure when sample sheet contains two index
MAPPING_BASESPACE_SAMPLE_SHEET_TWO_INDEX = [('SampleID','Unique_Sample_ID'),('Name','Sample_Name'), ('Project','Sample_Project'),('Well', 'Sample_Well'),
                ('Index1Name','I7_Index_ID'), ('Index1Sequence','index' ),('Index2Name','I5_Index_ID'),('Index2Sequence','index2') ]

# Sections to check in the IEM file created by user
SECTIONS_IN_IEM_SAMPLE_SHEET = ['[Header]', '[Reads]', '[Settings]', '[Data]']



##### HEADINGS VALUES

## Heading for pending Library Preparation state

HEADING_FOR_SAMPLES_TO_DEFINE_PROTOCOL = ['Sample Name', 'Molecule Code ID', 'Library Preparation Protocol']

HEADING_FOR_LIBRARY_PREPARATION_STATE = ['Sample extraction date', 'Sample', 'Molecule Code ID', 'Molecule Extraction Date',
                                    'Used Protocol', 'UserID']

#######HEADING_FOR_ADD_LIBRARY_PREPARATION = ['Molecule Code ID', 'Protocol', 'Extraction Date', 'To be included']
HEADING_FOR_ADD_LIBRARY_PREPARATION_PARAMETERS = ['Library Preparation Code ID', 'Sample Name', 'Protocol used', 'Add Library']
HEADING_FIX_FOR_ADDING_LIB_PARAMETERS = ['Sample Name', 'Library Preparation Code ID']
HEADING_FIX_FOR_ADDING_LIB_PROT_PARAMETERS = ['Sample Name','Library Preparation Code ID']
HEADING_FIX_FOR_ASSING_ADDITIONAL_KITS = ['Sample Name','Library Preparation Code ID']
HEADING_FOR_CREATION_LIBRARY_PREPARATION = ['Molecule Code ID', 'Protocol used', 'Single/Paired end', 'Length read']
HEADING_MAIN_DATA_SAMPLE_SHEET = ['Application', 'Instrument Type', 'Assay', 'Index Adapters', 'Reads', 'Adapter', 'Adapter 2']
HEADING_SUMMARY_DATA_SAMPLE_SHEET = ['Number of Samples', 'Projects Name', 'Users']


### Heading for display information on library Preparation definition
HEADING_FOR_LIBRARY_PREPARATION_DEFINITION = ['Library CodeID','Molecule CodeID ','Lib Preparation State', 'Protocol name',
                    'Project Name', 'I7 Index', 'I5 Index', 'Number of reused']

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

### Heading for displying additional kits used in library prepation
HEADING_FOR_DISPLAY_ADDITIONAL_KIT_LIBRARY_PREPARATION = ['Library Preparation Code ID', 'Additional Lot kit name', 'Commercial kit name', 'Lot number', 'Recorded Date']

### Heading for creating the pool for selected samples
HEADING_FOR_CREATING_RUN = ['Library CodeID', 'Sample Name', 'Pool Name', 'Sample Well','I7 Index', 'I7 Sequence', 'I5 Index', 'I5 Sequence', 'BaseSpace Library','Project Name', 'User Name']


### Heading for getting information when creating a new Run
HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_PAIREDEND = ['Unique_Sample_ID','Sample_Name','Sample_Plate','Sample_Well','Index_Plate_Well','I7_Index_ID','index','I5_Index_ID','index2','Sample_Project','Description']

HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_SINGLEREAD = ['Unique_Sample_ID','Sample_Name','Sample_Plate','Sample_Well','Index_Plate_Well','I7_Index_ID','index','Sample_Project','Description']


HEADING_FOR_STATISTICS_RUNS_BASIC_DATA = ['Run Name', 'Date sequencer start']

############### SUCCESSFUL_MESSAGES #######################################
SUCCESSFUL_RUN_STATE_CHANGE_FOR_RETRY = ['State of the Run has changed back to the state previous to error', 'Now the run is again in the process for updating information',
                'Check the run from time to time this run to verify that run is moving forward.']
SUCCESSFUL_REUSE_MOLECULE_EXTRACTION = ['Molecule Extraction has been set for reused ', 'Now is ready to be create a new library preparation']
SUCCESSFUL_REUSE_LIB_PREP = ['Library Preparation has been set for reused ', 'Now is ready to be assigned to a new pool']

########## ERROR MESSAGES  #########################
ERROR_INVALID_FILE_FORMAT = ['Invalid file format for the selected file', 'Select the valid file and submit it again']
ERROR_UNABLE_TO_DELETE_USER_FILE = 'Unable to delete user file form iSkyLIMS'
ERROR_SAMPLE_SHEET_CONTAINS_NOT_DEFINED_SAMPLES = ['Sample sheet cannot be uploaded because there are samples', 'which are not defined yet.']
ERROR_SAMPLES_INVALID_STATE_FOR_LIBRARY_PREPARATION = ['Sample sheet cannot be uploaded because there are samples', ' are in a state from which cannot accept index information data.']
ERROR_SAMPLES_INVALID_DUPLICATED_INDEXES = ['Sample sheet cannot be uploaded because there are samples', 'which have duplicated index']


ERROR_COLLECTION_INDEX_KIT_NOT_DEFINED = ['Sample sheet cannot be uploaded because collection Index Kit', 'is not defined']
ERROR_NO_COLLECTION_INDEX_KIT_ARE_DEFINED = ['There is not Collection Kit defined', 'Please add the ones that will be used in your unit']

ERROR_SAMPLE_SHEET_DOES_NOT_HAVE_COLLECTION_INDEX = ['Sample Sheet does not have Index Adapters field']
ERROR_SAMPLE_SHEET_INSTRUMENT_TYPE_NOT_INCLUDED = ['Sample Sheet does not have Instrument type']
ERROR_SAMPLE_SHEET_BOTH_INSTRUMENT_AND_INDEX_NOT_INCLUDED = ['Sample Sheet does not have Instrument type neither Index Adapters']
ERROR_SAMPLE_SHEET_USERS_ARE_NOT_DEFINED = ['Users in sample sheet are not defined']
ERROR_SAMPLE_SHEET_USER_IS_NOT_DEFINED = ['User in sample sheet is not defined']
ERROR_SAMPLE_SHEET_DOES_NOT_HAVE_DESCRIPTION_FIELD = ['Sample sheet does not have Description column ']
ERROR_SAMPLE_SHEET_WHEN_FETCHING_USERID_NAMES = ['Sample sheet does not have the rigth format on the Description column']

ERROR_USER_SAMPLE_SHEET_NO_LONGER_EXISTS = ['The Sample Sheet that you are uploaded does not longer exists', 'Upload again the sample sheet']

ERROR_EMPTY_VALUES = ['Your request cannot be recorded because ', 'it contains empty values']

ERROR_SAMPLE_PROJECT_ALREADY_EXISTS = ['Sample Project is already defined']
ERROR_SAMPLE_PROJECT_DOES_NOT_EXISTS = ['The Sample Project that you requested ', 'Does not exist']

ERROR_LIBRARY_PREPARATION_NOT_EXISTS = ['Library preparation are not defined', '']
ERROR_NOT_LIBRARY_PREPARATION_SELECTED = ['No Library preparation was selected','Select at least one library preparation sample to create the pool ']

ERROR_POOLS_WITH_NO_LIBRARY  = ['The selected Pools have not assigned to any Library Preparation.' ]

ERROR_NO_PROFILE_OR_CENTER_FOR_USER = ['Unable to save your request.', 'Update your Profile/ Center first']

ERROR_RUN_IN_WRONG_STATE = ['Unable to accept your Run definition', 'Because your Run is in state']

ERROR_UNABLE_SAVE_REQUEST = ['Unable to save your request.']

ERROR_INVALID_PARAMETERS_WHEN_REUSING_LIB_PREP = ['Invalid data when trying to reuse library preparation', '']

ERROR_NO_SAMPLE_FOUND = ['No sample found which  match your  conditions ','']

ERROR_TOO_SHORT_INDEX_BASE_SEQUENCE = ['Index Sequence must contains at least 6  caracters' ,'' ]

ERROR_TOO_SHORT_INDEX_LIBRAY_NAME = ['Index Library Name  must contains at least 5  caracters' ,'' ]

ERROR_NO_COLLECTION_INDEX_FOUND = ['There are no recorded information for the Collection index kit' ]

ERROR_INVALID_SEQUENCE_CHARACTERS = ['Invalid characters in Index sequene', '']

###### ERROR TEXT FOR SEACHING #############################################
ERROR_NO_MATCHES_FOR_RUN_SEARCH = ['There is not any run that matches your input conditions']
ERROR_NO_MATCHES_FOR_PROJECT_SEARCH = ['There is not any project that mathes your input conditions']
ERROR_NO_MATCHES_FOR_LIBRARY_STATISTICS = ['There is not any Index Library Kit  that mathes your input conditions']
ERROR_NO_MATCHES_FOR_USER_LOT_KIT = ['There is not any User Lot Kit  that matches your input conditions']
ERROR_NO_MATCHES_FOR_INPUT_CONDITIONS =['There is not any match for your input conditions ']
ERROR_NO_USER_LOT_KIT_DEFINED =['No User Lot Kit are defined']
ERROR_INVALID_FORMAT_FOR_DATES = ['Invalid date format. Use the format  (DD-MM-YYYY)']
ERROR_USER_NAME_TOO_SHORT = ['User name must contains at least 5 characters']
ERROR_NO_MATCHES_FOR_SEQUENCER_STATS = ['There is not any run that where using the sequencer']
ERROR_MANY_USER_MATCHES_FOR_INPUT_CONDITIONS =['There are many user names that matches your request']
ERROR_WRONG_SAMBA_CONFIGURATION_SETTINGS = ['Unsuccessful configuration settings for Samba connection']
ERROR_UNABLE_TO_SAVE_SAMBA_CONFIGURATION_SETTINGS = ['Unable to save the Samba configuration file ', 'check if folder iSkyLIMS_wetlab has write permision for apache user']
ERROR_UNABLE_TO_SAVE_EMAIL_CONFIGURATION_SETTINGS = ['Unable to save the email configuration file ', 'check if folder iSkyLIMS_wetlab has write permision for apache user']

############### HEADING FOR PROJECT DATA VISUALIZATION #####################
HEADING_FOR_PROJECT_DATES = ['Project Recorder date', 'Project date']




############### FIELD NAME TO COLLECT FROM RunParameter FILE #####################
FIELDS_TO_COLLECT_FROM_RUN_INFO_FILE = ['RunID','ExperimentName','RTAVersion','Chemistry','RunStartDate','RunManagementType','ApplicationVersion','NumTilesPerSwath',
                    'SystemSuiteVersion', 'LibraryID', 'AnalysisWorkflowType','PlannedRead1Cycles','PlannedRead2Cycles','PlannedIndex1ReadCycles','PlannedIndex2ReadCycles' ]

################ EMAIL CONFIGURATION FIELDS ###############################
EMAIL_CONFIGURATION_FIELDS = ['USER_NAME', 'USER_EMAIL', 'SENT_EMAIL_ON_ERROR']
EMAIL_CONFIGURATION_FILE_HEADING = '############# EMAIL CONFIGURATION FILE ########\n#DO NOT MODIFY MANUALLY THIS FILE\n#VALUES WILL BE MODIFIED WHEN USING THE CONFIGURATION FORM\n'
EMAIL_CONFIGURATION_FILE_END = '########## END EMAIL CONFIGURATION FILE'


################ SAMBA CONFIGURATION FIELDS ###############################
SAMBA_CONFIGURATION_FIELDS = ['SAMBA_USER_ID', 'SAMBA_USER_PASSWORD', 'SAMBA_SHARED_FOLDER_NAME', 'SAMBA_APPLICATION_FOLDER_NAME', 'SAMBA_REMOTE_SERVER_NAME',
                'SAMBA_NTLM_USED', 'SAMBA_IP_SERVER', 'SAMBA_HOST_NAME', 'SAMBA_PORT_SERVER', 'IS_DIRECT_TCP', 'SAMBA_DOMAIN']

SAMBA_CONFIGURATION_FILE_HEADING = '############# SAMBA CONFIGURATION FILE ########\n#DO NOT MODIFY MANUALLY THIS FILE\n#VALUES WILL BE MODIFIED WHEN USING THE CONFIGURATION FORM\n'
SAMBA_CONFIGURATION_FILE_END = '########## END SAMBA CONFIGURATION FILE'

######### PROJECT HEADING  #########################################
HEADING_SINGLE_PROJECT_FL_SUMMARY = ['Cluster (Raw)', 'Cluster (PF)', 'Yield (MBases)', 'Number of Samples']

HEADING_SINGLE_PROJECT_STATS_LANE = ['Lane', 'PF Clusters', '% of the lane','% Perfect barcode', '% One mismatch barcode','Yield (Mbases)','% >= Q30 bases', 'Mean Quality Score']

HEADING_SINGLE_PROJECT_SAMPLES = ['Sample','Barcode','PF Clusters','Percent of Project', 'Yield (Mbases)','% >= Q30 bases', 'Mean Quality Score']
