# Allow to import the configuration samba files from configuration folder
import sys

sys.path.append

##############################################################
####
WETLAB_MANAGER = "WetlabManager"
####
##############################################################

########################################################
# SETTINGS FOR CLEANUP RUNS AND NOT VALID FILES ##########

RETENTION_TIME = "7"  # in days

##############################################################

# #############  Define in proyect names can be the same in different Runs ######
# PROJECTS_ALLOWED_IN_MULTIPLE_RUNS = 'False' ## ('True'/'False')
##############################################################

# ############# FOLDER SETTINGS ###############################
# # Directory settings for processing the run data files ######
# Relative path from settings.BASE_DIR
LOG_DIRECTORY = "logs/"
# Relative path from settings.MEDIA_ROOT
RUN_TEMP_DIRECTORY_RECORDED = "wetlab/tmp/recorded/"
RUN_TEMP_DIRECTORY = "wetlab/tmp"
RUN_TEMP_DIRECTORY_PROCESSING = "wetlab/tmp/processing"
RUN_IMAGES_DIRECTORY = "wetlab/images_plot"
RUN_SAMPLE_SHEET_DIRECTORY = "wetlab/SampleSheets/"
# SAMPLE_SHEET_CREATED_ON_LAB = 'wetlab/SampleSheetsFromLab'
TEMPLATE_FILES_DIRECTORY = "wetlab/templates"

# Directory to store the imported user sampleSheets
LIBRARY_PREPARATION_SAMPLE_SHEET_DIRECTORY = "wetlab/SampleSheets4LibPrep/"

##############################################################

##############################################
#                 LOG NAMES
#############################
LOG_NAME_FIND_UPDATE_RUNS = "update_run"

###############################
#              CONFIG FILE LOG NAME
###############################
LOGGING_CONFIG_FILE = "logging_config.ini"

##############################################################

#  FILE NAME CONTAINING THE PROCESSED RUNS ###################
PROCESSED_RUN_FILE = "processed_run_file"
##############################################################

# ############ ILLUMINA OUTPUT FILES ##########################
RUN_PARAMETER_FILE = "RunParameters.xml"
# RUN_PARAMETER_MISEQ = 'runParameters.xml'
RUN_INFO = "RunInfo.xml"
RUN_COMPLETION_XML_FILE = "RunCompletionStatus.xml"
RUN_COMPLETION_TXT_FILE = "SequenceComplete.txt"
SAMPLE_SHEET = "samplesheet.csv"
# # sample sheet to be copied on the remote folder
COPY_SAMPLE_SHEET_TO_REMOTE = True  # boolean constant True if NextSeq
# sample sheet needs to be copied to remote server
RUN_LOG_FOLDER = "Logs"

STATS_FILE_PATH = "Data/Intensities/BaseCalls/Stats"


CONVERSION_STATS_FILE = "ConversionStats.xml"
DEMULTIPLEXION_STATS_FILE = "DemultiplexingStats.xml"
##############################################################

PLATFORM_WAY_TO_CHECK_RUN_COMPLETION = [
    ["NextSeq", "xml_file"],
    ["MiSeq", "xml_file"],
    ["NovaSeq", "txt_file"],
]

# ########### VALUE TAG FOR XML FILES #########################
COMPLETION_TAG = "CompletionStatus"
COMPLETION_SUCCESS = ["CompletedAsPlanned", "SuccessfullyCompleted"]
EXPERIMENT_NAME_TAG = "ExperimentName"
APPLICATION_NAME_TAG = "ApplicationName"
NUMBER_CYCLES_TAG = "NumCycles"
RUN_INFO_READ_TAG = "RunInfoRead"
NUMBER_TAG = "Number"

##############################################################
RUN_METRIC_GRAPHIC_COMMANDS = [
    "plot_by_cycle  ",
    "plot_by_lane  ",
    "plot_flowcell  ",
    "plot_qscore_histogram  ",
    "plot_qscore_heatmap  ",
    "plot_sample_qc  ",
]
##############################################################

######################################
# FIELD TAG NAMES TO COLLECT FROM RunParameter FILE
#####################################
FIELDS_TO_COLLECT_FROM_RUN_INFO_FILE = [
    "RunID",
    "ExperimentName",
    "RTAVersion",
    "Chemistry",
    "RunStartDate",
    "RunManagementType",
    "SystemSuiteVersion",
    "LibraryID",
    "AnalysisWorkflowType",
    "PlannedRead1Cycles",
    "PlannedRead2Cycles",
    "PlannedIndex1ReadCycles",
    "PlannedIndex2ReadCycles",
]
SETUP_TAG = "Setup"
FIELDS_TO_FETCH_FROM_SETUP_TAG = [
    "NumLanes",
    "ApplicationName",
    "ApplicationVersion",
    "NumTilesPerSwath",
]
READ_NUMBER_OF_CYCLES = [
    "PlannedRead1Cycles",
    "PlannedIndex1ReadCycles",
    "PlannedIndex2ReadCycles",
    "PlannedRead2Cycles",
]
# NOVASEQ 6000
FIELDS_NOVASEQ_TO_FETCH_TAG = [
    "NumLanes",
    "Application",
    "ApplicationVersion",
    "NumTilesPerSwath",
]


##############################################################

# ########### DEFAULT VALUES FOR MISEQ SAMPLE SHEET  ##########
DEFAULT_LIBRARY_KIT = "Unknown"
DEFAULT_CENTER = "CNM"
##############################################################

# ########### RUN METRIC FOLDERS AND FILES ####################
INTEROP_PATH = "/opt/interop/bin/"
RUN_METRIC_FOLDER = "InterOp"
PLOT_EXTENSION = ".png"

##############################################################

# ########  FOLDER SETTINGS FOR PROCESSING LIBRARY KITS #######
# # Relative path from settings.BASE_DIR
LIBRARY_KITS_DIRECTORY = "wetlab/collection_index_kits/"
COLLECTION_INDEX_KITS_DIRECTORY = "wetlab/collection_index_kits/"
# # Maximum file size allowed for the index library kits (in bytes)
LIBRARY_MAXIMUM_SIZE = "3145728"
# # Configuration for index library file
COLLECTION_INDEX_HEADING = [
    "[Version]",
    "[Name]",
    "[PlateExtension]",
    "[Settings]",
    "[I7]",
]
##############################################################

MIGRATION_DIRECTORY_FILES = "wetlab/BaseSpaceMigrationFiles/"

##############################################################

SAMPLE_SHEET_MISEQ_VERSION_4_TEMPLATE_NAME = (
    "miseq_version_4_sample_sheet_one_adapters_template.csv"
)
SAMPLE_SHEET_MISEQ_VERSION_5_TEMPLATE_NAME = (
    "miseq_version_5_sample_sheet_one_adapters_template.csv"
)
SAMPLE_SHEET_NEXTSEQ_VERSION_5_TEMPLATE_NAME = (
    "nextseq_version_5_sample_sheet_template.csv"
)

# #### CONFIGURATION FOR ADDING KITS FOR LIBRARY PREPARATION ########
HEADING_ADDING_COMMERCIAL_KITS_TO_PROTOCOL = [
    "Given name ",
    "Order",
    "Used",
    "Commercial Kit Name",
    "Description",
]
HEADING_FOR_MODIFYING_ADDITIONAL_KITS = [
    "Old name",
    "New name",
    "Order",
    "Used",
    "Commercial Kit Name",
    "Description",
]

HEADING_FOR_PENDING_PROCESS_SAMPLES = [
    "Sample name",
    "Laboratory",
    "Sample type",
    "Sample project",
    "User name",
]

# CONFIGURATION FOR SAMPLE SHEET CONVERSION #
# FILE TO BASESPACE FORMAT     #
# column names when sample sheet has only one index

HEADING_FOR_SAMPLE_SHEET_ONE_INDEX = [
    "Unique_Sample_ID",
    "Sample_Name",
    "Sample_Plate",
    "Sample_Well",
    "Index_Plate_Well",
    "I7_Index_ID",
    "index",
    "Sample_Project",
    "Description",
]

HEADING_FOR_SAMPLE_SHEET_TWO_INDEX = [
    "Unique_Sample_ID",
    "Sample_Name",
    "Sample_Plate",
    "Sample_Well",
    "Index_Plate_Well",
    "I7_Index_ID",
    "index",
    "I5_Index_ID",
    "index2",
    "Sample_Project",
    "Description",
]


# ######## MAPPING COLUMNS IN SAMPLE SHEET FROM USER TO DATABASE   #############
MAP_USER_SAMPLE_SHEET_TO_DATABASE_NEXTSEQ_SINGLE_READ = [
    ("Sample_ID", "userSampleID"),
    ("Sample_Name", "sample_name"),
    ("Sample_Plate", "samplePlate"),
    ("Sample_Well", "sampleWell"),
    ("I7_Index_ID", "i7IndexID"),
    ("index", "i7Index"),
    ("Sample_Project", "projectInSampleSheet"),
    ("Description", "userInSampleSheet"),
]

MAP_USER_SAMPLE_SHEET_TO_DATABASE_NEXTSEQ_PAIRED_END = [
    ("Sample_ID", "userSampleID"),
    ("Sample_Name", "sample_name"),
    ("Sample_Plate", "samplePlate"),
    ("Sample_Well", "sampleWell"),
    ("I7_Index_ID", "i7IndexID"),
    ("index", "i7Index"),
    ("I5_Index_ID", "i5IndexID"),
    ("index2", "i5Index"),
    ("Sample_Project", "projectInSampleSheet"),
    ("Description", "userInSampleSheet"),
]

MAP_USER_SAMPLE_SHEET_TO_DATABASE_MISEQ_SINGLE_READ_VERSION_5 = [
    ("Sample_ID", "userSampleID"),
    ("Sample_Name", "sample_name"),
    ("Sample_Plate", "samplePlate"),
    ("Sample_Well", "sampleWell"),
    ("I7_Index_ID", "i7IndexID"),
    ("index", "i7Index"),
    ("Sample_Project", "projectInSampleSheet"),
    ("Description", "userInSampleSheet"),
]

MAP_USER_SAMPLE_SHEET_TO_DATABASE_MISEQ_PAiRED_END_VERSION_5 = [
    ("Sample_ID", "userSampleID"),
    ("Sample_Name", "sample_name"),
    ("Sample_Plate", "samplePlate"),
    ("Sample_Well", "sampleWell"),
    ("I7_Index_ID", "i7IndexID"),
    ("index", "i7Index"),
    ("I5_Index_ID", "i5IndexID"),
    ("index2", "i5Index"),
    ("Sample_Project", "projectInSampleSheet"),
    ("Description", "userInSampleSheet"),
]

MAP_USER_SAMPLE_SHEET_TO_DATABASE_MISEQ_SINGLE_READ_VERSION_4 = [
    ("Sample_ID", "userSampleID"),
    ("Sample_Name", "sample_name"),
    ("Sample_Plate", "samplePlate"),
    ("Sample_Well", "sampleWell"),
    ("I7_Index_ID", "i7IndexID"),
    ("index", "i7Index"),
    ("Sample_Project", "projectInSampleSheet"),
    ("Description", "userInSampleSheet"),
]

MAP_USER_SAMPLE_SHEET_TO_DATABASE_MISEQ_PAiRED_END_VERSION_4 = [
    ("Sample_ID", "userSampleID"),
    ("Sample_Name", "sample_name"),
    ("Sample_Plate", "samplePlate"),
    ("Sample_Well", "sampleWell"),
    ("I7_Index_ID", "i7IndexID"),
    ("index", "i7Index"),
    ("I5_Index_ID", "i5IndexID"),
    ("index2", "i5Index"),
    ("Sample_Project", "projectInSampleSheet"),
    ("Description", "userInSampleSheet"),
]


MAP_USER_SAMPLE_SHEET_TO_DATABASE_ALL_PLATFORMS = [
    ("Sample_ID", "userSampleID"),
    ("Sample_Name", "sample_name"),
    ("Sample_Plate", "samplePlate"),
    ("Sample_Well", "sampleWell"),
    ("I7_Index_ID", "i7IndexID"),
    ("index", "i7Index"),
    ("I5_Index_ID", "i5IndexID"),
    ("index2", "i5Index"),
    ("Index_Plate_Well", "indexPlateWell"),
    ("Manifest", "manifest"),
    ("GenomeFolder", "genomeFolder"),
    ("Sample_Project", "projectInSampleSheet"),
    ("Description", "userInSampleSheet"),
]
# ######## MAPPING OPTIONAL COLUMNS THAT COULD BE IN SAMPLE SHEET FROM USER TO DATABASE   #############
MAP_USER_SAMPLE_SHEET_ADDITIONAL_FIELDS_FROM_TYPE_OF_SECUENCER = [
    ("Index_Plate_Well", "indexPlateWell"),
    ("Manifest", "manifest"),
    ("GenomeFolder", "genomeFolder"),
]

# Sections to check in the IEM file created by user
SECTIONS_IN_IEM_SAMPLE_SHEET = ["[Header]", "[Reads]", "[Settings]", "[Data]"]

FIELDS_IN_SAMPLE_SHEET_HEADER_IEM_VERSION_5 = [
    "Date",
    "Experiment Name",
    "Workflow",
    "Application",
    "Instrument Type",
    "Assay",
    "Index Adapters",
    "Chemistry",
    "Description",
]

# #### HEADINGS VALUES

# # Heading for pending Library Preparation state

HEADING_FOR_SAMPLES_TO_DEFINE_PROTOCOL = [
    "Sample Name",
    "Molecule Code ID",
    "Library Preparation Protocol",
]

HEADING_FOR_LIBRARY_PREPARATION_STATE = [
    "Sample extraction date",
    "Sample",
    "Molecule Code ID",
    "Molecule Extraction Date",
    "Used Protocol",
    "UserID",
]

# ######HEADING_FOR_ADD_LIBRARY_PREPARATION = ['Molecule Code ID', 'Protocol', 'Extraction Date', 'To be included']
HEADING_FOR_ADD_LIBRARY_PREPARATION_PARAMETERS = [
    "Library Preparation Code ID",
    "Sample Name",
    "Protocol used",
    "Add Library",
]
HEADING_FIX_FOR_ADDING_LIB_PARAMETERS = ["Sample Name", "Library Preparation Code ID"]
HEADING_FIX_FOR_ADDING_LIB_PROT_PARAMETERS = [
    "Sample Name",
    "Library Preparation Code ID",
]
HEADING_FIX_FOR_ASSING_ADDITIONAL_KITS = ["Sample Name", "Library Preparation Code ID"]
HEADING_FOR_CREATION_LIBRARY_PREPARATION = [
    "Molecule Code ID",
    "Protocol used",
    "Single/Paired end",
    "Length read",
]
HEADING_MAIN_DATA_SAMPLE_SHEET = [
    "Application",
    "Instrument Type",
    "Assay",
    "Index Adapters",
    "Reads",
    "Adapter",
    "Adapter 2",
]
HEADING_SUMMARY_DATA_SAMPLE_SHEET = ["Number of Samples", "Projects Name", "Users"]


# ## Heading for display information on library Preparation definition
HEADING_FOR_LIBRARY_PREPARATION_DEFINITION = [
    "Library CodeID",
    "Molecule CodeID ",
    "Lib Preparation State",
    "Protocol name",
    "Project Name",
    "I7 Index",
    "I5 Index",
    "Number of reused",
]

# ## Heading for display pool with the samples belongs to
HEADING_FOR_DISPLAY_SAMPLES_IN_POOL = [
    "Library CodeID",
    "Sample Name",
    "User Name",
    "Collection Index",
    "I7 Index",
    "I5 Index",
    "Include in Pool",
]

# ## Heading for display pool with the samples belongs to
HEADING_FOR_DISPLAY_CREATED_POOL = [
    "Pool Name",
    "Pool Code",
    "Number of Samples in Pool",
]

HEADING_FOR_DISPLAY_LIB_PREP_IN_POOL = ["Library Code", "Sample Name", "User Name"]

# ## Heading for display samples to be selected for Run
HEADING_FOR_SELECTING_POOLS = [
    "Pool Name",
    "Pool Code",
    "Number of samples",
    "Include in the Run",
]

HEADING_FOR_INCOMPLETED_SELECTION_POOLS = [
    "Pool Name",
    "Pool Code",
    "Number of Samples",
]

# ## Heading for display pool information when showing sample information
HEADING_FOR_DISPLAY_POOL_INFORMATION_IN_SAMPLE_INFO = [
    "Library Code",
    "Pool Name",
    "Pool Code",
    "Run Name",
]

# ## Heading for displying additional kits used in library prepation
HEADING_FOR_DISPLAY_ADDITIONAL_KIT_LIBRARY_PREPARATION = [
    "Library Preparation Code ID",
    "Additional Lot kit name",
    "Commercial kit name",
    "Lot number",
    "Recorded Date",
]

# ## Heading for displying  kits used in run preparation
HEADING_FOR_DISPLAY_KIT_IN_RUN_PREPARATION = [
    "Library Preparation Code ID",
    "Lot NUMBER",
    "Commercial kit name",
    "Lot exp. date",
    "Run recorded Date",
]

# ## Heading for creating the pool for selected samples
HEADING_FOR_CREATING_RUN = [
    "Library CodeID",
    "Sample Name",
    "Pool Name",
    "Sample Well",
    "I7 Index",
    "I7 Sequence",
    "I5 Index",
    "I5 Sequence",
    "BaseSpace Library",
    "Project Name",
    "User Name",
]


# ## Heading for getting information when creating a new Run
HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_NEXTSEQ_PAIRED_END = [
    "Sample_ID",
    "Sample_Name",
    "Sample_Plate",
    "Sample_Well",
    "Index_Plate_Well",
    "I7_Index_ID",
    "index",
    "I5_Index_ID",
    "index2",
    "Sample_Project",
    "Description",
]
HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_NEXTSEQ_SINGLE_READ = [
    "Sample_ID",
    "Sample_Name",
    "Sample_Plate",
    "Sample_Well",
    "Index_Plate_Well",
    "I7_Index_ID",
    "index",
    "Sample_Project",
    "Description",
]

HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_MISEQ_PAIRED_END_VERSION_4 = [
    "Sample_ID",
    "Sample_Name",
    "Sample_Plate",
    "Sample_Well",
    "I7_Index_ID",
    "index",
    "I5_Index_ID",
    "index2",
    "GenomeFolder",
    "Sample_Project",
    "Description",
]
HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_MISEQ_SINGLE_READ_VERSION_4 = [
    "Sample_ID",
    "Sample_Name",
    "Sample_Plate",
    "Sample_Well",
    "I7_Index_ID",
    "index",
    "GenomeFolder",
    "Sample_Project",
    "Description",
]

HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_MISEQ_PAIRED_END_VERSION_5 = [
    "Sample_ID",
    "Sample_Name",
    "Sample_Plate",
    "Sample_Well",
    "I7_Index_ID",
    "index",
    "I5_Index_ID",
    "index2",
    "Manifest",
    "GenomeFolder",
    "Sample_Project",
    "Description",
]
HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_MISEQ_SINGLE_READ_VERSION_5 = [
    "Sample_ID",
    "Sample_Name",
    "Sample_Plate",
    "Sample_Well",
    "I7_Index_ID",
    "index",
    "Manifest",
    "GenomeFolder",
    "Sample_Project",
    "Description",
]


HEADING_FOR_STATISTICS_RUNS_BASIC_DATA = ["Run Name", "Date sequencer start"]

HEADING_STATISTICS_FOR_RESEARCHER_SAMPLE = [
    "Samples",
    "Project name",
    "Run name",
    "Platform",
]
HEADING_STATISTICS_FOR_TIME_RUN = ["Run name", "Run state", "Sequencer", "Run date"]
HEADING_STATISTICS_FOR_TIME_SAMPLE = [
    "Sample name",
    "Researcher",
    "Project name",
    "Run name",
    "Barcode",
]
HEADING_STATISTICS_FOR_SEQUENCER_RUNS = [
    "Run name",
    "Run state",
    "Sequencer",
    "Run date",
]
NUMBER_OF_VALUES_TO_FETCH_FROM_RESEARCHER = 5

# ############## SUCCESSFUL_MESSAGES #######################################
SUCCESSFUL_RUN_STATE_CHANGE_FOR_RETRY = [
    "State of the Run has changed back to the state previous to error",
    "Now the run is again in the process for updating information",
    "Check the run from time to time this run to verify that run is moving forward.",
]
SUCCESSFUL_REUSE_MOLECULE_EXTRACTION = [
    "Molecule Extraction has been set for reused ",
    "Now is ready to be create a new library preparation",
]
SUCCESSFUL_REUSE_LIB_PREP = [
    "Library Preparation has been set for reused ",
    "Now is ready to be assigned to a new pool",
]

# ######### ERROR MESSAGES  #########################

ERROR_USER_NOT_WETLAB_MANAGER = "You do not have enough privileges to see this page "
ERROR_INVALID_FILE_FORMAT = "Invalid file format for the selected file"
ERROR_UNABLE_TO_DELETE_USER_FILE = "Unable to delete user file form iSkyLIMS"
ERROR_SAMPLE_SHEET_CONTAINS_NOT_DEFINED_SAMPLES = (
    "Sample sheet with not defined samples"
)
ERROR_SAMPLE_SHEET_CONTAINS_NOT_DEFINED_SAMPLES_WITH_PROTOCOL = (
    "Sample sheet with not prefix protocol in samples"
)
ERROR_SAMPLES_INVALID_STATE_FOR_LIBRARY_PREPARATION = (
    "Sample sheet with samples which are not in library preparation state."
)
ERROR_SAMPLES_INVALID_DUPLICATED_INDEXES = (
    "Sample sheet have samples with duplicated index"
)


ERROR_COLLECTION_INDEX_KIT_NOT_DEFINED = [
    "Sample sheet cannot be uploaded because collection Index Kit",
    "is not defined",
]
ERROR_NO_COLLECTION_INDEX_KIT_ARE_DEFINED = [
    "There is not Collection Kit defined",
    "Please add the ones that will be used in your unit",
]

ERROR_SAMPLE_SHEET_DOES_NOT_HAVE_COLLECTION_INDEX = [
    "Sample Sheet does not have Index Adapters field"
]
ERROR_SAMPLE_SHEET_INSTRUMENT_TYPE_NOT_INCLUDED = [
    "Sample Sheet does not have Instrument type"
]
ERROR_SAMPLE_SHEET_BOTH_INSTRUMENT_AND_INDEX_NOT_INCLUDED = [
    "Sample Sheet does not have Instrument type neither Index Adapters"
]
# ERROR_SAMPLE_SHEET_USERS_ARE_NOT_DEFINED = ["Users in sample sheet are not defined"]
ERROR_SAMPLE_SHEET_USER_IS_NOT_DEFINED = ["User in sample sheet is not defined"]
ERROR_SAMPLE_SHEET_DOES_NOT_HAVE_DESCRIPTION_FIELD = [
    "Sample sheet does not have Description column "
]
ERROR_SAMPLE_SHEET_WHEN_FETCHING_USERID_NAMES = [
    "Sample sheet does not have the rigth format on the Description column"
]
ERROR_SAMPLE_SHEET_USER_ARE_NOT_DEFINED = [
    "Sample sheet has users which are not defined"
]
ERROR_USER_SAMPLE_SHEET_NO_LONGER_EXISTS = [
    "The Sample Sheet that you are uploaded does not longer exists",
    "Upload again the sample sheet",
]


ERROR_EMPTY_VALUES = [
    "Your request cannot be recorded because ",
    "it contains empty values",
]

ERROR_SAMPLE_PROJECT_ALREADY_EXISTS = ["Sample Project is already defined"]
ERROR_SAMPLE_PROJECT_DOES_NOT_EXISTS = "Sample Project requested does not exist"

ERROR_LIBRARY_PREPARATION_NOT_EXISTS = ["Library preparation are not defined", ""]
ERROR_NOT_LIBRARY_PREPARATION_SELECTED = [
    "No Library preparation was selected",
    "Select at least one library preparation sample to create the pool ",
]

ERROR_POOLS_WITH_NO_LIBRARY = [
    "The selected Pools have not assigned to any Library Preparation."
]

ERROR_NO_PROFILE_OR_CENTER_FOR_USER = [
    "Unable to save your request.",
    "Update your Profile/ Center first",
]

ERROR_RUN_IN_WRONG_STATE = [
    "Unable to accept your Run definition",
    "Because the Run is in state",
]
ERROR_RUN_DOES_NOT_EXIST = "Requested run does not exist"
ERROR_UNABLE_SAVE_REQUEST = ["Unable to save your request."]

ERROR_INVALID_PARAMETERS_WHEN_REUSING_LIB_PREP = [
    "Invalid data when trying to reuse library preparation",
    "",
]

ERROR_NO_SAMPLE_FOUND = "No sample found which  match your  conditions "

ERROR_TOO_SHORT_INDEX_BASE_SEQUENCE = [
    "Index Sequence must contains at least 6  caracters",
    "",
]

ERROR_TOO_SHORT_INDEX_LIBRAY_NAME = [
    "Index Library Name  must contains at least 5  caracters",
    "",
]

ERROR_INVALID_SEQUENCE_CHARACTERS = ["Invalid characters in Index sequene", ""]

ERROR_NO_POOL_WAS_SELECTED_IN_FORM = [
    "There was not selected any pool to create the new Run"
]

ERROR_RUN_NAME_ALREADY_DEFINED = [
    "Run name ",
    "is already defined",
    "Write other Run Name on the field bellow",
]

ERROR_RUN_NAME_CREATED_ALREADY = " Run already created"

ERROR_DIFFERENT_ADAPTERS_USED_IN_POOL = [
    "Different Adapters ",
    "Were used in the pools",
]

ERROR_DUPLICATED_INDEXES_FOUND_IN_DIFFERENT_POOLS = [
    "Found duplicated index",
    "in the pools",
]

# ERROR TEXT FOR SEACHING #############################################
ERROR_INVALID_FORMAT_FOR_DATES = "Invalid date format. Use the format  (DD-MM-YYYY)"


ERROR_MANY_USER_MATCHES_FOR_INPUT_CONDITIONS = [
    "There are many user names that matches your request"
]
ERROR_NO_MATCHES_FOR_INPUT_CONDITIONS = (
    "There is not any match for your input conditions"
)
ERROR_NO_MATCHES_FOR_LIBRARY_STATISTICS = (
    "There is not any Index Library Kit that mathes your input conditions"
)
ERROR_NO_MATCHES_FOR_PROJECT_SEARCH = (
    "There is not any project that mathes your input conditions"
)
ERROR_NO_MATCHES_FOR_RUN_SEARCH = (
    "There are not matches that fits your input conditions"
)
ERROR_NO_MATCHES_FOR_SEQUENCER_STATS = [
    "There is not any run that where using the sequencer"
]
ERROR_NO_MATCHES_FOR_USER_LOT_KIT = [
    "There is not any User Lot Kit  that matches your input conditions"
]
ERROR_NO_SAMPLES_SELECTED = "They were not selected any Sample on your request"

ERROR_NOT_RUNS_FOUND_IN_SELECTED_PERIOD = (
    "There are not runs for the selected period of time"
)
ERROR_NOT_SAMPLES_FOR_USER_FOUND_BECAUSE_OF_CONFIGURATION_SETTINGS = [
    "No results. This could because the DESCRIPTION_IN_SAMPLE_SHEET_MUST_HAVE_USERNAME setting is set fo FALSE"
]
ERROR_NO_USER_LOT_KIT_DEFINED = ["No User Lot Kit are defined"]
ERROR_UNABLE_TO_SAVE_SAMBA_CONFIGURATION_SETTINGS = [
    "Unable to save the Samba configuration file ",
    "check if folder wetlab has write permision for apache user",
]
ERROR_UNABLE_TO_SAVE_EMAIL_CONFIGURATION_SETTINGS = [
    "Unable to save the email configuration file ",
    "check if folder wetlab has write permision for apache user",
]
ERROR_USER_DOES_NOT_HAVE_ANY_SAMPLE = [
    "User does not have any sample to perform statistics"
]
ERROR_USER_NAME_TOO_SHORT = "User name must contains at least 5 characters"
ERROR_WRONG_SAMBA_CONFIGURATION_SETTINGS = (
    "Unsuccessful configuration settings for Samba connection"
)
ERROR_WRONG_SAMBA_FOLDER_SETTINGS = (
    "Unsuccessful configuration. Samba folder not reachable"
)
ERROR_WRONG_SAMBA_AUTHENTICATION_SETTINGS = (
    "Unsuccessful configuration. Invalid credential"
)

# ########################  Sequencer errors #####################################
ERROR_SEQUENCER_ALREADY_DEFINED = [
    "Unable to save the Sequencer, because it already exists"
]

ERROR_SEQUENCER_CONFIGURATION_ALREADY_DEFINED = [
    "Unable to save the Sequencer configuration, because it already exists"
]

ERROR_NOT_ALLOWED_REPEATED_PROJECTS = (
    "Configuration settigs are set that same project can not be in different runs"
)

# ########################  Configuration test errors #####################################
ERROR_NOT_FOLDER_RUN_TEST_WAS_FOUND = [
    "Unable to run the configuration test",
    "Run test folder was found on remote server",
]
ERROR_NO_RUN_TEST_WAS_CREATED = [
    "Unable to continue with configuration testing",
    "Run test was not created during the execution",
]

# ############### Batch file errors ##########################################
ERROR_MESSAGE_UPLOAD_FILE_NOT_EXISTS = ["Upload file does not exist"]

# ######################## API ERRORS ########################################
ERROR_API_NO_SAMPLE_DEFINED = "There are no samples defined yet"
ERROR_API_SAMPLE_STATE_VALUE_IS_NOT_DEFINED = "The sample state is not defined"
ERROR_API_NO_SAMPLE_PROJECT_DEFINED = "The sample project is not defined"
ERROR_API_NO_SAMPLE_PROJECT_FIELD_DEFINED = "The sample project field is not defined"

# ############## HEADING FOR PROJECT DATA VISUALIZATION #####################
HEADING_FOR_PROJECT_DATES = ["Project Recorder date", "Project date"]


# ############### EMAIL CONFIGURATION FIELDS ###############################
EMAIL_CONFIGURATION_FIELDS = [
    "EMAIL_HOST",
    "EMAIL_PORT",
    "USER_PASSWORD",
    "USER_NAME",
    "USER_EMAIL",
    "SENT_EMAIL_ON_ERROR",
    "USE_TLS",
]


# ############### SAMBA CONFIGURATION FIELDS ###############################
SAMBA_CONFIGURATION_FIELDS = [
    "user_id",
    "user_password",
    "shared_folder_name",
    "samba_folder_name",
    "remote_server_name",
    "ntlm_used",
    "ip_server",
    "host_name",
    "port_server",
    "is_direct_tcp",
    "domain",
]

# ######## PROJECT HEADING  #########################################
HEADING_SINGLE_PROJECT_FL_SUMMARY = [
    "Cluster (Raw)",
    "Cluster (PF)",
    "Yield (MBases)",
    "Number of Samples",
]

HEADING_SINGLE_PROJECT_STATS_LANE = [
    "Lane",
    "PF Clusters",
    "% of the lane",
    "% Perfect barcode",
    "% One mismatch barcode",
    "Yield (Mbases)",
    "% >= Q30 bases",
    "Mean Quality Score",
]

HEADING_SINGLE_PROJECT_SAMPLES = [
    "Sample",
    "Barcode",
    "PF Clusters",
    "Percent of Project",
    "Yield (Mbases)",
    "% >= Q30 bases",
    "Mean Quality Score",
]


# #######   Sequencer data   #########################
EMPTY_FIELDS_IN_SEQUENCER = [
    "platformID",
    "sequencerDescription",
    "sequencerLocation",
    "sequencerSerialNumber",
    "sequencerOperationStart",
]


# ####### Sample comparation #########################
HEADING_COMPARATION_SAMPLE_LIST = [
    "Run Name",
    "User name",
    "Sample Nane",
    "Finish date",
    "Select Sample",
]
HEADING_COMPARATION_SAMPLE_INFORMATION = [
    "Sample Name",
    "Barcode Name",
    "Run Name",
    "Number of samples",
    "pfClusters",
    "% in Project",
    "yielddMB",
    "Quality >Q30",
    "Mean Quality",
]
