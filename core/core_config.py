# ################ SAMPLE SETTINGS ##############################
# ## Headings used when recording information
HEADING_FOR_RECORD_SAMPLES = [
    "Patient Code ID",
    "Sample Name",
    "Lab requested",
    "Type of Sample",
    "Species",
    "Project/Service",
    "Date sample reception",
    "Collection Sample Date",
    "Sample Storage Location",
    "Only recorded",
]

HEADING_FOR_OPTIONAL_FIELD_SAMPLES = [
    "Patient Code ID",
    "Lab requested",
    "Species",
    "Date sample reception",
    "Collection Sample Date",
    "Sample Storage Location",
    "Only recorded",
]

HEADING_FOR_DISPLAY_RECORDED_SAMPLES = [
    "Unique Sample ID",
    "Sample CodeID",
    "Sample Name",
    "Date for entry in Lab",
    "Type of Sample",
]

HEADING_FOR_COMPLETION_SAMPLES_PRE_DEFINED = [
    "Date sample extraction",
    "Sample CodeID",
    "Sample name",
    "Select sample",
]

MAPPING_SAMPLE_FORM_TO_DDBB = [
    ("Patient Code ID", "p_code_id"),
    ("Sample Name", "sampleName"),
    ("Lab requested", "labRequest"),
    ("Type of Sample", "sampleType"),
    ("Species", "species"),
    ("Project/Service", "project_service"),
    ("Date sample extraction", "sampleEntryDate"),
    ("Collection Sample Date", "collectionSampleDate"),
    ("Sample Storage Location", "sampleLocation"),
    ("Only recorded", "onlyRecorded"),
]

HEADING_FOR_MOLECULE_PROTOCOL_DEFINITION = [
    "Sample ID",
    "Sample Type",
    "Molecule type",
    "Type of Extraction",
    "Extraction date",
    "Protocol to be used",
]
HEADING_FOR_MOLECULE_ADDING_PARAMETERS = ["Molecule Code ID", "Lot Commercial Kit"]


# ########### Headings to confirm the sucessful recorded
HEADING_CONFIRM_MOLECULE_RECORDED = ["Molecule Code ID", "Used Protocol"]

# ## Heading values when showing pending samples
"""
HEADING_FOR_DEFINED_SAMPLES_STATE = [
    "Sample extraction date",
    "Sample Code ID",
    "Sample",
    "To be included",
]
"""
HEADING_FOR_EXTRACTED_MOLECULES_STATE = [
    "Sample extraction date",
    "Sample",
    "Molecule Code ID",
    "Molecule Extraction Date",
    "Used Protocol",
    "Select sample",
]

HEADING_FOR_PENDING_MOLECULES = [
    "Sample",
    "Molecule Code ID",
    "Molecule Extraction Date",
    "Used Protocol",
    "Select Molecule",
]

HEADING_FOR_DEFINED_SAMPLES_STATE_WETLAB_MANAGER = [
    "Sample extraction date",
    "Sample",
    "Sample Code ID",
    "UserID",
]
HEADING_FOR_EXTRACTED_MOLECULES_STATE_WETLAB_MANAGER = [
    "Sample extraction date",
    "Sample",
    "Molecule Code ID",
    "Molecule Extraction Date",
    "Used Protocol",
    "UserID",
]

# ## Heading for display information on sample definition
HEADING_FOR_SAMPLE_DEFINITION = [
    "Sample Name",
    "Sample CodeID",
    "Sample State",
    "Recorded Date",
    "Collection Sample Date",
    "Date sample reception",
    "Sample Type",
    "Species",
    "Number of reused",
    "User",
]
# ## Heading for display information on molecule definition
HEADING_FOR_MOLECULE_DEFINITION = [
    "Molecule CodeID",
    "Molecule State",
    "Extraction Date",
    "Extraction Type",
    "Molecule Type",
    "Used for",
    "Used Protocol",
    "Number of reused",
]

HEADING_FOR_SELECTING_MOLECULE_USE = [
    "Sample Name",
    "Molecule CodeID",
    "Molecule use for",
]

# ################ PROTOCOL PARAMETER SETTINGS ##############################
# ## Headings used when defining the custom protocol parameters
HEADING_FOR_DEFINING_PROTOCOL_PARAMETERS = [
    "Parameter name",
    "Order",
    "Used",
    "Parameter Type",
    "Option Values",
    "Min Value",
    "Max Value",
    "Description",
]
HEADING_FOR_MODIFY_PROTOCOL_FIELDS = [
    "Old field name",
    "New field name",
    "Order",
    "Used",
    "Parameter Type",
    "Option Values",
    "Description",
]

HEADING_FOR_SAMPLE_PROJECT_FIELDS = [
    "Field name",
    "Order",
    "Used",
    "Searchable",
    "Field type",
    "Option Values",
    "Description",
    "Classification",
]

HEADING_FOR_MODIFY_SAMPLE_PROJECT_FIELDS = [
    "Field name",
    "Change field name",
    "Order",
    "Used",
    "Searchable",
    "Field type",
    "Option Values",
    "Description",
]

# ################ PROJECT  SETTINGS #######################################
# ## Headings used when defining the custom project fields
HEADING_FOR_DEFINING_PROJECT_FIELDS = ["Field name", "Order", "Used", "Description"]


# #############################################################
# ################# KITS Heading ##############################
# #############################################################

HEADING_FOR_NEW_SAVED_COMMERCIAL_PROTOCOL_KIT = [
    "Commercial Kit Name",
    "Provider",
    "Protocol used for the kit",
]

HEADING_FOR_NEW_SAVED_COMMERCIAL_PLATFORM_KIT = [
    "Commercial Kit Name",
    "Provider",
    "Platform used for kit",
]

HEADING_FOR_COMMERCIAL_PROTOCOL_KIT_BASIC_DATA = [
    "Protocols used for Commercial Kit",
    "Provider",
    "Cat Number",
]

HEADING_FOR_COMMERCIAL_PLATFORM_KIT_BASIC_DATA = [
    "Platform used for Commercial Kit",
    "Provider",
    "Cat Number",
]

HEADING_FOR_COMMERCIAL_KIT_BASIC_DATA = [
    "Name",
    "Provider",
    "Protocols used for Commercial Kit",
    "Platform",
    "Cat Number",
]

HEADING_FOR_LOT_USER_COMMERCIAL_KIT_BASIC_DATA = [
    "Commercial Kit",
    "Lot Barcode",
    "Expiration Date",
]

HEADING_FOR_USER_LOT_KIT_INVENTORY = [
    "Kit used for",
    "Lot number",
    "Recorded by",
    "Number of use",
    "Used last date",
    "Expiration date",
]

HEADING_FOR_USER_LOT_SEARCH_RESULTS = [
    "Lot Number",
    "Commercial Kit",
    "Expiration Date",
    "Protocol",
    "Platform",
]

HEADING_FOR_DISPLAY_IN_SAMPLE_INFO_USER_KIT_DATA = [
    "Molecule Code ID",
    "Lot number",
    "Commercial kit name",
    "Expiration Date",
]


# ########  FOLDER SETTINGS TO STORE COLLECTION KITS #######
# # Relative path from settings.MEDIA_ROOT
COLLECTION_INDEX_KITS_DIRECTORY = "core/collection_index_kits/"

# # Configuration for index library file
COLLECTION_INDEX_HEADING = [
    "[Version]",
    "[Name]",
    "[PlateExtension]",
    "[Settings]",
    "[I7]",
]


# ##################### FIELDS NAME USED IN USER FORM  ########################################

SAMPLE_PROJECT_MAIN_DATA = [
    "Sample Project Name",
    "Project Manager",
    "Contact data",
    "Description",
]

FORM_PROJECT_CREATION = [
    "projectName",
    "projectManager",
    "projectContact",
    "projectDescription",
]


# #################### DEFAULT INFORMATION  ################
VALUE_NOT_PROVIDED = "Not Provided"

SUCCESSFUL_JSON_SCHEMA = ["Json schema was loaded successfully"]
# #################### ERROR MESSAGES  #####################
ERROR_TYPE_OF_SAMPLE_EXISTS = ["Type of sample is already recorded"]
ERROR_TYPE_OF_SAMPLE_ID_DOES_NOT_EXISTS = (
    "The type of sample that you request does not exist"
)
ERROR_MOLECULE_USE_FOR_EXISTS = ["Molecule use has been already recorded"]

ERROR_SPECIES_ALREADY_DEFINED = ["Species name is already defined"]
ERROR_LABORATORY_REQUEST_ALREADY_DEFINED = [
    "Laboratory/Institution request is already defined"
]
ERROR_MOLECULE_TYPE_ALREADY_DEFINED = ["Molecule Type is already defined"]
ERROR_PROTOCOL_TYPE_ALREADY_DEFINED = ["Protocol Type is already defined"]

ERROR_STATE_ALREADY_DEFINED = ["State is already defined"]
ERROR_CITY_ALREADY_DEFINED = ["City is already defined"]

ERROR_SAMPLE_NOT_FOUND = "Sample was not found"

ERROR_SAMPLE_ALREADY_DEFINED = "Samples already defined"
ERROR_SAMPLE_INCOMPLETED = " Samples with imcompleted data"

# ######################  Batch file ###############################################
ERROR_MESSAGE_FOR_EMPTY_SAMPLE_BATCH_FILE = [
    "The uploaded sample batch file does not have any sample",
    "Upload a valid batch file",
]
ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_LAB_REQUESTED = [
    "The Laboratory",
    "from where the samples are received is not defined",
]
ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_DEFINED_LAB_REQUESTED = [
    "No Laboratory is defined yet",
    "Check documentation to define the Laboratory",
]
ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_EMPTY_VALUE = ["Batch file contains empty values"]
ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_INVALID_FORMAT = [
    "Batch file does not have the correct format"
]
ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_DEFINED_TYPE_OF_SAMPLES = [
    "No Type of Samples are defined yet",
    "Check documentation to define them",
]
ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_SAMPLE_TYPE = [
    "The Type of Sample",
    "is not defined",
]
ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_DEFINED_SPECIES = [
    "No Species are defined yet",
    "Check documentation to define them",
]
ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_SPECIES = ["The specie", "is not defined"]
ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_DEFINED_SAMPLE_PROJECTS = [
    "No Sample Projects are defined yet",
    "Check documentation to define them",
]
ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_SAMPLE_PROJECTS = [
    "The Sample Project",
    "is not defined",
]
ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NOT_SAME_SAMPLE_PROTOCOL = [
    "The batch file must have the same type of samples and same project assigned"
]
ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NOT_SAMPLE_FIELD_DEFINED = [
    "Required project field",
    "is not included in batch file",
]
ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_MOLECULE_NOT_DEFINED = [
    "The mandatory field for molecule",
    "must be included in the batch file",
]
ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_DEFINED_MOLECULE_TYPES = [
    "No molecule types are defined yet",
    "Check documentation to define them",
]
ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_MOLECULE_TYPE = [
    "The molecule Type",
    "is not defined",
]
ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_DEFINED_PROTOCOL = [
    "No protocols name are defined yet",
    "Check documentation to define them",
]
ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_MOLECULE_PROTOCOL_NAME = [
    "The protocol name",
    "is not defined",
]
ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_SAMPLE_PROJECT_NO_DEFINED = [
    "Project/Service is not defined",
]
ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_TOO_MANY_PROJECTS = [
    "Only one project must be defined per sample batch file"
]
ERROR_MESSAGE_INVALID_JSON_SCHEMA = [
    "Upload schema cannot used because it contains errors "
]
ERROR_MESSAGE_PROPERTY_NOT_FOUND_IN_SCHEMA = ["Property was not found in schema"]
ERROR_FIELD_NOT_EXIST_IN_SCHEMA = ["Field does not exists in schema"]
ERROR_MESSAGE_FOR_REPEATED_SAMPLE_BATCH_FILE = [
    "The sample(s)",
    "of the excel file is/are repeated"
]
ERROR_MESSAGE_FOR_REPEATED_SAMPLE_BATCH_DATABASE = [
    "The sample",
    "of the excel file is already repeated on the Data Base"
]
ERROR_MESSAGE_FOR_EMPTY_SAMPLE_BATCH_FILE = ["At least one sample name is empty"]
]