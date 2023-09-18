# ################ SAMPLE SETTINGS ##############################

HEADING_BATCH = [
    "patient_core",
    "sample_name",
    "lab_request",
    "sample_type",
    "species",
    "sample_project",
    "sample_entry_date",
    "collection_sample_date",
    "sample_location",
    "only_recorded",
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
HEADING_FOR_PENDING_MOLECULES = [
    "Sample",
    "Molecule Code ID",
    "Molecule Extraction Date",
    "Used Protocol",
    "Select Molecule",
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
    "Classification",
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

ERROR_SAMPLE_ALREADY_DEFINED = ["Sample", "already exist in the database"]

# ######################  Batch file ###############################################
ERROR_EMPTY = [
    "The uploaded table or batch file does not have any sample. Upload a valid batch file"
]
ERROR_NO_LAB_REQUESTED = [
    "The Laboratory",
    "from where the samples are received is not defined",
]
ERROR_NO_DEFINED_LAB_REQUESTED = [
    "No Laboratory is defined yet. Check documentation to define the Laboratory"
]
ERROR_BATCH_INVALID_HEADER = ["The following columns don't have correct format:"]
ERROR_NO_DEFINED_TYPE_OF_SAMPLES = [
    "No Type of Samples are defined yet. Check documentation to define them"
]
ERROR_EMPTY_SAMPLE_TYPE = ["The Type of sample is empty and it is mandatory"]
ERROR_NO_SAMPLE_TYPE = [
    "The Type of sample",
    "is not defined in the database.",
]
ERROR_NO_DEFINED_SPECIES = [
    "No Species are defined yet. Check documentation to define them"
]
ERROR_NO_SPECIES = ["The specie", "is not defined"]
ERROR_NO_DEFINED_SAMPLE_PROJECTS = [
    "No Sample Projects are defined yet. Check documentation to define them"
]
ERROR_EMPTY_PROJECT = ["At least one project value is empty. If there is no project associated write 'None' as the project." ]
ERROR_NO_SAMPLE_PROJECTS = [
    "The Sample Project",
    "is not defined in the database",
]
ERROR_NOT_SAME_SAMPLE_PROTOCOL = [
    "The batch file must have the same type of samples",
    "Same Sample project",
    "and the same Protocol Name",
]
ERROR_PROJECT_FIELD_NOT_DEFINED = [
    "Required sample project field",
    "is not included in batch file",
]
ERROR_MOLECULE_FIELD_NOT_DEFINED = [
    "The mandatory field for molecule",
    "must be included in the batch file",
]
ERROR_NO_DEFINED_MOLECULE_TYPES = [
    "No molecule types are defined yet",
    "Check documentation to define them",
]
ERROR_NO_MOLECULE_TYPE = [
    "The molecule Type",
    "is not defined",
]
ERROR_NO_DEFINED_PROTOCOL = [
    "No protocols are defined yet",
    "Check documentation to define them",
]
ERROR_NO_MOLECULE_PROTOCOL_NAME = [
    "The protocol name",
    "is not defined",
]
ERROR_SAMPLE_PROJECT_NO_DEFINED = [
    "Project/Service is not defined",
]
ERROR_TOO_MANY_PROJECTS = ["Only one project must be defined per sample batch file"]
ERROR_INVALID_JSON_SCHEMA = ["Upload schema cannot used because it contains errors "]
ERROR_FIELD_NOT_EXIST_IN_SCHEMA = ["Field does not exists in schema"]
ERROR_EMPTY_SAMPLE_NAME = ["Sample name in line", "is empty"]
ERROR_REPEATED_SAMPLE_NAME = ["Sample name in line", "is repeated in the table"]
ERROR_ONLY_RECORDED_FIELD = ["Only recorded field must be True, False or empty"]
ERROR_DATE_FORMAT_FIELD = ["Date must have date format. For example YYYY-MM-DD"]
ERROR_MISSING_MANDATORY = [
    "The following columns are empty and are mandatory for the sample type:"
]
ERROR_PROJECT_FIELD_NOTSTRING = ["Project field", "muest be a normal word or sentence"]
ERROR_PROJECT_FIELD_NODATE = ["Project field", "must have date format. For example YYYY-MM-DD"]
ERROR_PROJECT_FIELD_NOOPTION = ["Project field", "only has the following options:"]
ERROR_PROJECT_FIELD_EMPTY = ["Project field", "is empty"]