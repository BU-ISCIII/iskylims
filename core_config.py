



################# SAMPLE SETTINGS ##############################
### Headings used when recording information
HEADING_FOR_RECORD_SAMPLES = [ 'Patient Code ID', 'Sample Name', 'Lab requested', 'Type of Sample', 'Species', 'Project/Service', 'Date sample reception', 'Sample Storage Location']

HEADING_FOR_OPTIONAL_FIELD_SAMPLES = [ 'Patient Code ID', 'Lab requested', 'Species', 'Date sample reception', 'Sample Storage Location']

HEADING_FOR_DISPLAY_RECORDED_SAMPLES = ['Unique Sample ID', 'Sample CodeID', 'Sample Name', 'Date for entry in Lab', 'Type of Sample']

HEADING_FOR_COMPLETION_SAMPLES_PRE_DEFINED = ['Date sample extraction' , 'Sample CodeID', 'Sample name']

MAPPING_SAMPLE_FORM_TO_DDBB = [('Patient Code ID','p_code_id'), ('Sample Name', 'sampleName'), ('Lab requested', 'samplesOrigin'),
                ('Type of Sample','sampleType'), ('Species', 'species'),('Project/Service', 'project_service'),('Date sample extraction', 'sampleEntryDate'),
                ('Sample Storage Location', 'sampleLocation'), ()]

HEADING_FOR_MOLECULE_PROTOCOL_DEFINITION = ['Sample ID','Sample Type', 'Molecule type', 'Type of Extraction', 'Extraction date', 'Protocol to be used']
HEADING_FOR_MOLECULE_ADDING_PARAMETERS = ['Molecule Code ID', 'Lot Commercial Kit']


#########
### Headings to confirm the sucessful recorded
HEADING_CONFIRM_MOLECULE_RECORDED = ['Molecule Code ID', 'Used Protocol']

### Heading values when showing pending samples
HEADING_FOR_DEFINED_SAMPLES_STATE = ['Sample extraction date', 'Sample Code ID', 'Sample','To be included']

HEADING_FOR_EXTRACTED_MOLECULES_STATE = ['Sample extraction date', 'Sample', 'Molecule Code ID', 'Molecule Extraction Date',
                                    'Used Protocol', 'Select sample']

HEADING_FOR_USER_PENDING_MOLECULES = ['Sample', 'Molecule Code ID', 'Molecule Extraction Date', 'Used Protocol', 'Select Molecule']

HEADING_FOR_DEFINED_SAMPLES_STATE_WETLAB_MANAGER = ['Sample extraction date', 'Sample', 'Sample Code ID', 'UserID']
HEADING_FOR_EXTRACTED_MOLECULES_STATE_WETLAB_MANAGER = ['Sample extraction date', 'Sample', 'Molecule Code ID', 'Molecule Extraction Date',
                                    'Used Protocol', 'UserID']

### Heading for display information on sample definition
HEADING_FOR_SAMPLE_DEFINITION = ['Sample Name', 'Sample CodeID','Sample State', 'Recorded Date', 'Sample Type', 'Species', 'Number of reused', 'User']
### Heading for display information on molecule definition
HEADING_FOR_MOLECULE_DEFINITION = ['Molecule CodeID', 'Molecule State','Extraction Date', 'Extraction Type', 'Molecule Type', 'Used for', 'Used Protocol', 'Number of reused']

HEADING_FOR_SELECTING_MOLECULE_USE = ['Sample Name','Molecule CodeID', 'Molecule use for']


################# PROTOCOL PARAMETER SETTINGS ##############################
### Headings used when defining the custom protocol parameters
HEADING_FOR_DEFINING_PROTOCOL_PARAMETERS = ['Parameter name', 'Order', 'Used', 'Parameter Type' ,'Option Values' ,'Min Value', 'Max Value', 'Description']

HEADING_FOR_SAMPLE_PROJECT_FIELDS = ['Field name', 'Order', 'Used', 'Searchable','Field type', 'Option Values', 'Description']

HEADING_FOR_MODIFY_SAMPLE_PROJECT_FIELDS = ['Old field name','New field name', 'Order', 'Used', 'Searchable','Field type', 'Option Values', 'Description']

################# PROJECT  SETTINGS #######################################
### Headings used when defining the custom project fields
HEADING_FOR_DEFINING_PROJECT_FIELDS = ['Field name', 'Order', 'Used', 'Description']


##############################################################
################## KITS Heading ##############################
##############################################################

HEADING_FOR_NEW_SAVED_COMMERCIAL_PROTOCOL_KIT = ['Commercial Kit Name', 'Provider',  'Protocol used for the kit']

HEADING_FOR_NEW_SAVED_COMMERCIAL_PLATFORM_KIT = ['Commercial Kit Name', 'Provider', 'Platform used for kit']

HEADING_FOR_COMMERCIAL_PROTOCOL_KIT_BASIC_DATA = ['Protocols used for Commercial Kit', 'Provider',  'Cat Number']

HEADING_FOR_COMMERCIAL_PLATFORM_KIT_BASIC_DATA = ['Platform used for Commercial Kit', 'Provider',  'Cat Number']

HEADING_FOR_COMMERCIAL_KIT_BASIC_DATA = ['Name', 'Provider', 'Protocols used for Commercial Kit', 'Platform', 'Cat Number']

HEADING_FOR_LOT_USER_COMMERCIAL_KIT_BASIC_DATA = ['Commercial Kit', 'Lot Barcode', 'Expiration Date']

HEADING_FOR_USER_LOT_INVENTORY = ['Lot Number' , 'Expiration date', 'Number of use' , 'Run Out']

HEADING_FOR_RUNOUT_USER_LOT_INVENTORY = ['Lot Number' , 'Expiration date', 'Number of use']

HEADING_FOR_USER_LOT_SEARCH_RESULTS = ['Lot Number', 'Commercial Kit', 'Expiration Date', 'Protocol', 'Platform']

HEADING_FOR_DISPLAY_IN_SAMPLE_INFO_USER_KIT_DATA = ['Molecule Code ID','Lot number', 'Commercial kit name', 'Expiration Date']


#########  FOLDER SETTINGS TO STORE COLLECTION KITS #######
## Relative path from settings.MEDIA_ROOT
COLLECTION_INDEX_KITS_DIRECTORY = 'core/collection_index_kits/'

## Configuration for index library file
COLLECTION_INDEX_HEADING = ['[Version]','[Name]', '[PlateExtension]','[Settings]', '[I7]']



###################### FIELDS NAME USED IN USER FORM  ########################################

SAMPLE_PROJECT_MAIN_DATA = ['Sample Project Name', 'Project Manager', 'Contact data', 'Description']

FORM_PROJECT_CREATION = ['projectName', 'projectManager', 'projectContact', 'projectDescription']


##################### DEFAULT INFORMATION  ################
VALUE_NOT_PROVIDED = 'Not Provided'


##################### ERROR MESSAGES  #####################
ERROR_TYPE_OF_SAMPLE_EXISTS = ['Type of sample is already recorded',]
ERROR_TYPE_OF_SAMPLE_ID_DOES_NOT_EXISTS = ['The type of sample that you request does not exist']
ERROR_MOLECULE_USE_FOR_EXISTS = ['Molecule use has been already recorded']

ERROR_SPECIES_ALREADY_DEFINED = ['Species name is already defined']
ERROR_LABORATORY_REQUEST_ALREADY_DEFINED = ['Laboratory/Institution request is already defined']
ERROR_MOLECULE_TYPE_ALREADY_DEFINED = ['Molecule Type is already defined']
ERROR_PROTOCOL_TYPE_ALREADY_DEFINED = ['Protocol Type is already defined']
