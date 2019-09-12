



################# SAMPLE SETTINGS ##############################
### Headings used when recording information
HEADING_FOR_RECORD_SAMPLES = [ 'Sample Name', 'Laboratory', 'Type of Sample', 'Species', 'Sample Extraction date']

HEADING_FOR_DISPLAY_RECORDED_SAMPLES = ['Unique Sample ID', 'Sample CodeID', 'Sample Name', 'Sample Extraction date', 'Type of Sample']

HEADING_FOR_MOLECULE_PROTOCOL_DEFINITION = ['Sample ID', 'Molecule type', 'Type of Extraction', 'Extraction date', 'Protocol type',
                'Protocol to be used']
HEADING_FOR_MOLECULE_ADDING_PARAMETERS = ['Molecule Code ID', 'Used Protocol']

#########
### Headings to confirm the sucessful recorded
HEADING_CONFIRM_MOLECULE_RECORDED = ['Molecule Code ID', 'Used Protocol']

### Heading values when showing pending samples
HEADING_FOR_DEFINED_SAMPLES_STATE = ['Sample extraction date', 'Sample Code ID', 'Sample','To be included']
HEADING_FOR_EXTRACTED_MOLECULES_STATE = ['Sample extraction date', 'Sample', 'Molecule Code ID', 'Molecule Extraction Date',
                                    'Used Protocol', 'To be included']


### Heading for display information on sample definition
HEADING_FOR_SAMPLE_DEFINITION = ['Sample Name', 'Sample CodeID','Sample State', 'Recorded Date', 'Sample Type', 'Species', 'Number of reused', 'User']
### Heading for display information on molecule definition
HEADING_FOR_MOLECULE_DEFINITION = ['Molecule CodeID', 'Molecule State','Extraction Date', 'Extraction Type', 'Molecule Type', 'Used Protocol', 'Number of reused']



################# PROTOCOL PARAMETER SETTINGS ##############################
### Headings used when defining the custom protocol parameters
HEADING_FOR_DEFINING_PROTOCOL_PARAMETERS = ['Parameter name', 'Order', 'Used', 'Min Value', 'Max Value', 'Description']



#########  FOLDER SETTINGS TO STORE COLLECTION KITS #######
## Relative path from settings.MEDIA_ROOT
COLLECTION_INDEX_KITS_DIRECTORY = 'core/collection_index_kits/'

## Configuration for index library file
COLLECTION_INDEX_HEADING = ['[Version]','[Name]', '[PlateExtension]','[Settings]', '[I7]']
