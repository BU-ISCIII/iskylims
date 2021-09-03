import pandas as pd
from iSkyLIMS_core.core_config import *
from iSkyLIMS_core.models import *


def check_samples_belongs_to_same_type_and_molecule_protocol(sample_batch_data):
    '''
    Description:
        The Function check if only one type of sample is in data and if the mandatory parameters
        for the type of sample are included.
        In case that are more than 1 type of sample it checks one by one the mandatory parameters
    Input:

    Return:
        True or False
    '''
    sample_type = sample_batch_data['Type of Sample'].unique().tolist()
    sample_project = sample_batch_data['Project/Service'].unique().tolist()
    protocol = sample_batch_data['ProtocolName'].unique().tolist()
    if len(sample_type) != 1 or len(sample_project) != 1 or len(protocol) != 1:
        return False
    return True

def check_defined_option_values_in_samples(sample_batch_df, package):
    '''
    Description:
        The Function checks if values (related to new sample creation) inside data frame are already defined in database.
        It also checks if the mandatory parameters in the type of sample are in dataframe
    Input:
        sample_batch_df     # sample data in dataframe
        package     # name of the apps that request the checking
    Constants:
        ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_LAB_REQUESTED
        ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_DEFINED_TYPE_OF_SAMPLES
        ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_DEFINED_SPECIES
        ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_DEFINED_SAMPLE_PROJECTS
        ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_SAMPLE_PROJECTS
        ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NOT_SAMPLE_FIELD_DEFINED
    Return:
        OK or error code
    '''
    # Check if option values are already defined
    unique_lab_request = sample_batch_df['Lab requested'].unique().tolist()
    if not  LabRequest.objects.all().exists():
        return ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_DEFINED_LAB_REQUESTED
    lab_requested_values = list( LabRequest.objects.all().values_list('labNameCoding', flat=True))
    for lab_request in unique_lab_request:
        if lab_request not in lab_requested_values:
            error_cause = ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_LAB_REQUESTED.copy()
            error_cause.insert(1,lab_request)
            return error_cause
    if not SampleType.objects.filter(apps_name = package).exists():
        return ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_DEFINED_TYPE_OF_SAMPLES
    sample_type_values = list(SampleType.objects.filter(apps_name = package).values_list('sampleType', flat=True))
    unique_sample_types = sample_batch_df['Type of Sample'].unique().tolist()
    for sample_type in unique_sample_types:
        if sample_type not in sample_type_values:
            error_cause = ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_SAMPLE_TYPE.copy()
            error_cause.insert(1,sample_type)
            return error_cause
    if not Species.objects.filter(apps_name = package).exists():
        return ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_DEFINED_SPECIES
    species_values = list(Species.objects.filter(apps_name = package).values_list('speciesName', flat=True))
    unique_species = sample_batch_df['Species'].unique().tolist()
    for specie in unique_species:
        if specie not in species_values:
            error_cause = ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_SPECIES.copy()
            error_cause.insert(1,specie)
            return error_cause
    # check if sample projects are defined
    if not SampleProjects.objects.filter(apps_name = package).exists():
        return ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_DEFINED_SAMPLE_PROJECTS
    sample_project_values = list(SampleProjects.objects.filter(apps_name = package).values_list('sampleProjectName', flat=True))
    unique_sample_projects = sample_batch_df['Project/Service'].unique().tolist()
    for sample_project in unique_sample_projects:
        if sample_project not in sample_project_values:
            error_cause = ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_SAMPLE_PROJECTS.copy()
            error_cause.insert(1,sample_project)
            return error_cause

    # check if additional sample Project parameters are include in bathc file
    if SampleProjectsFields.objects.filter(sampleProjects_id__sampleProjectName__exact = sample_project, sampleProjects_id__apps_name = package).exists():
        sample_project_fields_objs = SampleProjectsFields.objects.filter(sampleProjects_id__sampleProjectName__exact = sample_project, sampleProjects_id__apps_name = package)
        for sample_project_fields_obj in sample_project_fields_objs:
            if not sample_project_fields_obj.get_field_name() in sample_batch_df:
                error_cause = ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NOT_SAMPLE_FIELD_DEFINED.copy()
                error_cause.insert(1,sample_project_fields_obj.get_field_name())
                return error_cause

    return 'OK'

def check_molecule_has_same_data_type(sample_batch_df, package):
    '''
    Description:
        The Function checks if values (related to molecule creation) inside data frame are already defined in database.
    Input:
        sample_batch_df     # sample data in dataframe
        package     # name of the apps that request the checking
    Constants:
        ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_DEFINED_MOLECULE_TYPES

        ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_MOLECULE_PROTOCOL_NAME
    '''
    mandatory_fields = ['Molecule Type', 'Type of Extraction', 'Extraction date', 'ProtocolName'  ]
    for m_field in mandatory_fields:
        if not m_field in sample_batch_df :
            error_cause = ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_MOLECULE_NOT_DEFINED.copy()
            error_cause.insert(1, m_field)
            return error_cause
    if not  MoleculeType.objects.filter(apps_name__exact = package).exists():
        return ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_DEFINED_MOLECULE_TYPES
    molecule_values = list(MoleculeType.objects.filter(apps_name__exact = package).values_list('moleculeType', flat=True))
    unique_molecule_types = sample_batch_df['Molecule Type'].unique().tolist()
    for molecule_type in unique_molecule_types:
        if molecule_type not in molecule_values:
            error_cause = ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_MOLECULE_TYPE.copy()
            error_cause.insert(1,molecule_type)
            return error_cause
    if not Protocols.objects.filter(type__molecule__moleculeType__exact = molecule_type, type__molecule__apps_name__exact = package ).exists():
        return ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_DEFINED_PROTOCOL
    protocol_values = list(Protocols.objects.filter(type__molecule__moleculeType__exact = molecule_type, type__molecule__apps_name__exact = package).values_list('name', flat=True))
    unique_protocol_names = sample_batch_df['ProtocolName'].unique().tolist()
    for protocol_name in unique_protocol_names:
        if protocol_name not in protocol_values:
            error_cause = ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_MOLECULE_PROTOCOL_NAME.copy()
            error_cause.insert(1,protocol_name)
            return error_cause

    return 'OK'

def read_batch_sample_file(batch_file):
    '''
    Description:
        The Function read the batch file and return a panda data frame with the information
    Constants:
        ERROR_MESSAGE_FOR_EMPTY_SAMPLE_BATCH_FILE
    Return:
        sample_batch_data
    '''

    sample_batch_df = pd.read_excel(batch_file, sheet_name=0)
    sample_batch_heading = sample_batch_df.columns.values
    num_rows, num_cols = sample_batch_df.shape
    if num_rows == 0 :
        sample_data = {}
        sample_data['ERROR'] = ERROR_MESSAGE_FOR_EMPTY_SAMPLE_BATCH_FILE
        return sample_data
    '''
    for index, row_data in sample_batch_df.iterrows():
        sample_data = {}
        for index_col in range(len(sample_batch_heading)):
            sample_data[sample_batch_heading[index_col]] = row_data[sample_batch_heading[index_col]]
        #import pdb; pdb.set_trace()
        sample_batch_data.append(sample_data)
    #import pdb; pdb.set_trace()
    '''
    return sample_batch_df


def valid_sample_batch_file (sample_batch_df, package):
    '''
    Description:
        The Function check if all parameters required in the samples are included in file
    Input:
        sample_batch_df     # sample data in dataframe
        package     # package name "iSkyLIMS_wetlab"
    Functions:
        check_samples_belongs_to_same_type_and_molecule_protocol # located at this file
        check_defined_option_values_in_samples # located at this file
        check_record_samples_optional_parameters # located at this file
         # located at this file
    Constants:
        ERROR_MESSAGE_FOR_EMPTY_SAMPLE_BATCH_FILE
        ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NOT_SAME_SAMPLE_PROTOCOL
    Return:
        sample_batch_data
    '''

    # Check if dataframe has empty values
    if sample_batch_df.isnull().values.any() :
        return ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_EMPTY_VALUE
    if not check_samples_belongs_to_same_type_and_molecule_protocol(sample_batch_df):
        return ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NOT_SAME_SAMPLE_PROTOCOL
    check_opt_values = check_defined_option_values_in_samples(sample_batch_df, package)
    if check_opt_values != 'OK':
        return ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NOT_SAME_SAMPLE_PROTOCOL
    # check molecule columns data
    check_molecule_par = check_molecule_has_same_data_type(sample_batch_df, package)
    if check_molecule_par != 'OK':
        return check_molecule_par
    return 'OK'


def save_samples_in_batch_file (sample_batch_df):
    '''
    Description:
        The Function save the sample and the molecule information in database
    Input:
        sample_batch_df     # sample data in dataframe

    '''

    for sample_data in sample_batch_df:
        new_sample =  create_sample_from_file(sample_data)
        new_molecule = create_molecule_from_file(new_sample,sample_data)
        create_molecule_parameter_DNA_from_file(new_molecule, sample_data,new_sample)
        create_sample_project_fields_value(new_sample,sample_data)
