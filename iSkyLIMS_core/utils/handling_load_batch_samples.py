import pandas as pd
from iSkyLIMS_core.core_config import *
from iSkyLIMS_core.models import *


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

    # check if mandatory fields in type of sample are included

    import pdb; pdb.set_trace()

    return 'OK'

def check_defined_option_values_in_samples_projects(sample_batch_df, package):
    '''
    Description:
        The Function checks if parameters and values (related to new sample project) inside data frame are already defined in database.

    Input:
        sample_batch_df     # sample data in dataframe
        package     # name of the apps that request the checking
    Constants:
    '''
    return

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
        check_record_samples_optional_parameters # located at this file
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
    check_opt_par = check_record_samples_optional_parameters(sample_batch_df)
    if check_opt_par != 'OK':
        return check_opt_par
    return
