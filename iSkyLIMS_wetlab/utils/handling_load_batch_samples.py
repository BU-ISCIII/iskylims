import pandas as pd
from iSkyLIMS_wetlab.wetlab_config import *
from iSkyLIMS_core.models import *


def check_defined_option_values_in_samples(sample_batch_df, package):
    '''
    Description:
        The Function checks if all values in the data frame are already defined in database.
        It also checks if the mandatory parameters in the type of sample are in dataframe
    Input:
        sample_batch_df     # sample data in dataframe
        package     # package name "iSkyLIMS_wetlab"
    Constants:
        ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_NO_LAB_REQUESTED
    Return:
        sample_batch_df
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
    # check if mandatory fields in type of sample are included
    
    import pdb; pdb.set_trace()

    return

def check_record_samples_optional_parameters(sample_batch_data):
    '''
    Description:
        The Function check if the only one type of sample is in data and if the mandatory parameters
        for the type of sample are included.
        In case that are more than 1 type of sample it checks one by one the mandatory parameters
    Constants:
        ERROR_MESSAGE_FOR_EMPTY_SAMPLE_BATCH_FILE
    Return:
        sample_batch_data
    '''

    return

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
    Constants:
        ERROR_MESSAGE_FOR_EMPTY_SAMPLE_BATCH_FILE
    Return:
        sample_batch_data
    '''

    # Check if dataframe has empty values
    if sample_batch_df.isnull().values.any() :
        error_cause = ERROR_MESSAGE_FOR_SAMPLE_BATCH_FILE_EMPTY_VALUE
        return error_cause
    check_opt_values = check_defined_option_values_in_samples(sample_batch_df, package)
    if check_opt_values != 'OK':
        return check_opt_values
    check_record_samples_optional_parameters(sample_batch_df)
    return
