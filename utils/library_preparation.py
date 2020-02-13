import json
from iSkyLIMS_core.models import Samples, MoleculePreparation, Protocols
from iSkyLIMS_core.utils.handling_commercial_kits import *
from iSkyLIMS_core.utils.handling_samples import  get_sample_obj_from_sample_name
from iSkyLIMS_wetlab.models import *
from iSkyLIMS_wetlab.wetlab_config import *
from iSkyLIMS_wetlab.utils.sample_sheet_utils import *
from iSkyLIMS_wetlab.utils.collection_index_functions import check_collection_index_exists
from ..fusioncharts.fusioncharts import FusionCharts
from .stats_graphics import *



def analyze_input_param_values(request):
    if  'lib_prep_in_list' in request.POST:
        lib_prep_ids = request.POST.getlist('lib_prep_id')
        if len('lib_prep_in_list') == 0:
            lib_prep_ids = list(request.POST['lib_prep_id'])
    else:
        lib_prep_ids = request.POST['lib_prep_id'].split(',')
    headings = request.POST['heading_in_excel'].split(',')
    json_data = json.loads(request.POST['protocol_data'])
    fixed_heading_length = len(HEADING_FIX_FOR_ADDING_LIB_PARAMETERS)
    parameters_length = len(headings)
    stored_params = []
    for i in range(len(lib_prep_ids)):
        library_prep_obj = LibraryPreparation.objects.get(pk = lib_prep_ids[i])

        for p_index in range(fixed_heading_length, parameters_length):
            lib_parameter_value ={}
            lib_parameter_value['parameter_id'] = ProtocolParameters.objects.get(protocol_id = library_prep_obj.protocol_id,
                                parameterName__exact = headings[p_index])
            lib_parameter_value['library_id'] = library_prep_obj
            lib_parameter_value['parameterValue'] = json_data[i] [p_index]

            new_parameters_data = LibParameterValue.objects.create_library_parameter_value (lib_parameter_value)
        kit_index = HEADING_FIX_FOR_ADDING_LIB_PARAMETERS.index('Lot Regents Kit used')
        library_prep_obj.set_reagent_user_kit(json_data[i] [kit_index])
        stored_params.append([library_prep_obj.get_sample_name(), library_prep_obj.get_lib_prep_code()])
        library_prep_obj.set_state('Completed')
        sample_obj = library_prep_obj.get_sample_obj ()
        # Update the sample state to "Create Pool"
        sample_obj.set_state('Pool Preparation')
    #import pdb; pdb.set_trace()
    return stored_params

def check_samples_for_library_preparation():
    '''
    Description:
        The function checks if there are samples in run was success
    Return:

    '''
    if Samples.objects.filter(sampleState__sampleStateName__exact = 'Library preparation').exists():
        return True
    else:
        return False


def extract_sample_data (s_data):
    headings = s_data['headings']
    sample_list = []
    #columns = ['Sample_ID','Sample_Name','Sample_Plate','Sample_Well','Index_Plate_Well','I7_Index_ID','index','I5_Index_ID','index2','Sample_Project']
    for sample_row in s_data['samples']:
        '''
        if Samples.objects.filter(sampleName__exact = sample_row.index('Sample_Name'), sampleState__sampleStateName = 'Add Library Preparation' ).exists():
            sample_obj = Samples.objects.filter(sampleName__exact = sample_row.index('Sample_Name'), sampleState__sampleStateName = 'Add Library Preparation')
        '''
        lib_prep_data = {}
        for column in wetlab_config.MAP_USER_SAMPLE_SHEET_TO_DATABASE :
            if column[0] in headings:
                lib_prep_data[column[1]] = sample_row[headings.index(column[0])]
            else:
                lib_prep_data[column[1]] = ''
        sample_list.append(lib_prep_data)

    return sample_list


def extract_user_sample_sheet_data(file_in):
    '''
    Description:
        The function checks if the run was success
    Input:
        file_in    # csv file from IEM
    Functions:
        store_user_input_file  # located at utils/sample_sheet_utils.py
        valid_user_iem_file    # located at utils/sample_sheet_utils.py
        get_sample_sheet_data  # located at utils/sample_sheet_utils.py
    Constant:
        ERROR_INVALID_FILE_FORMAT
        ERROR_UNABLE_TO_DELETE_USER_FILE
    Return
        data with the extracted information from sample sheet or error if file is invalid
    '''
    data = {}

    data['full_path_file'], data['file_name'] = store_user_input_file (file_in)

    if not valid_user_iem_file(data['full_path_file']):
        data['ERROR'] = ERROR_INVALID_FILE_FORMAT
        if not delete_stored_file(full_path_file):
            data['ERROR'].append(ERROR_UNABLE_TO_DELETE_USER_FILE)
        return data
    data.update(get_sample_sheet_data(data['full_path_file']))
    return data


def valid_samples_for_lib_preparation(samples):
    '''
    Description:
        The function checks if samples are defined are they are in Library Preparation state ,
    Input:
        samples  # contain a list with sample names
    Functions:
        get_sample_obj_from_sample_name  # located at iSkyLIMS_core/handling_samples.py
    Constant:
        ERROR_SAMPLE_SHEET_CONTAINS_NOT_DEFINED_SAMPLES
        ERROR_SAMPLES_INVALID_STATE_FOR_LIBRARY_PREPARATION
    Return
        error or True
    '''
    sample_objs = []
    invalid_samples = []
    for sample in samples:
        s_obj = get_sample_obj_from_sample_name(sample)
        if not s_obj :
            invalid_samples.append(sample)
        else:
            sample_objs.append(s_obj)

    if len(invalid_samples) > 0:
        error = {}
        error_message = ERROR_SAMPLE_SHEET_CONTAINS_NOT_DEFINED_SAMPLES.copy()
        error_message.insert(1,' , '.join(invalid_samples))
        error['ERROR'] = error_message
        return error
    for sample_obj in sample_objs:
        if not 'Library Preparation' == sample_obj.get_sample_state():
            invalid_samples.append(sample_obj.get_sample_name())
    if len(invalid_samples) > 0:
        error = {}
        error_message = ERROR_SAMPLES_INVALID_STATE_FOR_LIBRARY_PREPARATION.copy()
        error_message.insert(1,' , '.join(invalid_samples))
        error['ERROR'] = error_message
        return error
    return True

def validate_sample_sheet_data (input_data ):
    '''
    Description:
        The function checks if samples are defined are they are in Library Preparation state ,
        if collection index is already defined and no duplication index exists.

    Input:
        input_data  # contain a dictionary with sample names, heading to match the values and
                    a list with all  information samples
    Functions:
        valid_samples_for_lib_preparation  # located at this file
        delete_stored_file              # located at this file
        find_duplicate_index             # located at this file
        check_collection_index_exists   #  located at iSkyLIMS_wetlab/utils/collection_index_functions.py
    Constant:
        ERROR_INVALID_FILE_FORMAT
        ERROR_UNABLE_TO_DELETE_USER_FILE
    Return
        data with the extracted information from sample sheet or error if file is invalid
    '''
    # check that samples are defined and in the right state
    valid_samples = valid_samples_for_lib_preparation (input_data['samples'])
    if 'ERROR'in valid_samples:
        delete_stored_file(input_data['full_path_file'])
        return valid_samples
    # check if sample sheet has duplicate index
    duplicate_index = find_duplicate_index(input_data['sample_data'], input_data['heading'] )
    if 'ERROR' in duplicate_index:
        delete_stored_file(input_data['full_path_file'])
        return duplicate_index
    # check if collection index is defined
    if not check_collection_index_exists (input_data['index_adapter']):
        error = {}
        error_message = ERROR_COLLECTION_INDEX_KIT_NOT_DEFINED.copy()
        error_message.insert(1, input_data['index_adapter'])
        error['ERROR'] = error_message
        delete_stored_file(input_data['full_path_file'])
        return error
    return True

def find_duplicate_index (sample_row_data, heading):
    '''
    Description:
        The function get the sample rows from sample sheet to check if the samples have duplicated
        index.
    Input:
        sample_row_data  # contains the row sample list
        heading     # contains the list of heading to get the index colums
    Constant:
        ERROR_SAMPLES_INVALID_DUPLICATED_INDEXES
    Return
        False if not duplication found. error message if exists duplicted index
    '''
    index_values = {}
    duplicated_index_sample = []
    i5_index = False
    i7_index = heading.index('I7_Index_ID')
    if 'I5_Index_ID' in heading:
        i5_index = heading.index('I5_Index_ID')
    sample_name_index = heading.index('Sample_Name')
    for sample_row in sample_row_data:
        if i5_index :
            indexes_in_sample = str(sample_row[i7_index] + '_' + sample_row[i5_index])
        else:
            indexes_in_sample = sample_row[i7_index]
        if  indexes_in_sample not in index_values:
            index_values[indexes_in_sample] = []
        else:
            duplicated_index_sample.append(sample_row[sample_name_index])
        index_values [sample_row[sample_name_index]].append([indexes_in_sample])

    if len(duplicated_index_sample) == 0:
        return False
    else:
        error = {}
        error_message = ERROR_SAMPLES_INVALID_DUPLICATED_INDEXES.copy()
        error_message.insert(1,' , '.join(duplicated_index_sample))
        error['ERROR'] = error_message
        return error



def get_protocols_for_library_preparation ():
    protocol_list = []
    if Protocols.objects.filter(type__protocol_type__exact ='Library Preparation').exists():
        protocols = Protocols.objects.filter(type__protocol_type__exact = 'Library Preparation')
        for protocol in protocols:
            protocol_list.append(protocol.get_name())
    return protocol_list


def get_all_library_information(sample_id):
    library_information = {}
    if LibraryPreparation.objects.filter(sample_id__pk__exact = sample_id).exists():
        library_information['library_definition_heading'] = HEADING_FOR_LIBRARY_PREPARATION_DEFINITION
        library_information['library_definition'] = []
        library_information['pool_information'] = []
        library_preparation_items = LibraryPreparation.objects.filter(sample_id__pk__exact = sample_id).exclude(libPrepState__libPrepState__exact = 'Created for Reuse')
        library_information['lib_prep_param_value'] = []
        for library_item in library_preparation_items:
            library_information['library_definition'].append(library_item.get_info_for_display())
            protocol_used_obj = library_item.get_protocol_obj()
            if ProtocolParameters.objects.filter(protocol_id = protocol_used_obj).exists():
                parameter_names = ProtocolParameters.objects.filter(protocol_id = protocol_used_obj).order_by('parameterOrder')
                lib_prep_param_heading = ['Lib Preparation CodeID']
                lib_prep_param_value = [library_item.get_lib_prep_code()]
                for p_name in parameter_names:
                    lib_prep_param_heading.append(p_name.get_parameter_name())
                    if LibParameterValue.objects.filter(library_id = library_item).exists():
                        lib_prep_param_value.append(LibParameterValue.objects.get(library_id = library_item, parameter_id = p_name).get_parameter_information())
                library_information['lib_prep_param_value'].append(lib_prep_param_value)
                library_information['lib_prep_param_heading'] = lib_prep_param_heading

            if library_item.pools.all().exists() :
                pools = library_item.pools.all()
                lib_prep_code_id = library_item.get_lib_prep_code()
                for pool in pools:
                    pool_name = pool.get_pool_name()
                    pool_code = pool.get_pool_code_id()
                    library_information['pool_information'].append([lib_prep_code_id, pool_name,pool_code, ''])

        if library_information['pool_information']:
            library_information['pool_heading'] = HEADING_FOR_DISPLAY_POOL_INFORMATION_IN_SAMPLE_INFO
    return library_information

def get_lib_prep_to_add_parameters():
    '''
    Description:
        The function will return a list with samples which are needs to add library preparation parameters
    Input:

    Variables:
        library_prep_information # Dictionary with the heading and the molecule information
    Return:
        lib_prep_parameters.
    '''
    lib_prep_parameters = {}
    lib_prep_parameters['length'] = 0
    if LibraryPreparation.objects.filter(libPrepState__libPrepState__exact = 'Defined').exists():
        samples = LibraryPreparation.objects.filter(libPrepState__libPrepState__exact = 'Defined')
        sample_info = []
        for sample in samples:
            lib_prep_info = []
            lib_prep_info.append(sample.get_lib_prep_code())
            lib_prep_info.append(sample.get_sample_name())
            lib_prep_info.append(sample.get_protocol_used())
            lib_prep_info.append(sample.get_id())
            sample_info.append(lib_prep_info)
        lib_prep_parameters['lib_prep_info'] = sample_info
        lib_prep_parameters['lib_prep_heading'] = HEADING_FOR_ADD_LIBRARY_PREPARATION_PARAMETERS
        lib_prep_parameters['length'] = len(sample_info)
    return lib_prep_parameters



def get_samples_in_lib_prep_state ():
    '''
    Description:
        The function will return a list with samples which are in add_library_preparation state.
        Include the ones that are requested to reprocess
    Return:
        samples_in_lib_prep.
    '''

    samples_in_lib_prep = {}
    lib_prep_data = []
    if Samples.objects.filter(sampleState__sampleStateName__exact = 'Library preparation').exists():
        samples_obj = Samples.objects.filter(sampleState__sampleStateName__exact =  'Library preparation').order_by('sampleUser').order_by('sampleEntryDate')

        for sample in samples_obj :
            if LibraryPreparation.objects.filter(sample_id = sample, libPrepState__libPrepState__exact = 'Defined').exists():
                continue
            if not LibraryPreparation.objects.filter(sample_id = sample).exists() or LibraryPreparation.objects.filter(sample_id = sample, libPrepState__libPrepState__exact = 'Created for Reuse').exists() :
                sample_information = sample.get_info_in_defined_state()
                sample_information.append(sample.get_register_user())
                molecule = MoleculePreparation.objects.filter(sample = sample,state__moleculeStateName = 'Completed').last()
                molecule_data = molecule.get_molecule_information()
                lib_prep_data.append(sample_information + molecule_data)

        samples_in_lib_prep['library_information'] = lib_prep_data
        samples_in_lib_prep['lib_prep_heading'] = HEADING_FOR_LIBRARY_PREPARATION_STATE
        samples_in_lib_prep['length'] = len(lib_prep_data)

        return samples_in_lib_prep
    else :
        samples_in_lib_prep['length'] = 0
        return samples_in_lib_prep

    #if LibraryPreparation.objects.filter(libPrepState__libPrepState__exact = 'Defined').exists():
    #lib_preparations = libraryPreparation.objects.filter(libPrepState__libPrepState__exact = 'Defined')

def pending_samples_for_grafic(pending):
    number_of_pending = {}
    number_of_pending ['DEFINED'] = pending['defined']['length']
    number_of_pending ['EXTRACTED MOLECULE'] = pending['extract_molecule']['length']
    number_of_pending ['LIBRARY PREPARATION'] = pending['create_library_preparation']['length']

    data_source = graphic_3D_pie('Number of Pending Samples', '', '', '', 'fint',number_of_pending)
    graphic_pending_samples = FusionCharts("pie3d", "ex1" , "430", "450", "chart-1", "json", data_source)
    return graphic_pending_samples


def store_library_preparation_samples(sample_sheet_data, form_data,user):
    '''
    Description:
        The function will get the user data in the form and the extracted sample information from sample sheet.
        Then store the libraryPreparation data for each sample and update the sample state to "library Preparation"
    Return:
        True .
    '''

    return


    protocol = request.POST['lib_protocols']
    single_paired = request.POST['singlePairedEnd']
    read_length = request.POST['readlength']
    user_sample_sheet_data = {}
    stored_lib_prep = {}
    stored_lib_prep['data'] = []

    if CollectionIndexKit.objects.filter(collectionIndexName__exact = index_adapters).exists():
        collection_index_kit_id = CollectionIndexKit.objects.get(collectionIndexName__exact = index_adapters)
    else:
        collection_index_kit_id = None
    register_user_obj = User.objects.get(username__exact = request.user.username)
    user_sample_sheet_data['registerUser'] = register_user_obj
    protocol_obj = Protocols.objects.get(name__exact = protocol)
    user_sample_sheet_data['collectionIndexKit_id'] = collection_index_kit_id

    user_sample_sheet_data['sampleSheet'] = file_name
    user_sample_sheet_data['assay'] = assay
    user_sample_sheet_data['adapter1'] = adapter1
    user_sample_sheet_data['adapter2'] = adapter2
    new_user_s_sheet_obj = libPreparationUserSampleSheet.objects.create_lib_prep_user_sample_sheet(user_sample_sheet_data)



    return
