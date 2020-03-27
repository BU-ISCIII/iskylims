import json
from iSkyLIMS_core.models import Samples, MoleculePreparation, Protocols
from iSkyLIMS_core.utils.handling_commercial_kits import *
from iSkyLIMS_core.utils.handling_protocols import *
from iSkyLIMS_core.utils.handling_samples import  get_sample_obj_from_sample_name
from iSkyLIMS_wetlab.models import *
from iSkyLIMS_wetlab.wetlab_config import *
from iSkyLIMS_wetlab.utils.sample_sheet_utils import *
from iSkyLIMS_wetlab.utils.collection_index_functions import check_collection_index_exists
from ..fusioncharts.fusioncharts import FusionCharts
from .stats_graphics import *

def check_empty_fields (data):
    '''
    Description:
        The function check if row_data contains empty values.
    Input:
        data:       # data to be checked
    Return:
        False is all fields contain data. True if any of them are empty
    '''

    for row_data in data :
        for field in row_data :
            if field == '':
                return True
    return False

def analyze_input_param_values(form_data):
    '''
    Description:
        The function get the user input  for the library preparation parameters and store them
        in database.
    Input:

    Constant:
        HEADING_FIX_FOR_ADDING_LIB_PARAMETERS
        ERROR_EMPTY_VALUES
    Return:
        ERROR message if some of the data are missing, or a list of recorded library prepartion obj
    '''
    if  'lib_prep_in_list' in form_data:
        lib_prep_ids = form_data.getlist('lib_prep_ids')
        if len('lib_prep_in_list') == 0:
            lib_prep_ids = list(form_data['lib_prep_ids'])
    else:
        lib_prep_ids = form_data['lib_prep_ids'].split(',')
    lib_prep_code_ids = form_data['lib_prep_code_ids'].split(',')
    headings = form_data['heading_in_excel'].split(',')
    json_data = json.loads(form_data['protocol_data'])
    fixed_heading_length = len(HEADING_FIX_FOR_ADDING_LIB_PARAMETERS)
    parameters_length = len(headings)

    if check_empty_fields(json_data) :
        stored_params = {}
        stored_params['ERROR'] = ERROR_EMPTY_VALUES
        return stored_params

    stored_params = []
    for row_index in range(len(json_data)):
        right_id = lib_prep_ids[lib_prep_code_ids.index(json_data[row_index][1])]

        library_prep_obj = get_lib_prep_obj_from_id(right_id)

        for p_index in range(fixed_heading_length, parameters_length):
            lib_parameter_value ={}
            lib_parameter_value['parameter_id'] = ProtocolParameters.objects.get(protocol_id__exact = form_data['protocol_id'],
                                parameterName__exact = headings[p_index])
            lib_parameter_value['library_id'] = library_prep_obj
            lib_parameter_value['parameterValue'] = json_data[row_index][p_index]

            new_parameters_data = LibParameterValue.objects.create_library_parameter_value (lib_parameter_value)

        kit_index = HEADING_FIX_FOR_ADDING_LIB_PARAMETERS.index('Lot Regents Kit used')
        library_prep_obj.set_reagent_user_kit(json_data[row_index] [kit_index])
        stored_params.append([library_prep_obj.get_sample_name(), library_prep_obj.get_lib_prep_code()])

        library_prep_obj.set_state('Completed')
        sample_obj = library_prep_obj.get_sample_obj ()
        # Update the sample state to "Create Pool"
        sample_obj.set_state('Pool Preparation')
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
    '''

    BORRAR
    '''
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
    Variables:
        sample_objs  # list of the sample objects
    Return
        error mesage or sample_objs
    '''
    sample_objs = []
    invalid_samples = []
    for sample in samples:
        s_obj = get_sample_obj_from_sample_name(sample)
        if not s_obj :
            invalid_samples.append(sample)
        else:
            sample_objs.append(s_obj)
    import pdb; pdb.set_trace()
    if len(invalid_samples) > 0:
        error = {}
        error_message = ERROR_SAMPLE_SHEET_CONTAINS_NOT_DEFINED_SAMPLES.copy()
        error_message.insert(1,' , '.join(invalid_samples))
        error['ERROR'] = error_message
        return error

    for sample_obj in sample_objs:
        # mark as invalid sample when it is not in Library preparation state
        if not 'Library preparation' == sample_obj.get_sample_state() :
            invalid_samples.append(sample_obj.get_sample_name())
        # mark as invalid when library prepation object is alreade created and it is not in Created for reused state
        #elif  LibraryPreparation.objects.filter(sample_id = sample_obj).exists() and not LibraryPreparation.objects.filter(sample_id = sample_obj ,libPrepState__libPrepState__exact = 'Created for Reuse').exists():
        elif  LibraryPreparation.objects.filter(sample_id = sample_obj , libPrepState__libPrepState__exact = 'Defined').exists():
                invalid_samples.append(sample_obj.get_sample_name())

    if len(invalid_samples) > 0:
        error = {}
        error_message = ERROR_SAMPLES_INVALID_STATE_FOR_LIBRARY_PREPARATION.copy()
        error_message.insert(1,' , '.join(invalid_samples))
        error['ERROR'] = error_message
        return error
    return sample_objs

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
        sample data objects if all checks are valid or ERROR if file is invalid
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
    # check if collection index kit is defined

    if not check_collection_index_exists (input_data['index_adapter']):
        error = {}
        error_message = ERROR_COLLECTION_INDEX_KIT_NOT_DEFINED.copy()
        error_message.insert(1, input_data['index_adapter'])
        error['ERROR'] = error_message
        delete_stored_file(input_data['full_path_file'])
        return error
    return valid_samples

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

        if not sample_row[sample_name_index] in index_values :
            index_values [sample_row[sample_name_index]] = []
        index_values [sample_row[sample_name_index]].append([indexes_in_sample])

    if len(duplicated_index_sample) == 0:
        return 'False'
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
    '''
    Description:
        The function get the library preparation information for sample.
        It return a dictionary with heading and lib_prep_data which is a list
        having index values, parameter heading, parameter values, library preparation id
        and library preparation codeID
    Input:
        sample_id  # sample id
    Constants:
        HEADING_FOR_LIBRARY_PREPARATION_DEFINITION
        HEADING_FOR_DISPLAY_POOL_INFORMATION_IN_SAMPLE_INFO
    Return
        library_information
    '''
    library_information = {}
    if LibraryPreparation.objects.filter(sample_id__pk__exact = sample_id).exists():
        library_information['library_definition_heading'] = HEADING_FOR_LIBRARY_PREPARATION_DEFINITION
        library_information['library_definition'] = []
        library_information['pool_information'] = []
        library_preparation_items = LibraryPreparation.objects.filter(sample_id__pk__exact = sample_id).exclude(libPrepState__libPrepState__exact = 'Created for Reuse')
        library_information['lib_prep_param_value'] = []
        library_information['lib_prep_data'] = []
        for library_item in library_preparation_items:

            lib_prep_data = []
            lib_prep_data.append(library_item.get_info_for_display())
            protocol_used_obj = library_item.get_protocol_obj()
            if ProtocolParameters.objects.filter(protocol_id = protocol_used_obj).exists():
                parameter_names = ProtocolParameters.objects.filter(protocol_id = protocol_used_obj).order_by('parameterOrder')
                lib_prep_param_heading = ['Lib Preparation CodeID']
                lib_prep_param_value = [library_item.get_lib_prep_code()]
                for p_name in parameter_names:
                    lib_prep_param_heading.append(p_name.get_parameter_name())
                    if LibParameterValue.objects.filter(library_id = library_item).exists():
                        lib_prep_param_value.append(LibParameterValue.objects.get(library_id = library_item, parameter_id = p_name).get_parameter_information())
                lib_prep_data.append(lib_prep_param_heading)
                lib_prep_data.append(lib_prep_param_value)
            else:

                lib_prep_data.append('')
                lib_prep_data.append('')
            lib_prep_data.append(library_item.get_id())
            lib_prep_data.append(library_item.get_lib_prep_code())
            library_information['lib_prep_data'].append(lib_prep_data)


            if library_item.pools.all().exists() :
                pools = library_item.pools.all()
                lib_prep_code_id = library_item.get_lib_prep_code()
                for pool in pools:
                    pool_name = pool.get_pool_name()
                    pool_code = pool.get_pool_code_id()
                    run_name= pool.get_run_name()
                    library_information['pool_information'].append([lib_prep_code_id, pool_name,pool_code, run_name])

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

def store_library_preparation_sample_sheet(sample_sheet_data, user) :
    '''
    Description:
        The function will get the extracted data from sample sheet.
        Then store the libraryPreparation data for each sample and update the sample state to "library Preparation"
    Input:
        sample_sheet_data   # extracted data from sample sheet in dictionary format
        user        # user object
    Return:
        new_user_s_sheet_obj .
    '''
    sample_sheet_data['user'] = user
    new_user_s_sheet_obj = libPreparationUserSampleSheet.objects.create_lib_prep_user_sample_sheet(sample_sheet_data)

    return new_user_s_sheet_obj

def get_library_code_and_unique_id (sample_obj):
    '''
    Description:
        The function will find out the latest library preparation uniqueID", increment the value
        and will return the updated value to use
    Input:
        sample_obj        # sample object
    Variables:

    Return:
        uniqueID .
    '''
    if LibraryPreparation.objects.filter(sample_id = sample_obj, libPrepState__libPrepState__exact = 'Created for Reuse').exists():
        lib_prep_obj = LibraryPreparation.objects.get(sample_id = sample_obj, libPrepState__libPrepState__exact = 'Created for Reuse')
        molecule_obj = lib_prep_obj.get_molecule_obj()
        last_lib_prep_for_molecule = LibraryPreparation.objects.filter(sample_id = sample_obj, molecule_id = molecule_obj).exclude(libPrepState__libPrepState__exact = 'Created for Reuse').last()
        if last_lib_prep_for_molecule :
            last_lib_prep_code_id = last_lib_prep_for_molecule.get_lib_prep_code()
            split_code = re.search('(.*_)(\d+)$',last_lib_prep_code_id)
            index_val = int(split_code.group(2))
            new_index = str(index_val +1).zfill(2)
            lib_prep_code_id = split_code.group(1) + new_index
            s_uniqueID = sample_obj.get_unique_sample_id()
            #unique_s_id_split = lib_prep_obj.get_unique_id().split('-')
            # count the number that library preparation was used on the same sample
            lib_prep_times = str(LibraryPreparation.objects.filter(sample_id = sample_obj).count())
            #inc_value = int(unique_s_id_split[-1]) + 1
            #unique_s_id_split[-1] = str(inc_value)
            uniqueID = s_uniqueID + '-' + lib_prep_times
        else:
            lib_prep_code_id = molecule_obj.get_molecule_code_id() + '_LIB_01'
            split_code = lib_prep_code_id.split('_')
            uniqueID = sample_obj.get_unique_sample_id() +'-' + split_code[-3][1:] + '-' + split_code[-1]
    else:
        molecule_obj = MoleculePreparation.objects.filter(sample = sample_obj).last()
        lib_prep_code_id = molecule_obj.get_molecule_code_id() + '_LIB_01'

        #split_code = lib_prep_code_id.split('_')
        #uniqueID = sample_obj.get_unique_sample_id() +'-' + split_code[-3][1:] + '-' + split_code[-1]
        uniqueID = sample_obj.get_unique_sample_id() +'-1'
    return lib_prep_code_id, uniqueID

def store_library_preparation_samples(sample_sheet_data, user, protocol , user_sample_sheet_obj):
    '''
    Description:
        The function will get the sample names, extracted data from sample sheet, index, .
        Then store the libraryPreparation data for each sample and update the sample state to "library Preparation"
    Input:
        sample_sheet_data   # extracted data from sample sheet in dictionary format
        user        # user object
        protocol    # protocol name to be used for these library preparation samples
        user_sample_sheet_obj # user_sample_sheet object for assigning to each library preparation
    Constamt:
        MAP_USER_SAMPLE_SHEET_TO_DATABASE_ONE_INDEX
        MAP_USER_SAMPLE_SHEET_TO_DATABASE_TWO_INDEX
    Functions:
        get_library_unique_id # located at this file
        get_sample_obj_from_sample_name  # located at iSkyLIMS_core/utils/handling_samples.py
        get_protocol_parameters   # located at iSkyLIMS_core/utils/handling_protocols.py
    Variables:
        stored_lib_prep     # dictionary to get data to create the library preparation object
    Return:
        stored_lib_prep .Including the lib_prep_code_id and lib_prep_id
    '''
    stored_lib_prep = {}
    lib_prep_id = []
    protocol_obj = Protocols.objects.get(name__exact = protocol)

    if 'I5_Index_ID' in  sample_sheet_data['heading'] :
        single_paired = 'Paired End'
        mapping = MAP_USER_SAMPLE_SHEET_TO_DATABASE_TWO_INDEX
    else:
        single_paired = 'Single Reads'
        mapping = MAP_USER_SAMPLE_SHEET_TO_DATABASE_ONE_INDEX
    for lib_prep_sample_data in sample_sheet_data['sample_data'] :
        stored_lib_prep = {}
        stored_lib_prep['protocol_obj'] = protocol_obj
        stored_lib_prep['single_paired'] = single_paired
        stored_lib_prep['read_length'] = sample_sheet_data['reads'][0]
        for item in mapping :
            stored_lib_prep[item[1]] =  lib_prep_sample_data[sample_sheet_data['heading'].index(item[0])]
        # if Single reads then set the index 5 to empty
        if not 'I5_Index_ID' in sample_sheet_data['heading'] :
            stored_lib_prep['i5IndexID'] = ''
            stored_lib_prep['i5Index'] = ''
        stored_lib_prep['protocol_id'] = protocol_obj
        stored_lib_prep['user_sample_sheet'] = user_sample_sheet_obj

        sample_obj = get_sample_obj_from_sample_name(stored_lib_prep['sample_name'])
        # get the latest molecule extraction to assing it by default to the library preparation
        molecule_obj = MoleculePreparation.objects.filter(sample = sample_obj).last()
        stored_lib_prep['lib_prep_code_id'], stored_lib_prep['uniqueID'] = get_library_code_and_unique_id(sample_obj)

        # update the library preparation when it was already created for reuse
        if LibraryPreparation.objects.filter(sample_id = sample_obj, libPrepState__libPrepState__exact = 'Created for Reuse').exists():
            lib_prep_obj = LibraryPreparation.objects.get(sample_id = sample_obj, libPrepState__libPrepState__exact = 'Created for Reuse')
            new_library_preparation = lib_prep_obj.update_lib_preparation_info_in_reuse_state(stored_lib_prep)
        else:
            stored_lib_prep['molecule_obj'] = molecule_obj
            new_library_preparation = LibraryPreparation.objects.create_lib_preparation(stored_lib_prep)

        lib_prep_id.append(new_library_preparation.get_id())
        #lib_prep_code_id.append(new_library_preparation.get_lib_prep_code())

        import pdb; pdb.set_trace()

    #stored_lib_prep['lib_prep_id'] = ','.join(lib_prep_id)
    #stored_lib_prep['lib_prep_code_id'] = ','.join(lib_prep_code_id)
    return lib_prep_id

def get_library_preparation_heading_for_samples (lib_prep_ids , protocol):
    '''
    Description:
        The function gets the information to display the library preparation to assing the protocol parameter values sample names, extracted data from sample sheet, index, .
        Then store the libraryPreparation data for each sample and update the sample state to "library Preparation"
    Input:
        lib_prep_ids   # Library preparation ids
        user        # user object
        protocol    # protocol name to get the protocol parameters
    Constamt:
        HEADING_FIX_FOR_ADDING_LIB_PARAMETERS

    Functions:
        get_lib_prep_obj_from_id # located at this file
        get_lot_commercial_kits  # located at iSkyLIMS_core/utils/handling_commercial_kits.py
        get_protocol_parameters   # located at iSkyLIMS_core/utils/handling_protocols.py
    Variables:
        stored_lib_prep     # dictionary to get data to create the library preparation object
    Return:
        stored_lib_prep .Including the lib_prep_code_id and lib_prep_id
    '''
    stored_lib_prep_data = {}
    stored_lib_prep_data['data'] = []
    valid_lib_prep_ids = []
    lib_prep_code_ids =  []
    user_list = []
    protocol_obj = Protocols.objects.get(name__exact = protocol)
    parameter_heading = get_protocol_parameters(protocol_obj)
    length_heading = len(HEADING_FIX_FOR_ADDING_LIB_PARAMETERS) + len (parameter_heading)

    for lib_prep in lib_prep_ids :
        lib_prep_obj = get_lib_prep_obj_from_id (lib_prep)
        if lib_prep_obj == 'None':
            continue

        lib_prep_code = lib_prep_obj.get_lib_prep_code()
        data = ['']*length_heading
        data[0] = lib_prep_obj.get_sample_name()
        data[1] = lib_prep_code

        stored_lib_prep_data['data'].append(data)
        valid_lib_prep_ids.append(lib_prep)
        lib_prep_code_ids.append(lib_prep_code)
        user_obj = lib_prep_obj.get_user_obj()
        if not user_obj in user_list:
            user_list.append(user_obj)
    # collect the reagents kits from the user in the sample_sheet

    reagents_kits = []
    for user_obj in user_list:
        reagents_kits += get_lot_commercial_kits(user_obj, protocol_obj)
    # get only unique regents Kits
    unique_reagents_kits = list(set(reagents_kits))
    stored_lib_prep_data['heading'] = HEADING_FIX_FOR_ADDING_LIB_PARAMETERS
    stored_lib_prep_data['param_heading'] = parameter_heading
    stored_lib_prep_data['lib_prep_ids'] = ','.join(valid_lib_prep_ids)
    stored_lib_prep_data['lib_prep_code_ids'] = ','.join(lib_prep_code_ids)
    stored_lib_prep_data['heading_in_excel'] = ','.join(HEADING_FIX_FOR_ADDING_LIB_PARAMETERS + parameter_heading)
    stored_lib_prep_data['protocol_id'] = protocol_obj.get_protocol_id()
    stored_lib_prep_data['reagents_kits'] = unique_reagents_kits

    return stored_lib_prep_data


def get_lib_prep_obj_from_id (library_preparation_id):
    '''
    Description:
        The function gets the library preparation id and it returns the object instance
    Input:

    Return:
        library_preparation_obj or None if not match
    '''

    if LibraryPreparation.objects.filter(pk__exact = library_preparation_id).exists():
        library_preparation_obj = LibraryPreparation.objects.get(pk__exact = library_preparation_id)
        return library_preparation_obj
    else:
        return 'None'


def update_batch_lib_prep_sample_state(lib_prep_ids,  sample_state):
    '''
    Description:
        The function set the sample state having as input the list of library preparation ids
    Input:
        lib_prep_ids        # list of library preparation ids
        sample_state        # state to be set
    Return:
        None
    '''
    for lib_id in lib_prep_ids:
        lib_obj = LibraryPreparation.objects.get(pk__exact = lib_id)
        sample_obj = lib_obj.get_sample_obj().set_state(sample_state)

    return
