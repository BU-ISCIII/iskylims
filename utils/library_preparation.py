import json
from django.contrib.auth.models import User
from iSkyLIMS_core.models import Samples, MoleculePreparation, Protocols
from iSkyLIMS_core.utils.handling_commercial_kits import *
from iSkyLIMS_core.utils.handling_protocols import *
from iSkyLIMS_core.utils.handling_samples import  get_sample_obj_from_sample_name, get_sample_obj_from_id
from iSkyLIMS_wetlab.models import *
from iSkyLIMS_wetlab.wetlab_config import *
from iSkyLIMS_wetlab.utils.sample_sheet_utils import *
from iSkyLIMS_wetlab.utils.collection_index_functions import check_collection_index_exists , get_list_of_collection_kits
from iSkyLIMS_wetlab.utils.handling_sequencers import *
from ..fusioncharts.fusioncharts import FusionCharts
from .stats_graphics import *
from Bio.Seq import Seq
from django.contrib.auth.models import User




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

def check_users_exists(user_list):
    '''
    Description:
        The function check if users are defined on database

    Retrun:
        all_valid if all users are defined. List of users that are not in database
    '''
    user_not_found = []
    for user in user_list :
        if User.objects.filter(username__exact = user).exists():
            continue
        user_not_found.append(user)
    if user_not_found :
        return user_not_found
    return 'all_valid'


def analyze_and_store_input_param_values(form_data):
    '''
    Description:
        The function get the user input  for the library preparation parameters and store them
        in database.
    Input:
        form_data   # User form
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
    protocol_id_obj = LibraryPreparation.objects.filter(pk__exact = lib_prep_ids[0]).last().get_protocol_obj()
    if AdditionaKitsLibraryPreparation.objects.filter(protocol_id =  protocol_id_obj).exists():
        additional_kits = True
    else:
        additional_kits = False
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

        '''
        ### Moving index library data to a dedicated form
        kit_index = HEADING_FIX_FOR_ADDING_LIB_PARAMETERS.index('Lot Regents Kit used')
        library_prep_obj.set_reagent_user_kit(json_data[row_index] [kit_index])
        '''
        stored_params.append([library_prep_obj.get_sample_name(), library_prep_obj.get_lib_prep_code()])
        if additional_kits:
            library_prep_obj.set_state('Updated parameters')
        else:
            library_prep_obj.set_state('Updated additional kits')

        #sample_obj = library_prep_obj.get_sample_obj ()
        # Update the sample state to "Create Pool"
        #sample_obj.set_state('Pool Preparation')

    return stored_params

def create_library_preparation_instance(samples_data, user):
    '''
    Description:
        The function create a new library preparation instance in database with
        user, moleculeID, sampleID, protocolID
        protocols.
    Input:
        samples_data :  Contains  moleculeID, sampleID, protocolID for each sample
    Return:
        library_preparation_objs
    '''
    library_preparation_objs = []
    for key, values in samples_data.items():
        lib_prep_data = {}
        lib_prep_data['sample_id'] = values[0]
        lib_prep_data['molecule_id'] = values[1]
        lib_prep_data['protocol_obj'] = Protocols.objects.filter(type__protocol_type__exact ='Library Preparation', name__exact = values[2]).last()
        lib_prep_data['registerUser'] = user
        lib_prep_data['user_sampleID'] = get_sample_obj_from_id(lib_prep_data['sample_id']).get_sample_code()
        lib_prep_data['lib_prep_code_id'], lib_prep_data['uniqueID'] = get_library_code_and_unique_id(lib_prep_data['sample_id'])
        library_preparation_objs.append(LibraryPreparation.objects.create_lib_preparation(lib_prep_data))

    return library_preparation_objs

def extract_protocol_library_preparation_form(form_data):
    '''
    Description:
        The function get the user input to assign sample to library preparation
        protocols.
    Input:
        form_data :     User input form
    Return:
        extraction_data
    '''
    samples_ids = form_data['samplesID'].split(',')
    samples_names = form_data['samplesNames'].split(',')
    molecules_ids = form_data['moleculesID'].split(',')
    json_data = json.loads(form_data['protocol_data'])
    extraction_data = {}

    for row_index in range(len(json_data)):
        if json_data[row_index][2] == '':
            continue
        right_index = samples_names.index(json_data[row_index][0])
        extraction_data[json_data[row_index][0]] = [samples_ids[right_index],molecules_ids[right_index], json_data[row_index][2]]
    return extraction_data


def get_protocol_parameters_for_library_preparation(library_preparation_objs):
    '''
    Description:
        The function get the protocols parameters for the list of library preparation.
        In case that library preparations do not have the same protocol, only the
        ones that matches with the protocol of the first library preparation are
        considered
    Functions:
        get_protocol_parameters_and_type # located at iSkyLIMS_core.utils.handling_protocols.py
        get_protocol_parameters  # located at iSkyLIMS_core.utils.handling_protocols.py
    Constant:
        HEADING_FIX_FOR_ADDING_LIB_PROT_PARAMETERS

    Return:
        lib_prep_same_prot_parameters
    '''
    protocol_considered = ''
    lib_prep_same_prot_parameters = {}
    lib_prep_same_prot_parameters['data'] = []
    samples_names = []
    lib_prep_ids = []
    lib_prep_code_ids = []
    for library_preparation_obj in library_preparation_objs:
        lib_prep_protocol = library_preparation_obj.get_protocol_used()
        if protocol_considered == '':
            protocol_considered = lib_prep_protocol
            # get protocol parameters
            parameters_heading = get_protocol_parameters(library_preparation_obj.get_protocol_obj())
            lib_prep_same_prot_parameters['protocol_parameters_heading_type'] = get_protocol_parameters_and_type(library_preparation_obj.get_protocol_obj())
            lib_prep_same_prot_parameters['protocol_used'] = protocol_considered
            lib_prep_same_prot_parameters['protocol_id'] = library_preparation_obj.get_protocol_id()
            lib_prep_same_prot_parameters['fix_heading'] = HEADING_FIX_FOR_ADDING_LIB_PROT_PARAMETERS
        if protocol_considered == lib_prep_protocol:
            data = ['']* (len(lib_prep_same_prot_parameters['protocol_parameters_heading_type'])+ len(HEADING_FIX_FOR_ADDING_LIB_PROT_PARAMETERS))
            sample_name = library_preparation_obj.get_sample_name()
            data[0] = sample_name
            data[1] = library_preparation_obj.get_lib_prep_code()
            samples_names.append(sample_name)
            lib_prep_ids.append(library_preparation_obj.get_lib_prep_id())
            lib_prep_code_ids.append(library_preparation_obj.get_lib_prep_code())
            lib_prep_same_prot_parameters['data'].append(data)
    lib_prep_same_prot_parameters['samples_names'] = ','.join(samples_names)
    lib_prep_same_prot_parameters['lib_prep_ids'] = ','.join(lib_prep_ids)
    lib_prep_same_prot_parameters['lib_prep_code_ids'] = ','.join(lib_prep_code_ids)
    lib_prep_same_prot_parameters['heading_in_excel'] = ','.join( HEADING_FIX_FOR_ADDING_LIB_PROT_PARAMETERS + parameters_heading)

    return lib_prep_same_prot_parameters

def get_samples_for_library_preparation():
    '''
    Description:
        The function checks if there are samples that are in library preparation state.
        samples are split according to the library preparation state in:
        - Not defined.
        - Defined. (only Protocol was defined)
        - Updated parameters
        - Updated additional Index
    Constant:
        HEADING_FOR_SAMPLES_TO_DEFINE_PROTOCOL
    Functions:
        configuration_sequencer_exists     # located at iSkyLIMS_wetlab.utils.handling_sequencers
        get_configuration_sequencers       # located at iSkyLIMS_wetlab.utils.handling_sequencers
    Return:
        samples_in_lib_prep
    '''
    samples_in_lib_prep = {}
    samples_in_lib_prep['avail_samples'] = {}
    samples_in_lib_prep['avail_samples']['data'] = []

    samples_id = []
    molecules_id = []
    samples_names = []

    if Samples.objects.filter(sampleState__sampleStateName__exact = 'Library preparation').exists():
        # data = ['']* len(HEADING_FOR_SAMPLES_TO_DEFINE_PROTOCOL)
        samples_in_lib_prep['avail_samples']['heading'] = HEADING_FOR_SAMPLES_TO_DEFINE_PROTOCOL
        samples_objs = Samples.objects.filter(sampleState__sampleStateName__exact = 'Library preparation')
        for samples_obj in samples_objs:
            if LibraryPreparation.objects.filter(sample_id = samples_obj).exists():
                library_preparation_obj = LibraryPreparation.objects.filter(sample_id = samples_obj).last()
                lib_prep_obj_state = library_preparation_obj.get_state()

                if lib_prep_obj_state == 'Defined':
                    # get the library preparations that need to add parameters
                    if not 'lib_prep_defined' in samples_in_lib_prep:
                        samples_in_lib_prep['lib_prep_defined'] = {}
                    protocol_name = library_preparation_obj.get_protocol_used()
                    if not protocol_name in samples_in_lib_prep['lib_prep_defined']:
                        samples_in_lib_prep['lib_prep_defined'][protocol_name] = []
                    lib_prep_data = []
                    lib_prep_data.append(library_preparation_obj.get_sample_name())
                    lib_prep_data.append(library_preparation_obj.get_lib_prep_code())
                    lib_prep_data.append(library_preparation_obj.get_lib_prep_id())
                    samples_in_lib_prep['lib_prep_defined'][protocol_name].append(lib_prep_data)

                elif lib_prep_obj_state == 'Updated parameters':
                    # get the library preparations that need to add kits
                    if not 'lib_prep_updated_param' in samples_in_lib_prep:
                        samples_in_lib_prep['lib_prep_updated_param'] = {}
                    protocol_name = library_preparation_obj.get_protocol_used()
                    if not protocol_name in samples_in_lib_prep['lib_prep_updated_param']:
                        samples_in_lib_prep['lib_prep_updated_param'][protocol_name] = []

                    lib_prep_param_data = []
                    lib_prep_param_data.append(library_preparation_obj.get_sample_name())
                    lib_prep_param_data.append(library_preparation_obj.get_lib_prep_code())
                    lib_prep_param_data.append(library_preparation_obj.get_lib_prep_id())
                    samples_in_lib_prep['lib_prep_updated_param'][protocol_name].append(lib_prep_param_data)
                elif lib_prep_obj_state == 'Updated additional kits':
                    samples_in_lib_prep['display_sample_sheet'] = True

            else:

                data = ['']* len(HEADING_FOR_SAMPLES_TO_DEFINE_PROTOCOL)
                sample_name = samples_obj.get_sample_name()
                data[0] = sample_name
                if MoleculePreparation.objects.filter(sample = samples_obj, state__moleculeStateName__exact = 'Completed',usedForMassiveSequencing = True).exists():
                    molecule_obj = MoleculePreparation.objects.filter(sample = samples_obj, state__moleculeStateName__exact = 'Completed',usedForMassiveSequencing = True).last()
                    data[1] = molecule_obj.get_molecule_code_id()
                    samples_in_lib_prep['avail_samples']['data'].append(data)
                    samples_id.append(samples_obj.get_sample_id())
                    samples_names.append(sample_name)
                    molecules_id.append(molecule_obj.get_molecule_id())
        # Get the information for sample sheet form
        if  'display_sample_sheet' in samples_in_lib_prep :
            if configuration_sequencer_exists():
                import pdb; pdb.set_trace()
                samples_in_lib_prep['configuration_platform'], samples_in_lib_prep['configuration_platform_option'] = get_configuration_sequencers()
            if SequencingConfiguration.objects.all().exists():
                platforms_used  = get_platform_name_of_defined_sequencers()
                seq_conf_objs = SequencingConfiguration.objects.filter(platformID__platformName__in = platforms_used).order_by('platformID')
                samples_in_lib_prep['configuration_platform'] = []
                samples_in_lib_prep['configuration_platform_option'] = []
                conf_data = {}
                for seq_conf_obj in seq_conf_objs:
                    platform_name = seq_conf_obj.get_platform_name()
                    if platform_name not in conf_data:
                        samples_in_lib_prep['configuration_platform'].append(platform_name)
                        conf_data[platform_name] = []
                    conf_data[platform_name].append(seq_conf_obj.get_configuration_name())
                for platform in samples_in_lib_prep['configuration_platform']:
                    samples_in_lib_prep['configuration_platform_option'].append([platform, conf_data[platform]])
            #import pdb; pdb.set_trace()
        samples_in_lib_prep['avail_samples']['lib_prep_protocols'] = get_protocols_for_library_preparation()
        samples_in_lib_prep['avail_samples']['samplesID'] = ','.join(samples_id)
        samples_in_lib_prep['avail_samples']['samplesNames'] = ','.join(samples_names)
        samples_in_lib_prep['avail_samples']['moleculesID'] = ','.join(molecules_id)

    else:
        samples_in_lib_prep ['no_samples'] = 'No samples'
    return samples_in_lib_prep



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
        get_userid_in_user_iem_file # located at utils/sample_sheet_utils.py
        read_user_iem_file     # located at utils/sample_sheet_utils.py
    Constant:
        ERROR_UNABLE_TO_DELETE_USER_FILE
        ERROR_SAMPLE_SHEET_DOES_NOT_HAVE_DESCRIPTION_FIELD
    Return
        data with the file name and the content of the sample sheet.
        data['Error'] if file is invalid
    '''
    data = {}

    data['full_path_file'], data['file_name'] = store_user_input_file (file_in)
    file_read = read_user_iem_file(data['full_path_file'])
    if not valid_user_iem_file(file_read):
        data['ERROR'] = ERROR_INVALID_FILE_FORMAT
    else:
        data['userid_names'] =  get_userid_in_user_iem_file(file_read)
        if  'ERROR' in data['userid_names'] :
            data['ERROR'] = ERROR_SAMPLE_SHEET_DOES_NOT_HAVE_DESCRIPTION_FIELD
        else:
            data.update(get_sample_sheet_data(file_read))
    if 'ERROR' in data:
        if not delete_stored_file(data['full_path_file']):
            data['ERROR'].append(ERROR_UNABLE_TO_DELETE_USER_FILE)
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
        # mark as invalid when library prepation object is already created and it is not in "Updated additional kits" state
        #elif  LibraryPreparation.objects.filter(sample_id = sample_obj).exists() and not LibraryPreparation.objects.filter(sample_id = sample_obj ,libPrepState__libPrepState__exact = 'Created for Reuse').exists():
        elif  not LibraryPreparation.objects.filter(sample_id = sample_obj , libPrepState__libPrepState__exact = 'Updated additional kits').exists():
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
        ERROR_SAMPLE_SHEET_DOES_NOT_HAVE_COLLECTION_INDEX
        ERROR_SAMPLE_SHEET_BOTH_INSTRUMENT_AND_INDEX_NOT_INCLUDED
        ERROR_SAMPLE_SHEET_BOTH_INSTRUMENT_AND_INDEX_NOT_INCLUDED
        ERROR_SAMPLE_SHEET_USERS_ARE_NOT_DEFINED
        ERROR_SAMPLE_SHEET_USER_IS_NOT_DEFINED
    Return
        sample data objects if all checks are valid or ERROR if file is invalid
    '''
    # check that samples are defined and in the right state
    error = {}
    # check for older IEM versions which do no have index adapter/ Instrument type in sample sheet
    if input_data['index_adapter'] == '' and not 'instrument' in input_data:
        error['ERROR'] = ERROR_SAMPLE_SHEET_BOTH_INSTRUMENT_AND_INDEX_NOT_INCLUDED
        error['detail_error'] = 'no_index_no_instrument'
    elif  input_data['index_adapter'] == '':
        error['ERROR'] = ERROR_SAMPLE_SHEET_DOES_NOT_HAVE_COLLECTION_INDEX
        error['detail_error'] = 'no_index'
    elif not 'instrument' in input_data :
        error['ERROR'] = ERROR_SAMPLE_SHEET_BOTH_INSTRUMENT_AND_INDEX_NOT_INCLUDED
        error['detail_error'] = 'no_instrument'
    else:
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
            error_message = ERROR_COLLECTION_INDEX_KIT_NOT_DEFINED.copy()
            error_message.insert(1, input_data['index_adapter'])
            error['ERROR'] = error_message

    if 'ERROR' in error :
        delete_stored_file(input_data['full_path_file'])
        return error
    else:
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


def get_data_for_library_preparation_in_defined():
    '''
    Description:
        The function get the basic data for Library preparation which are in defined state
    Return
        lib_prep_data
    '''

    lib_prep_data = []
    if LibraryPreparation.objects.filter(libPrepState__libPrepState__exact = 'Defined').exists():
        libs_preps_defined = LibraryPreparation.objects.filter(libPrepState__libPrepState__exact = 'Defined').order_by('libPrepCodeID')
        for lib_prep in libs_preps_defined :
            lib_prep_data.append(lib_prep.get_basic_data())
    return lib_prep_data


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
            lib_prep_data.append(library_item.get_molecule_code_id())
            lib_prep_data.append(library_item.get_sample_id())
            library_information['lib_prep_data'].append(lib_prep_data)

            if library_item.pools.all().exists() :
                pools = library_item.pools.all()
                lib_prep_code_id = library_item.get_lib_prep_code()
                for pool in pools:
                    pool_name = pool.get_pool_name()
                    pool_code = pool.get_pool_code_id()
                    run_name= pool.get_run_name()
                    library_information['pool_information'].append([lib_prep_code_id, pool_name,pool_code, run_name, library_item.get_id()])

        if library_information['pool_information']:
            library_information['pool_heading'] = HEADING_FOR_DISPLAY_POOL_INFORMATION_IN_SAMPLE_INFO

    return library_information

def get_lib_prep_to_add_parameters():
    '''
    Description:
        The function will return a list with samples which are needs to add library preparation parameters
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


def get_protocol_from_library_id (library_prep_id):
    '''
    Description:
        The function will return a list with samples which are needs to add library preparation parameters
    Input:
        library_prep_id # id to get the protocol name
    Return:
        protocol name or empty if library id does not exists.
    '''
    if LibraryPreparation.objects.filter(pk__exact = library_prep_id).exists():
         return LibraryPreparation.objects.get(pk__exact = library_prep_id).get_protocol_used()
    return ''

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
    states_excluded = ['Completed', 'Reused pool']
    if Samples.objects.filter(sampleState__sampleStateName__exact = 'Library preparation').exists():
        samples_obj = Samples.objects.filter(sampleState__sampleStateName__exact =  'Library preparation').order_by('sampleUser').order_by('sampleEntryDate')

        for sample in samples_obj :
            if (not LibraryPreparation.objects.filter(sample_id = sample).exists()) or (LibraryPreparation.objects.filter(sample_id = sample).exclude(libPrepState__libPrepState__in = states_excluded).exists()):
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

def find_index_sequence_collection_values_kit(sequence):
    '''
    Description:
        The function will try to find the sequence by looking on the I7 index and if not matched
        check on the I5 first on the forward and then on the reverse sequence.
    Input:
        sequence        : Sequence for searching
    Return:
        index_found and the sequence
    '''

    if CollectionIndexValues.objects.filter(i_7_seq__icontains =sequence).exists():
        return ['I7', sequence]
    if CollectionIndexValues.objects.filter(i_5_seq__icontains =sequence).exists():
        return ['I5', sequence]
    rev_sequence = str(Seq(sequence).reverse_complement())
    if CollectionIndexValues.objects.filter(i_5_seq__icontains =rev_sequence).exists():
        return ['I5', rev_sequence]
    return 'None', sequence


def store_library_preparation_index(form_data):
    '''
    Description:
        The function will fetch the indexes defined in the confirmed sample sheet
        and updated the library preparation sample with index information
    Input:
        form_data   # data included in the form
    Constant:
        ERROR_USER_SAMPLE_SHEET_NO_LONGER_EXISTS
        ERROR_LIBRARY_PREPARATION_NOT_EXISTS
    Return:
        store_result .
    '''
    store_result = {}
    unable_store_lib_prep = []
    json_data = json.loads(form_data['index_data'])
    heading = form_data['heading_excel'].split(',')
    if 'I5_Index_ID' in heading :
        single_paired = 'Paired End'
        mapping = MAP_USER_SAMPLE_SHEET_TO_DATABASE_TWO_INDEX
    else:
        single_paired = 'Single Reads'
        mapping = MAP_USER_SAMPLE_SHEET_TO_DATABASE_ONE_INDEX
    sample_name_index = heading.index('Sample_Name')
    if not libPreparationUserSampleSheet.objects.filter(pk__exact = form_data['libPrepUserSampleSheetId']).exists():
        store_result['ERROR'] = ERROR_USER_SAMPLE_SHEET_NO_LONGER_EXISTS
        return store_result
    user_sample_sheet_obj = libPreparationUserSampleSheet.objects.get(pk__exact = form_data['libPrepUserSampleSheetId'])

    for row_index in range(len(json_data)):

        lib_prep_data = {}
        sample_name = json_data[row_index][sample_name_index]

        if LibraryPreparation.objects.filter(sample_id__sampleName__exact = sample_name, libPrepState__libPrepState__exact = 'Updated additional kits').exists():
            lib_prep_obj = LibraryPreparation.objects.filter(sample_id__sampleName__exact = sample_name, libPrepState__libPrepState__exact = 'Updated additional kits').last()

            for item in mapping :

                lib_prep_data[item[1]] = json_data[row_index][heading.index(item[0])]
            for item in MAP_USER_SAMPLE_SHEET_ADDITIONAL_FIELDS_FROM_TYPE_OF_SECUENCER :
                try:
                    lib_prep_data[item[1]] =  json_data[row_index][heading.index(item[0])]
                except:
                    lib_prep_data[item[1]] = None

            # if Single reads then set the index 5 to empty
            if not 'I5_Index_ID' in heading :
                lib_prep_data['i5IndexID'] = ''
                lib_prep_data['i5Index'] = ''
            lib_prep_data['user_sample_sheet'] = user_sample_sheet_obj

            lib_prep_obj.update_library_preparation_with_indexes(lib_prep_data)
            # Update library preparation and sample state
            lib_prep_obj.set_state('Completed')
            lib_prep_obj.get_sample_obj().set_state('Pool Preparation')


        else:
            #### ERROR #####
            unable_store_lib_prep.append(sample_name)
            continue

    if len(unable_store_lib_prep) > 0 :
        store_result['ERROR'] = ERROR_LIBRARY_PREPARATION_NOT_EXISTS
        store_result['ERROR'].append(unable_store_lib_prep)
    else:
        user_sample_sheet_obj.update_confirm_used(True)
        store_result['Successful'] = True
    return store_result


def store_library_preparation_sample_sheet(sample_sheet_data, user, platform, configuration) :
    '''
    Description:
        The function will get the extracted data from sample sheet.
        Then store the libraryPreparation data for each sample and update the sample state to "library Preparation"
    Input:
        sample_sheet_data   # extracted data from sample sheet in dictionary format
        user        # user object
        platform    # platform used in the sample sheet
        configuration # configuration used in the sample sheet
    Return:
        new_user_s_sheet_obj .
    '''
    sample_sheet_data['user'] = user
    sample_sheet_data['platform'] = platform
    sample_sheet_data['configuration'] = configuration
    new_user_s_sheet_obj = libPreparationUserSampleSheet.objects.create_lib_prep_user_sample_sheet(sample_sheet_data)

    return new_user_s_sheet_obj

def get_library_code_and_unique_id (sample_id):
    '''
    Description:
        The function will find out the latest library preparation uniqueID", increment the value
        and will return the updated value to use
    Input:
        sample_id        # id of the sample
    Variables:

    Return:
        uniqueID .
    '''
    sample_obj = get_sample_obj_from_id(sample_id)
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

#############################################
# Posiblemente haya que borrarlo/modificarlo
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
        MAP_USER_SAMPLE_SHEET_TO_DATABASE_TWO_INDEX
        MAP_USER_SAMPLE_SHEET_TO_DATABASE_ONE_INDEX
        MAP_USER_SAMPLE_SHEET_ADDITIONAL_FIELDS_FROM_TYPE_OF_SECUENCER
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
        for item in MAP_USER_SAMPLE_SHEET_ADDITIONAL_FIELDS_FROM_TYPE_OF_SECUENCER :
            try:
                stored_lib_prep[item[1]] =  lib_prep_sample_data[sample_sheet_data['heading'].index(item[0])]
            except:
                stored_lib_prep[item[1]] = None

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



    #stored_lib_prep['lib_prep_id'] = ','.join(lib_prep_id)
    #stored_lib_prep['lib_prep_code_id'] = ','.join(lib_prep_code_id)
    return lib_prep_id

def get_user_for_sample_sheet():
    '''
    Descripion:
        The function collect the user_id defined in iSkyLIMS
    Return:
        user_list
    '''
    user_list = []
    user_objs = User.objects.all().order_by('username')
    for user_obj in user_objs:
        user_list.append(user_obj.username)
    return user_list

def format_sample_sheet_to_display_in_form (sample_sheet_data):
    '''
    Description:
        The function gets the information to display the library preparation to assing the protocol parameter values sample names, extracted data from sample sheet, index, .
        Then store the libraryPreparation data for each sample and update the sample state to "library Preparation"
    Input:
        lib_prep_ids   # Library preparation ids
        user        # user object
        user_sample_sheet_obj # user_sample_sheet object for assigning to each library preparation
    Constamt:
        HEADING_MAIN_DATA_SAMPLE_SHEET
        HEADING_SUMMARY_DATA_SAMPLE_SHEET
    Return:
        display_data
    '''
    display_data = {}
    main_data_values = []
    main_data_heading = HEADING_MAIN_DATA_SAMPLE_SHEET.copy()
    extract_values =['application', 'instrument', 'assay', 'index_adapter','reads', 'adapter1', 'adapter2']
    if '' == sample_sheet_data['adapter2']:
        extract_values.pop()
        main_data_heading.pop()
        display_data['adapter2'] = False
    else:
          display_data['adapter2'] = True

    display_data['sample_data'] = sample_sheet_data['sample_data']
    display_data['heading'] = sample_sheet_data['heading']
    main_values = []
    for item in extract_values:
        main_values.append(sample_sheet_data[item])
    summary_values = []
    summary_values.append(len(sample_sheet_data['samples']))
    summary_values.append(sample_sheet_data['proyects'])
    summary_values.append(sample_sheet_data['userid_names'])
    display_data['main_data'] = list(zip(main_data_heading, main_values))
    display_data['summary_data'] = list(zip(HEADING_SUMMARY_DATA_SAMPLE_SHEET, summary_values))
    display_data['heading_excel'] = ','.join(sample_sheet_data['heading'])
    if len(sample_sheet_data['userid_names']) == 0:
        display_data['no_user_defined'] = True

    return display_data



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
    '''

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


def update_library_preparation_for_reuse(sample):
    '''
    Description:
        The function step the reuse of the library preparation
    Input:
        sample_list        # list of samples for updating the reuse value

    Return:
        None
    '''
    if LibraryPreparation.objects.filter(sample_id__sampleName__exact = sample).exists():
        lib_prep_obj = LibraryPreparation.objects.filter(sample_id__sampleName__exact = sample).last()
        lib_prep_obj.set_increase_reuse()
    return
