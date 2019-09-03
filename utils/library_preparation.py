from iSkyLIMS_core.models import Samples, MoleculePreparation, Protocols
from iSkyLIMS_wetlab.models import *
from iSkyLIMS_wetlab.wetlab_config import *

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
    import pdb; pdb.set_trace()


    return sample_list


def get_protocol_lib ():
    protocol_list = []
    if Protocols.objects.filter(type__protocol_type__exact ='Library Preparation').exists():
        protocols = Protocols.objects.filter(type__protocol_type__exact = 'Library Preparation')
        for protocol in protocols:
            protocol_list.append(protocol.get_name())
    return protocol_list


def get_samples_add_lib_prep_parameters():
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
    if libraryPreparation.objects.filter(state__exact = 'Recorded').exists():
        samples = libraryPreparation.objects.filter(state__exact = 'Recorded')
        sample_info = []
        for sample in samples:
            lib_prep_info = []
            lib_prep_info.append(sample.get_lib_prep_code())
            lib_prep_info.append(sample.get_protocol_used())
            lib_prep_info.append(sample.get_id())
            sample_info.append(lib_prep_info)
        lib_prep_parameters['sample_info'] = sample_info
        lib_prep_parameters['heading'] = HEADING_FOR_ADD_LIBRARY_PREPARATION_PARAMETERS
    return lib_prep_parameters


def get_samples_in_add_library_preparation_state ():
    '''
    Description:
        The function will return a list with samples which are in add_library_preparation state.
    Input:

    Variables:
        library_prep_information # Dictionary with the heading and the molecule information
    Return:
        library_prep_information.
    '''
    library_prep_information = {}
    molecule_information = []
    if Samples.objects.filter(sampleState__sampleStateName__exact = 'add Library preparation').exists():
        samples_obj = Samples.objects.filter(sampleState__sampleStateName__exact =  'add Library preparation')
        for sample_obj in samples_obj:
            molecules = MoleculePreparation.objects.filter(sample = sample_obj)

            for molecule in molecules:
                m_information = []

                if molecule.get_state() == 'Completed':
                    m_information.append(molecule.get_molecule_code_id())
                    m_information.append(molecule.get_protocol())
                    m_information.append(molecule.get_extraction_date())
                    m_information.append(molecule.get_id())
                    molecule_information.append(m_information)
        library_prep_information['lib_protocols'] = get_protocol_lib()
        library_prep_information['molecule_information'] = molecule_information
        library_prep_information['molecule_heading'] = HEADING_FOR_ADD_LIBRARY_PREPARATION

    return library_prep_information

def get_user_reagents_kits(user_obj):
    user_reagents_kits = []
    if ReagentsUserCommercialKits.objects.filter(registerUser = user_obj).exists():
        available_kits =  ReagentsUserCommercialKits.objects.filter(registerUser = user_obj)
        for kit in available_kits:
            user_reagents_kits.append(kit.get_nick_name())
    return user_reagents_kits
