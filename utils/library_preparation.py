from iSkyLIMS_core.models import Samples, MoleculePreparation
from iSkyLIMS_wetlab.models import *
from iSkyLIMS_wetlab.wetlab_config import *

def extract_sample_data (s_data,user_sample_sheet_obj):
    headings = s_data['headings']
    columns = ['Sample_ID','Sample_Name','Sample_Plate','Sample_Well','Index_Plate_Well','I7_Index_ID','index','I5_Index_ID','index2','Sample_Project']
    for sample_row in s_data['samples']:
        '''
        if Samples.objects.filter(sampleName__exact = sample_row.index('Sample_Name'), sampleState__sampleStateName = 'Add Library Preparation' ).exists():
            sample_obj = Samples.objects.filter(sampleName__exact = sample_row.index('Sample_Name'), sampleState__sampleStateName = 'Add Library Preparation')
        '''
        lib_prep_data = {}
        for column in wetlab_config.MAP_USER_SAMPLE_SHEET_TO_DATABASE[0] :
            import pdb; pdb.set_trace()
            if column in headings:
                lib_prep_data[column[1]] = sample_row.index(column)
            else:
                lib_prep_data[column] = ''




    return


def get_protocol_lib ():
    protocol_list = []
    if ProtocolLibrary.objects.all().exists():
        protocols = ProtocolLibrary.objects.all()
        for protocol in protocols:
            protocol_list.append(protocol.get_name())
    return protocol_list

def get_samples_in_add_library_preparation_state ():
    '''sampleName__exact =
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
