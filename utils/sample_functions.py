from iSkyLIMS_wetlab.models import *

def sample_already_defined(sample_name):
    if SamplesInProject.objects.filter(sampleName__exact = sample_name, sampleState__sampleStateName = 'Defined').exists():
        return True
    else:
        return False

def get_data_from_excel_form (spreadsheet, n_columns):
    '''
    Description:
        The function split the spreadsheet data. Every row in the spreadsheet is stored
        in a list. Return a array with evey item list.
    Input:
        spreadsheet     # contains all data
        n_columns       # number of columns used to
    Variables:
        spread_data # array containing row data list
    Return:
        d_list # data splited in rows
    '''
    split_data = spreadsheet.split(',')
    d_list = []
    for i in range(0,len(split_data),n_columns):
        d_list.append(split_data[i:i+n_columns])
    return d_list

def get_laboratory ():
    '''
    Description:
        The function will return the laboratories defined in database.
    Input:
        none
    Variables:
        laboratories # list containing all laboratory names
    Return:
        laboratories.
    '''
    laboratories = []
    if Laboratory.objects.filter().exists():
        laboratory_names = Laboratory.objects.all()

        for laboratory in laboratory_names:
            laboratories.append(laboratory.get_name())
    return laboratories

def get_protocols_name ():
    '''
    Description:
        The function will return the name of the protols defined in database.
    Input:
        none
    Variables:
        prot_names # list containing all protocol names
    Return:
        prot_names.
    '''
    prot_names = []
    if ProtocolInLab.objects.filter().exists():
        protocols = ProtocolInLab.objects.all()

        for protocol in protocols:
            prot_names.append(protocol.get_name())
    return prot_names

def get_sample_type ():
    '''
    Description:
        The function will return the type of samples defined in database.
    Input:
        none
    Variables:
        sample_type_names # list containing all sample types names
    Return:
        sample_type_names.
    '''
    sample_type_names = []
    if SampleType.objects.filter().exists():
        sample_types = SampleType.objects.all()

        for sample in sample_types:
            sample_type_names.append(sample.get_name())
    return sample_type_names

def get_species ():
    '''
    Description:
        The function will return the name of the species defined in database.
    Input:
        none
    Variables:
        species_names # list containing all species names
    Return:
        species_names.
    '''
    species_names = []
    if Species.objects.filter().exists():
        all_species = Species.objects.all()

        for species in all_species:
            species_names.append(species.get_name())
    return species_names

def increase_unique_value (old_unique_number):

    split_value = old_unique_number.split('-')
    number = int(split_value[1]) +1
    letter = split_value[0]

    if number > 9999 :
        number = 1
        index_letter = list(split_value[0])
        if index_letter[2] == 'Z':
            if index_letter[1] == 'Z':
                index_letter[0] = chr(ord(index_letter[0])+1)
                index_letter[1] = 'A'
                index_letter[2] = 'A'
            else:
                index_letter[1] = chr(ord(index_letter[1])+1)
                split_index_letter[2] = 'A'

            index_letter = ''.join(split_index_letter)
        else:
            index_letter[2] = chr(ord(index_letter[2])+1)

        letter = ''.join(index_letter)

    number_str = str(number)
    number_str = number_str.zfill(4)
    return str(letter + '-' + number_str)
