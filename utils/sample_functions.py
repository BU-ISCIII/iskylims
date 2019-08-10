from iSkyLIMS_wetlab.models import *

def build_record_sample_form () :
    sample_information = {}
    sample_information['species'] = get_species()
    sample_information['laboratory'] = get_laboratory()
    sample_information['sampleType'] = get_sample_type()
    return sample_information


def sample_already_defined(sample_name):
    if SamplesInProject.objects.filter(sampleName__exact = sample_name, sampleState__sampleStateName = 'Defined').exists():
        return True
    else:
        return False

def get_available_lib_kit (register_user):
    lib_kits = []
    if ReagentsCommercialKits.objects.filter(registerUser__username__exact = register_user).exists():
        kits = ReagentsCommercialKits.objects.filter(registerUser__username__exact = register_user)
        for kit in kits:
            lib_kits.append(kit.get_name())

    return lib_kits

def get_data_from_excel_form (spreadsheet, n_columns):
    '''
    Description:        The function split the spreadsheet data. Every row in the spreadsheet is stored
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


def get_nucleic_accid_kits(na_type, register_user) :
    '''
    Description:
        The function will return the Nucleic Accid Kits defined in database.
    Input:
        na_type     # type of nucleic accid (DNA/RNA)
    Variables:
        laboratories # list containing all laboratory names
    Return:
    '''
    nucleic_kits = []
    if NucleotidesComercialKits.objects.filter(naType__exact = na_type, registerUser__username__exact = register_user).exists():
        kits = NucleotidesComercialKits.objects.filter(naType__exact = na_type, registerUser__username__exact = register_user)
        for kit in kits :
            nucleic_kits.append(kit.get_name()+ '_'+ kit.get_chipLot())

    return nucleic_kits

def get_full_defined_protocols_name ():
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
            if ( NAProtocolParameters.objects.filter(protocol_id = protocol).exists()
                and LibraryProtocolParameters.objects.filter(protocol_id = protocol).exists()):
                prot_names.append(protocol.get_name())
    return prot_names

def get_samples_for_library_definition (register_user):
    samples_obj = {}
    if SamplesInProject.objects.filter(sampleState__sampleStateName__exact = 'Added DNA Info',
                    registerUser__username__exact = register_user).exists():
        sample_list = SamplesInProject.objects.filter(sampleState__sampleStateName__exact = 'Added DNA Info',
                        registerUser__username__exact = register_user).order_by('sampleProtocol')
        for sample in sample_list :
            protocol = str (sample.get_sample_protocol())
            if not protocol in samples_obj:
                samples_obj[protocol] = []
            samples_obj[protocol].append(sample)

    return samples_obj

def get_samples_for_na_definition (register_user):
    samples_obj = {}
    if SamplesInProject.objects.filter(sampleState__sampleStateName__exact = 'Defined',
                    registerUser__username__exact = register_user).exists():
        sample_list = SamplesInProject.objects.filter(sampleState__sampleStateName__exact = 'Defined',
                        registerUser__username__exact = register_user).order_by('sampleProtocol')
        for sample in sample_list :
            protocol_dna = str (sample.get_sample_protocol() + '_' + sample.get_sample_nucleic_accid())
            if not protocol_dna in samples_obj:
                samples_obj[protocol_dna] = []
            samples_obj[protocol_dna].append(sample)

    return samples_obj

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


def select_samples_4_dna(valid_samples) :
    '''
    Description:
        The function will get the first sample and it fetch the value of protocol
        and the Nucleic Accid value. It returns the sample_id of the all valid samples
        that match with the same values.
    Input:
        valid_samples   # list of the valid samples
    Variables:
        protocol # protocol of the first sample
        nucleic  # nucleic accid of the first sample
    Return:
        filter_samples.
    '''
    filter_samples = []

    for sample in valid_samples :
        if sample[3] == nucleic and sample[4] == protocol :
            filter_samples.append(str(SamplesInProject.objects.get(sampleCodeID__exact = sample[1]).pk))
    return ','.join(filter_samples)
