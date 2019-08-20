import json, re
from iSkyLIMS_core.core_config import *
from iSkyLIMS_core.models import *
from django.contrib.auth.models import User

def analize_input_samples (request):
    '''
    Description:
        The function will return the laboratories defined in database.
    Input:
        request
    Variables:
        laboratories # list containing all laboratory names
    Return:
        laboratories.
    '''
    excel_data = request.POST['table_data']
    na_json_data = json.loads(request.POST['table_data'])
    heading_in_excel = ['patientCodeName', 'laboratory', 'labSampleName', 'extractionDate',
                    'sampleName', 'sampleType', 'species']

    sample_recorded = {}
    not_valid_samples_same_user = []
    not_valid_samples_other_user = []
    valid_samples =[]
    sample_recorded['all_samples_valid'] = True
    samples_continue = []
    for row in na_json_data :
        sample_data = {}
        sample_name = row[heading_in_excel.index('sampleName')]
        if sample_name == '' :
            continue
        if not sample_already_defined (row[heading_in_excel.index('sampleName')]):
            for i in range(len(heading_in_excel)) :
                sample_data[heading_in_excel[i]] = row[i]

            sample_data['user'] = request.user.username
            lab_code = Laboratory.objects.get(labName__exact = row[heading_in_excel.index('laboratory')] ).get_lab_code()
            sample_data['sample_id'] = str(lab_code + '_' + sample_name)
            if not Samples.objects.exclude(uniqueSampleID__isnull = True).exists():
                sample_data['new_unique_value'] = 'AAA-0001'
            else:
                last_unique_value = Samples.objects.exclude(uniqueSampleID__isnull = True).last().uniqueSampleID
                sample_data['new_unique_value'] = increase_unique_value(last_unique_value)
            new_sample = Samples.objects.create_sample(sample_data)
            valid_samples.append(new_sample.get_sample_definition_information().split(';'))
            samples_continue.append(new_sample.get_sample_id())
        else: # get the invalid sample to displays information to user
            sample_recorded['all_samples_valid'] = False
            invalid_sample = Samples.objects.get(sampleName__exact = sample_name)
            if request.user.username == invalid_sample.get_register_user() :
                not_valid_samples_same_user.append(invalid_sample.get_sample_definition_information().split(';'))
            else:
                not_valid_samples_other_user.append([invalid_sample.get_sample_name(), invalid_sample.get_registered_sample() ])

        sample_recorded['valid_samples'] = valid_samples
        sample_recorded['not_valid_samples_same_user'] = not_valid_samples_same_user
        sample_recorded['not_valid_samples_other_user'] = not_valid_samples_other_user
        sample_recorded['heading'] = ['Sample Recorded Date', 'Sample Code ID','Sample Type']
        sample_recorded['heading_same_user'] = ['Sample Recorded Date', 'Sample Code ID', 'Sample Type',
                        'New DNA', 'New Library', 'Repetition']
        sample_recorded['heading_other_user'] = [ 'Sample Name', 'Registered Sample Date']
        if sample_recorded['all_samples_valid']:
            sample_recorded['samples_to_continue'] = ','.join(samples_continue)
    #import pdb; pdb.set_trace()
    return sample_recorded

def analize_input_molecules (request):
    '''
    Description:
        The function will return the molecules defined in database.
    Input:
        request
    Variables:
        laboratories # list containing all laboratory names
    Return:
        laboratories.
    '''
    #excel_data = request.POST['table_data']
    molecule_json_data = json.loads(request.POST['nucleic_data'])
    samples = request.POST['samples'].split(',')
    heading_in_excel = ['sampleID', 'molecule_type', 'type_extraction', 'extractionDate',
                    'protocol_type', 'protocol_used']

    molecule_recorded = {}
    showed_molecule =[]
    pending_molecule =[]
    prot_used_in_display = ''
    molecule_recorded['molecule_id'] = []
    molecule_recorded['data'] = []
    for row_index in range(len(molecule_json_data)) :
        molecule_data = {}

        sample_obj = Samples.objects.get(pk = int(samples[row_index]))
        #import pdb; pdb.set_trace()
        if MoleculePreparation.objects.filter(sample = sample_obj).exists():
            last_molecule_code = MoleculePreparation.objects.filter(sample = sample_obj).last().get_molecule_code_id()
            code_split = re.search(r'(.*_E)(\d+)$', last_molecule_code)
            number_code = int(code_split.group(2))
            number_code +=1
            molecule_code_id = code_split.group(1) + str(number_code)
        else:
            number_code = 1
            molecule_code_id = sample_obj.get_sample_code() + '_E1'
        protocol_used = molecule_json_data[row_index][heading_in_excel.index('protocol_used')]
        protocol_used_obj = Protocols.objects.get(name__exact = protocol_used)
        molecule_used_obj = MoleculeType.objects.get(moleculeType__exact = molecule_json_data[row_index][heading_in_excel.index('molecule_type')])

        molecule_data['protocolUsed'] =  protocol_used_obj
        molecule_data['sample'] = sample_obj
        molecule_data['moleculeUsed'] =  molecule_used_obj
        molecule_data['moleculeCodeId'] = molecule_code_id
        molecule_data['extractionType'] =  molecule_json_data[row_index][heading_in_excel.index('type_extraction')]
        molecule_data['moleculeExtractionDate'] = molecule_json_data[row_index][heading_in_excel.index('extractionDate')]
        molecule_data['numberOfReused'] = str(number_code - 1)

        #new_molecule = MoleculePreparation.object.create_molecule(molecule_data)

        if prot_used_in_display == '':
            if ProtocolParameters.objects.filter(protocol_id__exact = protocol_used_obj).exists():
                prot_used_in_display = protocol_used
                protocol_parameters = ProtocolParameters.objects.filter(protocol_id__exact = protocol_used_obj, parameterUsed = True).order_by('parameterOrder')
                parameter_list = []
                for parameter in protocol_parameters :
                    parameter_list.append(parameter.get_parameter_name())
                length_heading = len(HEADING_FOR_MOLECULE_ADDING_PARAMETERS + parameter_list)
                molecule_recorded['fix_heading'] = HEADING_FOR_MOLECULE_ADDING_PARAMETERS
                molecule_recorded['param_heading'] = parameter_list

        if protocol_used == prot_used_in_display :
            #showed_molecule.append(new_molecule.get_id())
            data = ['']*length_heading
            data[0] = molecule_code_id
            data[1] = protocol_used
            molecule_recorded['data'].append(data)
        else:
            pending_molecule.append(new_molecule.get_id())
        #import pdb; pdb.set_trace()
    molecule_recorded['molecule_id'] = ','.join(showed_molecule)
    molecule_recorded['pending_id'] = ','.join(pending_molecule)
    return molecule_recorded

def build_record_sample_form () :
    sample_information = {}
    sample_information['species'] = get_species()
    sample_information['laboratory'] = get_laboratory()
    sample_information['sampleType'] = get_sample_type()
    return sample_information

def get_extraction_kits(username) :
    '''
    Description:
        The function will return the kits that are associated to the user.
    Input:
        username
    Variables:
        laboratories # list containing all laboratory names
    Return:
        laboratories.
    '''
    kits = 1
    return


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

def get_defined_samples (register_user):
    '''
    Description:
        The function will return the samples created by the user which are in Defined state.
    Input:
        register_user
    Variables:
        laboratories # list containing all laboratory names
    Return:
        samples_list.
    '''
    defined_samples = {}
    #sample_ids = []
    defined_samples['heading'] = HEADING_FOR_DEFINED_SAMPLES
    sample_information = []
    if Samples.objects.filter(sampleState__sampleStateName__exact = 'Defined',
                    sampleUser__username__exact = register_user).exists():
        sample_list = Samples.objects.filter(sampleState__sampleStateName__exact = 'Defined',
                        sampleUser__username__exact = register_user).order_by('generated_at')
        for sample in sample_list :
            sample_information.append(sample.get_info_in_defined_state())
        defined_samples['sample_information'] = sample_information
    return defined_samples

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

def get_molules_type ():
    '''
    Description:
        The function will return the list of type of molecules defined.
    Input:
        none
    Variables:
        molecule_type # list containing all type of molecules
    Return:
        molecule_type.
    '''
    molecule_type = []
    types = MoleculeType.objects.all()
    for type in types:
        molecule_type.append(type.get_name())
    return molecule_type

def get_molecule_protocols ():
    '''
    Description:
        The function will return the protocols defined.
    Input:
        none
    Variables:
        protocols # dictionary containing all protocols using "key" as protocol type
    Return:
        protocols.
    '''
    protocol_types = []
    protocol_list = []
    protocols = {}
    p_types = ProtocolType.objects.filter(molecule__isnull = False)
    for p_type in p_types :
        protocol_types.append(p_type.get_name())
    for protocol_type in protocol_types:
        if protocol_type not in protocols :
            protocols[protocol_type] = []
        protocols_in_type = Protocols.objects.filter(type__protocol_type__exact = protocol_type)
        for p_in_type in protocols_in_type:
            protocol_name = p_in_type.get_name()
            protocols[protocol_type].append(protocol_name)
            protocol_list.append(protocol_name)

    return protocols, protocol_list

def increase_unique_value (old_unique_number):
    '''
    Description:
        The function will increase the sample unique value.
    Input:
        old_unique_number # contains the last unique value store in Database
    Return:
        The increased number.
    '''
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


def prepare_sample_input_table ():
    '''
    Description:    The function collect the species, laboratory, type of samples, and heading
                    used in the input table. Return a dictionary with collected information.
    Input:
    Variables:
        s_information # dictionary which collects all info
    Return:
        s_information #
    '''
    # get the choices to be included in the form
    s_information = build_record_sample_form()
    s_information['heading'] = HEADING_FOR_RECORD_SAMPLES
    s_information ['table_size']= len(HEADING_FOR_RECORD_SAMPLES)
    return s_information

def prepare_molecule_input_table (samples):
    '''
    Description:    The function .
    Input:
        samples     # list of the samples to be include in the table
    Variables:
        molecule_information # dictionary which collects all info
    Return:
        molecule_information #
    '''
    molecule_information = {}
    molecule_information['headings'] = HEADING_FOR_MOLECULE_INITIAL_SETTINGS
    molecule_information['data'] = []
    valid_samples = []
    for sample in samples :
        try:
            if  Samples.objects.filter(pk__exact = int(sample)).exists():
                valid_samples.append(int(sample))
        except:
            continue
        if len(valid_samples) == 0:
            molecule_information['ERROR'] = True
            return molecule_information
    sample_code_id = []
    molecule_information['data'] =[]
    for sample in valid_samples:

        #sample_code_id.append(Samples.objects.get(pk__exact = sample).get_sample_code())
        data = ['']*len(HEADING_FOR_MOLECULE_INITIAL_SETTINGS)
        data[0] = Samples.objects.get(pk__exact = sample).get_sample_code()
        molecule_information['data'].append(data)
    molecule_information['type_of_molecules'] = get_molules_type ()
    molecule_information['protocols_dict'],molecule_information['protocol_list']  = get_molecule_protocols()
    molecule_information['number_of_samples'] = len(valid_samples)
    molecule_information['table_length']  = len(HEADING_FOR_MOLECULE_INITIAL_SETTINGS)
    return molecule_information

def sample_already_defined(sample_name):
    if Samples.objects.filter(sampleName__exact = sample_name, sampleState__sampleStateName = 'Defined').exists():
        return True
    else:
        return False
