import json, re
from iSkyLIMS_core.core_config import *
from iSkyLIMS_core.models import *
from django.contrib.auth.models import User

def display_molecule_protocol_parameters (molecules):
    '''
    Description:
        The function return the quality parameters defined for the
        selected protocol.

    Input:
        request
    Variables:

    Return:
        laboratories.
    '''

    molecule_recorded = {}
    showed_molecule =[]
    molecule_recorded['data'] = []
    pending_molecule =[]
    prot_used_in_display = ''
    for molecule in molecules :
        if not MoleculePreparation.objects.filter(pk = int(molecule)).exists():
            continue
        molecule_obj = MoleculePreparation.objects.get(pk = int(molecule))
        protocol_used = molecule_obj.get_protocol()
        protocol_used_obj = Protocols.objects.get(name__exact = protocol_used)


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
        #import pdb; pdb.set_trace()
        if protocol_used == prot_used_in_display :
            showed_molecule.append(molecule)
            data = ['']*length_heading
            data[0] = molecule_obj.get_molecule_code_id()
            data[1] = protocol_used
            molecule_recorded['data'].append(data)
        else:
            pending_molecule.append(new_molecule.get_id())

    molecule_recorded['molecule_id'] = ','.join(showed_molecule)
    molecule_recorded['pending_id'] = ','.join(pending_molecule)
    molecule_recorded['heading_in_excel'] = ','.join(parameter_list)
    #import pdb; pdb.set_trace()
    return molecule_recorded





def add_molecule_protocol_parameters(request):
    '''
    Description:
        The function will store in database the molecule parameters.
        Return the list of the molecules updated
    Input:
        request
    Variables:
        laboratories # list containing all laboratory names
    Return:
        laboratories.
    '''

    molecule_parameter_value = {}
    molecule_updated_list = []
    molecule_json_data = json.loads(request.POST['parameters_data'])
    molecules = request.POST['molecules'].split(',')
    parameter_heading = request.POST['heading_in_excel'].split(',')
    parameters_length = len(molecule_json_data[0])
    fixed_heading_length = len(HEADING_FOR_MOLECULE_ADDING_PARAMETERS)
    protocol_used_obj = Protocols.objects.get(name__exact = molecule_json_data[0][1])
    for row_index in range(len(molecule_json_data)) :
        molecule_obj = MoleculePreparation.objects.get(pk = int(molecules[row_index]))
        molecule_obj.set_state('Completed')
        molecule_updated_list.append(molecule_obj.get_molecule_code_id())

        # Update sample state
        molecule_obj.get_sample_obj().set_state('Library preparation')

        for p_index in range(fixed_heading_length, parameters_length):
            molecule_parameter_value['moleculeParameter_id'] = ProtocolParameters.objects.get(protocol_id = protocol_used_obj,
                                parameterName__exact = parameter_heading[p_index - fixed_heading_length])
            molecule_parameter_value['molecule_id'] = molecule_obj
            molecule_parameter_value['parameterValue'] = molecule_json_data[row_index] [p_index]
            new_parameters_data = MoleculeParameterValue.objects.create_molecule_parameter_value (molecule_parameter_value)


    return molecule_updated_list


def analyze_input_samples (request):
    '''
    Description:
        The function will get the user information in the form.
        For new samples, are recorded on database and included in valid_samples Variable
        For alredy defined, no action is done on them and included in not_valid_samples.
    Input:
        request
    Variables:
         # list containing all laboratory names
    Return:
        valid_samples, not_valid_samples.
    '''

    na_json_data = json.loads(request.POST['table_data'])
    heading_in_form = HEADING_FOR_RECORD_SAMPLES

    sample_recorded = {}

    valid_samples =[]
    invalid_samples = []
    invalid_samples_id =[]
    sample_recorded['all_samples_valid'] = True
    samples_continue = []

    reg_user = request.user.username
    for row in na_json_data :
        sample_data = {}
        sample_name = str(row[heading_in_form.index('Sample Name')])
        if sample_name == '' :
            continue
        if not check_if_sample_already_defined (row[heading_in_form.index('Sample Name')], reg_user):
            for i in range(len(heading_in_form)) :
                sample_data[heading_in_form[i]] = row[i]

            sample_data['user'] = reg_user

            sample_data['sample_id'] = str(reg_user + '_' + sample_name)
            if not Samples.objects.exclude(uniqueSampleID__isnull = True).exists():
                sample_data['new_unique_value'] = 'AAA-0001'
            else:
                last_unique_value = Samples.objects.exclude(uniqueSampleID__isnull = True).last().uniqueSampleID
                sample_data['new_unique_value'] = increase_unique_value(last_unique_value)
            new_sample = Samples.objects.create_sample(sample_data)
            valid_samples.append(new_sample.get_sample_definition_information())
            samples_continue.append(new_sample.get_sample_id())
            #import pdb; pdb.set_trace()
        else: # get the invalid sample to displays information to user
            sample_recorded['all_samples_valid'] = False
            sample_id = Samples.objects.get(sampleName__exact = sample_name).get_sample_id()
            if not 'sample_id_for_action' in sample_recorded:
                # get the first no valid sample to ask user for new action on the sample
                sample_recorded['sample_data_for_action'] = Samples.objects.get(sampleName__exact = sample_name).get_sample_definition_information()
                sample_recorded['sample_id_for_action'] = sample_id
                invalid_samples.append(Samples.objects.get(sampleName__exact = sample_name).get_sample_definition_information())
            else:
                invalid_samples_id.append(sample_id)
                invalid_samples.append(Samples.objects.get(sampleName__exact = sample_name).get_sample_definition_information())
    sample_recorded['valid_samples'] = valid_samples
    sample_recorded['invalid_samples'] = invalid_samples
    sample_recorded['invalid_samples_id'] = ','.join(invalid_samples_id)
    sample_recorded['heading'] = HEADING_FOR_DISPLAY_RECORDED_SAMPLES

    if sample_recorded['all_samples_valid']:
        sample_recorded['samples_to_continue'] = ','.join(samples_continue)

    return sample_recorded

def analyze_input_molecules (request):
    '''
    Description:
        The function analyze the user data to assign samples to the molecule protocol.
        Molecule is created for the sample and sample state is updated to "Extracted molecule"

    Input:
        request
    Variables:

    Return:
        molecule_recorded.
    '''
    #excel_data = request.POST['table_data']
    molecule_json_data = json.loads(request.POST['molecule_data'])
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
        if not Samples.objects.filter(pk = int(samples[row_index])).exists():
            continue
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

        new_molecule = MoleculePreparation.objects.create_molecule(molecule_data)
        # Update Sample state to "Extracted molecule"
        sample_obj.set_state('Extracted Molecule')
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
            showed_molecule.append(new_molecule.get_id())
            data = ['']*length_heading
            data[0] = molecule_code_id
            data[1] = protocol_used
            molecule_recorded['data'].append(data)
        else:
            pending_molecule.append(new_molecule.get_id())
        #import pdb; pdb.set_trace()
    molecule_recorded['molecule_id'] = ','.join(showed_molecule)
    molecule_recorded['pending_id'] = ','.join(pending_molecule)
    molecule_recorded['heading_in_excel'] = ','.join(parameter_list)
    import pdb; pdb.set_trace()
    return molecule_recorded

def build_record_sample_form () :
    sample_information = {}
    sample_information['species'] = get_species()
    sample_information['laboratory'] = get_laboratory()
    sample_information['sampleType'] = get_sample_type()
    return sample_information

def check_if_sample_already_defined (sample_name,reg_user):
    if Samples.objects.filter(sampleName__exact = sample_name, sampleUser__username__exact = reg_user).exists():
        return True
    else:
        return False

def get_all_sample_information (sample_id ):
    sample_information = {}
    parameter_heading_values = []
    if not Samples.objects.filter(pk__exact = sample_id).exists():
        return 'Error'
    sample_obj = Samples.objects.get(pk__exact = sample_id)
    sample_information['sample_definition'] = sample_obj.get_info_for_display()
    sample_information['sample_definition_heading'] = HEADING_FOR_SAMPLE_DEFINITION
    # check if molecule information exists for the sample
    if MoleculePreparation.objects.filter(sample = sample_obj).exists():
        sample_information['molecule_definition_heading'] = HEADING_FOR_MOLECULE_DEFINITION
        molecules = MoleculePreparation.objects.filter(sample = sample_obj)
        sample_information['molecule_definition'] = []
        sample_information['molecule_parameter_values'] = []
        for molecule in molecules:
            sample_information['molecule_definition'].append(molecule.get_info_for_display())
            protocol_used_obj = molecule.get_protocol_obj()
            if ProtocolParameters.objects.filter(protocol_id = protocol_used_obj).exists():
                parameter_names = ProtocolParameters.objects.filter(protocol_id = protocol_used_obj).order_by('parameterOrder')
                molecule_param_heading = ['Molecule CodeID']
                mol_param_value = [molecule.get_molecule_code_id()]
                for p_name in parameter_names:
                    molecule_param_heading.append(p_name.get_parameter_name())
                    if MoleculeParameterValue.objects.filter(molecule_id = molecule).exists():
                        mol_param_value.append(MoleculeParameterValue.objects.get(molecule_id = molecule, moleculeParameter_id = p_name).get_param_value())
                parameter_heading_values.append([molecule_param_heading, mol_param_value ])
                #sample_information['molecule_parameter_values'].append(mol_param_value)
                #sample_information['molecule_parameter_heading'] = molecule_param_heading

        sample_information['parameter_heading_values'] = parameter_heading_values
    return sample_information

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
    defined_samples['heading'] = HEADING_FOR_DEFINED_SAMPLES_STATE
    sample_information = []
    if Samples.objects.filter(sampleState__sampleStateName__exact = 'Defined',
                    sampleUser__username__exact = register_user).exists():
        sample_list = Samples.objects.filter(sampleState__sampleStateName__exact = 'Defined',
                        sampleUser__username__exact = register_user).order_by('generated_at')
        for sample in sample_list :
            sample_information.append(sample.get_info_in_defined_state())
        defined_samples['sample_information'] = sample_information
    return defined_samples

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

def get_molecule_codeid_from_object(molecule_obj):
    return molecule_obj.get_molecule_code_id()

def get_molecule_obj_from_sample(sample_obj):
    if MoleculePreparation.objects.filter(sample = sample_obj).exists():
        molecules_obj = MoleculePreparation.objects.filter(sample = sample_obj)
        return molecules_obj
    else:
        return ''


##### For each state get samples
def get_samples_in_defined_state ():
    '''
    Description:
        The function will return a list with samples which are in defined state.
    Input:
        state  # string of the state to be matched
    Variables:
        sample_type_names # list containing all sample types names
    Return:
        sample_type_names.
    '''
    sample_information = []
    samples_in_state = {}
    if Samples.objects.filter(sampleState__sampleStateName__exact = 'Defined').exists():
        samples_obj = Samples.objects.filter(sampleState__sampleStateName__exact = 'Defined')
        for sample_obj in samples_obj:
            sample_information.append(sample_obj.get_info_in_defined_state())
        samples_in_state['sample_information'] = sample_information
        samples_in_state['sample_heading'] = HEADING_FOR_DEFINED_SAMPLES_STATE
        return samples_in_state

    else:
        return ''

def get_samples_in_extracted_molecule_state ():
    '''
    Description:
        The function will return a list with samples which are in extracted_molecule state,
        excluding the ones that molecule state are Completed.
    Input:

    Variables:
        molecule_state # Dictionnary with the heading and the molecule information
    Return:
        molecule_state.
    '''

    molecule_state = {}
    molecule_information = []
    if Samples.objects.filter(sampleState__sampleStateName__exact = 'Extract molecule').exists():
        samples_obj = Samples.objects.filter(sampleState__sampleStateName__exact =  'Extract molecule')
        for sample_obj in samples_obj:
            molecules = MoleculePreparation.objects.filter(sample = sample_obj)

            for molecule in molecules:
                if molecule.get_state() == 'Completed':
                    continue
                sample_information = []
                sample_information.append(sample_obj.get_extraction_date())
                sample_information.append(sample_obj.get_sample_name())
                molecule_data = molecule.get_molecule_information()
                molecule_information.append(sample_information + molecule_data)

        molecule_state['molecule_information'] = molecule_information
        molecule_state['molecule_heading'] = HEADING_FOR_EXTRACTED_MOLECULES_STATE
        return molecule_state

    else:
        return ''

##### End of getting samples by state

def get_sample_instance(sample_id, register_user):
    if Samples.objects.filter(sampleName__exact = sample_id, sampleUser__username__exact = register_user).exists():
        sample_obj = Samples.objects.get(sampleName__exact = sample_id, sampleUser__username__exact = register_user)
        return sample_obj
    return

def get_sample_obj_from_id(sample_id):
    sample_obj = Samples.objects.get(pk__exact = sample_id)
    return sample_obj

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

def get_sample_states ():
    '''
    Description:
        The function will return the sample states defined in database.
    Return:
        sample_states.
    '''
    sample_states = []
    if StatesForSample.objects.all().exists():
        states = StatesForSample.objects.all()
        for state in states:
            sample_states.append(state.get_sample_state())
    return sample_states

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

def record_molecules (request):
    '''
    Description:    The function store in database the new molecule and molecule_updated_list
                    the sample state to Extracted molecule.
    Input:
        request
    Variables:
        molecule_information # dictionary which collects all info
    Return:
        molecules_recorded with the list of the recorded molecules and the heading to
        display them
    '''
    molecule_json_data = json.loads(request.POST['molecule_data'])
    samples = request.POST['samples'].split(',')
    molecules_recorded = {}
    molecules_ids = []
    molecule_list = []
    heading_in_excel = ['sampleID', 'molecule_type', 'type_extraction', 'extractionDate',
                    'protocol_type', 'protocol_used']
    for row_index in range(len(molecule_json_data)) :
        molecule_data = {}
        if not Samples.objects.filter(pk = int(samples[row_index])).exists():
            continue
        sample_obj = Samples.objects.get(pk = int(samples[row_index]))

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
        #molecule_data['numberOfReused'] = str(number_code - 1)

        new_molecule = MoleculePreparation.objects.create_molecule(molecule_data)
        molecule_list.append([molecule_code_id, protocol_used])
        # Update Sample state to "Extracted molecule"
        sample_obj.set_state('Extract molecule')
        # Include index key to allow adding quality parameter data
        molecules_ids.append(str(new_molecule.pk))

    molecules_recorded['heading'] = HEADING_CONFIRM_MOLECULE_RECORDED
    molecules_recorded['molecule_list'] = molecule_list
    molecules_recorded['molecules'] = ','.join(molecules_ids)
    return molecules_recorded

def get_info_for_reprocess_samples(sample_ids, sample_in_action):
    sample_recorded = {}
    invalid_samples = []

    for sample_id in sample_ids :
        if sample_id == sample_in_action :
            sample_recorded['sample_data_for_action'] = Samples.objects.get(pk__exact = sample_id).get_sample_definition_information()
        invalid_samples.append(Samples.objects.get(pk__exact = sample_id).get_sample_definition_information())

    sample_recorded['invalid_samples'] = invalid_samples
    # sample_recorded['invalid_samples_id'] = ','.join(sample_ids)
    sample_recorded['heading'] = HEADING_FOR_DISPLAY_RECORDED_SAMPLES

    return sample_recorded

def get_table_record_molecule (samples):
    '''
    Description:    The function get the sampleID to create the molecule table.
    Input:
        samples     # list of the samples to be include in the table
    Variables:
        molecule_information # dictionary which collects all info
    Return:
        molecule_information #
    '''
    molecule_information = {}
    molecule_information['headings'] = HEADING_FOR_MOLECULE_PROTOCOL_DEFINITION
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
        data = ['']*len(HEADING_FOR_MOLECULE_PROTOCOL_DEFINITION)
        data[0] = Samples.objects.get(pk__exact = sample).get_sample_code()
        molecule_information['data'].append(data)
    molecule_information['type_of_molecules'] = get_molules_type ()
    molecule_information['protocols_dict'],molecule_information['protocol_list']  = get_molecule_protocols()
    molecule_information['number_of_samples'] = len(valid_samples)
    molecule_information['table_length']  = len(HEADING_FOR_MOLECULE_PROTOCOL_DEFINITION)
    return molecule_information



def search_samples(sample_name, user_name, sample_state, start_date, end_date ):
    sample_list = []

    if Samples.objects.all().exists():
        sample_founds = Samples.objects.all()
    else:
        return sample_list
    if user_name != '':
        user_name_obj = User.objects.get(username__exact = user_name)
        if not sample_founds.objects.filter(sampleUser = user_name_obj).exists():
            return sample_list
        else :
            sample_founds = sample_founds.objects.filter(sampleUser = user_name_obj)
    if sample_name != '' :

        if sample_founds.filter(sampleName__exact = sample_name).exists():
            sample_founds = sample_founds.filter(sampleName__exact = sample_name)
            if len(sample_founds) == 1:
                sample_list.append(sample_founds[0].pk)
                return sample_list

        elif sample_founds.filter(sampleName__icontains = sample_name).exists():
            sample_founds = sample_founds.filter(sampleName__icontains = sample_name)
        else:
            return sample_list


    for sample in sample_founds :

        sample_list.append(sample.get_info_for_searching())
    return sample_list

def update_molecule_reused(sample_id, molecule_code_id):
    sample_obj = Samples.objects.get(pk__exact = sample_id)
    try:
        molecule_obj = MoleculePreparation.objects.get(sample = sample_obj, moleculeCodeId__exact = molecule_code_id)
    except:
        return 'Error'
    molecule_obj.set_increase_reuse()
    sample_obj.set_state('Library Preparation')
    return molecule_obj

def update_sample_reused(reprocess_id):
    sample_obj = Samples.objects.get(pk__exact = reprocess_id)
    sample_obj.set_increase_reuse()
    sample_obj.set_state('Defined')

    return sample_obj
