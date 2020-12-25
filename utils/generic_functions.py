from django.conf import settings
from django.contrib.auth.models import User
from iSkyLIMS_core.models import *
from iSkyLIMS_core.core_config  import *

def get_installed_apps () :
    app_prefix = 'iSkyLIMS'
    core = 'core'
    apps_list = []
    #import pdb; pdb.set_trace()
    apps = list(apps for apps in settings.INSTALLED_APPS if (app_prefix in apps and not core in apps))
    for app in apps :
        apps_list.append([app,settings.APPS_NAMES[app]])
    return apps_list

def get_friend_list(user_name):
    friend_list = []
    user_groups = user_name.groups.values_list('name',flat=True)
    if len (user_groups) > 0 :
        for user in user_groups :

            if User.objects.filter(username__exact = user).exists():
                # friend_list.append(User.objects.get(username__exact = user).id)
                friend_list.append(User.objects.get(username__exact = user))

    friend_list.append(user_name)
    return friend_list

def get_inital_sample_settings_values(apps_name):
    '''
    Description:
        The function get the defined values for species, origin of the sample,
        type of molecules, and type of protocol.
    Input:
        apps_name    # application which the data belongs to
    Return:
        initial_data
    '''
    initial_data = {}
    if Species.objects.filter(apps_name__exact = apps_name).exists():
        species_objs = Species.objects.filter(apps_name__exact = apps_name)
        initial_data['species_data'] = []
        for species_obj in species_objs:
            initial_data['species_data'].append([species_obj.get_name(), species_obj.get_id()])
    if SamplesOrigin.objects.filter(apps_name__exact = apps_name).exists():
        samples_origin_objs = SamplesOrigin.objects.filter(apps_name__exact = apps_name)
        initial_data['samples_origin_data'] = []
        for samples_origin_obj in samples_origin_objs:
            initial_data['samples_origin_data'].append([samples_origin_obj.get_name(),samples_origin_obj.get_sample_origin_code(), samples_origin_obj.get_id()])
    if MoleculeType.objects.filter(apps_name__exact = apps_name).exists():
        initial_data['molecule_type_data'] = []
        molecule_type_objs = MoleculeType.objects.filter(apps_name__exact = apps_name)
        for molecule_type_obj in molecule_type_objs:
            initial_data['molecule_type_data'].append([molecule_type_obj.get_name() ,molecule_type_obj.get_id()])
    if ProtocolType.objects.filter(apps_name__exact = apps_name).exists():
        protocol_type_objs = ProtocolType.objects.filter(apps_name__exact = apps_name).order_by('molecule')
        initial_data['protocol_type_data'] = []
        for protocol_type_obj in protocol_type_objs:
            initial_data['protocol_type_data'].append([protocol_type_obj.get_name(), protocol_type_obj.get_id(), protocol_type_obj.get_molecule_type()])

    return initial_data

def save_inital_sample_setting_value (apps_name, data):
    '''
    Description:
        The function get the defined values for species, origin of the sample,
        type of molecules, and type of protocol.
    Input:
        apps_name    # application which the data belongs to
        data        # information to store on database
    Constants:
        ERROR_SPECIES_ALREADY_DEFINED
        ERROR_SAMPLES_ORIGIN_ALREADY_DEFINED
        ERROR_MOLECULE_TYPE_ALREADY_DEFINED
        ERROR_PROTOCOL_TYPE_ALREADY_DEFINED
    Return:
        setting_defined
    '''
    setting_defined = {}
    if 'species' in data:
        species_data = {}
        species_data['apps_name'] = apps_name
        species_data['name'] = data['species']
        if Species.objects.filter(speciesName__iexact = species_data['name'],  apps_name__exact = species_data['apps_name']).exists():
            setting_defined['ERROR'] = [ERROR_SPECIES_ALREADY_DEFINED, species_data['name']]
            return setting_defined

        new_specie_obj = Species.objects.create_new_specie(species_data)
        setting_defined['settings'] = 'Species'
        setting_defined['value'] = species_data['name']

    if 'samples_origin' in data:
        samples_origin_data = {}
        samples_origin_data['apps_name'] = apps_name
        samples_origin_data['originName'] = data['samples_origin'][0]
        samples_origin_data['originNameCoding'] = data['samples_origin'][1]
        samples_origin_data['location'] = data['samples_origin'][2]
        if SamplesOrigin.objects.filter(originNameCoding__iexact = samples_origin_data['originNameCoding'],  apps_name__exact = samples_origin_data['apps_name']).exists():
            setting_defined['ERROR'] = [ERROR_SAMPLES_ORIGIN_ALREADY_DEFINED, samples_origin_data['originNameCoding']]
            return setting_defined
        new_sample_origin_obj = SamplesOrigin.objects.create_samples_origin(samples_origin_data)
        setting_defined['settings'] = 'Samples Origin'
        setting_defined['value'] = samples_origin_data['originName']
    if 'molecule_type' in data:
        molecule_type_data = {}
        molecule_type_data['apps_name'] = apps_name
        molecule_type_data['moleculeType'] = data['molecule_type']
        if MoleculeType.objects.filter(moleculeType__iexact = molecule_type_data['moleculeType'],  apps_name__exact = molecule_type_data['apps_name']).exists():
            setting_defined['ERROR'] = [ERROR_MOLECULE_TYPE_ALREADY_DEFINED, molecule_type_data['moleculeType']]
            return setting_defined
        new_molecule_type_obj = MoleculeType.objects.create_molecule_type(molecule_type_data)
        setting_defined['settings'] = 'Molecule Type'
        setting_defined['value'] = molecule_type_data['moleculeType']

    if 'protocol_type' in data:
        protocol_type_data = {}
        protocol_type_data['apps_name'] = apps_name
        protocol_type_data['protocol_type'] = data['protocol_type'][0]
        if data['protocol_type'][1] == 'None':
            protocol_type_data['molecule'] = None
        else:
            protocol_type_data['molecule'] = data['protocol_type'][1]
        if ProtocolType.objects.filter(protocol_type__iexact = protocol_type_data['protocol_type'],  molecule__moleculeType__iexact = protocol_type_data['molecule'], apps_name__exact = protocol_type_data['apps_name']).exists():
            setting_defined['ERROR'] = [ERROR_PROTOCOL_TYPE_ALREADY_DEFINED, protocol_type_data['protocol_type']]
            return setting_defined
        new_protocol_type_obj = ProtocolType.objects.create_protocol_type(protocol_type_data)
        setting_defined['settings'] = 'Protocol Type'
        setting_defined['value'] = protocol_type_data['protocol_type']
    return setting_defined
