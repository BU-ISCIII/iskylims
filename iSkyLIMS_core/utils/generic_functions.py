import smtplib
from django.core.mail import send_mail
from django.conf import settings
from django.contrib.auth.models import User
from iSkyLIMS_core.models import *
from iSkyLIMS_core.core_config  import *

def get_installed_apps () :

    return settings.APPS_NAMES

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
    if Species.objects.filter(apps_name__exact=apps_name).exists():
        species_objs = Species.objects.filter(apps_name__exact=apps_name)
        initial_data['species_data'] = []
        for species_obj in species_objs:
            initial_data['species_data'].append([species_obj.get_name(), species_obj.get_id()])
    if LabRequest.objects.filter(apps_name__exact=apps_name).exists():
        lab_request_objs = LabRequest.objects.filter(apps_name__exact=apps_name)
        initial_data['lab_request_data'] = []
        for lab_request_obj in lab_request_objs:
            data = lab_request_obj.get_all_data()
            data.append(lab_request_obj.get_id())
            initial_data['lab_request_data'].append(data)
    if MoleculeType.objects.filter(apps_name__exact=apps_name).exists():
        initial_data['molecule_type_data'] = []
        molecule_type_objs = MoleculeType.objects.filter(apps_name__exact=apps_name)
        for molecule_type_obj in molecule_type_objs:
            initial_data['molecule_type_data'].append([molecule_type_obj.get_name(), molecule_type_obj.get_id()])
    if ProtocolType.objects.filter(apps_name__exact = apps_name).exists():
        protocol_type_objs = ProtocolType.objects.filter(apps_name__exact=apps_name).order_by('molecule')
        initial_data['protocol_type_data'] = []
        for protocol_type_obj in protocol_type_objs:
            initial_data['protocol_type_data'].append([protocol_type_obj.get_name(), protocol_type_obj.get_id(), protocol_type_obj.get_molecule_type()])
    if StateInCountry.objects.filter(apps_name__exact=apps_name).exists():
        state_objs = StateInCountry.objects.filter(apps_name__exact=apps_name).order_by('stateName')
        initial_data['states_data'] = []
        for state_obj in state_objs:
            initial_data['states_data'].append([state_obj.get_state_name(), state_obj.get_state_id()])
    if City.objects.filter(apps_name__exact=apps_name).exists():
        city_objs = City.objects.filter(apps_name__exact=apps_name).order_by('cityName')
        initial_data['cities_data'] = []
        for city_obj in city_objs:
            initial_data['cities_data'].append([city_obj.get_city_name(), city_obj.get_city_id()])


    return initial_data


def get_email_data():
    '''
    Description:
        Fetch the email configuration file
    Constant:
        EMAIL_HOST
        EMAIL_PORT
        EMAIL_HOST_USER
        EMAIL_HOST_PASSWORD
        EMAIL_USE_TLS
    Return:
        email_data
    '''
    email_data = {}
    email_data['EMAIL_HOST'] = settings.EMAIL_HOST
    email_data['EMAIL_PORT'] = settings.EMAIL_PORT
    email_data['USER_EMAIL'] = settings.EMAIL_HOST_USER
    email_data['USER_PASSWORD'] = settings. EMAIL_HOST_PASSWORD
    email_data['USE_TLS'] = settings.EMAIL_USE_TLS
    return email_data

def send_test_email(form_data):
    '''
    Description:
        Get the configuration data from the user form and send a test email
    Constant:
        EMAIL_HOST
        EMAIL_PORT
        EMAIL_HOST_USER
        EMAIL_HOST_PASSWORD
        EMAIL_USE_TLS
    Return:
        email_data
    '''
    settings.EMAIL_HOST = form_data['EMAIL_HOST']
    settings.EMAIL_PORT = form_data['EMAIL_PORT']
    settings.EMAIL_HOST_USER = form_data['USER_EMAIL']
    settings.EMAIL_HOST_PASSWORD = form_data['USER_PASSWORD']
    settings.EMAIL_USE_TLS = True if form_data['USE_TLS'] == 'True' else False

    from_user = form_data['EMAIL_ISKYLIMS']
    to_users = [form_data['test_email']]
    subject = 'testing email from iSlyLIMS'
    body_message = 'This is a email test to verify iSkyLIMS'
    try:
        send_mail (subject, body_message, from_user, to_users)
        return 'OK'
    except smtplib.SMTPException as e:
        return str(e)
    except Exception as e:
        return str(e)


def save_inital_sample_setting_value(apps_name, data):
    '''
    Description:
        The function get the defined values for species, origin of the sample,
        type of molecules, and type of protocol.
    Input:
        apps_name    # application which the data belongs to
        data        # information to store on database
    Constants:
        ERROR_SPECIES_ALREADY_DEFINED
        ERROR_LABORATORY_REQUEST_ALREADY_DEFINED
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

    if 'lab_request' in data:
        lab_request_data = {}
        lab_request_data['apps_name'] = apps_name
        lab_request_data['labName'] = data['lab_request']['labRequestName']
        lab_request_data['labNameCoding'] = data['lab_request']['labRequesCoding']
        lab_request_data['labUnit'] = data['lab_request']['department']
        lab_request_data['labContactName'] = data['lab_request']['contact']
        lab_request_data['labPhone'] = data['lab_request']['phone']
        lab_request_data['labEmail'] = data['lab_request']['email']
        lab_request_data['address'] = data['lab_request']['address']
        lab_request_data['city'] = data['lab_request']['city']
        if LabRequest.objects.filter(labNameCoding__iexact = lab_request_data['labNameCoding'],  apps_name__exact = lab_request_data['apps_name']).exists():
            setting_defined['ERROR'] = [ERROR_LABORATORY_REQUEST_ALREADY_DEFINED, lab_request_data['labNameCoding']]
            return setting_defined
        new_sample_origin_obj = LabRequest.objects.create_lab_request(lab_request_data)
        setting_defined['settings'] = 'Lab Request'
        setting_defined['value'] = lab_request_data['labName']
    if 'molecule_type' in data:
        molecule_type_data = {}
        molecule_type_data['apps_name'] = apps_name
        molecule_type_data['moleculeType'] = data['molecule_type']
        if MoleculeType.objects.filter(moleculeType__iexact = molecule_type_data['moleculeType'], apps_name__exact = molecule_type_data['apps_name']).exists():
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

    if 'state' in data:
        state_data = {}
        if StateInCountry.objects.filter(stateName__iexact=data['state'], apps_name__exact=apps_name).exists():
            setting_defined['ERROR'] = [ERROR_STATE_ALREADY_DEFINED, data['state']]
            return setting_defined
            state_data['appps_name'] = apps_name
            state_data['state'] = data['state']
        StateInCountry.objects.create_new_state(state_data)
        setting_defined['settings'] = 'State'
        setting_defined['value'] = data['state']

    if 'city' in data:
        city_data = {}
        if City.objects.filter(cityName__iexact=data['city']['cityName'], apps_name__exact=apps_name).exists():
            setting_defined['ERROR'] = [ERROR_CITY_ALREADY_DEFINED, data['city']['cityName']]
            return setting_defined
        city_data['apps_name'] = apps_name
        city_data['state'] = data['city']['state']
        city_data['cityName'] = data['city']['cityName']
        city_data['latitude'] = data['city']['latitude']
        city_data['longitude'] = data['city']['longitude']
        City.objects.create_new_city(city_data)
        setting_defined['settings'] = 'City'
        setting_defined['value'] = data['city']['cityName']
    return setting_defined
