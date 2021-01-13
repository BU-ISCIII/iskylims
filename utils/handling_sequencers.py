from iSkyLIMS_core.utils.handling_platforms import *
from iSkyLIMS_core.models import  SequencingConfiguration, SequencerInLab
from iSkyLIMS_wetlab.wetlab_config import ERROR_SEQUENCER_ALREADY_DEFINED, ERROR_SEQUENCER_CONFIGURATION_ALREADY_DEFINED


def configuration_sequencer_exists():
    '''
    Description:
        The function check if exists sequencer configuration
    Return:
        True if they are defined
    '''
    if SequencingConfiguration.objects.all().exists():
        return True
    return False

def define_new_sequencer(form_data):
    '''
    Description:
        The function get the input information from the user form and create a
        new sequencer.
        Return error message if sequencer name already exists
    Constants:
        ERROR_SEQUENCER_ALREADY_DEFINED
    Return:
        new_sequencer_name or error
    '''
    seq_data = {}
    if SequencerInLab.objects.filter(sequencerName__iexact = form_data['sequencerName']).exists():
        error = {}
        error ['ERROR'] = ERROR_SEQUENCER_ALREADY_DEFINED
        return  error
    fields_in_sequencer = ['platformID', 'sequencerName', 'sequencerDescription', 'sequencerLocation', 'sequencerSerialNumber',
            'sequencerOperationStart','sequencerNumberLanes' ]

    for item in fields_in_sequencer :
        seq_data[item] = form_data[item]
    new_sequencer_name = SequencerInLab.objects.create_sequencer_in_lab(seq_data).get_sequencer_name()
    return new_sequencer_name

def define_new_seq_configuration(form_data):
    '''
    Description:
        The function get the input information from the user form and create a
        new configuration for the sequencer.
        Return error message if sequencer configuration already exists
    Constants:
        ERROR_SEQUENCER_CONFIGURATION_ALREADY_DEFINED
    Return:
        new_configuration_data or error
    '''
    seq_configuration_data = {}
    seq_configuration_data['platformID'] = form_data['platformID']
    seq_configuration_data['configurationName'] = form_data['sequencerConfiguration']
    if SequencingConfiguration.objects.filter(configurationName__iexact = seq_configuration_data['configurationName']).exists():
        error = {}
        error ['ERROR'] = ERROR_SEQUENCER_CONFIGURATION_ALREADY_DEFINED
        return  error
    new_seq_configuration = SequencingConfiguration.objects.create_new_configuration(seq_configuration_data)
    new_configuration_data = [seq_configuration_data['configurationName'],new_seq_configuration.get_platform_name()]
    import pdb; pdb.set_trace()
    return new_configuration_data

def get_configuration_sequencers_data():
    '''
    Description:
        The function get sequencers defined per platform. Returns a dictionary with
        platform as key and the values as a list of sequencers
    Return:
        sequencer_configuration
    '''
    sequencer_configuration = {}
    sequencer_configuration['configuration_platform'] = []
    sequencer_configuration['configuration_platform_option'] = []
    def_platforms = get_platform_name_of_defined_sequencers()

    if len(def_platforms) == 0:
        return sequencer_configuration
    for platform_name , platform_id in def_platforms.items():
        seq_conf_objs = SequencingConfiguration.objects.filter(platformID__exact = platform_id).order_by('configurationName')
        sequencer_configuration['configuration_platform'].append(platform_name)
        conf_data =[]
        for seq_conf_obj in seq_conf_objs:
            conf_data.append(seq_conf_obj.get_configuration_name())
        sequencer_configuration['configuration_platform_option'].append([platform_name, conf_data])

    return sequencer_configuration

def get_defined_sequencers():
    '''
    Description:
        The function get sequencers defined per platform. Returns a dictionary with
        platform as key and the values as a list of sequencers
    Return:
        sequencer_names
    '''
    sequencer_names = {}
    if SequencerInLab.objects.all().exists():
        sequencer_objs = SequencerInLab.objects.all().order_by('platformID')
        for sequencer_obj in sequencer_objs:
            seq_platform = sequencer_obj.get_sequencing_platform_name()
            if seq_platform not in sequencer_names:
                sequencer_names[seq_platform] = []
            sequencer_names[seq_platform].append(sequencer_obj.get_sequencer_name())
    return sequencer_names


def get_platform_name_of_defined_sequencers ():
    '''
    Description:
        The function get a dictionary with the platform names of the defined sequencers
    Return:
        platforms
    '''
    platforms = {}
    if SequencerInLab.objects.all().exists():
        sequencer_objs = SequencerInLab.objects.all()
        for sequencer_obj in sequencer_objs:
            try:
                platform_name = sequencer_obj.get_sequencing_platform_name()
                if platform_name not in platforms :
                    platforms[platform_name]= sequencer_obj.get_sequencing_platform_id()
            except:
                continue
    return platforms



def get_list_sequencer_configuration ():
    '''
    Description:
        The function get the list of sequencer configuration, and the platform
        defined

    Return
        sequencer_configuration
    '''
    sequencer_configuration = {}
    sequencer_configuration['platforms_used']  = get_platform_name_of_defined_sequencers()
    if SequencingConfiguration.objects.all().exists():
        sequencer_data = {}
        seq_conf_objs = SequencingConfiguration.objects.all().order_by('platformID')
        for seq_conf_obj in seq_conf_objs:
            platform_name = seq_conf_obj.get_platform_name()
            if platform_name not in sequencer_data :
                sequencer_data[platform_name] = []
            sequencer_data[platform_name].append(seq_conf_obj.get_configuration_name())
        sequencer_configuration['sequencer_data'] = sequencer_data

    return sequencer_configuration



def get_platform_data ():
    '''
    Description:
        The function get the platform list and their id
    Functions:
        get_defined_platforms_and_ids  # located at iSkyLIMS_core.utils.handling_platforms
    Return
        sequencer_configuration
    '''
    return get_defined_platforms_and_ids('NGS')
