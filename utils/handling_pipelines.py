from .drylab_common_functions import *

def get_data_form_pipeline_actions():
    '''
    Description:
        The function collect the available services and the fields to present information
        in the form
    Functions:
        get_available_services  # located at drylab_common_functions
    Return:
        data_actions
    '''
    data_actions = {}
    data_actions ['available_services']= get_children_available_services ()
    data_actions['heading']= drylab_config.HEADING_ACTIONS_PIPELINES
    data_actions['available_actions'] = drylab_config.AVAILABLE_ACTIONS_IN_PIPELINE
    return data_actions
