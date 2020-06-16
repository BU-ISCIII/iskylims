import json
from .drylab_common_functions import *
from iSkyLIMS_drylab import drylab_config
from iSkyLIMS_drylab.models import *


def analyze_input_pipelines(request):
    '''
    Description:
        The function extract the information from the user form
    Return:
        pipeline_actions
    '''
    pipeline_actions = {}
    pipeline_actions['availableService_id'] = request.POST['service_id']
    pipeline_actions['userName'] = request.user
    pipeline_actions['pipelineName'] =request.POST['pipelineName']
    pipeline_actions['pipelineVersion'] = request.POST['pipelineVersion']
    pipeline_actions['pipelineStrFolder'] = request.POST['pipelineStrFolder']
    pipeline_actions['automatic'] = request.POST['automaticAction']
    pipeline_json_data = json.loads(request.POST['pipeline_data'])

    action_length = len(drylab_config.HEADING_ACTIONS_PIPELINES)
    pipeline_actions['actions'] = []
    pipeline_actions['parameters'] = []
    for row in pipeline_json_data :
        action_name = str(row[drylab_config.HEADING_ACTIONS_PIPELINES.index('Given name for action')])
        if action_name == '' :
            continue
        action_dict = {}
        param_dict = {}
        for i in range(len(drylab_config.HEADING_ACTIONS_PIPELINES)):
            action_dict[drylab_config.HEADING_ACTIONS_PIPELINES[i]] = row[i]
        pipeline_actions['actions'].append(action_dict)

        for i in range(len(drylab_config.HEADING_ACTIONS_PARAMETERS)):
            param_dict[drylab_config.HEADING_ACTIONS_PARAMETERS[i]] = row[action_length + i]
        pipeline_actions['parameters'].append(param_dict)
    pipeline_actions['data'] = pipeline_json_data

    return pipeline_actions

def check_if_pipelines_exists_for_service(service_id):
    '''
    Description:
        The function check if the service had already a pipeline
    Return:
        True or False
    '''
    if Pipelines.objects.filter(availableService__pk__exact = service_id).count() > 1:
        return True
    return False

def get_data_form_pipeline_actions():
    '''
    Description:
        The function collect the available services and the fields to present information
        in the form
    Functions:
        get_children_available_services  # located at drylab_common_functions
    Return:
        data_actions
    '''
    data_actions = {}
    data_actions ['available_services']= get_children_available_services ()
    data_actions['heading']= drylab_config.HEADING_ACTIONS_PIPELINES +  drylab_config.HEADING_ACTIONS_PARAMETERS
    data_actions['available_actions'] = drylab_config.AVAILABLE_ACTIONS_IN_PIPELINE
    return data_actions

def get_detail_pipeline_data(pipeline_id):
    '''
    Description:
        The function collect the all information for the pipeline
    Input:
        pipeline_id     # id of the pipeline
    Return:
        detail_pipelines_data
    '''
    detail_pipelines_data = {}
    if Pipelines.objects.filter(pk__exact = pipeline_id).exists():
        pipeline_obj = Pipelines.objects.get(pk__exact = pipeline_id)
        detail_pipelines_data['pipeline_name'] = pipeline_obj.get_pipeline_name()
        detail_pipelines_data['pipeline_basic'] = pipeline_obj.get_pipeline_basic()
        detail_pipelines_data['pipeline_basic_heading'] = drylab_config.DISPLAY_DETAIL_PIPELINE_BASIC_INFO
        detail_pipelines_data['pipeline_additional_data'] = zip( drylab_config.DISPLAY_DETAIL_PIPELINE_ADDITIONAL_INFO, pipeline_obj.get_pipeline_additional())
        detail_pipelines_data['actions_heading'] = drylab_config.HEADING_ACTIONS_PIPELINES + drylab_config.HEADING_ACTIONS_PARAMETERS
        pipeline_obj = get_pipeline_obj_from_id (pipeline_id)
        if ActionPipeline.objects.filter(pipeline = pipeline_obj).exists():
            action_objs = ActionPipeline.objects.filter(pipeline = pipeline_obj).order_by('order')
            detail_pipelines_data['actions'] = []
            for action in action_objs :
                parameter_obj = ParameterActionPipeline.objects.get(actionPipeline = action)
                detail_pipelines_data['actions'].append(action.get_action_pipeline_data() + parameter_obj.get_all_action_parameters())



    import pdb; pdb.set_trace()
    return detail_pipelines_data

def get_pipeline_data_to_display(pipeline_information):
    '''
    Description:
        The function collect data to display information for a new defined pipeline

    Return:
        pipeline_data
    '''
    pipeline_data = {}
    if check_if_pipelines_exists_for_service(pipeline_information['availableService_id']):
        return  get_pipelines_for_service(pipeline_information['availableService_id'])
    service_name = get_service_name_from_id (pipeline_information['availableService_id'])
    pipeline_data['one_pipeline'] = [service_name, pipeline_information['pipelineName'],pipeline_information['pipelineVersion']]
    pipeline_data['heading_one_pipeline']= drylab_config.DISPLAY_NEW_DEFINED_PIPELINE
    return pipeline_data

def get_pipelines_for_manage():
    '''
    Description:
        The function get all pipelines defined for the services
    Return:
        pipeline_data
    '''
    pipeline_data = {}
    if Pipelines.objects.all().exists():
        pipelines_objs = Pipelines.objects.all().order_by('availableService')
        pipeline_data['data'] = []
        for pipeline in pipelines_objs:
            pipeline_data['data'].append(pipeline.get_pipeline_info())
        pipeline_data['heading'] = drylab_config.HEADING_MANAGE_PIPELINES

    return pipeline_data

def get_pipelines_for_service(service_id):
    '''
    Description:
        The function get the service name from service id
    Return:
        service_name
    '''
    pipeline_data = {}
    pipeline_data['multiple_pipeline'] = []
    pipeline_objs = Pipelines.objects.filter(availableService__pk__exact = service_id)
    for pipeline_obj in pipeline_objs :
        pipeline_data['multiple_pipeline'].append(pipeline_obj.get_pipeline_info())
    pipeline_data['heading_multi_pipeline'] = drylab_config.DISPLAY_MULTYPLE_DEFINED_PIPELINE
    return pipeline_data

def get_pipeline_obj_from_id (pipeline_id):
    '''
    Description:
        The function get the pipeline object from the pipeline id
    Return:
        pipeline_obj
    '''
    return Pipelines.objects.get(pk__exact = pipeline_id)

def get_service_name_from_id(service_id):
    '''
    Description:
        The function get the service name from service id
    Return:
        service_name
    '''
    return AvailableService.objects.get(pk__exact = service_id).get_service_description()

def pipeline_version_exists(pipeline_name, pipeline_version):
    '''
    Description:
        The function check if pipeline name and version already exists
    Return:
        True if already exists
    '''
    if Pipelines.objects.filter(pipelineName__iexact = pipeline_name, pipelineVersion__iexact = pipeline_version).exists():
        return True
    return False


def store_pipeline_actions(pipeline_obj, pipeline_actions, pipeline_parameters):
    '''
    Description:
        The function store in database the actions and parameters for the pipeline

    Return:
        True if already exists
    '''
    for i in range(len(pipeline_actions)):
        new_pipeline_action = ActionPipeline.objects.create_pipeline_action(pipeline_obj, pipeline_actions[i])
        new_pipeline_parameters = ParameterActionPipeline.objects.create_pipeline_parameters( new_pipeline_action, pipeline_parameters[i])
    return
