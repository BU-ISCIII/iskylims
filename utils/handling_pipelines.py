import json
from .drylab_common_functions import *
from iSkyLIMS_drylab import drylab_config
from iSkyLIMS_drylab.models import *


def analyze_input_pipelines(request):
    '''
    Description:
        The function extract the information from the user form
    Return:
        pipeline_data
    '''
    pipeline_data = {}
    pipeline_data['availableService_id'] = request.POST['service_id']
    pipeline_data['userName'] = request.user
    pipeline_data['pipelineName'] =request.POST['pipelineName']
    pipeline_data['pipelineVersion'] = request.POST['pipelineVersion']
    pipeline_data['externalRequest'] = request.POST['externalRequest']
    pipeline_data['useRunFolder'] = request.POST['useRunFolder']
    pipeline_json_data = json.loads(request.POST['pipeline_data'])

    action_length = len(drylab_config.HEADING_ADDITIONAL_PIPELINE_PARAMETERS)

    additional_param_dict = {}
    for row in pipeline_json_data :
        if row[0] == '' :
            continue
        additional_param_dict[row[0]] = row[1]
    pipeline_data['additional_parameters'] = additional_param_dict
    return pipeline_data

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


def get_data_form_pipeline():
    '''
    Description:
        The function collect the available services and the fields to present information
        in the form
    Functions:
        get_children_available_services  # located at drylab_common_functions
    Return:
        additional_data
    '''
    additional_data = {}
    additional_data ['available_services']= get_children_available_services ()
    additional_data['heading']= drylab_config.HEADING_ADDITIONAL_PIPELINE_PARAMETERS

    return additional_data

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

    return detail_pipelines_data

def get_pipeline_data_to_display(pipeline_information):
    '''
    Description:
        The function collect data to display information for a new defined pipeline

    Return:
        pipeline_data
    '''
    pipeline_data = {}
    service_name = get_service_name_from_id (pipeline_information['availableService_id'])
    pipeline_data['pipeline_data'] = [service_name, pipeline_information['pipelineName'],pipeline_information['pipelineVersion']]
    pipeline_data['heading_pipeline']= drylab_config.DISPLAY_NEW_DEFINED_PIPELINE
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
        The function get the list of defined pipelines for a service
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

def send_required_preparation_pipeline_email(service_request_number):
    '''
    Description:
        The function sends an email to defined user in configuration with the service number.

    Input:
        service_request_number     # service number
    Constant:
        SUBJECT_SERVICE_ON_QUEUED
    Return:
        None
    '''
    subject = drylab_config.SUBJECT_SERVICE_ON_QUEUED.copy()
    subject.insert(1, service_request_number)

    body_preparation = list(map(lambda st: str.replace(st, 'SERVICE_NUMBER', service_request_number), drylab_config.BODY_SERVICE_ON_QUEUED))
    body_message = '\n'.join(body_preparation)

    from_user = drylab_config.USER_EMAIL
    to_users = [drylab_config.USER_EMAIL]

    send_mail (subject, body_message, from_user, to_users)
    return

def services_allow_external_data(service_obj, project_requested_objs):
    '''
    Description:
        The function check if services requires to populate in the  starting the pipeline.
        Then adds the service into the actionJobs table.
    Input:
        service_obj     # service instance to check if requires preparation
    Functions:
        get_run_folder_from_user_project  # API from iSkyLIMS_wetlab located at file wetlab_api
    Return:
        services_added
    '''
    services_added = []
    #
    all_services = service_obj.get_child_services()
    for service_id, service_name in all_services:
        if Pipelines.objects.filter(availableService__pk__exact = service_id, default = True).exists():
            pipeline_obj = Pipelines.objects.filter(availableService__pk__exact = service_id, default = True).last()
            if pipeline_obj.get_external_request() == 'True':
                preparation_data = {}
                preparation_data['pipeline'] = pipeline_obj
                preparation_data['availableService'] = AvailableService.objects.get(pk__exact = service_id)
                if pipeline_obj.get_used_run_folder() == 'True':
                    runID_folders = []
                    for project_request in project_requested_objs:
                        run_folder = get_run_folder_from_user_project(project_request.get_requested_project_id())
                        if run_folder :
                            if not run_folder in runID_folders :
                                runID_folders.append(run_folder)
                                preparation_data['pendingToSetFolder'] = False
                                preparation_data['folderData'] = run_folder
                                services_added.append(PipelineExternalDataJobs.objects.create_pipeline_external_data_job(preparation_data))
                        else:
                            preparation_data['pendingToSetFolder'] = True
                            preparation_data['folderData'] = None
                            services_added.append(PipelineExternalDataJobs.objects.create_pipeline_external_data_job(preparation_data))
                else:
                    preparation_data['pendingToSetFolder'] = False
                    preparation_data['folderData'] = ''
                    services_added.append(PreparationPipelineJobs.objects.create_preparation_pipeline_job(preparation_data))
    return services_added

def set_default_service_pipeline(new_pipeline):
    '''
    Description:
        The function check if an older version of the service pipeline exist .
        If true the old version is set unchecked for default uses
    Input:
        new_pipeline    # pipeline instance
    Return:
        None
    '''

    service_pipeline_obj = new_pipeline.get_pipleline_service_obj()
    if Pipelines.objects.filter(availableService = service_pipeline_obj, default__exact = True).exists():
        old_version_pipeline = Pipelines.objects.filter(availableService = service_pipeline_obj, default__exact = True).last()
        old_version_pipeline.remove_default_pipeline()
    new_pipeline.set_default_pipeline()
    return

def store_pipeline_actions(pipeline_obj, pipeline_parameters):
    '''
    Description:
        The function store in database the additional parameters for the pipeline

    Return:
        None
    '''
    for item, value in pipeline_parameters.items():
        param_pipeline = {}
        param_pipeline['parameterPipeline'] = pipeline_obj
        param_pipeline['parameterName'] = item
        param_pipeline['parameterValue'] = value
        new_pipeline_parameters = ParameterPipeline.objects.create_pipeline_parameters(param_pipeline)
    return
