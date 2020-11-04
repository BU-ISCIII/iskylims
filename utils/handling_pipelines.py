import json
from .drylab_common_functions import *
from iSkyLIMS_drylab import drylab_config
from iSkyLIMS_drylab.models import *
from iSkyLIMS_wetlab.utils.api.wetlab_api import get_run_folder_from_user_project


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
    pipeline_json_data = json.loads(request.POST['pipeline_data'])

    action_length = len(drylab_config.HEADING_PIPELINE_PARAMETERS)

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
    additional_data['heading']= drylab_config.HEADING_PIPELINE_PARAMETERS

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

        #pipeline_obj = get_pipeline_obj_from_id (pipeline_id)
        if ParameterPipeline.objects.filter(parameterPipeline = pipeline_obj).exists():
            parameter_objs = ParameterPipeline.objects.filter(parameterPipeline = pipeline_obj)
            detail_pipelines_data['parameter_heading'] = drylab_config.HEADING_PARAMETER_PIPELINE
            detail_pipelines_data['parameters'] = []
            for parameter_obj in parameter_objs :

                detail_pipelines_data['parameters'].append([parameter_obj.get_pipeline_parameter_name(), parameter_obj.get_pipeline_parameter_value()])

    return detail_pipelines_data

def get_defined_pipeline_data_to_display(pipeline_obj):
    '''
    Description:
        The function collect data to display information for a new defined pipeline
    Input:
        pipeline_obj    # pipeline obj to get the information
    Return:
        pipeline_data
    '''
    pipeline_data = {}

    service_name = pipeline_obj.get_pipleline_avail_service()
    pipeline_data['pipeline_data'] = [service_name, pipeline_obj.get_pipeline_name(), pipeline_obj.get_pipeline_version()]
    pipeline_data['heading_pipeline']= drylab_config.DISPLAY_NEW_DEFINED_PIPELINE
    pipeline_data['parameters'] = get_pipeline_parameters(pipeline_obj)
    if len(pipeline_data['parameters']) > 0:
        pipeline_data['heading_parameter_pipeline'] = drylab_config.HEADING_PARAMETER_PIPELINE
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

def get_pipeline_and_versions_for_available_service(available_service):
    '''
    Description:
        The function get the list of defined pipeline names and version for an
        available service
    Input:
        available_service
    Return:
        service_name
    '''
    pipeline_names = []
    if  Pipelines.objects.filter(availableService__pk__exact = available_service).exists():
        pipelines_objs = Pipelines.objects.filter(availableService__pk__exact = available_service).order_by('generated_at').reverse()
        for pipelines_obj in pipelines_objs:
            pipeline_names.append([pipelines_obj.pk, pipelines_obj.get_pipleline_avail_service(), pipelines_obj.get_pipeline_version()])
    return pipeline_names

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

def get_avail_service_name_from_id(service_id):
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
def get_pipeline_parameters(pipeline_obj):
    '''
    Description:
        The function return the parameters defined for the pipeline
    Input:
        pipeline_obj    # pipeline instance
    Return:
        pipeline_parameters
    '''
    pipeline_parameters = []
    if ParameterPipeline.objects.filter(parameterPipeline = pipeline_obj ).exists():
        parameter_objs = ParameterPipeline.objects.filter(parameterPipeline = pipeline_obj)
        for parameter_obj in parameter_objs:
            pipeline_parameters.append(parameter_obj.get_pipeline_parameters())
    return pipeline_parameters

def store_parameters_pipeline(pipeline_obj, parameters):
    '''
    Description:
        The function store the parameters pipeline
    Input:
        pipeline_obj    # pipeline instance
        paremeters      # dict paramters names and values
    Return:
        None
    '''
    import pdb; pdb.set_trace()
    for item, value in parameters.items():
        parameter_pipeline = {}
        parameter_pipeline['parameterPipeline'] = pipeline_obj
        parameter_pipeline['parameterName'] = item
        parameter_pipeline['parameterValue'] = value
        ParameterPipeline.objects.create_pipeline_parameters(parameter_pipeline)

    return None
