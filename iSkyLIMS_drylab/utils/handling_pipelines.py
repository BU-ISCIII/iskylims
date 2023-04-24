# Generic imports
import json, os

# Local imports
import iSkyLIMS_drylab.drylab_config
import iSkyLIMS_drylab.models
import iSkyLIMS_drylab.utils.handling_request_services


def analyze_input_pipelines(request):
    """
    Description:
        The function extract the information from the user form
    Return:
        pipeline_data
    """
    pipeline_data = {}
    additional_param_dict = {}

    pipeline_data["user_name"] = request.user
    pipeline_data["pipeline_name"] = request.POST["pipeline_name"]
    pipeline_data["pipeline_version"] = request.POST["pipeline_version"]
    pipeline_data["url"] = request.POST["urlLocation"]
    pipeline_data["description"] = request.POST["description"]
    pipeline_json_data = json.loads(request.POST["pipeline_data"])
    for row in pipeline_json_data:
        if row[0] == "":
            continue
        additional_param_dict[row[0]] = row[1]
    pipeline_data["additional_parameters"] = additional_param_dict
    return pipeline_data


def check_if_pipelines_exists_for_service(service_id):
    """
    Description:
        The function check if the service had already a pipeline
    Return:
        True or False
    """
    if (
        iSkyLIMS_drylab.models.Pipelines.objects.filter(
            availableService__pk__exact=service_id
        ).count()
        > 1
    ):
        return True
    return False


def get_all_defined_pipelines(only_in_used):
    """
    Description:
        The function return a list with the defined pipelines
    Input:
        only_in_used    # boolean varible ir True only pipelines in used are returned.
                        # if not all defined pipelines are returned back
    Return:
        defined_pipelines
    """
    defined_pipelines = []
    pipeline_objs = False
    if only_in_used:
        if iSkyLIMS_drylab.models.Pipelines.objects.filter(
            pipeline_in_use=True
        ).exists():
            pipeline_objs = (
                iSkyLIMS_drylab.models.Pipelines.objects.filter(pipeline_in_use=True)
                .order_by("pipeline_name")
                .order_by("pipeline_version")
            )
    else:
        if iSkyLIMS_drylab.models.Pipelines.objects.all().exists():
            pipeline_objs = (
                iSkyLIMS_drylab.models.Pipelines.objects.all()
                .order_by("pipeline_name")
                .order_by("pipeline_version")
            )
    if pipeline_objs:
        for pipeline_obj in pipeline_objs:
            defined_pipelines.append(
                [
                    pipeline_obj.get_pipeline_name(),
                    pipeline_obj.get_pipeline_version(),
                    pipeline_obj.get_pipeline_id(),
                ]
            )
    return defined_pipelines


def get_data_form_pipeline():
    """
    Description:
        The function collect the available services and the fields to present information
        in the form
    Return:
        additional_data
    """
    additional_data = {}
    # additional_data ['available_services']= get_children_available_services ()
    additional_data[
        "heading"
    ] = iSkyLIMS_drylab.drylab_config.HEADING_PARAMETER_PIPELINE

    return additional_data


def get_detail_pipeline_data(pipeline_id):
    """
    Description:
        The function collect the all information for the pipeline
    Input:
        pipeline_id     # id of the pipeline
    Constants:
        HEADING_SERVICES_IN_PIPELINE
        DISPLAY_DETAIL_PIPELINE_BASIC_INFO
        DISPLAY_DETAIL_PIPELINE_ADDITIONAL_INFO
        HEADING_PARAMETER_PIPELINE
    Return:
        detail_pipelines_data
    """
    detail_pipelines_data = {}
    if iSkyLIMS_drylab.models.Pipelines.objects.filter(pk__exact=pipeline_id).exists():
        pipeline_obj = iSkyLIMS_drylab.models.Pipelines.objects.get(
            pk__exact=pipeline_id
        )
        detail_pipelines_data["pipeline_name"] = pipeline_obj.get_pipeline_name()
        detail_pipelines_data["pipeline_basic"] = pipeline_obj.get_pipeline_basic()
        detail_pipelines_data[
            "pipeline_basic_heading"
        ] = iSkyLIMS_drylab.drylab_config.DISPLAY_DETAIL_PIPELINE_BASIC_INFO
        detail_pipelines_data["pipeline_additional_data"] = zip(
            iSkyLIMS_drylab.drylab_config.DISPLAY_DETAIL_PIPELINE_ADDITIONAL_INFO,
            pipeline_obj.get_pipeline_additional(),
        )

        # pipeline_obj = get_pipeline_obj_from_id (pipeline_id)
        if iSkyLIMS_drylab.models.ParameterPipeline.objects.filter(
            parameter_pipeline=pipeline_obj
        ).exists():
            parameter_objs = iSkyLIMS_drylab.models.ParameterPipeline.objects.filter(
                parameter_pipeline=pipeline_obj
            )
            detail_pipelines_data[
                "parameter_heading"
            ] = iSkyLIMS_drylab.drylab_config.HEADING_PARAMETER_PIPELINE
            detail_pipelines_data["parameters"] = []
            for parameter_obj in parameter_objs:
                detail_pipelines_data["parameters"].append(
                    [
                        parameter_obj.get_pipeline_parameter_name(),
                        parameter_obj.get_pipeline_parameter_type(),
                    ]
                )
        # get the services where the pipeline was used
        """
        req_serv_objs = iSkyLIMS_drylab.utils.handling_request_services.get_requested_services_obj_from_available_service(pipeline_obj.get_pipleline_avail_service_obj())
        if req_serv_objs != None :
            detail_pipelines_data['services_using_pipeline'] = []
            detail_pipelines_data['services_using_pipeline_heading'] = drylab_config.HEADING_SERVICES_IN_PIPELINE
            for req_serv_obj in req_serv_objs:
                detail_pipelines_data['services_using_pipeline'].append([req_serv_obj.get_service_id() ,req_serv_obj.get_service_request_number(),
                req_serv_obj.get_service_creation_time(), req_serv_obj.get_service_requested_user(), req_serv_obj.get_service_state()])
        """
    return detail_pipelines_data


def get_defined_pipeline_data_to_display(pipeline_obj):
    """
    Description:
        The function collect data to display information for a new defined pipeline
    Input:
        pipeline_obj    # pipeline obj to get the information
    Return:
        pipeline_data
    """
    pipeline_data = {}

    # service_name = pipeline_obj.get_pipleline_avail_service()
    pipeline_data["pipeline_data"] = [
        pipeline_obj.get_pipeline_name(),
        pipeline_obj.get_pipeline_version(),
        pipeline_obj.get_pipeline_description(),
    ]
    pipeline_data[
        "heading_pipeline"
    ] = iSkyLIMS_drylab.drylab_config.DISPLAY_NEW_DEFINED_PIPELINE
    return pipeline_data


def get_pipelines_for_manage():
    """
    Description:
        The function get all pipelines defined for the services
    Return:
        pipeline_data
    """
    pipeline_data = {}
    if iSkyLIMS_drylab.models.Pipelines.objects.all().exists():
        pipelines_objs = iSkyLIMS_drylab.models.Pipelines.objects.filter(
            pipeline_in_use__exact=True
        ).order_by("pipeline_name")
        pipeline_data["data"] = []
        for pipeline in pipelines_objs:
            pipeline_data["data"].append(pipeline.get_pipeline_info())
        pipeline_data[
            "heading"
        ] = iSkyLIMS_drylab.drylab_config.HEADING_MANAGE_PIPELINES

    return pipeline_data


def get_pipelines_for_resolution(resolution_obj):
    """
    Description:
        The function get the list of defined pipelines for the resolution
    Input:
        resolution_obj      # resolution object
    Return:
        pipeline_data
    """
    pipeline_data = {}
    pipeline_data["pipelines"] = []
    pipeline_objs = resolution_obj.resolution_pipelines.all()
    for pipeline_obj in pipeline_objs:
        pipeline_data["pipelines"].append(pipeline_obj.get_pipeline_info())
    if len(pipeline_data["pipelines"]) > 0:
        pipeline_data[
            "heading_pipelines"
        ] = iSkyLIMS_drylab.drylab_config.DISPLAY_PIPELINES_USED_IN_RESOLUTION
    return pipeline_data


def get_pipeline_obj_from_id(pipeline_id):
    """
    Description:
        The function get the pipeline object from the pipeline id
    Return:
        pipeline_obj
    """
    return iSkyLIMS_drylab.models.Pipelines.objects.get(pk__exact=pipeline_id)


def get_avail_service_name_from_id(service_id):
    """
    Description:
        The function get the service name from service id
    Return:
        service_name
    """
    return iSkyLIMS_drylab.models.AvailableService.objects.get(
        pk__exact=service_id
    ).get_service_description()


def pipeline_version_exists(pipeline_name, pipeline_version):
    """
    Description:
        The function check if pipeline name and version already exists
    Return:
        True if already exists
    """
    if iSkyLIMS_drylab.models.Pipelines.objects.filter(
        pipeline_name__iexact=pipeline_name, pipeline_version__iexact=pipeline_version
    ).exists():
        return True
    return False


def get_pipeline_parameters(pipeline_obj):
    """
    Description:
        The function return the parameters defined for the pipeline
    Input:
        pipeline_obj    # pipeline instance
    Return:
        pipeline_parameters
    """
    pipeline_parameters = []
    if iSkyLIMS_drylab.models.ParameterPipeline.objects.filter(
        parameter_pipeline=pipeline_obj
    ).exists():
        parameter_objs = iSkyLIMS_drylab.models.ParameterPipeline.objects.filter(
            parameter_pipeline=pipeline_obj
        )
        for parameter_obj in parameter_objs:
            pipeline_parameters.append(parameter_obj.get_pipeline_parameters())
    return pipeline_parameters


def store_parameters_pipeline(pipeline_obj, parameters):
    """
    Description:
        The function store the parameters pipeline
    Input:
        pipeline_obj    # pipeline instance
        paremeters      # dict paramters names and parameter type
    Return:
        None
    """

    for item, value in parameters.items():
        parameter_pipeline = {}
        parameter_pipeline["parameter_pipeline"] = pipeline_obj
        parameter_pipeline["parameterName"] = item
        parameter_pipeline["parameterType"] = value
        iSkyLIMS_drylab.models.ParameterPipeline.objects.create_pipeline_parameters(
            parameter_pipeline
        )

    return None
