# Generic imports
from datetime import date, datetime
import json
from smtplib import SMTPException
import django.core.mail
import django_utils.models
import django.contrib.auth.models

# Local imports
import iSkyLIMS_drylab.models
import iSkyLIMS_drylab.drylab_config
import iSkyLIMS_drylab.utils.handling_request_services
import iSkyLIMS_drylab.utils.handling_pipelines


def add_pipelines_to_resolution(resolution_obj, pipeline_ids):
    """
    Description:
        The function add pipelines to the resolution
    Input:
        resolution_obj  # resoution obj
        pipeline_ids    # ids of the pipelines
    Return:
        None
    """
    for pipeline_id in pipeline_ids:
        resolution_obj.resolution_pipelines.add(
            iSkyLIMS_drylab.utils.handling_pipelines.get_pipeline_obj_from_id(
                pipeline_id
            )
        )

    return None


def allow_to_service_update_state(resolution_obj, new_state):
    """
    Description:
        The function check if all partial resolutions are handled all requested
        services
    Input:
        resolution_obj  # resoution obj
        new_state       # 'in_progress' or 'delivered'
    Return:
        True or False
    """
    service_obj = resolution_obj.get_service_obj()
    if (
        iSkyLIMS_drylab.models.Resolution.objects.filter(
            resolution_serviceID=service_obj
        )
        .exclude(resolution_state__state_value__exact="Queued")
        .exists()
    ):
        all_resolution_objs = iSkyLIMS_drylab.models.Resolution.objects.filter(
            resolution_serviceID=service_obj
        ).exclude(resolution_state__state_value__exact="Queued")
        avail_services_handled = []
        for all_resolution_obj in all_resolution_objs:
            if new_state == "delivered":
                if not (
                    all_resolution_obj.get_state() == "Delivered"
                    or all_resolution_obj.get_state() == "Cancelled"
                    or all_resolution_obj == resolution_obj
                ):
                    continue
            resolution_handle_list = all_resolution_obj.get_available_services()
            for item in resolution_handle_list:
                avail_services_handled.append(item)
        if len(set(avail_services_handled)) == len(service_obj.get_child_services()):
            return True
    return False


def get_assign_resolution_full_number(service_id, acronymName):
    """
    Description:
        The function get the resolution full number if resolution already exists.
        Build the resolution  full number if it is the first resolution for the service
    Input:
        service_id # contains the service id
        acronymName # acronym name given to the service
    Functions:
        get_service_obj_from_id   # located at this file
    Return:
        resolution_full_number
    """
    service_obj = (
        iSkyLIMS_drylab.utils.handling_request_services.get_service_obj_from_id(
            service_id
        )
    )
    if iSkyLIMS_drylab.models.Resolution.objects.filter(
        resolution_serviceID=service_id
    ).exists():
        resolution_full_number = (
            iSkyLIMS_drylab.models.Resolution.objects.filter(
                resolution_serviceID=service_obj
            )
            .last()
            .get_resolution_number()
        )
    else:
        resolution_full_number = ""
        resolution_full_number += service_obj.get_service_request_number() + "_"
        resolution_full_number += str(date.today()).replace("-", "") + "_"
        resolution_full_number += acronymName + "_"
        resolution_full_number += service_obj.get_service_requested_user() + "_S"
    return resolution_full_number


def create_resolution_number(service_id):
    """
    Description:
        The function create the resolution number and step it if more than 1 resolution
        have been created for the service.
    Input:
        service_id # contains the service id
    Functions:
        get_service_obj_from_id   # located at this file
    Return:
        resolution_number
    """
    service_obj = (
        iSkyLIMS_drylab.utils.handling_request_services.get_service_obj_from_id(
            service_id
        )
    )
    service_request_number = service_obj.get_service_request_number()
    if iSkyLIMS_drylab.models.Resolution.objects.filter(
        resolution_serviceID=service_obj
    ).exists():
        resolution_count = iSkyLIMS_drylab.models.Resolution.objects.filter(
            resolution_serviceID=service_obj
        ).count()
        resolution_number = service_request_number + "." + str(resolution_count + 1)
    else:
        resolution_number = service_request_number + ".1"
    return resolution_number


def get_data_for_resolution(service_obj, resolution_obj):
    information, user, resolution_data = {}, {}, {}
    service_number, center = service_obj.get_service_information()

    resolution_info = resolution_obj.get_resolution_information()
    # get profile object
    user_id = service_obj.serviceUserId.id

    information["resolution_number"] = resolution_obj.get_resolution_number()
    information["requested_date"] = service_obj.get_service_creation_time()
    information["resolution_date"] = resolution_info[4]
    if resolution_obj.available_services.all() == ["None"]:
        information["nodes"] = service_obj.service_available_service.all()
    else:
        information["nodes"] = resolution_obj.available_services.all()
    user["name"] = service_obj.serviceUserId.first_name
    user["surname"] = service_obj.serviceUserId.last_name

    user["area"] = django_utils.models.Profile.objects.get(
        profileUserID=user_id
    ).profileArea
    user["center"] = django_utils.models.Profile.objects.get(
        profileUserID=user_id
    ).profileCenter
    user["position"] = django_utils.models.Profile.objects.get(
        profileUserID=user_id
    ).profilePosition
    user["phone"] = django_utils.models.Profile.objects.get(
        profileUserID=user_id
    ).profileExtension
    user["email"] = service_obj.get_user_email()
    information["user"] = user
    resolution_info_split = resolution_info[2].split("_")
    resolution_data["acronym"] = resolution_info_split[2]
    resolution_data["estimated_date"] = resolution_info[4]
    resolution_data["notes"] = resolution_info[7]
    resolution_data["decission"] = service_obj.get_service_state()
    information["service_data"] = service_obj.get_service_user_notes()

    resolution_data["folder"] = resolution_info[2]
    information["resolution_data"] = resolution_data

    return information


def get_add_resolution_data_form(form_data):
    """
    Description:
        The function extract the user form information and store it in a dictionary
    Input:
        form_data		# contains the user data form
    Constants:
        HEADING_ADDITIONAL_RESOLUTION_PARAMETERS
    Return:
        resolution_data_form
    """
    resolution_data_form = {}

    resolution_data_form["service_id"] = form_data["service_id"]
    resolution_data_form["resolutionEstimatedDate"] = datetime.strptime(
        form_data["resolutionEstimatedDate"], "%Y-%m-%d"
    ).date()
    resolution_data_form["acronymName"] = form_data["acronymName"]
    resolution_data_form["resolution_asigned_user"] = form_data[
        "resolution_asigned_user"
    ]
    resolution_data_form["serviceAccepted"] = form_data["serviceAccepted"]
    resolution_data_form["resolutionNotes"] = form_data["resolutionNotes"]

    if "pipeline_data" in form_data:
        resolution_data_form["pipeline_ids"] = []
        selected_pipelines_table = json.loads(form_data["pipeline_data"])
        for row in selected_pipelines_table:
            if row[-1]:
                resolution_data_form["pipeline_ids"].append(row[-2])
    if "select_available_services" in form_data:
        # remove the last 'comma' in the string
        resolution_data_form["select_available_services"] = form_data[
            "select_available_services"
        ][:-1].split(",")

    return resolution_data_form


def check_if_resolution_exists(resolution_id):
    """
    Description:
        The function check if the resolution id exists
    Input:
        resolution_id   # resolution id
    Return:
        True or False
    """
    if iSkyLIMS_drylab.models.Resolution.objects.filter(
        pk__exact=resolution_id
    ).exists():
        return True
    return False


def create_new_resolution(resolution_data_form):
    """
    Description:
        The function create a new resolution. If there is some selected services
        this services are added to the resolution object
    Input:
        resolution_data_form	# contains the formated user data form
    Functions:
        get_assign_resolution_full_number  # located at this file
        create_resolution_number  # located at this file
        get_service_obj_from_id  # located at iSkyLIMS_drylab.utils.handling_request_services
        store_resolution_additional_parameter  # located at this file
        get_available_service_obj_from_id # located at iSkyLIMS_drylab.utils.handling_request_services
    Return:
        new_resolution
    """
    service_obj = (
        iSkyLIMS_drylab.utils.handling_request_services.get_service_obj_from_id(
            resolution_data_form["service_id"]
        )
    )
    if iSkyLIMS_drylab.models.Resolution.objects.filter(
        resolution_serviceID=service_obj
    ).exists():
        resolution_data_form["resolution_full_number"] = (
            iSkyLIMS_drylab.models.Resolution.objects.filter(
                resolution_serviceID=service_obj
            )
            .last()
            .get_resolution_full_number()
        )
    else:
        resolution_data_form[
            "resolution_full_number"
        ] = get_assign_resolution_full_number(
            resolution_data_form["service_id"], resolution_data_form["acronymName"]
        )
    resolution_data_form["resolutionNumber"] = create_resolution_number(
        resolution_data_form["service_id"]
    )

    # service_request_number = service_obj.get_service_request_number()
    new_resolution = iSkyLIMS_drylab.models.Resolution.objects.create_resolution(
        resolution_data_form
    )
    if "select_available_services" not in resolution_data_form:
        # include all services in the resolution
        resolution_data_form["select_available_services"] = []
        avail_services = service_obj.get_child_services()
        for avail_service in avail_services:
            resolution_data_form["select_available_services"].append(avail_service[0])

    # Add selected available services to the new resolution
    for avail_sarvice in resolution_data_form["select_available_services"]:
        avail_service_obj = iSkyLIMS_drylab.utils.handling_request_services.get_available_service_obj_from_id(
            avail_sarvice
        )
        new_resolution.available_services.add(avail_service_obj)

    if "pipelines" in resolution_data_form:
        for pipeline in resolution_data_form["pipelines"]:
            if pipeline != "":
                pipeline_obj = (
                    iSkyLIMS_drylab.utils.handling_pipelines.get_pipeline_obj_from_id(
                        pipeline
                    )
                )
                new_resolution.servicePipelines.add(pipeline_obj)

    if "additional_parameters" in resolution_data_form:
        store_resolution_additional_parameter(
            resolution_data_form["additional_parameters"], new_resolution
        )
    if resolution_data_form["serviceAccepted"] == "Accepted":
        if "select_available_services" in resolution_data_form:
            if len(resolution_data_form["select_available_services"]) == len(
                service_obj.get_child_services()
            ):
                service_obj.update_service_state("queued")
            elif iSkyLIMS_drylab.models.Resolution.objects.filter(
                resolution_serviceID=service_obj
            ).exists():
                resolution_objs = iSkyLIMS_drylab.models.Resolution.objects.filter(
                    resolution_serviceID=service_obj
                )
                avail_services_handled = []
                for resolution_obj in resolution_objs:
                    resolution_handle_list = resolution_obj.get_available_services()
                    for item in resolution_handle_list:
                        avail_services_handled.append(item)
                if len(set(avail_services_handled)) == len(
                    service_obj.get_child_services()
                ):
                    service_obj.update_service_state("queued")
        else:
            service_obj.update_service_state("queued")
        service_obj.update_approved_date(date.today())
    else:
        service_obj.update_service_state("rejected")
        service_obj.update_service_rejected_date(date.today())

    return new_resolution


def get_resolution_obj_from_id(resolution_id):
    """
    Description:
        The function get the resolution obj from its id
    Input:
        resolution_id		# resolution id
    Return:
        resolution_obj
    """
    resolution_obj = None
    if iSkyLIMS_drylab.models.Resolution.objects.filter(
        pk__exact=resolution_id
    ).exists():
        resolution_obj = iSkyLIMS_drylab.models.Resolution.objects.filter(
            pk__exact=resolution_id
        ).last()
    return resolution_obj


def prepare_form_data_add_resolution(form_data):
    """
    Description:
        The function collect additional info and then save the form
    Input:
        form_data		# contains the user data form
    Return:
        resolution_form_data
    """
    resolution_form_data = {}
    selected_children_services = []
    # pipelines_data = []
    if "childrenServices" in form_data:
        list_of_ch_services = form_data.getlist("childrenServices")
    else:
        list_of_ch_services = False
    service_obj = (
        iSkyLIMS_drylab.utils.handling_request_services.get_service_obj_from_id(
            form_data["service_id"]
        )
    )
    resolution_form_data["service_number"] = service_obj.get_service_request_number()
    all_tree_services = service_obj.service_available_service.all()
    all_children_services = iSkyLIMS_drylab.utils.handling_request_services.get_available_children_services_and_id(
        all_tree_services
    )

    if list_of_ch_services:
        if len(list_of_ch_services) != len(all_children_services):
            for children in list_of_ch_services:
                avail_serv_obj = iSkyLIMS_drylab.utils.handling_request_services.get_available_service_obj_from_id(
                    children
                )
                selected_children_services.append(
                    [children, avail_serv_obj.get_service_description()]
                )
            resolution_form_data[
                "selected_avail_services_data"
            ] = selected_children_services

    if iSkyLIMS_drylab.models.Resolution.objects.filter(
        resolution_serviceID=service_obj
    ).exists():
        existing_resolution = iSkyLIMS_drylab.models.Resolution.objects.filter(
            resolution_serviceID=service_obj
        ).last()
        resolution_form_data[
            "resolution_full_number"
        ] = existing_resolution.get_resolution_full_number()
    users = django.contrib.auth.models.User.objects.filter(
        groups__name=iSkyLIMS_drylab.drylab_config.SERVICE_MANAGER
    )
    resolution_form_data["assigned_user"] = []
    for user in users:
        resolution_form_data["assigned_user"].append([user.pk, user.username])

    resolution_form_data["service_id"] = form_data["service_id"]

    # get available pipelines for services
    if len(selected_children_services) == 0:
        req_available_services_with_desc = all_children_services
    else:
        req_available_services_with_desc = selected_children_services  # resolution_form_data['selected_avail_services_data']
    req_available_services_id = []
    for req_service in req_available_services_with_desc:
        req_available_services_id.append(req_service[0])
    # data = get_active_pipeline_and_versions()
    # for avail_service in req_available_services_id:
    #    data = get_pipeline_and_versions_for_available_service(avail_service)
    #    if data :
    #        pipelines_data.append([ get_available_service_obj_from_id(avail_service).get_service_description() , data])
    resolution_form_data[
        "pipelines_data"
    ] = iSkyLIMS_drylab.utils.handling_pipelines.get_all_defined_pipelines(True)
    resolution_form_data[
        "pipelines_heading"
    ] = iSkyLIMS_drylab.drylab_config.HEADING_PIPELINES_SELECTION_IN_RESOLUTION

    return resolution_form_data


def send_resolution_creation_email(email_data):
    """
    Description:
        The function send the service email for resolution to user.
        Functions uses the send_email django core function to send the email
    Input:
        email_data      # Contains the information to include in the email
    Constant:
        SUBJECT_RESOLUTIONQUEUED
        USER_EMAIL
        EMAIL_FOR_NOTIFICATIONS

    Return:
        None
    """
    subject_tmp = iSkyLIMS_drylab.drylab_config.SUBJECT_RESOLUTION_QUEUED.copy()
    subject_tmp.insert(1, email_data["service_number"])
    subject = " ".join(subject_tmp)
    if email_data["status"] == "Accepted":
        date = email_data["date"].strftime("%d %B, %Y")
        body_preparation = list(
            map(
                lambda st: str.replace(
                    st, "SERVICE_NUMBER", email_data["service_number"]
                ),
                iSkyLIMS_drylab.drylab_config.BODY_RESOLUTION_ACCEPTED,
            )
        )
        body_preparation = list(
            map(
                lambda st: str.replace(st, "USER_NAME", email_data["user_name"]),
                body_preparation,
            )
        )
        body_preparation = list(
            map(
                lambda st: str.replace(st, "STATUS", email_data["status"]),
                body_preparation,
            )
        )
        body_preparation = list(
            map(lambda st: str.replace(st, "DATE", date), body_preparation)
        )
    else:
        body_preparation = list(
            map(
                lambda st: str.replace(
                    st, "SERVICE_NUMBER", email_data["service_number"]
                ),
                iSkyLIMS_drylab.drylab_config.BODY_RESOLUTION_REJECTED,
            )
        )
        body_preparation = list(
            map(
                lambda st: str.replace(st, "USER_NAME", email_data["user_name"]),
                body_preparation,
            )
        )
        body_preparation = list(
            map(
                lambda st: str.replace(st, "STATUS", email_data["status"]),
                body_preparation,
            )
        )
    body_message = "\n".join(body_preparation)
    notification_user = (
        iSkyLIMS_drylab.models.ConfigSetting.objects.filter(
            configuration_name__exact="EMAIL_FOR_NOTIFICATIONS"
        )
        .last()
        .get_configuration_value()
    )
    from_user = notification_user
    to_users = [
        email_data["user_email"],
        email_data["service_owner_email"],
        notification_user,
    ]
    try:
        django.core.mail.send_mail(subject, body_message, from_user, to_users)
    except SMTPException:
        pass
    return


def send_resolution_in_progress_email(email_data):
    """
    Description:
        The function send the service email for resolution in progress to user.
        Functions uses the send_email django core function to send the email
    Input:
        email_data      # Contains the information to include in the email
    Constant:
        SUBJECT_RESOLUTION_QUEUED
        BODY_RESOLUTION_IN_PROGRESS
        USER_EMAIL
    Return:
        None
    """
    subject_tmp = iSkyLIMS_drylab.drylab_config.SUBJECT_RESOLUTION_IN_PROGRESS.copy()
    subject_tmp.insert(1, email_data["resolution_number"])
    subject = " ".join(subject_tmp)
    body_preparation = list(
        map(
            lambda st: str.replace(
                st, "RESOLUTION_NUMBER", email_data["resolution_number"]
            ),
            iSkyLIMS_drylab.drylab_config.BODY_RESOLUTION_IN_PROGRESS,
        )
    )
    body_preparation = list(
        map(
            lambda st: str.replace(st, "USER_NAME", email_data["user_name"]),
            body_preparation,
        )
    )
    body_message = "\n".join(body_preparation)
    notification_user = (
        iSkyLIMS_drylab.models.ConfigSetting.objects.filter(
            configuration_name__exact="EMAIL_FOR_NOTIFICATIONS"
        )
        .last()
        .get_configuration_value()
    )
    from_user = notification_user
    to_users = [email_data["user_email"], notification_user]
    try:
        send_mail(subject, body_message, from_user, to_users)
    except SMTPException:
        pass
    return


def store_resolution_additional_parameter(additional_parameters, resolution_obj):
    """
    Description:
        The function store in database the additional resolution parameters.
    Input:
        additional_parameters       # Contains the list with the additional parameters
        resolution_obj              # resolution instance
    Constants:
        MAPPING_ADDITIONAL_RESOLUTION_PARAMETERS
    Return:
        None
    """

    for additional_parameter in additional_parameters:
        parameter = {}
        parameter["resolution"] = resolution_obj
        for (
            field
        ) in iSkyLIMS_drylab.drylab_config.MAPPING_ADDITIONAL_RESOLUTION_PARAMETERS:
            parameter[field[0]] = additional_parameter[field[1]]
        iSkyLIMS_drylab.models.ResolutionParameters.objects.create_resolution_parameters(
            parameter
        )
