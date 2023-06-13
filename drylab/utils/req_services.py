# Generic imports
import json

import django.conf
import django.contrib.auth.models
import django.core.mail

# Local imports
import core.fusioncharts.fusioncharts
import core.utils.samples
import django_utils.models
import drylab.api.serializers
import drylab.config
import drylab.models
import drylab.utils.common
import drylab.utils.graphics
import drylab.utils.multi_files

# API from Wetlab #
try:
    import wetlab.utils.api.wetlab_api

    wetlab_api_available = True
except Exception:
    wetlab_api_available = False


def add_files_to_service(file_ids, service_obj):
    """
    Description:
        The function update the upload service file with the service instance
    Input:
        file_ids 		# id of the files
        service_obj      # service instance
    Functions:
        update_upload_file_with_service   # located at drylab/utils/multi_files
    Return:
        True if service id exists
    """
    for file_id in file_ids:
        if file_id != "undefined":
            drylab.utils.multi_files.update_upload_file_with_service(
                file_id, service_obj
            )
    return


def save_sequencing_service(request):
    """
    Description:
        The function collect  the data in the user form and create a new instance
        for sequencing request
    Input:
        request      # user data form
    Functions:
        create_service_id
        increment_service_number
    Return:
        new_service		# recorded instance of the form
    """
    service_data = {}
    available_service_list = request.POST.getlist("requested_services")
    available_service_objs = []
    for av_service in available_service_list:
        available_service_objs.append(get_available_service_obj(av_service))

    if "user_id_request" in request.POST:
        if request.POST["user_id_request"] == "":
            request_user = request.user
        else:
            request_user = django.contrib.auth.models.User.objects.get(
                pk__exact=request.POST["user_id_request"]
            )
    else:
        request_user = request.user

    if django_utils.models.Profile.objects.filter(
        profile_user_id=request_user
    ).exists():
        try:
            service_data["service_center"] = (
                django_utils.models.Profile.objects.filter(profile_user_id=request_user)
                .last()
                .get_user_center_abbr()
            )
        except Exception:
            service_data["service_center"] = drylab.config.INTERNAL_SEQUENCING_UNIT

    service_data["service_notes"] = request.POST["description"]
    service_data["service_user_id"] = request_user
    service_data["service_request_int"] = drylab.utils.common.increment_service_number(
        request_user
    )

    service_data["service_request_number"] = drylab.utils.common.create_service_id(
        service_data["service_request_int"], request_user
    )

    # Save the new service
    new_service = drylab.models.Service.objects.create_service(service_data)

    # Save the many-to-many data for the form
    for av_service_obj in available_service_objs:
        new_service.service_available_service.add(av_service_obj)
    return new_service


def save_counseling_infrastructure_service(request):
    """
    Description:
        The function collect  the data in the user form and create a new instance
        for sequencing request
    Input:
        request      # user data form
    Return:
        new_service		# recorded instance of the form
    """
    service_data = {}
    available_service_list = request.POST.getlist("requested_services")
    available_service_objs = []
    for av_service in available_service_list:
        available_service_objs.append(get_available_service_obj(av_service))

    if "user_id_request" in request.POST:
        if request.POST["user_id_request"] == "":
            request_user = request.user
        else:
            request_user = django.contrib.auth.models.User.objects.get(
                pk__exact=request.POST["user_id_request"]
            )
    else:
        request_user = request.user

    if django_utils.models.Profile.objects.filter(
        profile_user_id=request_user
    ).exists():
        try:
            service_data["service_center"] = (
                django_utils.models.Profile.objects.filter(profile_user_id=request_user)
                .last()
                .profile_center.get_center_name()
            )
        except Exception:
            service_data["service_center"] = drylab.config.INTERNAL_SEQUENCING_UNIT

    service_data["service_notes"] = request.POST["description"]
    service_data["service_run_specs"] = ""
    service_data["service_sequencing_platform"] = ""
    service_data["serviceFileExt"] = ""
    service_data["service_run_specs"] = ""
    service_data["service_user_id"] = request.user
    service_data["service_request_int"] = drylab.utils.common.increment_service_number(
        request.user.id
    )
    service_data["service_request_number"] = drylab.utils.common.create_service_id(
        service_data["service_request_int"], request.user.id
    )
    # Save the new service
    new_service = drylab.models.Service.objects.create_service(service_data)
    # Save the many-to-many data for the form
    for av_service_obj in available_service_objs:
        new_service.service_available_service.add(av_service_obj)

    return new_service


def delete_samples_in_service(sample_list):
    """
    Description:
        The function delete the samples requested in service
    Input:
        sample_list  # list of samples to be deleted
    Return:
        deleted_sample_names
    """
    deleted_sample_names = []
    samples_in_services_objs = drylab.models.RequestedSamplesInServices.objects.filter(
        pk__in=sample_list
    )
    for samples_in_services_obj in samples_in_services_objs:
        deleted_sample_names.append(samples_in_services_obj.get_sample_name())
        samples_in_services_obj.delete()
    return deleted_sample_names


def get_children_services(all_tree_services):
    """
    Description:
        The function get the children available services from a query of service
    Input:
        all_tree_services  # queryset of available service
    Return:
        children_service
    """
    children_services = []
    for t_services in all_tree_services:
        if t_services.get_children():
            continue
        children_services.append([t_services.id, t_services.get_service_description()])
    return children_services


def get_available_service_obj(available_service_id):
    """
    Description:
        The function get the available service obj  from the id
    Input:
        available_service_id  # id of the available service
    Return:
        avail_service_obj
    """
    avail_service_obj = None
    if drylab.models.AvailableService.objects.filter(
        pk__exact=available_service_id
    ).exists():
        avail_service_obj = drylab.models.AvailableService.objects.filter(
            pk__exact=available_service_id
        ).last()
    return avail_service_obj


def get_available_service_states(add_internal_value=False):
    """Function that returns available services in the database

    Parameters
    ----------
    add_internal_value : bool, optional
        it defines if a tuple or only display string is returned
        , by default False

    Returns
    -------
    state_values : list
        If add_internal_value is set to True returns a list of tuples when for
        each item the first value is the internal string and the second the
        one to display
        If False each item only contains the string to display
    """
    if add_internal_value:
        # include a tuple were the first index is the internal value and the
        # second the one to display
        state_values = list(
            drylab.models.ServiceState.objects.all().values_list(
                "state_value", "state_display"
            )
        )
    else:
        # return only the display values in a list
        state_values = list(
            drylab.models.ServiceState.objects.all().values_list(
                "state_display", flat=True
            )
        )
    return state_values


def get_pending_services_info():
    """
    Description:
        The function get the pending service information
    Constants:
        MULTI_LEVEL_PIE_PENDING_MAIN_TEXT
    Functions:
        graphic_3D_pie
        graphic_multi_level_pie
    Return:
        pending_services_details
    """

    state_services = list(
        drylab.models.ServiceState.objects.filter(show_in_stats=True)
        .exclude(state_value__in=["recorded", "delivered"])
        .values_list("state_display", flat=True)
    )

    info_by_state = {}
    services_number = {}
    pending_services_per_unit = {}
    pending_services_details = {}
    pending_services_graphics = {}
    service_data = []

    service_objs = drylab.models.Service.objects.exclude(
        service_state__state_value__in=["delivered", "rejected", "cancelled"]
    ).order_by("-service_created_date")

    recorded_state = drylab.models.ServiceState.objects.get(
        state_value="recorded"
    ).get_state(to_display=True)
    info_by_state[recorded_state] = []
    services_number[recorded_state] = 0
    for service_obj in service_objs:
        if service_obj.get_state() == "recorded":
            services_number[service_obj.get_state(to_display=True)] += 1
            # fetch service id and name
            service_data = [service_obj.get_service_id()]
            service_data.append(service_obj.get_identifier())
            # fetch requested user
            service_data.append(service_obj.get_user_name())
            # get date creation service
            service_data.append(service_obj.get_creation_date())
            # set empty values when no resolution are defined for service
            service_data += ["--", "--", "--", "--"]
            info_by_state[service_obj.get_state(to_display=True)].append(service_data)

    for state in state_services:
        info_by_state[state] = []
        if drylab.models.Resolution.objects.filter(
            resolution_state__state_display=state
        ).exists():
            resolution_objs = drylab.models.Resolution.objects.filter(
                resolution_state__state_display=state
            )
            services_number[state] = len(resolution_objs)
            for resolution_obj in resolution_objs:
                # fetch service id and name
                service_obj = resolution_obj.get_service_obj()
                service_data = [service_obj.get_service_id()]
                service_data.append(service_obj.get_identifier())
                # fetch requested user
                service_data.append(service_obj.get_user_name())
                # get date creation service
                service_data.append(service_obj.get_creation_date())
                # fetch resolution info
                service_data.append(resolution_obj.get_resolution_number())
                service_data.append(resolution_obj.get_asigned_user())
                service_data.append(resolution_obj.get_on_queued_date())
                service_data.append(resolution_obj.get_resolution_estimated_date())
                info_by_state[state].append(service_data)

                # calculate the number of requests per center
                unit_req_serv = service_obj.get_user_center_name()
                if unit_req_serv not in pending_services_per_unit:
                    pending_services_per_unit[unit_req_serv] = {}
                if state not in pending_services_per_unit[unit_req_serv]:
                    pending_services_per_unit[unit_req_serv][state] = 0
                pending_services_per_unit[unit_req_serv][state] += 1

    pending_services_details["data"] = info_by_state

    data_source = drylab.utils.graphics.graphic_3D_pie(
        "Number of Pending Services", "", "", "", "fint", services_number
    )
    graphic_pending_services = core.fusioncharts.fusioncharts.FusionCharts(
        "pie3d", "ex1", "535", "435", "chart-1", "json", data_source
    )
    pending_services_graphics["all_state"] = graphic_pending_services.render()

    data_source = drylab.utils.graphics.graphic_multi_level_pie(
        "Pending Services per Unit",
        drylab.config.MULTI_LEVEL_PIE_PENDING_TEXT_IN_CHILD_SERVICE,
        drylab.config.MULTI_LEVEL_PIE_PENDING_MAIN_TEXT,
        "fint",
        drylab.config.COLORS_MULTI_LEVEL_PIE,
        pending_services_per_unit,
    )
    graphic_unit_pending_services = core.fusioncharts.fusioncharts.FusionCharts(
        "multilevelpie", "ex2", "535", "435", "chart-2", "json", data_source
    )
    pending_services_graphics[
        "graphic_pending_unit_services"
    ] = graphic_unit_pending_services.render()

    pending_services_details["graphics"] = pending_services_graphics
    return pending_services_details


def get_user_pending_services_info(user_name):
    """
    Description:
        The function get the services that are pending for a user
    Input:
        user_name 		# name of the user to fetch the service infomation
    Constants:
        HEADING_USER_PENDING_SERVICE_QUEUED
    Return:
        pending_services_details
    """
    user_pending_services_details = {}
    res_in_queued, res_in_progress = [], []
    if drylab.models.Resolution.objects.filter(
        resolution_assigned_user__username__exact=user_name,
        resolution_state__state_value__exact="queued",
    ).exists():
        resolution_recorded_objs = drylab.models.Resolution.objects.filter(
            resolution_assigned_user__username__exact=user_name,
            resolution_state__state_value__exact="queued",
        ).order_by("-resolution_service_id")
        for resolution_recorded_obj in resolution_recorded_objs:
            resolution_data = resolution_recorded_obj.get_info_pending_resolutions()
            del resolution_data[4]
            res_in_queued.append(resolution_data)
        user_pending_services_details["queued"] = res_in_queued
        user_pending_services_details[
            "heading_in_queued"
        ] = drylab.config.HEADING_USER_PENDING_SERVICE_QUEUED
    if drylab.models.Resolution.objects.filter(
        resolution_assigned_user__username__exact=user_name,
        resolution_state__state_value__exact="in_progress",
    ).exists():
        resolution_recorded_objs = drylab.models.Resolution.objects.filter(
            resolution_assigned_user__username__exact=user_name,
            resolution_state__state_value__exact="in_progress",
        ).order_by("-resolution_service_id")
        for resolution_recorded_obj in resolution_recorded_objs:
            resolution_data = resolution_recorded_obj.get_info_pending_resolutions()
            del resolution_data[4]
            res_in_progress.append(resolution_data)
        user_pending_services_details["in_progress"] = res_in_progress
        user_pending_services_details[
            "heading_in_progress"
        ] = drylab.config.HEADING_USER_PENDING_SERVICE_QUEUED
    return user_pending_services_details


def get_projects_in_requested_samples(service_obj):
    """
    Description:
        The function get the different projects that are involved in the requested samples
    Input:
        service_obj  # service instance
    Return:
        project_unique_list
    """
    project_unique_list = []
    if drylab.models.RequestedSamplesInServices.objects.filter(
        samples_in_service=service_obj
    ).exists():
        project_list = []
        req_sample_objs = drylab.models.RequestedSamplesInServices.objects.filter(
            samples_in_service=service_obj
        )
        for req_sample_obj in req_sample_objs:
            project_list.append(req_sample_obj.get_project_name())
        project_unique_list = list(set(project_list))
    return project_unique_list


def get_run_in_requested_samples(service_obj):
    """
    Description:
        The function get the different projects that are involved in the requested samples
    Input:
        service_obj  # service instance
    Return:
        run_unique_list
    """
    run_unique_list = []
    if drylab.models.RequestedSamplesInServices.objects.filter(
        samples_in_service=service_obj
    ).exists():
        run_list = []
        req_sample_objs = drylab.models.RequestedSamplesInServices.objects.filter(
            samples_in_service=service_obj
        )
        for req_sample_obj in req_sample_objs:
            run_list.append(req_sample_obj.get_run_name())
        run_unique_list = list(set(run_list))
    return run_unique_list


def get_service_data(request):
    """
    Description:
        The function get the information to display in the request sequencing service form
    Input:
        request      # user instance who request the service
    Constants:
        HEADING_SELECT_SAMPLE_IN_SERVICE
        HEADING_SELECT_ONLY_RECORDED_SAMPLE_IN_SERVICE
    Functions:
        get_only_recorded_samples_and_dates
        get_user_sharing_list
        get_defined_username_and_ids
        get_runs_projects_samples_and_dates
    Return:
        service_data
    """
    service_data = {}

    if drylab.utils.common.is_service_manager(request):
        service_data["users"] = drylab.utils.common.get_defined_username_and_ids()
    service_data["nodes"] = (
        drylab.models.AvailableService.objects.filter(
            avail_service_description__exact="Genomic data analysis"
        )
        .get_descendants(include_self=True)
        .exclude(service_in_use=False)
    )

    if wetlab_api_available:
        # get samples which have sequencing data in iSkyLIMS
        user_sharing_list = drylab.utils.common.get_user_sharing_list(request.user)

        service_data[
            "samples_data"
        ] = wetlab.utils.api.wetlab_api.get_runs_projects_samples_and_dates(
            user_sharing_list
        )

        if len(service_data["samples_data"]) > 0:
            service_data[
                "samples_heading"
            ] = drylab.config.HEADING_SELECT_SAMPLE_IN_SERVICE

    # get the samples that are only defined without sequencing data available from iSkyLIMS

    service_data[
        "sample_only_recorded"
    ] = core.utils.samples.get_only_recorded_samples_and_dates()
    if len(service_data["sample_only_recorded"]) > 0:
        service_data[
            "sample_only_recorded_heading"
        ] = drylab.config.HEADING_SELECT_ONLY_RECORDED_SAMPLE_IN_SERVICE

    return service_data


def get_counseling_service_data(request):
    """
    Description:
        The function get the information to display in the counseling service form
    Input:
        request_user      # user instance who request the service
    Return:
        service_data_info
    """
    service_data = {}
    if drylab.utils.common.is_service_manager(request):
        service_data["users"] = drylab.utils.common.get_defined_username_and_ids()

    service_data["nodes"] = drylab.models.AvailableService.objects.filter(
        avail_service_description__exact="Bioinformatics consulting and training"
    ).get_descendants(include_self=True)
    return service_data


def get_infrastructure_service_data(request):
    """
    Description:
        The function get the information to display in the infrastructure service form
    Input:
        request_user      # user instance who request the service
    Return:
        service_data_info
    """
    service_data = {}
    if drylab.utils.common.is_service_manager(request):
        service_data["users"] = drylab.utils.common.get_defined_username_and_ids()
    service_data["nodes"] = drylab.models.AvailableService.objects.filter(
        avail_service_description__exact="User support"
    ).get_descendants(include_self=True)
    return service_data


def send_service_confirmation_email(email_data):
    """
    Description:
        The function send the service email confirmation to user.
        Functions uses the send_email django core function to send the email
    Input:
        email_data      # Contains the information to include in the email
    Constant:
        SUBJECT_SERVICE_RECORDED
        BODY_SERVICE_RECORDED
        USER_EMAIL
        ERROR_UNABLE_TO_SEND_EMAIL
    Return:
        None
    """
    subject_tmp = drylab.config.SUBJECT_SERVICE_RECORDED.copy()
    subject_tmp.insert(1, email_data["service_number"])
    subject = " ".join(subject_tmp)
    body_preparation = list(
        map(
            lambda st: str.replace(st, "SERVICE_NUMBER", email_data["service_number"]),
            drylab.config.BODY_SERVICE_RECORDED,
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
        drylab.models.ConfigSetting.objects.filter(
            configuration_name__exact="EMAIL_FOR_NOTIFICATIONS"
        )
        .last()
        .get_configuration_value()
    )
    from_user = notification_user
    to_users = [email_data["user_email"], notification_user]
    try:
        django.core.mail.send_mail(subject, body_message, from_user, to_users)
    except Exception:
        return drylab.config.ERROR_UNABLE_TO_SEND_EMAIL
    return "OK"


def save_service_samples(form_data, new_service):
    """
    Description:
        The function get the samples that were selected and store them on database
    Input:
        form_data      # form with the internal and external samples
        new_service     # service obj
    Return:
        display_service
    """
    requested_sample_list = []
    heading = [
        "run_name",
        "run_id",
        "project_name",
        "project_id",
        "sample_name",
        "sample_id",
        "date",
        "sample_path",
    ]
    # get the internals samples
    if "samples_requested" in form_data:
        requested_services_table = json.loads(form_data["samples_requested"])
        for row in requested_services_table:
            if row[-1]:
                data = {}
                for i in range(len(heading)):
                    data[heading[i]] = row[i]
                data["samples_in_service"] = new_service
                data["only_recorded"] = False
                drylab.models.RequestedSamplesInServices.objects.create_request_sample(
                    data
                )
                requested_sample_list.append(data["sample_name"])
    # get external samples
    if "only_recorded_samples" in form_data:
        only_recorded_samples_service = json.loads(form_data["only_recorded_samples"])
        for row in only_recorded_samples_service:
            if not row[-1]:
                continue
            data = {}
            for item in heading:
                data[item] = None
            data["sample_name"] = row[0]
            data["project_name"] = row[1]
            data["sample_id"] = row[5]
            data["samples_in_service"] = new_service
            data["only_recorded"] = True
            drylab.models.RequestedSamplesInServices.objects.create_request_sample(data)
            requested_sample_list.append(data["sample_name"])

    return requested_sample_list
