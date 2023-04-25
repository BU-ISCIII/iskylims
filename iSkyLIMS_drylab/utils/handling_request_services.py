# Generic imports
import json
from datetime import datetime
import os
from django.conf import settings
from django.contrib.auth.models import User
from django.core.mail import send_mail
from django_utils.fusioncharts.fusioncharts import FusionCharts
import django_utils.models
import django.contrib.auth.models

# Local imports
import iSkyLIMS_drylab.drylab_config
import iSkyLIMS_drylab.models
import iSkyLIMS_drylab.utils
import iSkyLIMS_core.utils

#### API from Wetlab ######
try:
    from iSkyLIMS_wetlab.utils.api.wetlab_api import *

    wetlab_api_available = True
except:
    wetlab_api_available = False


def add_files_to_service(file_ids, service_obj):
    """
    Description:
        The function update the upload service file with the service instance
    Input:
        file_ids 		# id of the files
        service_obj      # service instance
    Functions:
        update_upload_file_with_service   # located at iSkyLIMS_drylab/utils/handling_multiple_files
    Return:
        True if service id exists
    """
    for file_id in file_ids:
        if file_id != "undefined":
            iSkyLIMS_drylab.utils.handling_multiple_files.update_upload_file_with_service(file_id, service_obj)
    return


def add_requested_samples_in_service(request):
    """
    Description:
        The get the list of samples to be added to the already requested service. it reused the
        same function prepare_form_data_request_service_sequencing, as used for getting the samples to
        request a new service
    Input:
        request

    Functions:
        prepare_form_data_request_service_sequencing		# located at this file
    Return:
        service_data_information
    """
    return prepare_form_data_request_service_sequencing(request)


def check_service_id_exists(service_id):
    """
    Description:
        The function check if service id exists
    Input:
        service_id      # id of the service
    Return:
        True if service id exists
    """
    if iSkyLIMS_drylab.models.Service.objects.filter(pk=service_id).exists():
        return True
    else:
        return False


def create_new_save_sequencing_service_request(request):
    """
    Description:
        The function collect  the data in the user form and create a new instance
        for sequencing request
    Input:
        request      # user data form
    Functions:
        create_service_id			# located at iSkyLIMS_drylab.utils.drylab_common_functions
        increment_service_number	# located at iSkyLIMS_drylab.utils.drylab_common_functions
    Return:
        new_service		# recorded instance of the form
    """
    service_data = {}
    available_service_list = request.POST.getlist("RequestedServices")
    available_service_objs = []
    for av_service in available_service_list:
        available_service_objs.append(get_available_service_obj_from_id(av_service))

    if "requestedForUserid" in request.POST:
        if request.POST["requestedForUserid"] == "":
            request_user = request.user
        else:
            request_user = django.contrib.auth.models.User.objects.get(
                pk__exact=request.POST["requestedForUserid"]
            )
    else:
        request_user = request.user

    if django_utils.models.Profile.objects.filter(profileUserID=request_user).exists():
        try:
            service_data["service_seq_center"] = (
                django_utils.models.Profile.objects.filter(profileUserID=request_user)
                .last()
                .profileCenter.get_center_name()
            )
        except:
            service_data["service_seq_center"] = iSkyLIMS_drylab.drylab_config.INTERNAL_SEQUENCING_UNIT

    service_data["service_notes"] = request.POST["description"]
    service_data["serviceUserId"] = request_user
    service_data["service_request_int"] = iSkyLIMS_drylab.utils.drylab_common_functions.increment_service_number(request_user)

    service_data["service_request_number"] = iSkyLIMS_drylab.utils.drylab_common_functions.create_service_id(
        service_data["service_request_int"], request_user
    )

    # Save the new service
    new_service = iSkyLIMS_drylab.models.Service.objects.create_service(service_data)

    # Save files
    """
    if "uploadfile" in request.FILES:
        path = os.path.join(drylab_config.USER_REQUESTED_SERVICE_FILE_DIRECTORY)
        files = request.FILES.getlist("uploadfile")
        for file in files:
            file_name, full_path_file_name = store_file_from_form(file, path)
            new_file = RequestedServiceFile.objects.create_request_service_file(
                file_data
            )
            # service_data['serviceFile'] = full_path_file_name
    """
    # Save the many-to-many data for the form
    for av_service_obj in available_service_objs:
        new_service.service_available_service.add(av_service_obj)
    return new_service


def create_new_save_counseling_infrastructure_service_request(request):
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
    available_service_list = request.POST.getlist("RequestedServices")
    available_service_objs = []
    for av_service in available_service_list:
        available_service_objs.append(get_available_service_obj_from_id(av_service))
    service_data["service_seq_center"] = iSkyLIMS_drylab.drylab_config.INTERNAL_SEQUENCING_UNIT
    service_data["service_notes"] = request.POST["description"]
    service_data["service_run_specs"] = ""
    service_data["service_sequencing_platform"] = ""
    service_data["serviceFileExt"] = ""
    service_data["service_run_specs"] = ""
    service_data["serviceUserId"] = request.user
    service_data["service_request_int"] = iSkyLIMS_drylab.utils.drylab_common_functions.increment_service_number(request.user.id)
    service_data["service_request_number"] = iSkyLIMS_drylab.utils.drylab_common_functions.create_service_id(
        service_data["service_request_int"], request.user.id
    )
    # Save the new service
    new_service = iSkyLIMS_drylab.models.Service.objects.create_service(service_data)
    # Save the many-to-many data for the form
    for av_service_obj in available_service_objs:
        new_service.service_available_service.add(av_service_obj)

    return new_service


def delete_requested_samples_in_service(sample_list):
    """
    Description:
        The function delete the samples requested in service
    Input:
        sample_list  # list of samples to be deleted
    Return:
        deleted_sample_names
    """
    deleted_sample_names = []
    samples_in_services_objs = iSkyLIMS_drylab.models.RequestedSamplesInServices.objects.filter(
        pk__in=sample_list
    )
    for samples_in_services_obj in samples_in_services_objs:
        deleted_sample_names.append(samples_in_services_obj.get_sample_name())
        samples_in_services_obj.delete()
    return deleted_sample_names


def get_available_children_services_and_id(all_tree_services):
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


def get_available_service_obj_from_id(available_service_id):
    """
    Description:
        The function get the available service obj  from the id
    Input:
        available_service_id  # id of the available service
    Return:
        avail_service_obj
    """
    avail_service_obj = None
    if iSkyLIMS_drylab.models.AvailableService.objects.filter(pk__exact=available_service_id).exists():
        avail_service_obj = iSkyLIMS_drylab.models.AvailableService.objects.filter(
            pk__exact=available_service_id
        ).last()
    return avail_service_obj


def get_data_for_service_confirmation(service_requested):
    """
    Description:
        The function get the data to confirm that service was created
    Input:
        service_requested # service instance
    Functions:
        get_uploaded_files_for_service   # located at iSkyLIMS_drylab.utils.handling_multiple_files
        get_projects_in_requested_samples # located at this file
    Return:
        information
    """
    information = {}
    user = {}
    service_data = {}
    service = iSkyLIMS_drylab.models.Service.objects.filter(
        service_request_number__exact=service_requested
    ).last()
    service_number, center = service.get_service_name_and_center()
    information["service_number"] = service_number
    information["requested_date"] = service.get_service_creation_time()
    information["nodes"] = service.service_available_service.all()
    user["name"] = service.serviceUserId.first_name
    user["surname"] = service.serviceUserId.last_name

    user_id = service.serviceUserId.id
    user["area"] = django_utils.models.Profile.objects.get(profileUserID=user_id).profileArea
    user["center"] = django_utils.models.Profile.objects.get(profileUserID=user_id).profileCenter
    user["phone"] = django_utils.models.Profile.objects.get(profileUserID=user_id).profileExtension
    user["position"] = django_utils.models.Profile.objects.get(profileUserID=user_id).profilePosition
    user["email"] = service.serviceUserId.email
    information["user"] = user
    service_data["projects"] = get_projects_in_requested_samples(service)
    # service_data['platform'] = platform
    # service_data['run_specifications'] = run_specs
    service_data["center"] = center
    service_data["notes"] = service.get_service_user_notes()
    files = iSkyLIMS_drylab.utils.handling_multiple_files.get_uploaded_files_for_service(service)
    if len(files) > 0:
        service_data["file"] = files
    else:
        service_data["file"] = ["Not provided"]
    information["service_data"] = service_data

    return information

    # def create_service_pdf_file (service_request_number, absolute_url):
    """
    Description:
        The function collect the information to create the pdf file
    Input:
        request # contains the session information
    Functions:
        get_data_for_service_confirmation   # located at this file
        create_pdf                          # located at iSkyLIMS_drylab.utils.drylab_common_functions
    Constants:
        OUTPUT_DIR_SERVICE_REQUEST_PDF
    Return:
        pdf_file which contains the full path and name of the pdf file
    """
    """
    information_to_include = get_data_for_service_confirmation(service_request_number)
    pdf_file_name = service_request_number + '.pdf'
    full_path_pdf_file = create_pdf(absolute_url, information_to_include, drylab_config.REQUESTED_CONFIRMATION_SERVICE, pdf_file_name,  drylab_config.OUTPUT_DIR_SERVICE_REQUEST_PDF)
    pdf_file = full_path_pdf_file.replace(settings.BASE_DIR,'')
    return pdf_file
    """


def get_pending_services_information():
    """
    Description:
        The function get the pending service information
    Constants:
        MULTI_LEVEL_PIE_PENDING_MAIN_TEXT
    Functions:
        graphic_3D_pie				# located at utils/graphics.py
        graphic_multi_level_pie		# located at utils/graphics.py
    Return:
        pending_services_details
    """
    # pending_services_details = {}
    stat_services = list(
        iSkyLIMS_drylab.models.ServiceState.objects.filter(show_in_stats=True).values_list(
            "state_display", flat=True
        )
    )
    # fileds_to_show = ["service_approved_date"]
    """
        variable_column = 'service_approved_date'
        search_type = 'contains'
        filter = variable_column + '__' + search_type
        info=members.filter(**{ filter: search_string })
    """
    info_for_stat_service = {}
    number_of_services = {}
    pending_services_per_unit = {}
    pending_services_details = {}
    pending_services_graphics = {}
    for stat_service in stat_services:
        info_for_stat_service[stat_service] = []
        service_objs = iSkyLIMS_drylab.models.Service.objects.filter(
            service_state__state_display__exact=stat_service
        ).order_by("-service_created_date")
        number_of_services[stat_service] = len(service_objs)
        for service_obj in service_objs:
            # fetch service id and name
            service_data = [service_obj.get_service_id()]
            service_data.append(service_obj.get_service_request_number())
            # fetch requested user
            service_data.append(service_obj.get_service_requested_user())
            # get date creation service
            service_data.append(service_obj.get_service_creation_time())
            # fetch latest resolution for the service
            if iSkyLIMS_drylab.models.Resolution.objects.filter(resolution_serviceID=service_obj).exists():
                resolution_obj = iSkyLIMS_drylab.models.Resolution.objects.filter(
                    resolution_serviceID=service_obj
                ).last()
                service_data.append(resolution_obj.get_resolution_number())
                service_data.append(resolution_obj.get_resolution_handler_user())
                service_data.append(resolution_obj.get_resolution_on_queued_date())
                service_data.append(resolution_obj.get_resolution_estimated_date())
            else:
                # set empty values when no resolution are defined for service
                service_data += ["--", "--", "--", "--"]
            info_for_stat_service[stat_service].append(service_data)
            # calculate the number of services per center
            unit_req_serv = service_obj.get_service_request_center_name()
            if unit_req_serv not in pending_services_per_unit:
                pending_services_per_unit[unit_req_serv] = {}
            if stat_service not in pending_services_per_unit[unit_req_serv]:
                pending_services_per_unit[unit_req_serv][stat_service] = 0
            pending_services_per_unit[unit_req_serv][stat_service] += 1

        pending_services_details[stat_service] = info_for_stat_service
    pending_services_details["data"] = info_for_stat_service

    data_source = iSkyLIMS_drylab.utils.graphics.graphic_3D_pie(
        "Number of Pending Services", "", "", "", "fint", number_of_services
    )
    graphic_pending_services = FusionCharts(
        "pie3d", "ex1", "740", "600", "chart-1", "json", data_source
    )
    pending_services_graphics["all_state"] = graphic_pending_services.render()

    data_source = iSkyLIMS_drylab.utils.graphics.graphic_multi_level_pie(
        "Pending Services per Unit",
        iSkyLIMS_drylab.drylab_config.MULTI_LEVEL_PIE_PENDING_TEXT_IN_CHILD_SERVICE,
        iSkyLIMS_drylab.drylab_config.MULTI_LEVEL_PIE_PENDING_MAIN_TEXT,
        "fint",
        iSkyLIMS_drylab.drylab_config.COLORS_MULTI_LEVEL_PIE,
        pending_services_per_unit,
    )
    graphic_unit_pending_services = FusionCharts(
        "multilevelpie", "ex2", "740", "600", "chart-2", "json", data_source
    )
    pending_services_graphics[
        "graphic_pending_unit_services"
    ] = graphic_unit_pending_services.render()
    pending_services_details["graphics"] = pending_services_graphics
    return pending_services_details


def get_user_pending_services_information(user_name):
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
    if iSkyLIMS_drylab.models.Resolution.objects.filter(
        resolution_asigned_user__username__exact=user_name,
        resolution_state__state_value__exact="Recorded",
    ).exists():
        resolution_recorded_objs = iSkyLIMS_drylab.models.Resolution.objects.filter(
            resolution_asigned_user__username__exact=user_name,
            resolution_state__state_value__exact="Recorded",
        ).order_by("-resolution_serviceID")
        for resolution_recorded_obj in resolution_recorded_objs:
            resolution_data = (
                resolution_recorded_obj.get_information_for_pending_resolutions()
            )
            del resolution_data[4]
            res_in_queued.append(resolution_data)
        user_pending_services_details["queued"] = res_in_queued
        user_pending_services_details[
            "heading_in_queued"
        ] = iSkyLIMS_drylab.drylab_config.HEADING_USER_PENDING_SERVICE_QUEUED
    if iSkyLIMS_drylab.models.Resolution.objects.filter(
        resolution_asigned_user__username__exact=user_name,
        resolution_state__state_value__exact="In Progress",
    ).exists():
        resolution_recorded_objs = iSkyLIMS_drylab.models.Resolution.objects.filter(
            resolution_asigned_user__username__exact=user_name,
            resolution_state__state_value__exact="In Progress",
        ).order_by("-resolution_serviceID")
        for resolution_recorded_obj in resolution_recorded_objs:
            resolution_data = (
                resolution_recorded_obj.get_information_for_pending_resolutions()
            )
            del resolution_data[4]
            res_in_progress.append(resolution_data)
        user_pending_services_details["in_progress"] = res_in_progress
        user_pending_services_details[
            "heading_in_progress"
        ] = iSkyLIMS_drylab.drylab_config.HEADING_USER_PENDING_SERVICE_QUEUED
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
    if iSkyLIMS_drylab.models.RequestedSamplesInServices.objects.filter(
        samples_in_service=service_obj
    ).exists():
        project_list = []
        req_sample_objs = iSkyLIMS_drylab.models.RequestedSamplesInServices.objects.filter(
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
    if iSkyLIMS_drylab.models.RequestedSamplesInServices.objects.filter(
        samples_in_service=service_obj
    ).exists():
        run_list = []
        req_sample_objs = iSkyLIMS_drylab.models.RequestedSamplesInServices.objects.filter(
            samples_in_service=service_obj
        )
        for req_sample_obj in req_sample_objs:
            run_list.append(req_sample_obj.get_run_name())
        run_unique_list = list(set(run_list))
    return run_unique_list


def get_service_obj_from_id(service_id):
    """
    Description:
        The function get the  service obj  from the id
    Input:
        service_id  # id of the  service
    Return:
        service_obj
    """
    service_obj = None
    if iSkyLIMS_drylab.models.Service.objects.filter(pk__exact=service_id).exists():
        service_obj = iSkyLIMS_drylab.models.Service.objects.filter(pk__exact=service_id).last()
    return service_obj


def get_requested_services_obj_from_available_service(avail_service_obj):
    """
    Description:
        The function get the requested service instances from the avilable
        service object
    Input:
        avail_service_obj  # instance of the  available service
    Return:
        request_service_objs
    """
    request_service_objs = None
    if iSkyLIMS_drylab.models.Service.objects.filter(service_available_service=avail_service_obj).exists():
        request_service_objs = iSkyLIMS_drylab.models.Service.objects.filter(
            service_available_service=avail_service_obj
        )
    return request_service_objs


def get_service_information(service_id, service_manager):
    """
    Description:
        The function get the  service information, which includes the requested services
        resolutions, deliveries, samples and the allowed actions on the service
    Input:
        service_id  # id of the  service
        service_manager 	# Boolean variable to know if user is service manager or not
    Return:
        display_service_details
    """
    service_obj = get_service_obj_from_id(service_id)
    display_service_details = {}

    # text_for_dates = ['Service Date Creation', 'Approval Service Date', 'Rejected Service Date']
    service_dates = []
    display_service_details["service_name"] = service_obj.get_service_request_number()
    display_service_details["service_id"] = service_id
    # get the list of samples
    if iSkyLIMS_drylab.models.RequestedSamplesInServices.objects.filter(
        samples_in_service=service_obj, only_recorded_sample__exact=False
    ).exists():
        samples_in_service = iSkyLIMS_drylab.models.RequestedSamplesInServices.objects.filter(
            samples_in_service=service_obj, only_recorded_sample__exact=False
        )
        display_service_details["samples_sequenced"] = []
        for sample in samples_in_service:
            display_service_details["samples_sequenced"].append(
                [
                    sample.get_sample_id(),
                    sample.get_sample_name(),
                    sample.get_project_name(),
                    sample.get_run_name(),
                    sample.get_requested_sample_id(),
                ]
            )
    if iSkyLIMS_drylab.models.RequestedSamplesInServices.objects.filter(
        samples_in_service=service_obj, only_recorded_sample__exact=True
    ).exists():
        samples_in_service = iSkyLIMS_drylab.models.RequestedSamplesInServices.objects.filter(
            samples_in_service=service_obj, only_recorded_sample__exact=True
        )
        display_service_details["only_recorded_samples"] = []
        for sample in samples_in_service:
            display_service_details["only_recorded_samples"].append(
                [sample.get_sample_name(), "--", sample.get_requested_sample_id()]
            )

    display_service_details["user_name"] = service_obj.get_service_requested_user()
    user_input_files = iSkyLIMS_drylab.utils.handling_multiple_files.get_uploaded_files_and_file_name_for_service(service_obj)
    if user_input_files:
        display_service_details["file"] = []
        for input_file in user_input_files:
            display_service_details["file"].append(
                [os.path.join(settings.MEDIA_URL, input_file[0]), input_file[1]]
            )
    display_service_details["state"] = service_obj.get_service_state(True)
    display_service_details["service_notes"] = service_obj.get_service_user_notes()
    display_service_details["service_dates"] = zip(
        iSkyLIMS_drylab.drylab_config.HEADING_SERVICE_DATES, service_obj.get_service_dates()
    )
    # get the proposal for the delivery date for the last resolution
    if iSkyLIMS_drylab.models.Resolution.objects.filter(resolution_serviceID=service_obj).exists():
        last_resolution = iSkyLIMS_drylab.models.Resolution.objects.filter(
            resolution_serviceID=service_obj
        ).last()
        display_service_details[
            "resolution_folder"
        ] = last_resolution.get_resolution_full_number()
        resolution_estimated_date = last_resolution.get_resolution_estimated_date()

    # get all services
    display_service_details["nodes"] = service_obj.service_available_service.all()
    display_service_details[
        "children_services"
    ] = get_available_children_services_and_id(display_service_details["nodes"])

    # adding actions fields
    if service_manager:
        display_service_details["service_manager"] = True
        if (
            service_obj.serviceStatus != "rejected"
            or service_obj.serviceStatus != "archived"
        ):
            if iSkyLIMS_drylab.models.Resolution.objects.filter(resolution_serviceID=service_obj).exists():
                ## get informtaion from the defined Resolutions
                resolution_objs = iSkyLIMS_drylab.models.Resolution.objects.filter(
                    resolution_serviceID=service_obj
                )
                display_service_details["resolution_for_progress"] = []
                display_service_details["resolution_for_delivery"] = []
                # display_service_details['resolution_delivered'] = []
                available_services_ids = []
                for resolution_obj in resolution_objs:
                    if resolution_obj.get_state() == "Queued":
                        req_available_services = resolution_obj.get_available_services()
                        if req_available_services != ["None"]:
                            req_available_service_ids = (
                                resolution_obj.get_available_services_ids()
                            )
                            for req_available_service in req_available_service_ids:
                                available_services_ids.append(req_available_service)
                            display_service_details["resolution_for_progress"].append(
                                [
                                    resolution_obj.get_resolution_id(),
                                    resolution_obj.get_resolution_number(),
                                    req_available_services,
                                ]
                            )
                        else:
                            display_service_details["resolution_for_progress"].append(
                                [
                                    resolution_obj.get_resolution_id(),
                                    resolution_obj.get_resolution_number(),
                                    [""],
                                ]
                            )
                    elif resolution_obj.get_state() == "In Progress":
                        req_available_services = resolution_obj.get_available_services()
                        if req_available_services != ["None"]:
                            req_available_service_ids = (
                                resolution_obj.get_available_services_ids()
                            )
                            for req_available_service in req_available_service_ids:
                                available_services_ids.append(req_available_service)
                            display_service_details["resolution_for_delivery"].append(
                                [
                                    resolution_obj.get_resolution_id(),
                                    resolution_obj.get_resolution_number(),
                                    req_available_services,
                                ]
                            )
                        else:
                            display_service_details["resolution_for_delivery"].append(
                                [
                                    resolution_obj.get_resolution_id(),
                                    resolution_obj.get_resolution_number(),
                                    [""],
                                ]
                            )
                    elif resolution_obj.get_state() == "Delivery":
                        display_service_details["resolution_delivered"] = True
                        """
                        # allow to display all requested service because issue #116
                        delivered_services = (
                            resolution_obj.get_available_services_and_ids()
                        )

                        if not "None" in delivered_services:
                            for delivered_service in delivered_services:
                                if (
                                    delivered_service
                                    in display_service_details["children_services"]
                                ):
                                    display_service_details["children_services"].remove(
                                        delivered_service
                                    )
                        """
                if len(available_services_ids) < len(
                    display_service_details["children_services"]
                ):
                    # if len(available_services_ids) > 0 and (len(available_services_ids) < len(display_service_details['children_services'])):
                    display_service_details["add_resolution_action"] = service_id
                    display_service_details["multiple_services"] = True
                    display_service_details["pending_to_add_resolution"] = []
                    for ch_service in display_service_details["children_services"]:
                        if ch_service[0] not in available_services_ids:
                            # available_service_obj = get_available_service_obj_from_id(ch_service[0])
                            display_service_details["pending_to_add_resolution"].append(
                                ch_service
                            )
                if "resolution_delivered" in display_service_details:
                    display_service_details[
                        "all_requested_services"
                    ] = get_available_children_services_and_id(
                        display_service_details["nodes"]
                    )
            else:
                display_service_details["add_resolution_action"] = service_id
                if len(display_service_details["children_services"]) > 1:
                    display_service_details["multiple_services"] = True
                    # if Resolution.objects.filter(resolution_serviceID = service_obj).exists():
                    #    resolutions = Resolution.objects.filter(resolution_serviceID = service_obj)
                    #    for resolution in resolutions:
                    #        pass
                    # else:
                display_service_details["first_resolution"] = True

        if service_obj.get_service_state(None) == "queued":
            resolution_id = (
                iSkyLIMS_drylab.models.Resolution.objects.filter(resolution_serviceID=service_obj).last().id
            )
            display_service_details["add_in_progress_action"] = resolution_id
        if service_obj.get_service_state(None) == "in_progress":
            if iSkyLIMS_drylab.models.Resolution.objects.filter(resolution_serviceID=service_obj).exists():
                resolution_id = (
                    iSkyLIMS_drylab.models.Resolution.objects.filter(resolution_serviceID=service_obj)
                    .last()
                    .id
                )
                display_service_details["add_delivery_action"] = resolution_id

    if iSkyLIMS_drylab.models.Resolution.objects.filter(resolution_serviceID=service_obj).exists():
        resolution_heading = iSkyLIMS_drylab.drylab_config.HEADING_FOR_RESOLUTION_INFORMATION
        resolution_objs = iSkyLIMS_drylab.models.Resolution.objects.filter(
            resolution_serviceID=service_obj
        ).order_by("resolution_state")
        resolution_info = []
        for resolution_obj in resolution_objs:
            resolution_info.append(
                [
                    resolution_obj.get_resolution_number(),
                    list(
                        zip(
                            resolution_heading,
                            resolution_obj.get_resolution_information(),
                        )
                    ),
                ]
            )
        display_service_details["resolutions"] = resolution_info

        delivery_info = []
        for resolution_obj in resolution_objs:
            if iSkyLIMS_drylab.models.Delivery.objects.filter(delivery_resolutionID=resolution_obj).exists():
                delivery = iSkyLIMS_drylab.models.Delivery.objects.filter(
                    delivery_resolutionID=resolution_obj
                ).last()
                delivery_info.append([delivery.get_delivery_information()])
                display_service_details["delivery"] = delivery_info

        display_service_details["piplelines_data"] = {}
        for resolution_obj in resolution_objs:
            pipelines_objs = resolution_obj.resolution_pipelines.all()
            if pipelines_objs:
                for pipelines_obj in pipelines_objs:
                    pipeline_name = pipelines_obj.get_pipeline_name()
                    if not pipeline_name in display_service_details["piplelines_data"]:
                        display_service_details["piplelines_data"][pipeline_name] = []
                    display_service_details["piplelines_data"][pipeline_name].append(
                        [
                            pipelines_obj.get_pipeline_id(),
                            pipeline_name,
                            pipelines_obj.get_pipeline_version(),
                            resolution_obj.get_resolution_number(),
                        ]
                    )

        if len(display_service_details["piplelines_data"]) > 0:
            display_service_details[
                "pipelines_heading"
            ] = iSkyLIMS_drylab.drylab_config.HEADING_PIPELINES_USED_IN_RESOLUTIONS
    created_date = service_obj.get_service_creation_time_no_format()
    delivery_date = service_obj.get_service_delivery_time_no_format()
    dates = []
    if iSkyLIMS_drylab.models.Resolution.objects.filter(resolution_serviceID=service_obj).exists():
        resolution_obj = iSkyLIMS_drylab.models.Resolution.objects.filter(
            resolution_serviceID=service_obj
        ).first()
        in_progress_date = resolution_obj.get_resolution_in_progress_date_no_format()
        if in_progress_date != None:
            time_in_queue = (in_progress_date - created_date).days
            dates.append(["Time in Queue", time_in_queue])
            if delivery_date != None:
                execution_time = (delivery_date - in_progress_date).days
                dates.append(["Execution time", execution_time])
    display_service_details["calculation_dates"] = dates
    # Allow that service could set on hold if state is other than delivery
    if display_service_details["state"].lower() != "delivered":
        display_service_details["allowed_waiting_info"] = "allowed"
    return display_service_details


def prepare_form_data_request_service_sequencing(request):
    """
    Description:
        The function get the information to display in the request sequencing service form
    Input:
        request      # user instance who request the service
    Constants:
        HEADING_SELECT_SAMPLE_IN_SERVICE
        HEADING_SELECT_ONLY_RECORDED_SAMPLE_IN_SERVICE
    Functions:
        get_only_recorded_samples_and_dates		# located at iSkyLIMS_core.utils.handling_samples
        get_user_sharing_list					# located at iSkyLIMS_drylab.utils.drylab_common_functions
        get_defined_username_and_ids   			# located at iSkyLIMS_drylab.utils.drylab_common_functions
        get_runs_projects_samples_and_dates		# located at iSkyLIMS_wetlab.utils.api.wetlab_api
    Return:
        service_data_information
    """
    service_data_information = {}
    # get requestiong sequencing data
    """
    if SequencingPlatform.objects.all().exists():
        service_data_information['platform'] = []
        platform_objs = SequencingPlatform.objects.all()
        for platform_obj in platform_objs:
            service_data_information['platform'].append([platform_obj.get_platform_id(), platform_obj.get_platform_name()])
    if FileExt.objects.all().exists():
        service_data_information['file_extension'] =[]
        file_ext_objs = FileExt.objects.all()
        for file_ext_obj in file_ext_objs:
            service_data_information['file_extension'].append([file_ext_obj.get_file_extension_id(),file_ext_obj.get_file_extension()])
    """

    if iSkyLIMS_drylab.utils.drylab_common_functions.is_service_manager(request):
        service_data_information["users"] = iSkyLIMS_drylab.utils.drylab_common_functions.get_defined_username_and_ids()
    service_data_information["nodes"] = (
        iSkyLIMS_drylab.models.AvailableService.objects.filter(
            avail_service_description__exact="Genomic data analysis"
        )
        .get_descendants(include_self=True)
        .exclude(service_in_use=False)
    )

    if wetlab_api_available:
        ## get samples which have sequencing data in iSkyLIMS
        user_sharing_list = iSkyLIMS_drylab.utils.drylab_common_functions.get_user_sharing_list(request.user)

        service_data_information["samples_data"] = get_runs_projects_samples_and_dates(
            user_sharing_list
        )

        if len(service_data_information["samples_data"]) > 0:
            service_data_information[
                "samples_heading"
            ] = iSkyLIMS_drylab.drylab_config.HEADING_SELECT_SAMPLE_IN_SERVICE

    ## get the samples that are only defined without sequencing data available from iSkyLIMS

    service_data_information[
        "sample_only_recorded"
    ] = iSkyLIMS_core.utils.handling_samples.get_only_recorded_samples_and_dates()
    if len(service_data_information["sample_only_recorded"]) > 0:
        service_data_information[
            "sample_only_recorded_heading"
        ] = iSkyLIMS_drylab.drylab_config.HEADING_SELECT_ONLY_RECORDED_SAMPLE_IN_SERVICE

    return service_data_information


def prepare_form_data_request_counseling_service():
    """
    Description:
        The function get the information to display in the counseling service form
    Input:
        request_user      # user instance who request the service
    Return:
        service_data_information
    """
    service_data_information = {}
    service_data_information["nodes"] = iSkyLIMS_drylab.models.AvailableService.objects.filter(
        avail_service_description__exact="Bioinformatics consulting and training"
    ).get_descendants(include_self=True)
    return service_data_information


def prepare_form_data_request_infrastructure_service():
    """
    Description:
        The function get the information to display in the infrastructure service form
    Input:
        request_user      # user instance who request the service
    Return:
        service_data_information
    """
    service_data_information = {}
    service_data_information["nodes"] = iSkyLIMS_drylab.models.AvailableService.objects.filter(
        avail_service_description__exact="User support"
    ).get_descendants(include_self=True)
    return service_data_information


def send_service_creation_confirmation_email(email_data):
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
    subject_tmp = iSkyLIMS_drylab.drylab_config.SUBJECT_SERVICE_RECORDED.copy()
    subject_tmp.insert(1, email_data["service_number"])
    subject = " ".join(subject_tmp)
    body_preparation = list(
        map(
            lambda st: str.replace(st, "SERVICE_NUMBER", email_data["service_number"]),
            iSkyLIMS_drylab.drylab_config.BODY_SERVICE_RECORDED,
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
    except:
        return iSkyLIMS_drylab.drylab_config.ERROR_UNABLE_TO_SEND_EMAIL
    return "OK"


def set_service_waiting_for_user(service_id):
    """
    Description:
        Function set the service to waiting info for user
    Input:
        service_id
    Output:
        Return service_name or None if no service_id is defined
    """
    service_obj = get_service_obj_from_id(service_id)
    if service_obj is not None:
        service_obj.update_service_state("waiting_information")
        return service_obj.get_service_request_number()
    return None


def stored_samples_for_sequencing_request_service(form_data, new_service):
    """
    Description:
        The function get the samples that were selected and store them on database
    Input:
        form_data      # form with the internal and external samples
        new_service     # service obj
    Functions:
        get_user_projects	# API from iSkyLIMS_wetlab located at file wetlab_api
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
                req_samp_obj = iSkyLIMS_drylab.models.RequestedSamplesInServices.objects.create_request_sample(
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
            ext_samp_obj = iSkyLIMS_drylab.models.RequestedSamplesInServices.objects.create_request_sample(
                data
            )
            requested_sample_list.append(data["sample_name"])

    return requested_sample_list
