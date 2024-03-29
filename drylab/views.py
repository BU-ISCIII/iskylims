# -*- coding: utf-8 -*-
# Generic imports
import json
import os
from datetime import date, datetime

import django.contrib.auth.models
from django.conf import settings
from django.contrib.auth.decorators import login_required
from django.core.files.storage import FileSystemStorage
from django.db.models import Prefetch
from django.http import HttpResponse
from django.shortcuts import redirect, render

# Local imports
import core.fusioncharts.fusioncharts
import core.utils.common
import django_utils.models
import drylab.config
import drylab.models
import drylab.utils.common
import drylab.utils.deliveries
import drylab.utils.graphics
import drylab.utils.multi_files
import drylab.utils.pipelines
import drylab.utils.req_services
import drylab.utils.resolutions
import drylab.utils.stats
import drylab.utils.test_conf


@login_required
def index(request):
    service_list = {}
    if drylab.models.Service.objects.filter(
        service_state__state_value__exact="recorded"
    ).exists():
        r_service_objs = drylab.models.Service.objects.filter(
            service_state__state_value__exact="recorded"
        ).order_by("service_created_date")
        service_list["recorded"] = []

        for r_service_obj in r_service_objs:
            s_info = []
            s_info.append(r_service_obj.get_identifier())
            s_info.append(r_service_obj.get_user_name())
            service_list["recorded"].append(s_info)

    if (
        drylab.models.Service.objects.all()
        .exclude(service_state__state_value__exact="delivered")
        .exclude(service_approved_date=None)
        .exists()
    ):
        ongoing_services_objs = (
            drylab.models.Service.objects.all()
            .exclude(service_state__state_value__exact="delivered")
            .exclude(service_approved_date=None)
            .order_by("service_approved_date")
        )
        service_list["ongoing"] = []

        for ongoing_services_obj in ongoing_services_objs:
            s_info = []
            s_info.append(ongoing_services_obj.get_identifier())
            s_info.append(ongoing_services_obj.get_delivery_date())
            service_list["ongoing"].append(s_info)
    org_name = drylab.utils.common.get_configuration_from_database("ORGANIZATION_NAME")

    return render(
        request,
        "drylab/index.html",
        {"service_list": service_list, "organization_name": org_name},
    )


@login_required
def configuration_email(request):
    if request.user.username != "admin":
        return redirect("/wetlab")
    email_conf_data = core.utils.common.get_email_data()
    email_conf_data[
        "EMAIL_ISKYLIMS"
    ] = drylab.utils.common.get_configuration_from_database("EMAIL_FOR_NOTIFICATIONS")
    if request.method == "POST" and (request.POST["action"] == "emailconfiguration"):
        result_email = core.utils.common.send_test_email(request.POST)
        if result_email != "OK":
            email_conf_data = core.utils.common.get_email_data()
            email_conf_data["EMAIL_ISKYLIMS"] = request.POST["EMAIL_ISKYLIMS"]
            email_conf_data["test_email"] = request.POST["test_email"]
            return render(
                request,
                "drylab/configuration_email.html",
                {"ERROR": result_email, "email_conf_data": email_conf_data},
            )
        drylab.utils.common.save_database_configuration_value(
            "EMAIL_FOR_NOTIFICATIONS", request.POST["EMAIL_ISKYLIMS"]
        )
        return render(
            request,
            "drylab/configuration_email.html",
            {"succesful_settings": True},
        )
    return render(
        request,
        "drylab/configuration_email.html",
        {"email_conf_data": email_conf_data},
    )


@login_required
def request_seq_service(request):
    if request.POST and request.FILES:
        if "file" in request.FILES:
            data = drylab.utils.multi_files.get_and_save_service_file(request)
            response = drylab.utils.multi_files.JSONResponse(
                data, mimetype="application/json"
            )
            response["Content-Disposition"] = "inline; filename=files.json"
            return response
        else:
            service_data_info = drylab.utils.req_services.get_service_data(request)
            return render(
                request,
                "drylab/request_seq_service.html",
                {"service_data_info": service_data_info},
            )

    if request.method == "POST" and request.POST["sub_action"] == "create_service":
        # check that at some services have been requested
        if len(request.POST.getlist("requested_services")) == 0:
            service_data_info = drylab.utils.req_services.get_service_data(request)
            error_message = drylab.config.ERROR_NO_SERVICES_ARE_SELECTED
            return render(
                request,
                "drylab/request_seq_service.html",
                {
                    "service_data_info": service_data_info,
                    "error_message": error_message,
                },
            )

        new_service = drylab.utils.req_services.save_sequencing_service(request)
        sample_stored = drylab.utils.req_services.save_service_samples(
            request.POST, new_service
        )
        if "files" in request.POST:
            drylab.utils.req_services.add_files_to_service(
                request.POST.getlist("files"), new_service
            )
        # Send mail to user and drylab notification email
        email_data = {}
        if "user_id_request" in request.POST and request.POST["user_id_request"] != "":
            user_obj = django.contrib.auth.models.User.objects.filter(
                pk__exact=request.POST["user_id_request"]
            ).last()
            email_data["user_name"] = user_obj.username
            email_data["user_email"] = user_obj.email
        else:
            email_data["user_email"] = request.user.email
            email_data["user_name"] = request.user.username
        email_data["service_number"] = new_service.get_identifier()
        email_result = drylab.utils.req_services.send_service_confirmation_email(
            email_data
        )

        confirmation_result = {}
        service_request_number = new_service.get_identifier()
        confirmation_result["text"] = list(
            map(
                lambda st: str.replace(st, "SERVICE_NUMBER", service_request_number),
                drylab.config.CONFIRMATION_TEXT_MESSAGE,
            )
        )

        if len(sample_stored) > 0:
            confirmation_result["samples"] = sample_stored

        if email_result != "OK":
            return render(
                request,
                "drylab/request_seq_service.html",
                {
                    "confirmation_result": confirmation_result,
                    "error_message": email_result,
                },
            )
        else:
            return render(
                request,
                "drylab/request_seq_service.html",
                {"confirmation_result": confirmation_result},
            )
    else:
        service_data_info = drylab.utils.req_services.get_service_data(request)
        return render(
            request,
            "drylab/request_seq_service.html",
            {"service_data_info": service_data_info},
        )


@login_required
def counseling_request(request):
    if request.method == "POST" and request.POST["sub_action"] == "create_service":
        # check that at some services have been requested
        if len(request.POST.getlist("requested_services")) == 0:
            service_data_info = drylab.utils.req_services.get_counseling_service_data(
                request
            )
            error_message = drylab.config.ERROR_NO_SERVICES_ARE_SELECTED
            return render(
                request,
                "drylab/request_counseling_service.html",
                {
                    "service_data_info": service_data_info,
                    "error_message": error_message,
                },
            )

        new_service = drylab.utils.req_services.save_counseling_infrastructure_service(
            request
        )

        email_data = {}
        email_data["user_email"] = request.user.email
        email_data["user_name"] = request.user.username
        email_data["service_number"] = new_service.get_identifier()
        drylab.utils.req_services.send_service_confirmation_email(email_data)
        confirmation_result = {}

        service_request_number = new_service.get_identifier()
        confirmation_result["text"] = list(
            map(
                lambda st: str.replace(st, "SERVICE_NUMBER", service_request_number),
                drylab.config.CONFIRMATION_TEXT_MESSAGE,
            )
        )
        return render(
            request,
            "drylab/request_counseling_service.html",
            {"confirmation_result": confirmation_result},
        )

    else:
        service_data_info = drylab.utils.req_services.get_counseling_service_data(
            request
        )
        return render(
            request,
            "drylab/request_counseling_service.html",
            {"service_data_info": service_data_info},
        )


@login_required
def infrastructure_request(request):
    if request.method == "POST" and request.POST["sub_action"] == "create_service":
        # check that at some services have been requested
        if len(request.POST.getlist("requested_services")) == 0:
            service_data_info = (
                drylab.utils.req_services.get_infrastructure_service_data(request)
            )
            error_message = drylab.config.ERROR_NO_SERVICES_ARE_SELECTED
            return render(
                request,
                "drylab/request_infrastructure.html",
                {
                    "service_data_info": service_data_info,
                    "error_message": error_message,
                },
            )

        new_service = drylab.utils.req_services.save_counseling_infrastructure_service(
            request
        )

        email_data = {}
        email_data["user_email"] = request.user.email
        email_data["user_name"] = request.user.username
        email_data["service_number"] = new_service.get_identifier()
        drylab.utils.req_services.send_service_confirmation_email(email_data)
        confirmation_result = {}

        service_request_number = new_service.get_identifier()
        confirmation_result["text"] = list(
            map(
                lambda st: str.replace(st, "SERVICE_NUMBER", service_request_number),
                drylab.config.CONFIRMATION_TEXT_MESSAGE,
            )
        )

        return render(
            request,
            "drylab/request_infrastructure.html",
            {"confirmation_result": confirmation_result},
        )
    else:
        service_data_info = drylab.utils.req_services.get_infrastructure_service_data(
            request
        )
        return render(
            request,
            "drylab/request_infrastructure.html",
            {"service_data_info": service_data_info},
        )


@login_required
def add_samples_service(request):
    if request.user.is_authenticated:
        if not drylab.utils.common.is_service_manager(request):
            return render(
                request,
                "drylab/add_samples_service.html",
                {"ERROR": [drylab.config.ERROR_USER_NOT_ALLOWED]},
            )
    else:
        # redirect to login webpage
        return redirect("/accounts/login")

    if request.method == "POST" and request.POST["action"] == "add_samples_service":
        if not drylab.models.Service.objects.filter(
            service_request_number__exact=request.POST["service_id"]
        ).exists():
            return render(
                request,
                "drylab/add_samples_service.html",
                {"ERROR": ["The service that you are trying to get does not exist "]},
            )
        service_obj = drylab.utils.common.get_service_obj(
            request.POST["service_id"], input="id"
        )
        samples_added = {}
        samples_added["samples"] = drylab.utils.req_services.save_service_samples(
            request.POST, service_obj
        )
        samples_added["service_id"] = service_obj.get_service_id()
        return render(
            request,
            "drylab/add_samples_service.html",
            {"samples_added": samples_added},
        )
    else:
        service_data_info = drylab.utils.req_services.get_service_data(request)
        service_data_info["service_id"] = request.POST["service_id"]
        return render(
            request,
            "drylab/add_samples_service.html",
            {"service_data_info": service_data_info},
        )


@login_required
def delete_samples_service(request):
    if request.user.is_authenticated:
        if not drylab.utils.common.is_service_manager(request):
            return render(
                request,
                "drylab/delete_samples_service.html",
                {"ERROR": [drylab.config.ERROR_USER_NOT_ALLOWED]},
            )
    else:
        # redirect to login webpage
        return redirect("/accounts/login")

    if request.method == "POST" and request.POST["action"] == "delete_samples_service":
        if not drylab.models.Service.objects.filter(
            service_request_number__exact=request.POST["service_id"]
        ).exists():
            return render(
                request,
                "drylab/delete_samples_service.html",
                {"ERROR": ["The service that you are trying to get does not exist "]},
            )

        if "samples_delete" not in request.POST:
            return redirect(
                "/drylab/display-service="
                + str(
                    drylab.utils.common.get_service_obj(
                        request.POST["service_id"], input="id"
                    ).get_service_id()
                )
            )

        samples_id = []
        samples_id = [
            sample[0]
            for sample in json.loads(request.POST["samples_delete"])
            if sample[4]
        ]

        deleted_samples = drylab.utils.req_services.delete_samples_in_service(
            samples_id
        )
        service_data = {
            "service_id": drylab.utils.common.get_service_obj(
                request.POST["service_id"], input="id"
            ).get_service_id(),
            "service_name": request.POST["service_id"],
        }

        return render(
            request,
            "drylab/delete_samples_service.html",
            {"deleted_samples": deleted_samples, "service_data": service_data},
        )
    else:
        return redirect(
            "/drylab/display-service="
            + str(
                drylab.utils.common.get_service_obj(
                    request.POST["service_id"], input="id"
                ).get_service_id()
            )
        )


@login_required
def display_service(request, service_id):
    if not request.user.is_authenticated:
        # redirect to login webpage
        return redirect("/accounts/login")

    if drylab.models.Service.objects.filter(pk=service_id).exists():
        service_manager = drylab.utils.common.is_service_manager(request)

        display_service_details = {}
        service_obj = (
            drylab.models.Service.objects.prefetch_related(
                Prefetch(
                    "resolutions",
                    queryset=drylab.models.Resolution.objects.all(),
                    to_attr="filtered_resolutions",
                )
            )
            .filter(pk=service_id)
            .last()
        )

        display_service_details = drylab.api.serializers.ServiceSerializer(
            service_obj, context={"output_label": True}
        ).data

        user_input_files = drylab.utils.multi_files.get_uploaded_files(service_obj)
        service_files = {}
        if user_input_files:
            for input_file in user_input_files:
                service_files[input_file[1]] = os.path.join(
                    settings.MEDIA_URL, input_file[0]
                )

        available_services = service_obj.service_available_service.all()
        return render(
            request,
            "drylab/display_service.html",
            {
                "display_service": display_service_details,
                "service_manager": service_manager,
                "service_files": service_files,
                "available_services": available_services,
            },
        )
    else:
        return render(
            request,
            "drylab/display_service.html",
            {
                "ERROR": [
                    "The service that you are trying to get does not exist ",
                    "Contact with your administrator .",
                ]
            },
        )


@login_required
def search_service(request):
    if not request.user.is_authenticated:
        # redirect to login webpage
        return redirect("/accounts/login")
    services_search_list = {}

    center_list_abbr = []
    center_availables = django_utils.models.Center.objects.all().order_by("center_abbr")

    for center in center_availables:
        center_list_abbr.append(center.center_abbr)
    services_search_list["centers"] = center_list_abbr
    services_search_list[
        "states"
    ] = drylab.utils.req_services.get_available_service_states(True)

    if "wetlab" in settings.INSTALLED_APPS:
        services_search_list["wetlab_app"] = True

    if not drylab.utils.common.is_service_manager(request):
        services_search_list["username"] = request.user.username

    if request.method == "POST" and request.POST["action"] == "search_service":
        service_id = request.POST["service_id"]
        service_state = request.POST["service_state"]
        start_date = request.POST["start_date"]
        end_date = request.POST["end_date"]
        center = request.POST["service_center"]

        service_user = request.POST["service_user"]
        assigned_user = request.POST["bioinfo_user"]
        # Data related to samples,and run from wetlab application
        sample_name = (
            request.POST["sample_name"] if "sample_name" in request.POST else ""
        )
        project_name = (
            request.POST["project_name"] if "project_name" in request.POST else ""
        )
        run_name = request.POST["run_name"] if "run_name" in request.POST else ""
        # check if all fields in form are empty
        if (
            service_id == ""
            and service_state == ""
            and start_date == ""
            and end_date == ""
            and center == ""
            and sample_name == ""
            and service_user == ""
            and project_name == ""
            and run_name == ""
            and assigned_user == ""
        ):
            return render(
                request,
                "drylab/search_service.html",
                {"services_search_list": services_search_list},
            )

        # check the right format of start and end date
        if (
            start_date and not drylab.utils.common.check_valid_date_format(start_date)
        ) or (
            end_date != "" and not drylab.utils.common.check_valid_date_format(end_date)
        ):
            error_message = drylab.config.ERROR_INCORRECT_FORMAT_DATE
            return render(
                request,
                "drylab/search_service.html",
                {"services_search_list": services_search_list, "ERROR": error_message},
            )

        if service_id != "":
            # check if the requested service in the form matches exactly with the existing service in DB
            if drylab.models.Service.objects.filter(
                service_request_number__exact=service_id
            ).exists():
                services_found = drylab.models.Service.objects.get(
                    service_request_number__exact=service_id
                )
                redirect_page = "/drylab/display-service=" + str(services_found.id)
                return redirect(redirect_page)
            if drylab.models.Service.objects.filter(
                service_request_number__icontains=service_id
            ).exists():
                services_found = drylab.models.Service.objects.filter(
                    service_request_number__icontains=service_id
                )
            else:
                return render(
                    request,
                    "drylab/search_service.html",
                    {
                        "ERROR": [
                            "No matches have been found for the service number ",
                            service_id,
                        ]
                    },
                )
        else:
            services_found = drylab.models.Service.objects.prefetch_related(
                "service_user_id", "service_state"
            ).all()

        if service_state != "":
            services_found = services_found.filter(
                service_state__state_value__exact=service_state
            )
        if start_date != "" and end_date != "":
            services_found = services_found.filter(
                service_created_date__range=(start_date, end_date)
            )
        if start_date != "" and end_date == "":
            services_found = services_found.filter(service_created_date__gte=start_date)
        if start_date == "" and end_date != "":
            services_found = services_found.filter(service_created_date__lte=end_date)
        if center != "":
            services_found = services_found.filter(
                service_request_number__icontains=center
            )
        if service_user != "":
            services_found = services_found.filter(
                service_user_id__username__iexact=service_user
            )
        if assigned_user != "":
            if not drylab.models.Resolution.objects.filter(
                resolution_assigned_user__username__icontains=service_user
            ).exists():
                error_message = (
                    drylab.config.ERROR_NO_MATCHES_FOUND_FOR_YOUR_SERVICE_SEARCH
                )
                return render(
                    request,
                    "drylab/search_service.html",
                    {
                        "services_search_list": services_search_list,
                        "ERROR": error_message,
                    },
                )
            services_handled_by = []
            services_handled_by_objs = drylab.models.Resolution.objects.filter(
                resolution_assigned_user__username__icontains=service_user
            )
            for services_handled_by_obj in services_handled_by_objs:
                services_handled_by.append(services_handled_by_obj.get_service_name())
            services_found = services_found.filter(
                service_request_number__in=services_handled_by
            )

        if project_name != "" or run_name != "" or sample_name != "":
            samples_in_services = drylab.models.RequestedSamplesInServices.objects.all()
            if project_name != "":
                samples_in_services = (
                    drylab.models.RequestedSamplesInServices.objects.filter(
                        project_name__icontains=project_name,
                        samples_in_service__in=services_found,
                    )
                )
            else:
                samples_in_services = (
                    drylab.models.RequestedSamplesInServices.objects.filter(
                        samples_in_service__in=services_found
                    )
                )
            if run_name != "":
                samples_in_services = samples_in_services.filter(
                    run_name__icontains=run_name
                )
            if sample_name != "":
                samples_in_services = samples_in_services.filter(
                    sample_name__icontains=sample_name
                )
            service_list = []
            for samples_in_service in samples_in_services:
                service_list.append(samples_in_service.samples_in_service.pk)
            services_found = services_found.filter(pk__in=service_list)
        if len(services_found) == 0:
            error_message = drylab.config.ERROR_NO_MATCHES_FOUND_FOR_YOUR_SERVICE_SEARCH
            return render(
                request,
                "drylab/search_service.html",
                {"services_search_list": services_search_list, "ERROR": error_message},
            )

        # If only 1 service mathes the user conditions, then get the user information
        if len(services_found) == 1:
            redirect_page = "/drylab/display-service=" + str(
                services_found[0].get_service_id()
            )
            return redirect(redirect_page)
        else:
            display_multiple_services = {}
            s_list = {}
            for service_item in services_found:
                data = []
                service_id = service_item.get_service_id()
                data.append(service_item.get_identifier())
                data.append(service_item.get_state(to_display=True))
                data.append(service_item.get_dates())
                data.append(service_item.get_user_center_name())
                data.append(
                    drylab.utils.req_services.get_projects_in_requested_samples(
                        service_item
                    )
                )
                s_list[service_id] = [data]
            display_multiple_services["s_list"] = s_list
            return render(
                request,
                "drylab/search_service.html",
                {"display_multiple_services": display_multiple_services},
            )

    return render(
        request,
        "drylab/search_service.html",
        {"services_search_list": services_search_list},
    )


@login_required
def pending_services(request):
    if request.user.is_authenticated:
        if not drylab.utils.common.is_service_manager(request):
            return render(
                request,
                "drylab/pending_services.html",
                {"ERROR": drylab.config.ERROR_USER_NOT_ALLOWED},
            )
    else:
        # redirect to login webpage
        return redirect("/accounts/login")

    pending_services_details = drylab.utils.req_services.get_pending_services_info()
    user_pending_services = drylab.utils.req_services.get_user_pending_services_info(
        request.user.username
    )
    return render(
        request,
        "drylab/pending_services.html",
        {
            "pending_services": pending_services_details,
            "user_pending_services": user_pending_services,
        },
    )


@login_required
def add_on_hold(request):
    if request.user.is_authenticated:
        try:
            groups = django.contrib.auth.models.Group.objects.get(
                name=drylab.config.SERVICE_MANAGER
            )
            if groups not in request.user.groups.all():
                return render(
                    request,
                    "drylab/add_on_hold.html",
                    {
                        "ERROR": [
                            "You do have the enough privileges to see this page ",
                            "Contact with your administrator .",
                        ]
                    },
                )
        except Exception:
            return render(
                request,
                "drylab/add_on_hold.html",
                {
                    "ERROR": [
                        "You do have the enough privileges to see this page ",
                        "Contact with your administrator .",
                    ]
                },
            )
    else:
        # redirect to login webpage
        return redirect("/accounts/login")

    if request.method == "POST" and request.POST["action"] == "add_on_hold":
        resolution_obj = drylab.utils.resolutions.get_resolution_obj(
            request.POST["resolution_id"], input="id"
        )
        if resolution_obj is not None:
            resolution_obj.update_state("on_hold")
            resolution_number = resolution_obj.get_identifier()

        service_obj = resolution_obj.get_service_obj()
        # update the service status and in_porgress date
        if drylab.utils.resolutions.check_allow_service_update(
            resolution_obj, "on_hold"
        ):
            service_obj = service_obj.update_state("on_hold")

        email_data = {}
        email_data["user_email"] = service_obj.get_user_email()
        email_data["user_name"] = service_obj.get_user_name()
        email_data["resolution_number"] = resolution_number
        drylab.utils.resolutions.send_resolution_in_progress_email(email_data)
        on_hold_resolution = {}
        on_hold_resolution["resolution_number"] = resolution_number

        return render(
            request,
            "drylab/add_on_hold.html",
            {"on_hold_resolution": on_hold_resolution},
        )
    else:
        return render(
            request,
            "drylab/add_on_hold.html",
            {"ERROR": ["There's been an unexpected error."]},
        )


@login_required
def add_resolution(request):
    if request.user.is_authenticated:
        if not drylab.utils.common.is_service_manager(request):
            return render(
                request,
                "django_utils/error_page.html",
                {"content": drylab.config.ERROR_USER_NOT_ALLOWED},
            )
    else:
        # redirect to login webpage
        return redirect("/accounts/login")

    if request.method != "POST" or "service_id" not in request.POST:
        return render(
            request,
            "drylab/add_resolution.html",
            {"ERROR": drylab.config.ERROR_SERVICE_ID_NOT_FOUND},
        )

    if request.method == "POST" and request.POST["action"] == "add_resolution":
        resolution_data_form = drylab.utils.resolutions.get_add_resolution_data(
            request.POST
        )
        new_resolution = drylab.utils.resolutions.create_new_resolution(
            resolution_data_form
        )
        if "pipeline_ids" in resolution_data_form:
            drylab.utils.resolutions.add_pipelines_to_resolution(
                new_resolution, resolution_data_form["pipeline_ids"]
            )

        # Send email
        email_data = {}
        email_data["user_email"] = request.user.email
        email_data["user_name"] = request.user.username
        email_data["service_number"] = new_resolution.get_identifier()
        email_data["status"] = resolution_data_form["service_accepted"]
        email_data["date"] = resolution_data_form["resolution_estimated_date"]
        # include the email for the user who requested the service
        email_data["service_owner_email"] = new_resolution.get_service_owner_email()
        drylab.utils.resolutions.send_resolution_creation_email(email_data)
        created_resolution = {}
        created_resolution["resolution_number"] = resolution_data_form[
            "resolution_number"
        ]
        # Display pipeline parameters

        return render(
            request,
            "drylab/add_resolution.html",
            {"created_resolution": created_resolution},
        )

    if request.method == "POST" and request.POST["action"] == "add_resolution_form":
        resolution_form_data = (
            drylab.utils.resolutions.prepare_form_data_add_resolution(request.POST)
        )
        return render(
            request,
            "drylab/add_resolution.html",
            {"resolution_form_data": resolution_form_data},
        )
    else:
        return render(
            request,
            "drylab/add_resolution.html",
            {"error": drylab.config.ERROR_SERVICE_ID_NOT_FOUND},
        )


@login_required
def add_in_progress(request):
    if request.user.is_authenticated:
        try:
            groups = django.contrib.auth.models.Group.objects.get(
                name=drylab.config.SERVICE_MANAGER
            )
            if groups not in request.user.groups.all():
                return render(
                    request,
                    "drylab/add_in_progress.html",
                    {
                        "ERROR": [
                            "You do have the enough privileges to see this page ",
                            "Contact with your administrator .",
                        ]
                    },
                )
        except Exception:
            return render(
                request,
                "drylab/add_in_progress.html",
                {
                    "ERROR": [
                        "You do have the enough privileges to see this page ",
                        "Contact with your administrator .",
                    ]
                },
            )
    else:
        # redirect to login webpage
        return redirect("/accounts/login")

    if request.method == "POST" and request.POST["action"] == "add_in_progress":
        resolution_id = request.POST["resolution_id"]
        if not drylab.utils.resolutions.check_if_resolution_exists(
            resolution_id, input="id"
        ):
            error_message = drylab.config.ERROR_RESOLUTION_DOES_NOT_EXISTS
            return render(
                request, "drylab/add_in_progress.html", {"ERROR": error_message}
            )

        resolution_obj = drylab.utils.resolutions.get_resolution_obj(
            resolution_id, input="id"
        )
        resolution_obj.update_to_in_progress()
        resolution_number = resolution_obj.get_resolution_number()
        service_obj = resolution_obj.get_service_obj()
        # check if services can change to "in progress"
        if drylab.utils.resolutions.check_allow_service_update(
            resolution_obj, "in_progress"
        ):
            # update the service status and in_porgress date
            service_obj = service_obj.update_state("in_progress")

        email_data = {}
        email_data["user_email"] = service_obj.get_user_email()
        email_data["user_name"] = service_obj.get_user_name()
        email_data["resolution_number"] = resolution_number
        drylab.utils.resolutions.send_resolution_in_progress_email(email_data)
        in_progress_resolution = {}
        in_progress_resolution["resolution_number"] = resolution_number
        return render(
            request,
            "drylab/add_in_progress.html",
            {"in_progress_resolution": in_progress_resolution},
        )

    error_message = drylab.config.ERROR_RESOLUTION_DOES_NOT_EXISTS
    return render(request, "drylab/add_in_progress.html", {"ERROR": error_message})


@login_required
def add_delivery(request):
    if request.user.is_authenticated:
        try:
            groups = django.contrib.auth.models.Group.objects.get(
                name=drylab.config.SERVICE_MANAGER
            )
            if groups not in request.user.groups.all():
                return render(
                    request,
                    "drylab/add_delivery.html",
                    {
                        "ERROR": [
                            "You do have the enough privileges to see this page ",
                            "Contact with your administrator .",
                        ]
                    },
                )
        except Exception:
            return render(
                request,
                "drylab/add_delivery.html",
                {
                    "ERROR": [
                        "You do have the enough privileges to see this page ",
                        "Contact with your administrator .",
                    ]
                },
            )
    else:
        # redirect to login webpage
        return redirect("/accounts/login")
    if request.method == "POST" and request.POST["action"] == "add_delivery":
        delivery_data = drylab.utils.deliveries.prepare_delivery_form(
            request.POST["resolution_id"]
        )

        return render(
            request,
            "drylab/add_delivery.html",
            {"delivery_data": delivery_data},
        )

    if request.method == "POST" and request.POST["action"] == "add_delivery_resolution":
        if (
            request.POST["startdate"] != ""
            and not drylab.utils.common.check_valid_date_format(
                request.POST["startdate"]
            )
        ) or (
            request.POST["enddate"] != ""
            and not drylab.utils.common.check_valid_date_format(request.POST["enddate"])
        ):
            delivery_data = drylab.utils.deliveries.prepare_delivery_form(
                request.POST["resolution_id"]
            )
            error_message = drylab.config.ERROR_INCORRECT_FORMAT_DATE
            return render(
                request,
                "drylab/add_delivery.html",
                {"delivery_data": delivery_data, "ERROR": error_message},
            )

        delivery_recorded = drylab.utils.deliveries.store_resolution_delivery(
            request.POST
        )
        resolution_obj = delivery_recorded["delivery_resolution_id"]
        if delivery_recorded is not None:
            email_data = {}
            email_data["user_email"] = request.user.email
            email_data["user_name"] = request.user.username
            email_data["resolution_number"] = delivery_recorded["resolution_number"]
            email_data["service_owner_email"] = resolution_obj.get_service_owner_email()
            drylab.utils.deliveries.send_delivery_service_email(email_data)
            if drylab.utils.resolutions.check_allow_service_update(
                resolution_obj, "delivered"
            ):
                service_obj = resolution_obj.get_service_obj()
                service_obj = service_obj.update_state("delivered")
                service_obj.update_delivered_date(date.today())
            return render(
                request,
                "drylab/add_delivery.html",
                {"delivery_recorded": delivery_recorded},
            )

    return render(
        request,
        "drylab/add_delivery.html",
        {
            "ERROR": [
                "The resolution that you are trying to upadate does not exists ",
                "Contact with your administrator .",
            ]
        },
    )


@login_required
def stats_by_user(request):
    if not drylab.utils.common.is_service_manager:
        return render(
            request,
            "drylab/stats_by_user.html",
            {"error_message": "You do have enough privileges to see this page"},
        )

    user_list = drylab.utils.common.get_users_requested_services()
    if request.method == "POST" and request.POST["action"] == "userStatistics":
        # validate the input data in the form
        user_id = request.POST["userID"]
        start_date = request.POST["start_date"]
        end_date = request.POST["end_date"]

        if start_date != "" and not drylab.utils.common.check_valid_date_format(
            start_date
        ):
            error_message = drylab.config.ERROR_INCORRECT_FORMAT_DATE
            return render(
                request,
                "drylab/stats_by_user.html",
                {"user_list": user_list, "ERROR": error_message},
            )
        if end_date != "":
            if not drylab.utils.common.check_valid_date_format(end_date):
                error_message = drylab.config.ERROR_INCORRECT_FORMAT_DATE
                return render(
                    request,
                    "drylab/stats_by_user.html",
                    {"user_list": user_list, "ERROR": error_message},
                )
        else:
            end_date = date.today().strftime("%Y-%m-%d")

        if not django.contrib.auth.models.User.objects.filter(
            pk__exact=user_id
        ).exists():
            error_message = drylab.config.ERROR_USER_NOT_DEFINED
            return render(
                request,
                "drylab/stats_by_user.html",
                {"user_list": user_list, "ERROR": error_message},
            )

        stats_info = drylab.utils.stats.create_statistics_by_user(
            user_id, start_date, end_date
        )

        if "ERROR" in stats_info:
            return render(
                request,
                "drylab/stats_by_user.html",
                {"user_list": user_list, "ERROR": stats_info["ERROR"]},
            )

        return render(request, "drylab/stats_by_user.html", {"stats_info": stats_info})
    else:
        return render(request, "drylab/stats_by_user.html", {"user_list": user_list})


@login_required
def stats_by_services_request(request):
    if request.user.is_authenticated:
        try:
            groups = django.contrib.auth.models.Group.objects.get(
                name=drylab.config.SERVICE_MANAGER
            )
            if groups not in request.user.groups.all():
                return render(
                    request,
                    "django_utils/error_page.html",
                    {
                        "content": [
                            "You do have the enough privileges to see this page ",
                            "Contact with your administrator .",
                        ]
                    },
                )
        except Exception:
            return render(
                request,
                "django_utils/error_page.html",
                {
                    "content": [
                        "You do have the enough privileges to see this page ",
                        "Contact with your administrator .",
                    ]
                },
            )
    else:
        # redirect to login webpage
        return redirect("/accounts/login")
    if request.method == "POST" and request.POST["action"] == "service_statistics":
        start_date = request.POST["start_date"]
        end_date = request.POST["end_date"]
        if start_date != "" and not drylab.utils.common.check_valid_date_format(
            start_date
        ):
            return render(request, "drylab/stats_services_time.html")
        if end_date != "":
            if not drylab.utils.common.check_valid_date_format(end_date):
                return render(request, "drylab/stats_services_time.html")
        else:
            end_date = date.today().strftime("%Y-%m-%d")

        start_date_format = datetime.strptime(start_date, "%Y-%m-%d")
        end_date_format = datetime.strptime(end_date, "%Y-%m-%d")

        if drylab.models.Service.objects.filter(
            service_created_date__range=(start_date, end_date)
        ).exists():
            services_found = drylab.models.Service.objects.filter(
                service_created_date__range=(start_date, end_date)
            ).order_by("-service_created_date")
            services_stats_info = {}
            # preparing stats for services request by users
            user_services = {}
            for service in services_found:
                user = service.get_user_name()
                if user in user_services:
                    user_services[user] += 1
                else:
                    user_services[user] = 1

            period_of_time_selected = str(
                " For the period between " + start_date + " and " + end_date
            )
            # creating the graphic for requested services
            data_source = drylab.utils.graphics.column_graphic_dict(
                "Requested Services by:",
                period_of_time_selected,
                "User names",
                "Number of Services",
                "fint",
                user_services,
            )
            graphic_requested_services = core.fusioncharts.fusioncharts.FusionCharts(
                "column3d", "ex1", "525", "350", "chart-1", "json", data_source
            )
            services_stats_info[
                "graphic_requested_services_per_user"
            ] = graphic_requested_services.render()

            # preparing stats for status of the services
            status_services = {}
            for service in services_found:
                status = service.get_state(to_display=True)
                if status in status_services:
                    status_services[status] += 1
                else:
                    status_services[status] = 1

            # creating the graphic for status services
            data_source = drylab.utils.graphics.graphic_3D_pie(
                "Status of Requested Services",
                period_of_time_selected,
                "",
                "",
                "fint",
                status_services,
            )
            graphic_status_requested_services = (
                core.fusioncharts.fusioncharts.FusionCharts(
                    "pie3d", "ex2", "525", "350", "chart-2", "json", data_source
                )
            )
            services_stats_info[
                "graphic_status_requested_services"
            ] = graphic_status_requested_services.render()

            # preparing stats for request by Area
            user_area_dict = {}
            for service in services_found:
                user_service_obj = service.get_user_obj()
                if django_utils.models.Profile.objects.filter(
                    profile_user_id=user_service_obj
                ).exists():
                    user_classification_area = django_utils.models.Profile.objects.get(
                        profile_user_id=user_service_obj
                    ).get_clasification_area()
                else:
                    user_classification_area = "No_user_area"

                if user_classification_area in user_area_dict:
                    user_area_dict[user_classification_area] += 1
                else:
                    user_area_dict[user_classification_area] = 1

            # creating the graphic for areas
            data_source = drylab.utils.graphics.column_graphic_dict(
                "Services requested per Area",
                period_of_time_selected,
                "Area ",
                "Number of Services",
                "fint",
                user_area_dict,
            )
            graphic_area_services = core.fusioncharts.fusioncharts.FusionCharts(
                "column3d", "ex3", "600", "350", "chart-3", "json", data_source
            )
            services_stats_info[
                "graphic_area_services"
            ] = graphic_area_services.render()

            # preparing stats for services request by Center
            user_center_dict = {}
            for service in services_found:
                user_service_obj = service.get_user_obj()
                if django_utils.models.Profile.objects.filter(
                    profile_user_id=user_service_obj
                ).exists():
                    user_center = django_utils.models.Profile.objects.get(
                        profile_user_id=user_service_obj
                    ).get_user_center_abbr()
                else:
                    user_center = "Not defined"
                if user_center in user_center_dict:
                    user_center_dict[user_center] += 1
                else:
                    user_center_dict[user_center] = 1
            # creating the graphic for areas
            data_source = drylab.utils.graphics.column_graphic_dict(
                "Services requested per Center",
                period_of_time_selected,
                "Center ",
                "Number of Services",
                "fint",
                user_center_dict,
            )
            graphic_center_services = core.fusioncharts.fusioncharts.FusionCharts(
                "column3d", "ex4", "600", "350", "chart-4", "json", data_source
            )
            services_stats_info[
                "graphic_center_services"
            ] = graphic_center_services.render()

            ################################################
            # Preparing the statistics per period of time
            ################################################
            # calculating the period to be displayed the graphic (per month o per year)
            delta_dates = (end_date_format - start_date_format).days
            if delta_dates > 366:
                period_year_month = "%Y"
            else:
                period_year_month = "%Y_%m"

            # Preparing the statistics for Center on period of time
            user_services_period = {}
            time_values_dict = {}
            for service in services_found:
                user_service_obj = service.get_user_obj()
                date_service = service.service_created_date.strftime(period_year_month)
                if django_utils.models.Profile.objects.filter(
                    profile_user_id=user_service_obj
                ).exists():
                    user_center = django_utils.models.Profile.objects.get(
                        profile_user_id=user_service_obj
                    ).get_user_center_abbr()
                else:
                    user_center = "Not defined"
                if date_service not in time_values_dict:
                    time_values_dict[date_service] = 1
                if user_center in user_services_period:
                    if date_service in user_services_period[user_center]:
                        user_services_period[user_center][date_service] += 1
                    else:
                        user_services_period[user_center][date_service] = 1
                else:
                    user_services_period[user_center] = {}
                    user_services_period[user_center][date_service] = 1
            time_values = []
            for key, values in sorted(time_values_dict.items()):
                time_values.append(key)
            # fill with zero for the centers that have no sevice during some period
            for center, value in user_services_period.items():
                for d_period in time_values:
                    if d_period not in user_services_period[center]:
                        user_services_period[center][d_period] = 0
            data_source = drylab.utils.graphics.column_graphic_per_time(
                "Services requested by center ",
                period_of_time_selected,
                "date",
                "number of services",
                time_values,
                user_services_period,
            )
            graphic_center_services_per_time = (
                core.fusioncharts.fusioncharts.FusionCharts(
                    "mscolumn3d", "ex5", "525", "350", "chart-5", "json", data_source
                )
            )
            services_stats_info[
                "graphic_center_services_per_time"
            ] = graphic_center_services_per_time.render()

            # Preparing the statistics for Area on period of time
            user_area_services_period = {}
            time_values_dict = {}
            for service in services_found:
                user_id = service.get_user_id()
                date_service = service.service_created_date.strftime(period_year_month)
                if django_utils.models.Profile.objects.filter(
                    profile_user_id=user_id
                ).exists():
                    user_area = django_utils.models.Profile.objects.get(
                        profile_user_id=user_id
                    ).profile_area
                else:
                    user_center = "Not defined"
                if date_service not in time_values_dict:
                    time_values_dict[date_service] = 1
                if user_area in user_area_services_period:
                    if date_service in user_area_services_period[user_area]:
                        user_area_services_period[user_area][date_service] += 1
                    else:
                        user_area_services_period[user_area][date_service] = 1
                else:
                    user_area_services_period[user_area] = {}
                    user_area_services_period[user_area][date_service] = 1
            time_values = []
            for key, values in sorted(time_values_dict.items()):
                time_values.append(key)
            # fill with zero for the centers that have no sevice during some period
            for area, value in user_area_services_period.items():
                for d_period in time_values:
                    if d_period not in user_area_services_period[area]:
                        user_area_services_period[area][d_period] = 0

            data_source = drylab.utils.graphics.column_graphic_per_time(
                "Services requested by Area ",
                period_of_time_selected,
                "date",
                "number of services",
                time_values,
                user_area_services_period,
            )
            graphic_area_services_per_time = (
                core.fusioncharts.fusioncharts.FusionCharts(
                    "mscolumn3d", "ex6", "525", "350", "chart-6", "json", data_source
                )
            )
            services_stats_info[
                "graphic_area_services_per_time"
            ] = graphic_area_services_per_time.render()

            services_stats_info["period_time"] = period_of_time_selected

            # statistics on Requested Level 2 Services

            service_dict = {}
            for service in services_found:
                service_request_list = service.service_available_service.filter(level=2)
                for service_requested in service_request_list:
                    service_name = service_requested.avail_service_description
                    if service_name in service_dict:
                        service_dict[service_name] += 1
                    else:
                        service_dict[service_name] = 1

            # creating the graphic for requested services
            data_source = drylab.utils.graphics.column_graphic_dict(
                "Requested Services:", "level 2 ", "", "", "fint", service_dict
            )
            graphic_req_l2_services = core.fusioncharts.fusioncharts.FusionCharts(
                "column3d", "ex7", "800", "375", "chart-7", "json", data_source
            )
            services_stats_info[
                "graphic_req_l2_services"
            ] = graphic_req_l2_services.render()

            # statistics on Requested Level 3 Services

            service_dict = {}
            for service in services_found:
                service_request_list = service.service_available_service.filter(level=3)
                for service_requested in service_request_list:
                    service_name = service_requested.avail_service_description
                    if service_name in service_dict:
                        service_dict[service_name] += 1
                    else:
                        service_dict[service_name] = 1

            # creating the graphic for requested services
            data_source = drylab.utils.graphics.column_graphic_dict(
                "Requested Services:", "level 3 ", "", "", "fint", service_dict
            )
            graphic_req_l3_services = core.fusioncharts.fusioncharts.FusionCharts(
                "column3d", "ex8", "800", "375", "chart-8", "json", data_source
            )
            services_stats_info[
                "graphic_req_l3_services"
            ] = graphic_req_l3_services.render()

            return render(
                request,
                "drylab/stats_services_time.html",
                {"services_stats_info": services_stats_info},
            )

        else:
            return render(
                request,
                "drylab/stats_services_time.html",
                {
                    "error_message": "There are no services created by "
                    + "For the time of period of between: "
                    + start_date
                    + " and "
                    + end_date
                },
            )
    else:
        # form = ByServicesRequest()
        return render(request, "drylab/stats_services_time.html")


@login_required
def configuration_test(request):
    # check user privileges
    if request.user.is_authenticated:
        if not request.user.is_staff or not request.user.is_superuser:
            return render(
                request,
                "django_utils/error_page.html",
                {
                    "content": [
                        "You do have the enough privileges to see this page ",
                        "Contact with your administrator .",
                    ]
                },
            )
    else:
        # redirect to login webpage
        return redirect("/accounts/login")

    if request.method == "POST" and request.POST["action"] == "basicTest":
        test_results = {}
        test_results["basic_checks_ok"] = "OK"
        config_file = os.path.join(settings.BASE_DIR, "drylab", "config.py")
        # check if access to databases are defined
        try:
            drylab.models.Service.objects.all()
            test_results["database_access"] = "OK"
        except Exception:
            test_results["database_access"] = "NOK"

        # check if available services are defined

        list_available_services = drylab.models.AvailableService.objects.all()

        if len(list_available_services) == 0:
            test_results["services"] = ("Available services", "NOK")
        else:
            test_results["services"] = ("Available services", "OK")
        test_results[
            "iSkyLIMS_settings"
        ] = drylab.utils.test_conf.get_iSkyLIMS_settings()
        test_results["config_file"] = drylab.utils.test_conf.get_config_file(
            config_file
        )
        for result in test_results:
            if test_results[result] == "NOK":
                test_results["basic_checks_ok"] = "NOK"
                break

        return render(
            request,
            "drylab/ConfigurationTest.html",
            {"test_results": test_results},
        )
    elif request.method == "POST" and request.POST["action"] == "resolutionTest":
        if "Delete" in request.POST:
            drylab.utils.test_conf.delete_test_service("SRVTEST-IIER001")
            return render(request, "drylab/ConfigurationTest.html")

        resolution_results = {}
        service_requested = "SRVTEST-IIER001"
        (
            resolution_results["CreateService"],
            result,
        ) = drylab.utils.test_conf.create_service_test(service_requested)

        if result == "NOK":
            resolution_results["create_service_ok"] = "NOK"
            return render(
                request,
                "drylab/ConfigurationTest.html",
                {"resolution_results": resolution_results},
            )

        else:
            resolution_results["create_service_ok"] = "OK"
            resolution_number = "SRVTEST-IIER001.1"
            resolution_results[
                "resolution_test"
            ] = drylab.utils.test_conf.create_resolution_test(
                resolution_number, service_requested
            )
            resolution_results["create_resolution_ok"] = "OK"

            resolution_results["completed_ok"] = "OK"
            for result in resolution_results["resolution_test"]:
                if result[1] == "NOK":
                    resolution_results["completed_ok"] = "NOK"
                break

            return render(
                request,
                "drylab/ConfigurationTest.html",
                {"resolution_results": resolution_results},
            )
    else:
        return render(request, "drylab/ConfigurationTest.html")


@login_required
def define_pipeline(request):
    if request.user.is_authenticated:
        if not drylab.utils.common.is_service_manager(request):
            return render(
                request,
                "django_utils/error_page.html",
                {"content": drylab.config.ERROR_USER_NOT_ALLOWED},
            )
    else:
        return redirect("/accounts/login")

    data_pipeline = drylab.utils.pipelines.get_data_form_pipeline()
    if request.method == "POST" and request.POST["action"] == "define_pipeline":
        pipeline_data_form = drylab.utils.pipelines.analyze_input_pipelines(request)
        if drylab.utils.pipelines.pipeline_version_exists(
            request.POST["pipeline_name"], request.POST["pipeline_version"]
        ):
            error_message = drylab.config.ERROR_PIPELINE_ALREADY_EXISTS
            data_pipeline.update(pipeline_data_form)
            return render(
                request,
                "drylab/define_pipeline.html",
                {"data_pipeline": data_pipeline, "error_message": error_message},
            )

        if request.FILES:
            fs = FileSystemStorage(
                location=os.path.join(
                    settings.MEDIA_ROOT,
                    drylab.config.PIPELINE_FILE_DIRECTORY,
                )
            )
            pipeline_data_form["file_name"] = fs.save(
                str(
                    pipeline_data_form["pipeline_name"]
                    + "_"
                    + pipeline_data_form["pipeline_version"]
                ),
                request.FILES["pipeline_file"],
            )
        else:
            pipeline_data_form["file_name"] = ""
        new_pipeline = drylab.models.Pipelines.objects.create_pipeline(
            pipeline_data_form
        )

        if "additional_parameters" in pipeline_data_form:
            drylab.utils.pipelines.store_parameters_pipeline(
                new_pipeline, pipeline_data_form["additional_parameters"]
            )
        defined_service_pipeline = (
            drylab.utils.pipelines.get_defined_pipeline_data_to_display(new_pipeline)
        )

        return render(
            request,
            "drylab/define_pipeline.html",
            {"defined_service_pipeline": defined_service_pipeline},
        )

    return render(
        request,
        "drylab/define_pipeline.html",
        {"data_pipeline": data_pipeline},
    )


@login_required
def manage_pipelines(request):
    if request.user.is_authenticated:
        if not drylab.utils.common.is_service_manager(request):
            return render(
                request,
                "django_utils/error_page.html",
                {"content": drylab.config.ERROR_USER_NOT_ALLOWED},
            )
    else:
        # redirect to login webpage
        return redirect("/accounts/login")
    pipelines_data = drylab.utils.pipelines.get_pipelines_for_manage()
    return render(
        request,
        "drylab/manage_pipelines.html",
        {"pipelines_data": pipelines_data},
    )


@login_required
def detail_pipeline(request, pipeline_id):
    if request.user.is_authenticated:
        if not drylab.utils.common.is_service_manager(request):
            return render(
                request,
                "django_utils/error_page.html",
                {"content": drylab.config.ERROR_USER_NOT_ALLOWED},
            )
    else:
        # redirect to login webpage
        return redirect("/accounts/login")
    detail_pipelines_data = drylab.utils.pipelines.get_detail_pipeline_data(pipeline_id)
    return render(
        request,
        "drylab/detail_pipeline.html",
        {"detail_pipelines_data": detail_pipelines_data},
    )


# Delete upload service file  #
def upload_service_file_delete(request, file_id):
    if request.method == "DELETE":
        if not drylab.utils.multi_files.check_if_file_is_linked_to_service(file_id):
            drylab.utils.multi_files.delete_service_file(file_id)
            response = drylab.utils.multi_files.JSONResponse(
                True, mimetype="application/json"
            )
            response["Content-Disposition"] = "inline; filename=files.json"

            return response
        return HttpResponse(
            content="Not allowed", status=403, content_type="application/json"
        )
    return HttpResponse(
        content="data not found", status=410, content_type="application/json"
    )
