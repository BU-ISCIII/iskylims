# Generic imports
import os
from datetime import date

from django.contrib.auth.models import User
from django.db import models
from django.utils.timezone import now as timezone_now
from django.utils.translation import gettext_lazy as _
from mptt.fields import TreeForeignKey, TreeManyToManyField
from mptt.models import MPTTModel

import core.models
import drylab.config

# STATUS_CHOICES will be deprecated from release 2.4.0 or higher
STATUS_CHOICES = (
    ("recorded", _("Recorded")),
    ("approved", _("Approved")),
    ("rejected", _("Rejected")),
    ("queued", _("Queued")),
    ("in_progress", _("In Progress")),
    ("delivered", _("Delivered")),
    ("archived", _("Archived")),
)


def service_files_upload(instance, filename):
    now = timezone_now()
    filename_base, filename_ext = os.path.splitext(filename)
    return "drylab/servicesRequest/%s_%s%s" % (
        now.strftime("%Y%m%d%H%M%S"),
        filename_base.lower(),
        filename_ext.lower(),
    )


class ServiceState(models.Model):
    state_value = models.CharField(max_length=50)
    state_display = models.CharField(max_length=80, null=True, blank=True)
    description = models.CharField(max_length=255, null=True, blank=True)
    show_in_stats = models.BooleanField(default=False)

    class Meta:
        db_table = "drylab_service_state"

    def __str__(self):
        return "%s" % (self.state_value)

    def get_state(self, to_display=None):
        if to_display:
            return "%s" % (self.state_display)
        else:
            return "%s" % (self.state_value)

    def get_description(self):
        return "%s" % (self.description)


class ResolutionStates(models.Model):
    state_value = models.CharField(max_length=50)
    state_display = models.CharField(max_length=80, null=True, blank=True)
    description = models.CharField(max_length=255, null=True, blank=True)

    class Meta:
        db_table = "drylab_resolution_states"

    def __str__(self):
        return "%s" % (self.state_value)

    def get_state(self, to_display=None):
        if to_display:
            return "%s" % (self.state_display)
        else:
            return "%s" % (self.state_value)

    def get_description(self):
        return "%s" % (self.description)


class AvailableService(MPTTModel):
    avail_service_description = models.CharField(
        _("Available services"), max_length=100
    )
    parent = TreeForeignKey("self", models.SET_NULL, null=True, blank=True)
    service_in_use = models.BooleanField(default=True)
    service_id = models.CharField(max_length=40, null=True, blank=True)
    description = models.CharField(max_length=200, null=True, blank=True)

    class Meta:
        db_table = "drylab_available_service"
        ordering = ["tree_id", "lft"]
        verbose_name = "AvailableService"
        verbose_name_plural = "AvailableServices"

    def __str__(self):
        return "%s" % (self.avail_service_description)

    def get_service_description(self):
        return "%s" % (self.avail_service_description)


class PipelinesManager(models.Manager):
    def create_pipeline(self, data):
        new_pipeline = self.create(
            user_name=data["user_name"],
            pipeline_name=data["pipeline_name"],
            pipeline_in_use=True,
            pipeline_version=data["pipeline_version"],
            pipeline_file=os.path.join(
                drylab.config.PIPELINE_FILE_DIRECTORY, data["filename"]
            ),
            pipeline_url=data["url"],
            pipeline_description=data["description"],
        )
        return new_pipeline


class Pipelines(models.Model):
    # availableService = models.ForeignKey(
    # 			AvailableService,
    # 			on_delete = models.CASCADE)
    user_name = models.ForeignKey(User, on_delete=models.CASCADE)
    pipeline_name = models.CharField(max_length=50)
    pipeline_version = models.CharField(max_length=10)
    pipeline_in_use = models.BooleanField(default=True)
    pipeline_file = models.FileField(
        upload_to=drylab.config.PIPELINE_FILE_DIRECTORY,
        null=True,
        blank=True,
    )
    pipeline_url = models.CharField(max_length=200, null=True, blank=True)
    pipeline_description = models.CharField(max_length=500, null=True, blank=True)
    generated_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        db_table = "drylab_pipelines"

    def __str__(self):
        return "%s_%s" % (self.pipeline_name, self.pipeline_version)

    def get_pipeline_name(self):
        return "%s" % (self.pipeline_name)

    def get_pipeline_id(self):
        return "%s" % (self.pk)

    def get_pipeline_version(self):
        return "%s" % (self.pipeline_version)

    def get_pipeline_additional(self):
        data = []
        data.append(self.user_name.username)
        data.append(self.pipeline_url)
        data.append(self.pipeline_file)
        data.append(self.pipeline_description)

        return data

    def get_pipeline_basic(self):
        data = []
        data.append(self.pipeline_name)
        data.append(self.pipeline_version)
        data.append(self.generated_at.strftime("%B %d, %Y"))
        data.append(self.pipeline_in_use)
        return data

    def get_pipeline_description(self):
        return "%s" % (self.pipeline_description)

    def get_pipeline_info(self):
        data = []
        data.append(self.user_name.username)
        data.append(self.pipeline_name)
        data.append(self.pipeline_version)
        data.append(self.generated_at.strftime("%d/ %B/ %Y"))
        data.append(self.pipeline_in_use)
        data.append(self.pk)
        return data

    objects = PipelinesManager()


class ParameterPipelineManager(models.Manager):
    def create_pipeline_parameters(self, pipeline_parameters):
        new_parameter_action_pipeline = self.create(
            parameter_pipeline=pipeline_parameters["parameter_pipeline"],
            parameter_name=pipeline_parameters["parameter_name"],
            parameter_type=pipeline_parameters["parameter_type"],
        )
        return new_parameter_action_pipeline


class ParameterPipeline(models.Model):
    parameter_pipeline = models.ForeignKey(Pipelines, on_delete=models.CASCADE)
    parameter_name = models.CharField(max_length=80)
    parameter_value = models.CharField(max_length=200, null=True, blank=True)
    parameter_type = models.CharField(max_length=20, null=True, blank=True)

    class Meta:
        db_table = "drylab_parameter_pipeline"

    def __str__(self):
        return "%s" % (self.parameter_name)

    def get_pipeline_parameter_name(self):
        return "%s" % (self.parameter_name)

    def get_pipeline_parameter_type(self):
        return "%s" % (self.parameter_type)

    def get_pipeline_parameter_value(self):
        return "%s" % (self.parameter_value)

    def get_pipeline_parameters(self):
        data = []
        data.append(self.parameter_name)
        data.append(self.parameter_value)
        return data

    objects = ParameterPipelineManager()


class ServiceManager(models.Manager):
    def create_service(self, data):
        service_sequencing_platform = None
        service_run_specs = ""
        service_state_obj = ServiceState.objects.filter(
            state_value__iexact="recorded"
        ).last()
        new_service = self.create(
            service_user_id=data["serviceUserId"],
            service_sequencing_platform=service_sequencing_platform,
            service_run_specs=service_run_specs,
            service_center=data["service_center"],
            service_request_number=data["service_request_number"],
            service_request_int=data["service_request_int"],
            service_state=service_state_obj,
            serviceStatus="Recorded",
            service_notes=data["service_notes"],
        )
        return new_service


class Service(models.Model):
    # User requesting service:
    # 'serviceUsername' refactored to 'serviceUserid' which shows better its true nature
    #  decision taken to change Foreign Key from 'Profile'  to 'User' until full develop of "user registration"
    service_user_id = models.ForeignKey(User, on_delete=models.CASCADE, null=True)
    service_sequencing_platform = models.ForeignKey(
        core.models.SequencingPlatform,
        on_delete=models.CASCADE,
        null=True,
        blank=True,
    )
    service_available_service = TreeManyToManyField(
        AvailableService, verbose_name=_("AvailableServices")
    )

    service_project_names = models.ManyToManyField(
        "wetlab.Projects", verbose_name=_("User's projects"), blank=True
    )
    service_state = models.ForeignKey(
        ServiceState,
        verbose_name="Service State",
        on_delete=models.CASCADE,
        null=True,
        blank=True,
    )
    service_center = models.CharField(
        _("Sequencing center"), max_length=50, blank=False, null=True
    )
    service_request_number = models.CharField(max_length=80, null=True)
    service_request_int = models.CharField(max_length=80, null=True)
    service_run_specs = models.CharField(
        _("Run specifications"), max_length=10, blank=True, null=True
    )
    # Not changed in refactorization because it will be deprecated from next release
    serviceStatus = models.CharField(
        _("Service status"), max_length=15, choices=STATUS_CHOICES
    )
    service_notes = models.TextField(
        _("Service Notes"), max_length=2048, null=True, blank=True
    )
    service_created_date = models.DateField(auto_now_add=True, null=True)

    service_approved_date = models.DateField(auto_now_add=False, null=True, blank=True)
    service_rejected_date = models.DateField(auto_now_add=False, null=True, blank=True)

    service_delivered_date = models.DateField(auto_now_add=False, null=True, blank=True)

    class Meta:
        db_table = "drylab_service"

    def __str__(self):
        return "%s" % (self.service_request_number)

    def get_service_name_and_center(self):
        return [self.service_request_number, self.service_center]

    def get_service_information(self):
        return [self.service_request_number, self.service_center]

    def get_service_id(self):
        return "%s" % self.pk

    def get_service_dates(self):
        service_dates = []
        service_dates.append(self.service_created_date.strftime("%d %B, %Y"))
        try:
            service_dates.append(self.service_approved_date.strftime("%d %B, %Y"))
        except (TypeError, AttributeError):
            service_dates.append("--")
        try:
            service_dates.append(self.service_rejected_date.strftime("%d %B, %Y"))
        except (TypeError, AttributeError):
            service_dates.append("--")
        try:
            service_dates.append(self.service_delivered_date.strftime("%d %B, %Y"))
        except (TypeError, AttributeError):
            service_dates.append("--")
        return service_dates

    def get_stats_information(self):
        stats_information = []
        stats_information.append(self.id)
        stats_information.append(self.service_request_number)
        stats_information.append(self.service_state.get_state(to_display=True))
        # stats_information.append(self.serviceStatus)

        stats_information.append(self.service_created_date.strftime("%d %B, %Y"))
        if self.service_approved_date is None:
            if self.service_rejected_date is None:
                stats_information.append("--")
            else:
                stats_information.append(
                    self.service_rejected_date.strftime("%d %B, %Y")
                )
        else:
            stats_information.append(self.service_approved_date.strftime("%d %B, %Y"))
        if self.service_delivered_date is None:
            stats_information.append("--")
        else:
            stats_information.append(self.service_delivered_date.strftime("%d %B, %Y"))

        return stats_information

    def get_service_creation_time(self):
        return self.service_created_date.strftime("%d %B, %Y")

    def get_service_creation_time_no_format(self):
        return self.service_created_date

    def get_service_delivery_time_no_format(self):
        return self.service_delivered_date

    def get_delivery_date(self):
        if self.service_delivered_date:
            return self.service_delivered_date.strftime("%d %B, %Y")
        else:
            return "Not yet defined"

    def get_service_request_integer(self):
        return "%s" % (self.service_request_int)

    def get_service_request_number(self):
        return "%s" % (self.service_request_number)

    def get_service_requested_user(self):
        if self.service_user_id is not None:
            return "%s" % (self.service_user_id.username)
        return "Not available"

    def get_user_service_obj(self):
        return self.service_user_id

    def get_service_request_center_abbr(self):
        return "%s" % (self.service_user_id.profile.profile_center.center_abbr)

    def get_service_state(self, display_type=False):
        return self.service_state.get_state(to_display=display_type)

    def get_service_request_center_name(self):
        return "%s" % (self.service_user_id.profile.profile_center.center_name)

    def get_service_user_notes(self):
        return "%s" % (self.service_notes)

    def get_time_to_delivery(self):
        if self.service_delivered_date == self.service_created_date:
            return 1
        else:
            number_days, _ = str(
                self.service_delivered_date - self.service_created_date
            ).split(",")
            number, _ = number_days.split(" ")
        return number

    def get_user_email(self):
        return "%s" % (self.service_user_id.email)

    def get_child_services(self):
        children_services = []
        all_services = self.service_available_service.all()
        for service in all_services:
            if not service.get_children().exists():
                children_services.append(
                    [service.pk, service.get_service_description()]
                )
        return children_services

    def update_approved_date(self, date):
        self.service_approved_date = date
        self.save()
        return self

    def update_service_state(self, state):
        state_obj = ServiceState.objects.filter(state_value__iexact=state).last()
        self.service_state = state_obj
        self.save()
        return self

    def update_service_approved_date(self, date):
        self.service_approved_date = date
        self.save()
        return self

    def update_service_delivered_date(self, date):
        self.service_delivered_date = date
        self.save()
        return self

    def update_service_rejected_date(self, date):
        self.service_rejected_date = date
        self.save()
        return self

    def update_sequencing_platform(self, data):
        self.service_sequencing_platform = data
        self.save()
        return self

    objects = ServiceManager()


class RequestedSamplesInServicesManager(models.Manager):
    def create_request_sample(self, data):
        new_req_sample_service = self.create(
            samples_in_service=data["samples_in_service"],
            run_name=data["run_name"],
            run_name_key=data["run_id"],
            project_name=data["project_name"],
            project_key=data["project_id"],
            sample_name=data["sample_name"],
            sample_key=data["sample_id"],
            only_recorded_sample=data["only_recorded"],
            sample_path=data["sample_path"],
        )
        return new_req_sample_service


class RequestedSamplesInServices(models.Model):
    samples_in_service = models.ForeignKey(
        Service, on_delete=models.CASCADE, related_name="samples"
    )

    sample_key = models.CharField(max_length=15, null=True, blank=True)
    sample_name = models.CharField(max_length=50, null=True, blank=True)
    sample_path = models.CharField(max_length=250, null=True, blank=True)
    run_name_key = models.CharField(max_length=15, null=True, blank=True)
    run_name = models.CharField(max_length=50, null=True, blank=True)
    project_key = models.CharField(max_length=15, null=True, blank=True)
    project_name = models.CharField(max_length=50, null=True, blank=True)
    only_recorded_sample = models.BooleanField(default=False)
    generated_at = models.DateField(auto_now_add=True)

    class Meta:
        db_table = "drylab_request_samples_in_services"

    def __str__(self):
        return "%s" % (self.sample_name)

    def get_requested_sample_id(self):
        return "%s" % (self.pk)

    def get_sample_name(self):
        return "%s" % (self.sample_name)

    def get_sample_id(self):
        return "%s" % (self.sample_key)

    def get_project_name(self):
        return "%s" % (self.project_name)

    def get_run_name(self):
        return "%s" % (self.run_name)

    objects = RequestedSamplesInServicesManager()


class UploadServiceFileManager(models.Manager):
    def create_upload_file(self, data):
        new_upload_file = self.create(
            upload_file=data["file"], upload_file_name=data["file_name"]
        )
        return new_upload_file


class UploadServiceFile(models.Model):
    upload_service = models.ForeignKey(
        Service, on_delete=models.CASCADE, null=True, blank=True
    )
    upload_file = models.FileField(
        upload_to=drylab.config.USER_REQUESTED_SERVICE_FILE_DIRECTORY
    )
    upload_file_name = models.CharField(null=True, blank=True, max_length=255)
    uploaded_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        db_table = "drylab_upload_service_file"

    def __str__(self):
        return "%s" % (self.upload_file)

    def get_upload_file_full_path_and_name(self):
        return "%s" % (self.upload_file)

    def get_upload_file_name(self):
        return "%s" % (self.upload_file_name)

    def get_upload_file_id(self):
        return "%s" % (self.pk)

    def update_service_id(self, service_obj):
        self.upload_service = service_obj
        self.save()
        return self

    objects = UploadServiceFileManager()


class ResolutionManager(models.Manager):
    def create_resolution(self, resolution_data):
        resolution_asigned_user = User.objects.get(
            pk__exact=resolution_data["resolution_asigned_user"]
        )
        resolution_service_id = Service.objects.get(
            pk__exact=resolution_data["service_id"]
        )
        state = ResolutionStates.objects.get(state_value__exact="Queued")
        new_resolution = self.create(
            resolution_service_id=resolution_service_id,
            resolution_asigned_user=resolution_asigned_user,
            resolution_number=resolution_data["resolutionNumber"],
            resolution_estimated_date=resolution_data["resolutionEstimatedDate"],
            resolution_queued_date=date.today(),
            resolution_notes=resolution_data["resolutionNotes"],
            resolution_full_number=resolution_data["resolution_full_number"],
            resolution_state=state,
        )
        return new_resolution


class Resolution(models.Model):
    resolution_service_id = models.ForeignKey(
        Service, on_delete=models.CASCADE, related_name="resolutions"
    )
    resolution_asigned_user = models.ForeignKey(
        User, related_name="groups+", on_delete=models.CASCADE, null=True, blank=True
    )
    resolution_state = models.ForeignKey(
        ResolutionStates, on_delete=models.CASCADE, null=True, blank=True
    )

    resolution_pipelines = models.ManyToManyField(Pipelines, blank=True)

    available_services = models.ManyToManyField(AvailableService, blank=True)
    resolution_number = models.CharField(
        _("Resolutions name"), max_length=255, null=True
    )
    resolution_estimated_date = models.DateField(
        _(" Estimated resolution date"), null=True, blank=False
    )
    resolution_date = models.DateField(
        _("Resolution date"), auto_now_add=True, blank=True
    )
    resolution_queued_date = models.DateField(auto_now_add=False, null=True, blank=True)
    resolution_in_progress_date = models.DateField(
        auto_now_add=False, null=True, blank=True
    )
    resolution_delivery_date = models.DateField(
        auto_now_add=False, null=True, blank=True
    )
    resolution_notes = models.TextField(
        _("Resolution notes"), max_length=1000, null=True, blank=True
    )
    resolution_full_number = models.CharField(
        _("Acronym Name"), max_length=255, null=True, blank=True
    )
    resolution_pdf_file = models.FileField(
        upload_to=drylab.config.RESOLUTION_FILES_DIRECTORY, null=True, blank=True
    )

    class Meta:
        db_table = "drylab_resolution"

    def __str__(self):
        return "%s" % (self.resolution_number)

    def get_resolution_number(self):
        return "%s" % self.resolution_number

    def get_service_request_number(self):
        return "%s" % self.resolution_service_id.get_service_request_number()

    def get_available_services(self):
        if self.available_services.all().exists():
            avail_services = self.available_services.all()
            service_list = []
            for avail_service in avail_services:
                service_list.append(avail_service.get_service_description())
            return service_list
        return ["None"]

    def get_available_services_ids(self):
        if self.available_services.all().exists():
            avail_services = self.available_services.all()
            service_ids_list = []
            for avail_service in avail_services:
                service_ids_list.append(avail_service.pk)
            return service_ids_list
        return ["None"]

    def get_available_services_and_ids(self):
        if self.available_services.all().exists():
            avail_services = self.available_services.all()
            service_list = []
            for avail_service in avail_services:
                service_list.append(
                    [avail_service.pk, avail_service.get_service_description()]
                )
            return service_list
        return ["None"]

    def get_resolution_id(self):
        return "%s" % (self.pk)

    def get_resolution_information(self):
        resolution_info = []
        resolution_info.append(self.get_available_services())
        if self.resolution_state is not None:
            resolution_info.append(self.resolution_state.get_state())
        else:
            resolution_info.append("Not assigned")
        resolution_info.append(self.resolution_full_number)
        if self.resolution_asigned_user is not None:
            resolution_info.append(self.resolution_asigned_user.username)
        else:
            resolution_info.append("Not assigned")
        if self.resolution_estimated_date is not None:
            resolution_info.append(self.resolution_estimated_date.strftime("%d %B, %Y"))
        else:
            resolution_info.append("Not defined yet")

        if self.resolution_queued_date is not None:
            resolution_info.append(self.resolution_queued_date.strftime("%d %B, %Y"))
        else:
            resolution_info.append("Not defined yet")
            
        if self.resolution_in_progress_date is None:
            resolution_info.append("--")
        else:
            resolution_info.append(
                self.resolution_in_progress_date.strftime("%d %B, %Y")
            )

        resolution_info.append(self.resolution_notes)
        resolution_info.append(self.resolution_pdf_file)

        return resolution_info

    def get_information_for_pending_resolutions(self):
        if self.resolution_queued_date is None:
            on_queued_date = "Not defined"
        else:
            on_queued_date = self.resolution_queued_date.strftime("%d %B, %Y")
        if self.resolution_estimated_date is None:
            on_estimated_date = "Not defined"
        else:
            on_estimated_date = self.resolution_estimated_date.strftime("%d %B, %Y")
        data = []
        data.append(self.resolution_service_id.get_service_id())
        data.append(self.resolution_service_id.get_service_request_number())
        data.append(self.resolution_number)
        data.append(self.resolution_full_number)
        if self.resolution_asigned_user is None:
            data.append("Not assigned yet")
        else:
            data.append(self.resolution_asigned_user.username)
        data.append(on_queued_date)
        data.append(on_estimated_date)
        return data

    def get_service_obj(self):
        return self.resolution_service_id

    def get_service_name(self):
        return "%s" % (self.resolution_service_id.get_service_request_number())

    def get_state(self):
        if self.resolution_state is not None:
            return "%s" % (self.resolution_state.get_state())
        else:
            return "Not assigned"

    def get_resolution_full_number(self):
        return "%s" % (self.resolution_full_number)

    def get_resolution_estimated_date(self):
        if self.resolution_estimated_date is not None:
            return "%s" % (self.resolution_estimated_date)
        return "--"

    def get_resolution_state(self):
        if self.resolution_state is not None:
            return "%s" % (self.resolution_state.get_resolution_state())
        else:
            return "Not assigned"

    def get_resolution_on_queued_date(self):
        if self.resolution_queued_date is not None:
            return self.resolution_queued_date.strftime("%d %B, %Y")
        else:
            return "--"

    def get_resolution_in_progress_date_no_format(self):
        return self.resolution_in_progress_date

    def get_resolution_request_center_abbr(self):
        return "%s" % (
            self.resolution_service_id.service_user_id.profile.profile_center.center_abbr
        )

    def get_resolution_handler_user(self):
        if self.resolution_asigned_user is not None:
            return self.resolution_asigned_user.username
        return "--"

    def get_service_owner_email(self):
        if self.resolution_service_id is not None:
            return self.resolution_service_id.get_user_email()
        return ""

    def update_resolution_in_progress_date(self):
        today = date.today()
        self.resolution_in_progress_date = today
        self.resolution_state = ResolutionStates.objects.get(
            state_value__exact="In Progress"
        )
        self.save()
        return self

    def update_resolution_in_delivered(self):
        today = date.today()
        self.resolution_delivery_date = today
        self.resolution_state = ResolutionStates.objects.get(
            state_value__exact="Delivery"
        )
        self.save()
        return self

    def update_resolution_file(self, file):
        self.resolution_pdf_file = file
        self.save()
        return self

    def update_resolution_state(self, state):
        self.resolution_state = ResolutionStates.objects.get(state_value__exact=state)
        self.save()
        return self

    objects = ResolutionManager()


class ResolutionParametersManager(models.Manager):
    def create_resolution_parameters(self, data):
        new_resolution_parameter = self.create(
            resolution=data["resolution"],
            resolution_parameter=data["resolution_parameter"],
            resolution_param_value=data["resolution_param_value"],
            resolution_param_notes=data["resolution_param_notes"],
        )
        return new_resolution_parameter


class ResolutionParameters(models.Model):
    resolution = models.ForeignKey(Resolution, on_delete=models.CASCADE)
    resolution_parameter = models.CharField(max_length=50)
    resolution_param_value = models.CharField(max_length=80)
    resolution_param_notes = models.CharField(max_length=200, null=True, blank=True)
    generated_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        db_table = "drylab_resolution_parameters"

    def __str__(self):
        return "%s" % (self.resolution_parameter)

    def get_resolution_parameter_name(self):
        return "%s" % (self.resolution_parameter)

    objects = ResolutionParametersManager()


class DeliveryManager(models.Manager):
    def create_delivery(self, delivery_data):
        new_delivery = self.create(
            delivery_resolution_id=delivery_data["delivery_resolutionID"],
            delivery_notes=delivery_data["deliveryNotes"],
            execution_start_date=delivery_data["executionStartDate"],
            execution_end_date=delivery_data["execution_end_date"],
            execution_time=delivery_data["execution_time"],
            permanent_used_space=delivery_data["permanent_used_space"],
            temporary_used_space=delivery_data["temporary_used_space"],
        )
        return new_delivery


class Delivery(models.Model):
    delivery_resolution_id = models.ForeignKey(
        Resolution, on_delete=models.CASCADE, related_name="delivery"
    )
    pipelines_in_delivery = models.ManyToManyField(Pipelines, blank=True)
    delivery_notes = models.TextField(max_length=1000, null=True, blank=True)
    execution_start_date = models.DateField(null=True, blank=True)
    execution_end_date = models.DateField(null=True, blank=True)
    execution_time = models.CharField(max_length=80, null=True, blank=True)
    permanent_used_space = models.CharField(max_length=80, null=True, blank=True)
    temporary_used_space = models.CharField(max_length=80, null=True, blank=True)

    class Meta:
        db_table = "drylab_delivery"

    def __str__(self):
        return "%s" % (self.delivery_resolution_id)

    def get_delivery_information(self):
        delivery_info = []
        delivery_info.append(self.delivery_resolution_id.resolution_number)
        delivery_info.append(self.delivery_resolution_id.resolution_delivery_date.strftime("%d %B, %Y"))
        delivery_info.append(self.delivery_notes)
        return delivery_info

    def get_resolution_obj(self):
        return self.delivery_resolution_id

    objects = DeliveryManager()


class ConfigSettingManager(models.Manager):
    def create_config_setting(self, configuration_name, configuration_value):
        new_config_settings = self.create(
            configuration_name=configuration_name,
            configuration_value=configuration_value,
        )
        return new_config_settings


class ConfigSetting(models.Model):
    configuration_name = models.CharField(max_length=80)
    configuration_value = models.CharField(max_length=255, null=True, blank=True)
    generated_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        db_table = "drylab_config_setting"

    def __str__(self):
        return "%s" % (self.configuration_name)

    def get_configuration_value(self):
        return "%s" % (self.configuration_value)

    def set_configuration_value(self, new_value):
        self.configuration_value = new_value
        self.save()
        return self

    objects = ConfigSettingManager()
