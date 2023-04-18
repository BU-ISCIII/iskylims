import os
from datetime import date, datetime
from django.db import models
from django import forms
from django.contrib.auth.models import User

from django.utils.translation import gettext_lazy as _
from mptt.models import MPTTModel
from mptt.fields import TreeForeignKey, TreeManyToManyField
from django.utils.timezone import now as timezone_now

# from iSkyLIMS_wetlab.models import RunProcess, Projects

from django_utils.models import Profile, Center
from django.contrib.auth.models import User

from iSkyLIMS_core.models import Samples, SequencingPlatform
from iSkyLIMS_drylab import drylab_config

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
    return "iSkyLIMS_drylab/servicesRequest/%s_%s%s" % (
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
    availServiceDescription = models.CharField(_("Available services"), max_length=100)
    parent = TreeForeignKey("self", models.SET_NULL, null=True, blank=True)
    inUse = models.BooleanField(default=True)
    serviceId = models.CharField(max_length=40, null=True, blank=True)
    description = models.CharField(max_length=200, null=True, blank=True)

    class Meta:
        db_table = "drylab_available_service"

    def __str__(self):
        return "%s" % (self.availServiceDescription)

    def get_service_description(self):
        return "%s" % (self.availServiceDescription)

    class Meta:
        ordering = ["tree_id", "lft"]
        verbose_name = "AvailableService"
        verbose_name_plural = "AvailableServices"


class PipelinesManager(models.Manager):
    def create_pipeline(self, data):
        new_pipeline = self.create(
            userName=data["userName"],
            pipelineName=data["pipelineName"],
            pipelineInUse=True,
            pipelineVersion=data["pipelineVersion"],
            pipelineFile=os.path.join(
                drylab_config.PIPELINE_FILE_DIRECTORY, data["filename"]
            ),
            pipelineUrl=data["url"],
            pipelineDescription=data["description"],
        )
        return new_pipeline


class Pipelines(models.Model):
    # availableService = models.ForeignKey(
    # 			AvailableService,
    # 			on_delete = models.CASCADE)
    userName = models.ForeignKey(User, on_delete=models.CASCADE)
    pipelineName = models.CharField(max_length=50)
    pipelineVersion = models.CharField(max_length=10)
    pipelineInUse = models.BooleanField(default=True)
    pipelineFile = models.FileField(
        upload_to=drylab_config.PIPELINE_FILE_DIRECTORY, null=True, blank=True
    )
    pipelineUrl = models.CharField(max_length=200, null=True, blank=True)
    pipelineDescription = models.CharField(max_length=500, null=True, blank=True)
    generated_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        db_table = "drylab_pipelines"

    def __str__(self):
        return "%s_%s" % (self.pipelineName, self.pipelineVersion)

    def get_pipeline_name(self):
        return "%s" % (self.pipelineName)

    def get_pipeline_id(self):
        return "%s" % (self.pk)

    def get_pipeline_version(self):
        return "%s" % (self.pipelineVersion)

    def get_pipeline_additional(self):
        data = []
        data.append(self.userName.username)
        data.append(self.pipelineUrl)
        data.append(self.pipelineFile)
        data.append(self.pipelineDescription)

        return data

    def get_pipeline_basic(self):
        data = []
        data.append(self.pipelineName)
        data.append(self.pipelineVersion)
        data.append(self.generated_at.strftime("%B %d, %Y"))
        data.append(self.pipelineInUse)
        # data.append(self.availableService.get_service_description())
        return data

    def get_pipeline_description(self):
        return "%s" % (self.pipelineDescription)

    def get_pipeline_info(self):
        data = []
        # data.append(self.availableService.get_service_description())
        data.append(self.userName.username)
        data.append(self.pipelineName)
        data.append(self.pipelineVersion)
        data.append(self.generated_at.strftime("%d/ %B/ %Y"))
        data.append(self.pipelineInUse)
        data.append(self.pk)
        return data

    objects = PipelinesManager()


class ParameterPipelineManager(models.Manager):
    def create_pipeline_parameters(self, pipeline_parameters):
        new_parameter_action_pipeline = self.create(
            parameterPipeline=pipeline_parameters["parameterPipeline"],
            parameterName=pipeline_parameters["parameterName"],
            parameterType=pipeline_parameters["parameterType"],
        )
        return new_parameter_action_pipeline


class ParameterPipeline(models.Model):
    parameterPipeline = models.ForeignKey(Pipelines, on_delete=models.CASCADE)
    parameterName = models.CharField(max_length=80)
    parameterValue = models.CharField(max_length=200, null=True, blank=True)
    parameterType = models.CharField(max_length=20, null=True, blank=True)

    class Meta:
        db_table = "drylab_parameter_pipeline"

    def __str__(self):
        return "%s" % (self.parameterName)

    def get_pipeline_parameter_name(self):
        return "%s" % (self.parameterName)

    def get_pipeline_parameter_type(self):
        return "%s" % (self.parameterType)

    def get_pipeline_parameter_value(self):
        return "%s" % (self.parameterValue)

    def get_pipeline_parameters(self):
        data = []
        data.append(self.parameterName)
        data.append(self.parameterValue)
        return data

    objects = ParameterPipelineManager()


class ServiceManager(models.Manager):
    def create_service(self, data):
        serviceFileExt = None
        serviceSequencingPlatform = None
        serviceRunSpecs = ""
        service_state_obj = ServiceState.objects.filter(
            state_value__iexact="recorded"
        ).last()
        new_service = self.create(
            serviceUserId=data["serviceUserId"],
            serviceSequencingPlatform=serviceSequencingPlatform,
            serviceRunSpecs=serviceRunSpecs,
            serviceSeqCenter=data["serviceSeqCenter"],
            serviceRequestNumber=data["serviceRequestNumber"],
            serviceRequestInt=data["serviceRequestInt"],
            service_state=service_state_obj,
            serviceStatus="Recorded",
            serviceNotes=data["serviceNotes"],
        )
        return new_service


class Service(models.Model):
    ## User requesting service:
    # 'serviceUsername' refactored to 'serviceUserid' which shows better its true nature
    #  decision taken to change Foreign Key from 'Profile'  to 'User' until full develop of "user registration"
    serviceUserId = models.ForeignKey(User, on_delete=models.CASCADE, null=True)
    serviceSequencingPlatform = models.ForeignKey(
        SequencingPlatform, on_delete=models.CASCADE, null=True, blank=True
    )
    serviceAvailableService = TreeManyToManyField(
        AvailableService, verbose_name=_("AvailableServices")
    )

    serviceProjectNames = models.ManyToManyField(
        "iSkyLIMS_wetlab.Projects", verbose_name=_("User's projects"), blank=True
    )
    service_state = models.ForeignKey(
        ServiceState,
        verbose_name="Service State",
        on_delete=models.CASCADE,
        null=True,
        blank=True,
    )
    serviceSeqCenter = models.CharField(
        _("Sequencing center"), max_length=50, blank=False, null=True
    )
    serviceRequestNumber = models.CharField(max_length=80, null=True)
    serviceRequestInt = models.CharField(max_length=80, null=True)
    serviceRunSpecs = models.CharField(
        _("Run specifications"), max_length=10, blank=True, null=True
    )
    # serviceFile=models.FileField(_("Service description file"),upload_to=service_files_upload, null=True,blank=True)
    serviceStatus = models.CharField(
        _("Service status"), max_length=15, choices=STATUS_CHOICES
    )
    serviceNotes = models.TextField(
        _("Service Notes"), max_length=2048, null=True, blank=True
    )
    serviceCreatedOnDate = models.DateField(auto_now_add=True, null=True)

    serviceOnApprovedDate = models.DateField(auto_now_add=False, null=True, blank=True)
    serviceOnRejectedDate = models.DateField(auto_now_add=False, null=True, blank=True)

    serviceOnDeliveredDate = models.DateField(auto_now_add=False, null=True, blank=True)

    class Meta:
        db_table = "drylab_service"

    def __str__(self):
        return "%s" % (self.serviceRequestNumber)

    def get_service_name_and_center(self):
        return [self.serviceRequestNumber, self.serviceSeqCenter]

    def get_service_id(self):
        return "%s" % self.pk

    def get_service_dates(self):
        service_dates = []
        service_dates.append(self.serviceCreatedOnDate.strftime("%d %B, %Y"))
        try:
            service_dates.append(self.serviceOnApprovedDate.strftime("%d %B, %Y"))
        except (TypeError, AttributeError):
            service_dates.append("--")
        try:
            service_dates.append(self.serviceOnRejectedDate.strftime("%d %B, %Y"))
        except (TypeError, AttributeError):
            service_dates.append("--")
        try:
            service_dates.append(self.serviceOnDeliveredDate.strftime("%d %B, %Y"))
        except (TypeError, AttributeError):
            service_dates.append("--")
        return service_dates

    def get_stats_information(self):
        stats_information = []
        stats_information.append(self.id)
        stats_information.append(self.serviceRequestNumber)
        stats_information.append(self.service_state.get_state(to_display=True))
        # stats_information.append(self.serviceStatus)

        stats_information.append(self.serviceCreatedOnDate.strftime("%d %B, %Y"))
        if self.serviceOnApprovedDate is None:
            if self.serviceOnRejectedDate is None:
                stats_information.append("--")
            else:
                stats_information.append(
                    self.serviceOnRejectedDate.strftime("%d %B, %Y")
                )
        else:
            stats_information.append(self.serviceOnApprovedDate.strftime("%d %B, %Y"))
        if self.serviceOnDeliveredDate is None:
            stats_information.append("--")
        else:
            stats_information.append(self.serviceOnDeliveredDate.strftime("%d %B, %Y"))

        return stats_information

    def get_service_creation_time(self):
        return self.serviceCreatedOnDate.strftime("%d %B, %Y")

    def get_service_creation_time_no_format(self):
        return self.serviceCreatedOnDate

    def get_service_delivery_time_no_format(self):
        return self.serviceOnDeliveredDate

    def get_delivery_date(self):
        if self.serviceOnDeliveredDate:
            return self.serviceOnDeliveredDate.strftime("%d %B, %Y")
        else:
            return "Not yet defined"

    def get_service_request_integer(self):
        return "%s" % (self.serviceRequestInt)

    def get_service_request_number(self):
        return "%s" % (self.serviceRequestNumber)

    def get_service_requested_user(self):
        if self.serviceUserId != None:
            return "%s" % (self.serviceUserId.username)
        return "Not available"

    def get_user_service_obj(self):
        return self.serviceUserId

    def get_service_request_center_abbr(self):
        return "%s" % (self.serviceUserId.profile.profileCenter.centerAbbr)

    def get_service_state(self, display_type):
        return self.service_state.get_state(to_display=display_type)

    def get_service_request_center_name(self):
        return "%s" % (self.serviceUserId.profile.profileCenter.centerName)

    def get_service_user_notes(self):
        return "%s" % (self.serviceNotes)

    def get_time_to_delivery(self):
        if self.serviceOnDeliveredDate == self.serviceCreatedOnDate:
            return 1
        else:
            number_days, time = str(
                self.serviceOnDeliveredDate - self.serviceCreatedOnDate
            ).split(",")
            number, string_value = number_days.split(" ")
        return number

    def get_user_email(self):
        return "%s" % (self.serviceUserId.email)

    def get_child_services(self):
        children_services = []
        all_services = self.serviceAvailableService.all()
        for service in all_services:
            if not service.get_children().exists():
                children_services.append(
                    [service.pk, service.get_service_description()]
                )
        return children_services

    def update_approved_date(self, date):
        self.serviceOnApprovedDate = date
        self.save()
        return self

    def update_service_state(self, state):
        state_obj = ServiceState.objects.filter(state_value__iexact=state).last()
        self.service_state = state_obj
        self.save()
        return self

    def update_service_approved_date(self, date):
        self.serviceOnApprovedDate = date
        self.save()
        return self

    def update_service_delivered_date(self, date):
        self.serviceOnDeliveredDate = date
        self.save()
        return self

    def update_service_rejected_date(self, date):
        self.serviceOnRejectedDate = date
        self.save()
        return self

    def update_sequencing_platform(self, data):
        self.serviceSequencingPlatform = data
        self.save()
        return self

    objects = ServiceManager()


class RequestedSamplesInServicesManager(models.Manager):
    def create_request_sample(self, data):
        new_req_sample_service = self.create(
            samplesInService=data["samplesInService"],
            runName=data["run_name"],
            runNameKey=data["run_id"],
            projectName=data["project_name"],
            projectKey=data["project_id"],
            sampleName=data["sample_name"],
            sampleKey=data["sample_id"],
            onlyRecordedSample=data["only_recorded"],
            samplePath=data["sample_path"],
        )
        return new_req_sample_service


class RequestedSamplesInServices(models.Model):
    samplesInService = models.ForeignKey(
        Service, on_delete=models.CASCADE, related_name="samples"
    )

    sampleKey = models.CharField(max_length=15, null=True, blank=True)
    sampleName = models.CharField(max_length=50, null=True, blank=True)
    samplePath = models.CharField(max_length=250, null=True, blank=True)
    runNameKey = models.CharField(max_length=15, null=True, blank=True)
    runName = models.CharField(max_length=50, null=True, blank=True)
    projectKey = models.CharField(max_length=15, null=True, blank=True)
    projectName = models.CharField(max_length=50, null=True, blank=True)
    onlyRecordedSample = models.BooleanField(default=False)
    generated_at = models.DateField(auto_now_add=True)

    class Meta:
        db_table = "drylab_request_samples_in_services"

    def __str__(self):
        return "%s" % (self.sampleName)

    def get_requested_sample_id(self):
        return "%s" % (self.pk)

    def get_sample_name(self):
        return "%s" % (self.sampleName)

    def get_sample_id(self):
        return "%s" % (self.sampleKey)

    def get_project_name(self):
        return "%s" % (self.projectName)

    def get_run_name(self):
        return "%s" % (self.runName)

    objects = RequestedSamplesInServicesManager()


class UploadServiceFileManager(models.Manager):
    def create_upload_file(self, data):
        new_upload_file = self.create(
            uploadFile=data["file"], uploadFileName=data["file_name"]
        )
        return new_upload_file


class UploadServiceFile(models.Model):
    uploadService = models.ForeignKey(
        Service, on_delete=models.CASCADE, null=True, blank=True
    )
    uploadFile = models.FileField(
        upload_to=drylab_config.USER_REQUESTED_SERVICE_FILE_DIRECTORY
    )
    uploadFileName = models.CharField(null=True, blank=True, max_length=255)
    uploaded_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        db_table = "drylab_upload_service_file"

    def __str__(self):
        return "%s" % (self.uploadFile)

    def get_uploadFile_full_path_and_name(self):
        return "%s" % (self.uploadFile)

    def get_uploadFile_name(self):
        return "%s" % (self.uploadFileName)

    def get_uploadFile_id(self):
        return "%s" % (self.pk)

    def update_service_id(self, service_obj):
        self.uploadService = service_obj
        self.save()
        return self

    objects = UploadServiceFileManager()


class ResolutionManager(models.Manager):
    def create_resolution(self, resolution_data):
        today = date.today()
        resolutionAsignedUser = User.objects.get(
            pk__exact=resolution_data["resolutionAsignedUser"]
        )
        resolutionServiceID = Service.objects.get(
            pk__exact=resolution_data["service_id"]
        )
        state = ResolutionStates.objects.get(state_value__exact="Queued")
        new_resolution = self.create(
            resolutionServiceID=resolutionServiceID,
            resolutionAsignedUser=resolutionAsignedUser,
            resolutionNumber=resolution_data["resolutionNumber"],
            resolutionEstimatedDate=resolution_data["resolutionEstimatedDate"],
            resolutionOnQueuedDate=date.today(),
            resolutionNotes=resolution_data["resolutionNotes"],
            resolutionFullNumber=resolution_data["resolutionFullNumber"],
            resolutionState=state,
        )
        return new_resolution


class Resolution(models.Model):
    resolutionServiceID = models.ForeignKey(
        Service, on_delete=models.CASCADE, related_name="resolutions"
    )
    resolutionAsignedUser = models.ForeignKey(
        User, related_name="groups+", on_delete=models.CASCADE, null=True, blank=True
    )
    resolutionState = models.ForeignKey(
        ResolutionStates, on_delete=models.CASCADE, null=True, blank=True
    )

    resolutionPipelines = models.ManyToManyField(Pipelines, blank=True)

    availableServices = models.ManyToManyField(AvailableService, blank=True)
    resolutionNumber = models.CharField(
        _("Resolutions name"), max_length=255, null=True
    )
    resolutionEstimatedDate = models.DateField(
        _(" Estimated resolution date"), null=True, blank=False
    )
    resolutionDate = models.DateField(
        _("Resolution date"), auto_now_add=True, blank=True
    )
    resolutionOnQueuedDate = models.DateField(auto_now_add=False, null=True, blank=True)
    resolutionOnInProgressDate = models.DateField(
        auto_now_add=False, null=True, blank=True
    )
    resolutionDeliveryDate = models.DateField(auto_now_add=False, null=True, blank=True)
    resolutionNotes = models.TextField(
        _("Resolution notes"), max_length=1000, null=True, blank=True
    )
    resolutionFullNumber = models.CharField(
        _("Acronym Name"), max_length=255, null=True, blank=True
    )
    resolutionPdfFile = models.FileField(
        upload_to=drylab_config.RESOLUTION_FILES_DIRECTORY, null=True, blank=True
    )

    class Meta:
        db_table = "drylab_resolution"

    def __str__(self):
        return "%s" % (self.resolutionNumber)

    def get_resolution_number(self):
        return "%s" % self.resolutionNumber

    def get_service_request_number(self):
        return "%s" % self.resolutionServiceID.get_service_request_number()

    def get_available_services(self):
        if self.availableServices.all().exists():
            avail_services = self.availableServices.all()
            service_list = []
            for avail_service in avail_services:
                service_list.append(avail_service.get_service_description())
            return service_list
        return ["None"]

    def get_available_services_ids(self):
        if self.availableServices.all().exists():
            avail_services = self.availableServices.all()
            service_ids_list = []
            for avail_service in avail_services:
                service_ids_list.append(avail_service.pk)
            return service_ids_list
        return ["None"]

    def get_available_services_and_ids(self):
        if self.availableServices.all().exists():
            avail_services = self.availableServices.all()
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
        if self.resolutionState != None:
            resolution_info.append(self.resolutionState.get_state())
        else:
            resolution_info.append("Not assigned")
        resolution_info.append(self.resolutionFullNumber)
        if self.resolutionAsignedUser != None:
            resolution_info.append(self.resolutionAsignedUser.username)
        else:
            resolution_info.append("Not assigned")
        if self.resolutionEstimatedDate is not None:
            resolution_info.append(self.resolutionEstimatedDate.strftime("%d %B, %Y"))
        else:
            resolution_info.append("Not defined yet")

        if self.resolutionOnQueuedDate is not None:
            resolution_info.append(self.resolutionOnQueuedDate.strftime("%d %B, %Y"))
        else:
            resolution_info.append("Not defined yet")

        if self.resolutionOnInProgressDate is None:
            resolution_info.append("--")
        else:
            resolution_info.append(
                self.resolutionOnInProgressDate.strftime("%d %B, %Y")
            )

        resolution_info.append(self.resolutionNotes)
        resolution_info.append(self.resolutionPdfFile)

        return resolution_info

    def get_information_for_pending_resolutions(self):
        if self.resolutionOnQueuedDate is None:
            on_queued_date = "Not defined"
        else:
            on_queued_date = self.resolutionOnQueuedDate.strftime("%d %B, %Y")
        if self.resolutionEstimatedDate is None:
            on_estimated_date = "Not defined"
        else:
            on_estimated_date = self.resolutionEstimatedDate.strftime("%d %B, %Y")
        data = []
        data.append(self.resolutionServiceID.get_service_id())
        data.append(self.resolutionServiceID.get_service_request_number())
        data.append(self.resolutionNumber)
        data.append(self.resolutionFullNumber)
        if self.resolutionAsignedUser is None:
            data.append("Not assigned yet")
        else:
            data.append(self.resolutionAsignedUser.username)
        data.append(on_queued_date)
        data.append(on_estimated_date)
        return data

    def get_service_obj(self):
        return self.resolutionServiceID

    def get_service_name(self):
        return "%s" % (self.resolutionServiceID.get_service_request_number())

    def get_state(self):
        if self.resolutionState != None:
            return "%s" % (self.resolutionState.get_state())
        else:
            return "Not assigned"

    def get_resolution_full_number(self):
        return "%s" % (self.resolutionFullNumber)

    def get_resolution_estimated_date(self):
        if self.resolutionEstimatedDate != None:
            return "%s" % (self.resolutionEstimatedDate)
        return "--"

    def get_resolution_on_queued_date(self):
        if self.resolutionOnQueuedDate is not None:
            return self.resolutionOnQueuedDate.strftime("%d %B, %Y")
        else:
            return "--"

    def get_resolution_in_progress_date_no_format(self):
        return self.resolutionOnInProgressDate

    def get_resolution_request_center_abbr(self):
        return "%s" % (
            self.resolutionServiceID.serviceUserId.profile.profileCenter.centerAbbr
        )

    def get_resolution_request_center_abbr(self):
        return "%s" % (
            self.resolutionServiceID.serviceUserId.profile.profileCenter.centerName
        )

    def get_resolution_handler_user(self):
        if self.resolutionAsignedUser is not None:
            return self.resolutionAsignedUser.username
        return "--"

    def get_service_owner_email(self):
        if self.resolutionServiceID != None:
            return self.resolutionServiceID.get_user_email()
        return ""

    def update_resolution_in_progress_date(self):
        today = date.today()
        self.resolutionOnInProgressDate = today
        self.resolutionState = ResolutionStates.objects.get(
            state_value__exact="In Progress"
        )
        self.save()
        return self

    def update_resolution_in_delivered(self):
        today = date.today()
        self.resolutionDeliveryDate = today
        self.resolutionState = ResolutionStates.objects.get(
            state_value__exact="Delivery"
        )
        self.save()
        return self

    def update_resolution_file(self, file):
        self.resolutionPdfFile = file
        self.save()
        return self

    def update_resolution_state(self, state):
        self.resolutionState = ResolutionStates.objects.get(state_value__exact=state)
        self.save()
        return self

    objects = ResolutionManager()


class ResolutionParametersManager(models.Manager):
    def create_resolution_parameters(self, data):
        new_resolution_parameter = self.create(
            resolution=data["resolution"],
            resolutionParameter=data["resolutionParameter"],
            resolutionParamValue=data["resolutionParamValue"],
            resolutionParamNotes=data["resolutionParamNotes"],
        )
        return new_resolution_parameter


class ResolutionParameters(models.Model):
    resolution = models.ForeignKey(Resolution, on_delete=models.CASCADE)
    resolutionParameter = models.CharField(max_length=50)
    resolutionParamValue = models.CharField(max_length=80)
    resolutionParamNotes = models.CharField(max_length=200, null=True, blank=True)
    generated_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        db_table = "drylab_resolution_parameters"

    def __str__(self):
        return "%s" % (self.resolutionParameter)

    def get_resolution_parameter_name(self):
        return "%s" % (self.resolutionParameter)

    objects = ResolutionParametersManager()


class DeliveryManager(models.Manager):
    def create_delivery(self, delivery_data):
        new_delivery = self.create(
            deliveryResolutionID=delivery_data["deliveryResolutionID"],
            deliveryNotes=delivery_data["deliveryNotes"],
            executionStartDate=delivery_data["executionStartDate"],
            executionEndDate=delivery_data["executionEndDate"],
            executionTime=delivery_data["executionTime"],
            permanentUsedSpace=delivery_data["permanentUsedSpace"],
            temporaryUsedSpace=delivery_data["temporaryUsedSpace"],
        )
        return new_delivery


class Delivery(models.Model):
    deliveryResolutionID = models.ForeignKey(
        Resolution, on_delete=models.CASCADE, related_name="delivery"
    )
    pipelinesInDelivery = models.ManyToManyField(Pipelines, blank=True)

    deliveryDate = models.DateField(auto_now_add=True, null=True, blank=True)

    deliveryNotes = models.TextField(max_length=1000, null=True, blank=True)
    executionStartDate = models.DateField(null=True, blank=True)
    executionEndDate = models.DateField(null=True, blank=True)
    executionTime = models.CharField(max_length=80, null=True, blank=True)
    permanentUsedSpace = models.CharField(max_length=80, null=True, blank=True)
    temporaryUsedSpace = models.CharField(max_length=80, null=True, blank=True)

    class Meta:
        db_table = "drylab_delivery"

    def __str__(self):
        return "%s" % (self.deliveryResolutionID)

    def get_delivery_information(self):
        delivery_info = []
        delivery_info.append(self.deliveryResolutionID.resolutionNumber)
        delivery_info.append(self.deliveryDate.strftime("%d %B, %Y"))
        delivery_info.append(self.deliveryNotes)
        return delivery_info

    def get_resolution_obj(self):
        return self.deliveryResolutionID

    objects = DeliveryManager()


class ConfigSettingManager(models.Manager):
    def create_config_setting(self, configuration_name, configuration_value):
        new_config_settings = self.create(
            configurationName=configuration_name, configurationValue=configuration_value
        )
        return new_config_settings


class ConfigSetting(models.Model):
    configurationName = models.CharField(max_length=80)
    configurationValue = models.CharField(max_length=255, null=True, blank=True)
    generated_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        db_table = "drylab_config_setting"

    def __str__(self):
        return "%s" % (self.configurationName)

    def get_configuration_value(self):
        return "%s" % (self.configurationValue)

    def set_configuration_value(self, new_value):
        self.configurationValue = new_value
        self.save()
        return self

    objects = ConfigSettingManager()

