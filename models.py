import datetime, os
from django.db import models
from django import forms
from django.contrib.auth.models import User
from django.utils.encoding import python_2_unicode_compatible
from django.utils.translation import ugettext_lazy as _
from mptt.models import MPTTModel
from mptt.fields import TreeForeignKey, TreeManyToManyField
from django.utils.timezone import now as timezone_now

#from iSkyLIMS_wetlab.models import RunProcess, Projects

from django_utils.models import Profile,Center
from django.contrib.auth.models import User

from iSkyLIMS_core.models import Samples


STATUS_CHOICES = (
			('recorded',_("Recorded")),
	   		('approved',_("Approved")),
			('rejected',_("Rejected")),
			('queued',_("Queued")),
			('in_progress',_('In progress')),
			('delivered',_("Delivered")),
			('archived',_("Archived"))
		)

def service_files_upload(instance,filename):
	now = timezone_now()
	filename_base,filename_ext = os.path.splitext(filename)
	return 'iSkyLIMS_drylab/servicesRequest/%s_%s%s' % (
			now.strftime("%Y%m%d%H%M%S"),
			filename_base.lower(),
			filename_ext.lower(),
	)

# Create your models here.
class FileExt(models.Model):
	fileExt=models.CharField(_("File extension"),max_length=10)
	def __str__ (self):
		return '%s' %(self.fileExt)

class Platform(models.Model):
	platformName=models.CharField(_("Sequencing platform"),max_length=20)

	def __str__ (self):
 		return '%s' %(self.platformName)

	def get_platform_name(self):
		return '%s'  %(self.platformName)

class Machines (models.Model) :
	platformID = models.ForeignKey(Platform ,on_delete=models.CASCADE)
	machineName = models.CharField(_("Machine Name"),max_length=255)
	machineDescription = models.CharField(_("Description"),max_length=255,null=True,blank=True)
	machineLocation = models.CharField(_("Location"),max_length=255,null=True,blank=True)
	machineProvider = models.CharField(_("Machine owner brand"),max_length=255,null=True,blank=True)
	machineSerialNumber = models.CharField(_("Serial Number"),max_length=255,null=True,blank=True)
	machineState =  models.CharField(_("Machine State"),max_length=50,null=True,blank=True)
	machineOperationStart = models.DateField(auto_now_add=False, null=True,blank=True)
	machineOperationEnd = models.DateField(auto_now_add=False, null=True,blank=True)
	machineNumberLanes = models.CharField("Number of Lanes", max_length= 5)

	def __str__ (self) :
		return '%s' %(self.machineName)


	def get_machine_name(self):
		return '%s'  %(self.machineName)

	def get_number_of_lanes(self):
		return '%s' %(self.machineNumberLanes)

class AvailableService(MPTTModel):
	availServiceDescription=models.CharField(_("Available services"),max_length=100)
	parent=TreeForeignKey('self',models.SET_NULL,null=True,blank=True)

	def __str__(self):
		return '%s' %(self.availServiceDescription)

	def get_service_description(self):
		return '%s' %(self.availServiceDescription)

	class Meta:
		ordering = ["tree_id","lft"]
		verbose_name = ("AvailableService")
		verbose_name_plural = ("AvailableServices")

class Service(models.Model):
	## User requesting service:
	# 'serviceUsername' refactored to 'serviceUserid' which shows better its true nature
	#  decision taken to change Foreign Key from 'Profile'  to 'User' until full develop of "user registration"
	serviceUserId = models.ForeignKey(
				User ,
				on_delete=models.CASCADE, null=True)
	servicePlatform = models.ForeignKey(
				Platform ,
				on_delete=models.CASCADE ,
				verbose_name=_("Sequencing platform"),blank=True,null=True)
	serviceFileExt = models.ForeignKey(
				FileExt ,
				on_delete=models.CASCADE ,
				verbose_name=_("File extension"),blank=True,null=True)
	serviceAvailableService = TreeManyToManyField(
				AvailableService,
				verbose_name=_("AvailableServices"))

	serviceProjectNames = models.ManyToManyField(
				'iSkyLIMS_wetlab.Projects',
				verbose_name=_("User's projects"),blank=True)

	serviceSeqCenter=models.CharField(_("Sequencing center"),max_length=50,blank=False,null=True)
	serviceRequestNumber=models.CharField(max_length=80, null=True)
	serviceRequestInt=models.CharField(max_length=80, null=True)
	serviceRunSpecs=models.CharField(_("Run specifications"),max_length=10,blank=True,null=True)
	serviceFile=models.FileField(_("Service description file"),upload_to=service_files_upload, null=True,blank=True)
	serviceStatus=models.CharField(_("Service status"),max_length=15,choices=STATUS_CHOICES)
	serviceNotes=models.TextField(_("Service Notes"),max_length=2048,null=True)
	serviceCreatedOnDate= models.DateField(auto_now_add=True,null=True)

	serviceOnApprovedDate = models.DateField(auto_now_add=False, null=True,blank=True)
	serviceOnRejectedDate = models.DateField(auto_now_add=False, null=True,blank=True)

	serviceOnDeliveredDate = models.DateField(auto_now_add=False, null=True)



	def __str__ (self):
		return '%s' %(self.serviceRequestNumber)

	def get_service_information (self):
		platform = str(self.servicePlatform)

		return '%s;%s;%s;%s'  %(self.serviceRequestNumber ,self.serviceRunSpecs, self.serviceSeqCenter, platform)

	def get_service_id (self):
		return '%s' %self.serviceUserId.pk

	def get_service_information_with_service_name (self):
		platform = str(self.servicePlatform)
		if Resolution.objects.filter(resolutionServiceID__exact = self).exists():
			resolutions = Resolution.objects.filter(resolutionServiceID__exact = self).last()
			folder_name = resolutions.resolutionFullNumber
			if folder_name is None:
				resolution_for_service = ''
			else:
				folder_name_split= folder_name.split('_')
				resolution_for_service= '_'.join( folder_name_split[2:-1])
			assigned_to = resolutions.resolutionAsignedUser
			if assigned_to is None:
				assigned_to = ''
		else:
			resolutions_for_service = ''
			assigned_to = ''
		if self.serviceOnApprovedDate is None:
			approved_date = 'Not defined'
		else:
			approved_date = self.serviceOnApprovedDate.strftime("%d %B, %Y")



		if resolutions.resolutionEstimatedDate is None:
			estimated_date = 'Not defined'
		else:
			estimated_date = resolutions.resolutionEstimatedDate.strftime("%d %B, %Y")

		return '%s;%s;%s;%s;%s;%s;%s;%s'  %(self.serviceRequestNumber ,resolution_for_service , assigned_to, approved_date,
										estimated_date, self.serviceRunSpecs, self.serviceSeqCenter, platform)

	def get_service_dates (self):
		service_dates =[]
		service_dates.append(self.serviceCreatedOnDate.strftime("%d %B, %Y"))
		if self.serviceOnApprovedDate is None:
			if self.serviceStatus == 'rejected':
				service_dates.append('--')
			else:
				service_dates.append('Approved Date not set')
		else:
			service_dates.append(self.serviceOnApprovedDate.strftime("%d %B, %Y"))
		if self.serviceOnRejectedDate is None:
			if self.serviceStatus != 'recorded':
				service_dates.append('--')
			else:
				service_dates.append('Rejected Date not set')
		else:
			service_dates.append(self.serviceOnRejectedDate.strftime("%d %B, %Y"))

		return service_dates

	def get_stats_information (self):

		stats_information =[]
		stats_information.append(self.id)
		stats_information.append(self.serviceRequestNumber)
		stats_information.append(self.serviceStatus)

		stats_information.append(self.serviceCreatedOnDate.strftime("%d %B, %Y"))
		if self.serviceOnApprovedDate is None:
			if self.serviceOnRejectedDate is None:
				stats_information.append('--')
			else:
				stats_information.append(self.serviceOnRejectedDate.strftime("%d %B, %Y"))
		else:
			stats_information.append(self.serviceOnApprovedDate.strftime("%d %B, %Y"))
		if self.serviceOnDeliveredDate is None:
			stats_information.append('--')
		else:
			stats_information.append(self.serviceOnDeliveredDate.strftime("%d %B, %Y"))

		return stats_information

	def get_service_creation_time (self):
		return self.serviceCreatedOnDate.strftime("%d %B, %Y")

	def get_service_request_number(self):
		return '%s' %(self.serviceRequestNumber)

	def get_time_to_delivery (self):
		if self.serviceOnDeliveredDate == self.serviceCreatedOnDate :
			return 1
		else:
			number_days, time = str(self.serviceOnDeliveredDate - self.serviceCreatedOnDate).split(',')
			number, string_value = number_days.split(' ')
		return number

	def get_child_services (self):
		children_services = []
		all_services = self.serviceAvailableService.all()
		for service in all_services :
			if not service.get_children().exists():
				children_services.append([service.pk, service.get_service_description()])
		return children_services

class RequestedProjectInServicesManager(models.Manager):
	def create_request_project_service(self, data):
		new_project_request =self.create( projectService = data['projectService'], externalProjectKey = data['externalProjectKey'],
					externalProjectName = data['externalProjectName'])
		return new_project_request

class RequestedProjectInServices (models.Model):
	projectService = models.ForeignKey(
					Service,
					on_delete = models.CASCADE)
	externalProjectKey = models.CharField(max_length = 5, null = True, blank = True)
	externalProjectName = models.CharField(max_length = 70, null = True, blank = True)
	generated_at = models.DateField(auto_now_add = True)

	def __str__ (self):
		return '%s' %(self.externalProjectName)

	def get_requested_project_name (self):
		return '%s' %(self.externalProjectName)

	def get_requested_project_id (self):
		return '%s' %(self.pk)

	objects = RequestedProjectInServicesManager()

class RequestedSamplesInServices (models.Model):
	samplesInService = models.ForeignKey(
					Service,
					on_delete = models.CASCADE)
	sample = models.ForeignKey(
					Samples,
					null = True, blank = True, on_delete = models.CASCADE)
	externalSampleKey = models.CharField(max_length = 5, null = True, blank = True)
	externalSampleName = models.CharField(max_length = 50, null = True, blank = True)
	externalSamplePath =  models.CharField(max_length = 100, null = True, blank = True)

	generated_at = models.DateField(auto_now_add = True)


class Resolution(models.Model):
	resolutionServiceID=models.ForeignKey(Service ,on_delete=models.CASCADE)
	resolutionNumber=models.CharField(_("Resolutions name"),max_length=255,null=True)
	#resolutionServiceSRV=models.CharField(_("Service identifier"),max_length=10)
	resolutionEstimatedDate=models.DateField(_(" Estimated resolution date"), null = True,blank=False)
	resolutionDate=models.DateField(_("Resolution date"),auto_now_add=True,blank=True)
	#resolutionDate=models.DateField(_("Resolution date"),auto_now_add=False,blank=True)
	resolutionOnQueuedDate = models.DateField(auto_now_add=False, null=True,blank=True)
	resolutionOnInProgressDate = models.DateField(auto_now_add=False, null=True,blank=True)
	resolutionNotes=models.TextField(_("Resolution notes"),max_length=1000, null=True)
	resolutionFullNumber = models.CharField(_("Acronym Name"),max_length=255,null=True,blank=True)
	resolutionAsignedUser = models.ForeignKey(User, related_name='groups+', on_delete=models.CASCADE, null=True,blank=True )
	#resolutionAsignedUser = models.ForeignKey(User ,on_delete=models.CASCADE, null=True)
	def __str__ (self):
		return '%s' %(self.resolutionNumber)
	def get_resolution_information (self):
		resolution_info =[]
		resolution_info.append(self.resolutionNumber)
		resolution_info.append(self.resolutionFullNumber)
		resolution_info.append(self.resolutionAsignedUser)
		if self.resolutionEstimatedDate is not None:
		    resolution_info.append(self.resolutionEstimatedDate.strftime("%d %B, %Y"))
		else:
		    resolution_info.append("Not defined yet")

		if self.resolutionOnQueuedDate is not None:
		    resolution_info.append(self.resolutionOnQueuedDate.strftime("%d %B, %Y"))
		else:
		    resolution_info.append("Not defined yet")

		if self.resolutionOnInProgressDate is None:
			resolution_info.append('--')
		else:
			resolution_info.append(self.resolutionOnInProgressDate.strftime("%d %B, %Y"))

		resolution_info.append(self.resolutionNotes)

		return resolution_info


class Delivery(models.Model):
	deliveryResolutionID=models.ForeignKey(Resolution ,on_delete=models.CASCADE )
	#deliveryResolutionID=models.OneToOneField(Resolution ,on_delete=models.CASCADE )
	#deliveryNumber=models.IntegerField(_("Number of deliveries"))
	#deliveryEstimatedDate=models.DateField(_("Delivery estimated date"))
	deliveryDate=models.DateField(_("Delivery date"),auto_now_add=True,null=True,blank=True)
	#deliveryDate=models.DateField(_("Delivery date"),auto_now_add=False,null=True,blank=True)
	deliveryNotes=models.TextField(_("Delivery notes"),max_length=1000, null=True)
	def __str__ (self):
		return '%s' %(self.deliveryResolutionID)
	def get_delivery_information (self):
		delivery_info = []
		delivery_info.append(self.deliveryResolutionID.resolutionNumber)
		delivery_info.append(self.deliveryDate.strftime("%d %B, %Y"))
		delivery_info.append(self.deliveryNotes)
		return delivery_info

class RunIDFolder (models.Model):
	servicerequest = models.ForeignKey(
				Service,
				on_delete = models.CASCADE)
	run_id_Folder = models.CharField(max_length = 50)

	def __str__ (self):
		return '%s' %(self.run_id_Folder)

class PipelinesManager (models.Manager):
	def create_pipeline(self, data):
		availableService = AvailableService.objects.get(pk__exact = data['availableService_id'])

		new_pipeline = self.create(availableService = availableService, userName = data['userName'],
				pipelineName = data['pipelineName'], pipelineInUse = True,
				pipelineVersion	= data['pipelineVersion'],
				externalRequest = data['externalRequest'], useRunFolder = data['useRunFolder'])
		return new_pipeline

class Pipelines(models.Model):
	availableService = models.ForeignKey(
				AvailableService,
				on_delete = models.CASCADE)
	userName = models.ForeignKey(
                User,
                on_delete=models.CASCADE)
	pipelineName = models.CharField(max_length = 50)
	pipelineVersion = models.CharField(max_length = 10)


	pipelineStrFolder = models.CharField(max_length = 20)
	useRunFolder = models.NullBooleanField(default = True, null = True)
	automatic = models.BooleanField(default = True)

	useRunFolder = models.NullBooleanField(default = True, null = True)


	externalRequest = models.BooleanField(default = True)

	default = models.BooleanField(default = False)
	pipelineInUse = models.BooleanField(default = True)
	generated_at = models.DateTimeField(auto_now_add = True)

	def __str__ (self):
		return '%s' %(self.pipelineName)

	'''
	def get_all_pipeline_data(self):
		data = []
		data.append(self.availableService.get_service_description())
		data.append(self.userName.username)
		data.append(self.pipelineName)
		data.append(self.pipelineVersion)
		data.append(self.generated_at.strftime("%B %d, %Y"))
		data.append(self.default)
		data.append(self.pipelineInUse)
		return data
	'''


	def get_pipeline_name (self):
		return '%s' %(self.pipelineName)

	def get_external_request(self):
		return '%s' %(self.externalRequest)

	def get_pipleline_service_obj(self):
		return self.availableService

	def get_pipeline_version(self):
		return '%s' %(self.pipelineVersion)

	def get_pipeline_additional(self):
		data = []
		data.append(self.userName.username)
		data.append(self.generated_at.strftime("%B %d, %Y"))

		data.append(self.default)
		data.append(self.pipelineInUse)
		data.append(self.externalRequest)
		return data

	def get_pipeline_basic (self):
		data = []
		data.append(self.pipelineName)
		data.append(self.pipelineVersion)
		data.append(self.availableService.get_service_description())
		return data

	def get_pipeline_info(self):
		data = []
		data.append(self.availableService.get_service_description())
		data.append(self.userName.username)
		data.append(self.pipelineName)
		data.append(self.pipelineVersion)
		data.append(self.generated_at.strftime("%B %d, %Y"))
		data.append(self.default)
		data.append(self.pipelineInUse)
		data.append(self.pk)
		return data

	def get_used_run_folder (self):
		return '%s' %(self.useRunFolder)

	def remove_default_pipeline(self):
		self.default = False
		self.save()
		return self

	def set_default_pipeline(self):
		self.default = True
		self.save()
		return self


	objects = PipelinesManager()


class ParameterPipelineManager(models.Manager):
	def create_pipeline_parameters(self, pipeline_parameters):
		new_parameter_action_pipeline = self.create(parameterPipeline = pipeline_parameters['parameterPipeline'],
					parameterName = pipeline_parameters['parameterName'], parameterValue = pipeline_parameters['parameterValue'])
		return new_parameter_action_pipeline

class ParameterPipeline (models.Model):
	parameterPipeline = models.ForeignKey(
				Pipelines,
				on_delete = models.CASCADE)
	parameterName = models.CharField(max_length = 80)
	parameterValue = models.CharField(max_length = 80)

	def __str__ (self):
		return '%s' %(self.parameterName)

	def get_pipeline_parameter_name (self):
		return '%s' %(self.parameterName)

	def get_pipeline_parameter_value (self):
		return '%s' %(self.parameterValue)

	def get_pipeline_parameters(self):
		data = []
		data.append(self.parameterName)
		data.append(self.parameterValue)
		return data

	objects = ParameterPipelineManager()

class JobStates (models.Model):
	jobStateName = models.CharField(max_length = 20)

	def __str__ (self):
		return '%s' %(self.jobStateName)

	def get_job_state (self):
		return '%s' %(self.jobStateName)


class PipelineExternalDataJobsManager(models.Manager):
	def create_pipeline_external_data_job (self, preparation_data):
		jobState = JobStates.objects.get(jobStateName__exact = 'Queued')
		new_preparation_pipeline = self.create(pipeline = preparation_data['pipeline'],
				availableService = preparation_data['availableService'], jobState =  jobState,
				pipelineName = preparation_data['pipelineName'],  pipelineVersion = preparation_data['pipelineVersion'])
		return new_preparation_pipeline

class PipelineExternalDataJobs (models.Model):
	pipeline = models.ForeignKey(
				Pipelines,
				on_delete = models.CASCADE)
	availableService = models.ForeignKey(
				AvailableService,
				on_delete = models.CASCADE)
	jobState =  models.ForeignKey(
				JobStates,
				on_delete = models.CASCADE)

	pipelineName = models.CharField(max_length = 50)
	pipelineVersion = models.CharField(max_length = 10)

	generated_at = models.DateTimeField(auto_now_add = True)
	lastRequestedTime = models.DateTimeField(auto_now_add = False, null = True, blank = True)

	def __str__ (self):
		return '%s_%s' %(self.pipeline, self.availableService)


	objects = PipelineExternalDataJobsManager()
