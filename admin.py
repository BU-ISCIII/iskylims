from django.contrib import admin
from iSkyLIMS_drylab.models import *
from django_mptt_admin.admin import DjangoMpttAdmin
from django.contrib.auth.admin import UserAdmin as BaseUserAdmin
from django.contrib.auth.models import User

class ResolutionStatesAdmin(admin.ModelAdmin):
	list_display = ['resolutionStateName']
class FileExtAdmin(admin.ModelAdmin):
	list_display = ('fileExt',)
'''
class PlatformAdmin(admin.ModelAdmin):
	list_display=('platformName',)

class MachinesAdmin (admin.ModelAdmin) :
	list_display=('machineName', 'platformID','machineDescription', 'machineLocation','machineProvider','machineSerialNumber', 'machineState','machineOperationStart','machineOperationEnd','machineNumberLanes')
'''
# 'serviceUsername' refactored to 'serviceUserid' which shows better its real nature
class ServiceAdmin(admin.ModelAdmin):
	list_display =('serviceRequestNumber', 'serviceSeqCenter', 'serviceUserId','serviceStatus','serviceCreatedOnDate','serviceOnDeliveredDate')
	list_filter = ['serviceCreatedOnDate']
	search_fields = ['serviceRequestNumber__icontains']

class RequestedSamplesInServicesAdmin(admin.ModelAdmin):
	list_display = ['samplesInService','sampleName', 'projectName', 'runName']

class UploadServiceFileAdmin(admin.ModelAdmin):
	list_display = ['uploadService','uploadFile' , 'uploadFileName']
	search_fields = ('uploadService__serviceRequestNumber__icontains','uploadFileName__icontains')


class AvailableServiceAdmin(DjangoMpttAdmin):
	list_display=['availServiceDescription']

class ResolutionAdmin(admin.ModelAdmin):
	list_display=('resolutionServiceID','resolutionNumber', 'resolutionState','resolutionDate','resolutionEstimatedDate','resolutionOnQueuedDate','resolutionOnInProgressDate','resolutionFullNumber','resolutionAsignedUser','resolutionNotes')

class ResolutionParametersAdmin(admin.ModelAdmin):
	list_display = ['resolution','resolutionParameter', 'resolutionParamValue','resolutionParamNotes']

class DeliveryAdmin(admin.ModelAdmin):
	list_display=['deliveryResolutionID','deliveryDate', 'executionStartDate', 'executionEndDate', 'executionTime', 'permanentUsedSpace', 'temporaryUsedSpace','deliveryNotes']

class PipelinesManager(admin.ModelAdmin):
	list_display = ['pipelineName','pipelineVersion','pipelineInUse','pipelineFile']

class ParameterPipelineManager(admin.ModelAdmin):
	list_display = ['parameterPipeline', 'parameterName', 'parameterValue']

class ConfigSettingAdmin(admin.ModelAdmin):
    list_display = ['configurationName', 'configurationValue']


'''
class JobStatesManager(admin.ModelAdmin):
	list_display = ['jobStateName']

class PipelineExternalDataJobsManager(admin.ModelAdmin):
	list_display = ['pipeline', 'availableService', 'serviceRequestNumber', 'folderData', 'jobState', 'lastRequestedTime', 'pipelineName', 'pipelineVersion']
'''
admin.site.register(ResolutionStates,ResolutionStatesAdmin)
admin.site.register(FileExt,FileExtAdmin)
#admin.site.register(Platform,PlatformAdmin)
# admin.site.register(Machines,MachinesAdmin)
admin.site.register(Service,ServiceAdmin)
admin.site.register(UploadServiceFile, UploadServiceFileAdmin)
admin.site.register(ResolutionParameters, ResolutionParametersAdmin)
#admin.site.register(RequestedProjectInServices, RequestedProjectInServicesAdmin)
admin.site.register(RequestedSamplesInServices, RequestedSamplesInServicesAdmin)
admin.site.register(AvailableService,AvailableServiceAdmin)
admin.site.register(Resolution,ResolutionAdmin)
admin.site.register(Delivery,DeliveryAdmin)

admin.site.register(Pipelines,PipelinesManager)
admin.site.register(ParameterPipeline,ParameterPipelineManager)

admin.site.register(ConfigSetting, ConfigSettingAdmin)

#admin.site.register(JobStates,JobStatesManager)
#admin.site.register(PipelineExternalDataJobs, PipelineExternalDataJobsManager)
