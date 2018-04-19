from django.contrib import admin
from drylab.models import *
from django_mptt_admin.admin import DjangoMpttAdmin
from django.contrib.auth.admin import UserAdmin as BaseUserAdmin
from django.contrib.auth.models import User

class FileExtAdmin(admin.ModelAdmin):
	list_display = ('fileExt',)

class PlatformAdmin(admin.ModelAdmin):
	list_display=('platformName',)


# 'serviceUsername' refactored to 'serviceUserid' which shows better its real nature
class ServiceAdmin(admin.ModelAdmin):
	list_display=('serviceRequestNumber','serviceUserId','serviceSeqCenter','servicePlatform','serviceRunSpecs','serviceFileExt','serviceFile','serviceStatus','serviceNotes','serviceCreatedOnDate','serviceOnApprovedDate','serviceOnRejectedDate','serviceOnDeliveredDate')

class AvailableServiceAdmin(DjangoMpttAdmin):
	list_display=('availServiceDescription',)

class ResolutionAdmin(admin.ModelAdmin):
	list_display=('resolutionServiceID','resolutionNumber','resolutionDate','resolutionEstimatedDate','resolutionOnQueuedDate','resolutionOnInProgressDate','resolutionFullNumber','resolutionAsignedUser','resolutionNotes')

class DeliveryAdmin(admin.ModelAdmin):
	list_display=('deliveryResolutionID','deliveryDate','deliveryNotes')

admin.site.register(FileExt,FileExtAdmin)
admin.site.register(Platform,PlatformAdmin)
admin.site.register(Service,ServiceAdmin)
admin.site.register(AvailableService,AvailableServiceAdmin)
admin.site.register(Resolution,ResolutionAdmin)
admin.site.register(Delivery,DeliveryAdmin)

