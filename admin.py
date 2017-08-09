from django.contrib import admin
from drylab.models import *
from django_mptt_admin.admin import DjangoMpttAdmin

class CenterAdmin(admin.ModelAdmin):                           
     list_display = ('centerName','centerAbbr') 
     #import pdb; pdb.set_trace()                                          

class FileExtAdmin(admin.ModelAdmin):
	list_display = ('fileExt',)

class PlatformAdmin(admin.ModelAdmin):
	list_display=('platformName',)

class ServiceAdmin(admin.ModelAdmin):
	list_display=('serviceName','servicePosition','serviceCenter','serviceArea','serviceExtension','serviceEmail','serviceSeqCenter','servicePlatform','serviceRunSpecs','serviceFileExt','serviceFile','serviceStatus','serviceNotes')

class AvailableServiceAdmin(DjangoMpttAdmin):
	list_display=('availServiceDescription',)

class ResolutionAdmin(admin.ModelAdmin):
	list_display=('resolutionServiceID','resolutionNumber','resolutionServiceSRV','resolutionDate')

class DeliveryAdmin(admin.ModelAdmin):
	list_display=('deliveryResolutionID','deliveryEstimatedDate','deliveryDate','deliveryNumber','deliveryNotes')



admin.site.register(Center,CenterAdmin)
admin.site.register(FileExt,FileExtAdmin)
admin.site.register(Platform,PlatformAdmin)
admin.site.register(Service,ServiceAdmin)
admin.site.register(AvailableService,AvailableServiceAdmin)
admin.site.register(Resolution,ResolutionAdmin)
admin.site.register(Delivery,DeliveryAdmin)
