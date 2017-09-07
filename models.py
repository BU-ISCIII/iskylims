import os
from django.db import models
from django import forms
from django.utils.encoding import python_2_unicode_compatible
from django.utils.translation import ugettext_lazy as _
from mptt.models import MPTTModel
from mptt.fields import TreeForeignKey, TreeManyToManyField
from django.utils.timezone import now as timezone_now
from wetlab.models import RunProcess
from utils.models import Profile,Center

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
	return 'drylab/servicesRequest/%s_%s%s' % (
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

class AvailableService(MPTTModel):
	availServiceDescription=models.CharField(_("Available services"),max_length=100)
	parent=TreeForeignKey('self',models.SET_NULL,null=True,blank=True) 

	def __str__(self):
		return '%s' %(self.availServiceDescription)
	
	class Meta: 
		ordering = ["tree_id","lft"]
		verbose_name = ("AvailableService")
		verbose_name_plural = ("AvailableServices")

class Service(models.Model):
	serviceUsername=models.ForeignKey(Profile)
	serviceSeqCenter=models.CharField(_("Sequencing center"),max_length=50,blank=False,null=True)
	serviceRunID=models.ForeignKey(RunProcess,null=True)
	servicePlatform=models.ForeignKey(Platform,verbose_name=_("Sequencing platform"),blank=False,null=True)
	serviceRunSpecs=models.CharField(_("Run specifications"),max_length=10,blank=False,null=True)
	serviceFileExt=models.ForeignKey(FileExt,verbose_name=_("File extension"),blank=False,null=True)
	serviceAvailableService=TreeManyToManyField(AvailableService,verbose_name=_("AvailableServices"))
	serviceFile=models.FileField(_("Service description file"),upload_to=service_files_upload)
	serviceStatus=models.CharField(_("Service status"),max_length=10,choices=STATUS_CHOICES)
	serviceNotes=models.TextField(_("Service Notes"),max_length=500)
	
	def __str__ (self):
		return '%s' %(self.serviceName)

class Resolution(models.Model):
	resolutionServiceID=models.ForeignKey(Service)
	resolutionNumber=models.IntegerField(_("Number of resolutions"))
	resolutionServiceSRV=models.CharField(_("Service identifier"),max_length=10)
	resolutionDate=models.DateField(_("Resolution date"),auto_now_add=True)


class Delivery(models.Model):
	deliveryResolutionID=models.ForeignKey(Resolution)
	deliveryNumber=models.IntegerField(_("Number of deliveries"))
	deliveryEstimatedDate=models.DateField(_("Delivery estimated date"))
	deliveryDate=models.DateField(_("Delivery date"),auto_now_add=True)
	deliveryNotes=models.TextField(_("Delivery notes"))

