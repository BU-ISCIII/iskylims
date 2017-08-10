from django.db import models
import datetime
from django.db import models
from django.utils import timezone
from django import forms
from django.utils.encoding import python_2_unicode_compatible
from django.utils.translation import ugettext_lazy as _
from mptt.models import MPTTModel
from mptt.fields import TreeForeignKey, TreeManyToManyField


# Create your models here.

class Center(models.Model):
	centerName=models.CharField(_("Center"),max_length=50)
	centerAbbr=models.CharField(_("Center"),max_length=25)

	def __str__ (self):
		return '%s' %(self.centerName)

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
		return self.availServiceDescription
	
	class Meta: 
		ordering = ["tree_id","lft"]
		verbose_name = ("AvailableService")
		verbose_name_plural = ("AvailableServices")

class Service(models.Model):
	serviceName=models.CharField(_("Complete name"),max_length=50)
	servicePosition=models.CharField(_("Position"),max_length=50)
	serviceCenter=models.ForeignKey(Center,verbose_name=_("Center"))
	serviceArea=models.CharField(_("Area / Unit"),max_length=50)
	serviceExtension=models.CharField(_("Phone extension"),max_length=5)
	serviceEmail=models.EmailField(_("email"),max_length=45)
	serviceSeqCenter=models.CharField(_("Sequencing center"),max_length=50)
	servicePlatform=models.ForeignKey(Platform,verbose_name=_("Sequencing platform"))
	serviceRunSpecs=models.CharField(_("Run specifications"),max_length=10)
	serviceFileExt=models.ForeignKey(FileExt,verbose_name=_("File extension"))
	serviceAvailableService=TreeManyToManyField(AvailableService,verbose_name=_("AvailableServices"))
	serviceFile=models.FileField(_("Service description file"),upload_to='documents/')
	serviceStatus=models.CharField(_("Service status"),max_length=10)
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

