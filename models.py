from django.db import models
import datetime
from django.db import models
from django.utils import timezone
from django import forms
from django.utils.encoding import python_2_unicode_compatible
from mptt.models import MPTTModel
from mptt.fields import TreeForeignKey, TreeManyToManyField


# Create your models here.

class Center(models.Model):
	centerName=models.CharField(max_length=50)
	centerAbbr=models.CharField(max_length=25)

	def __str__ (self):
		return '%s' %(self.centerName)

class FileExt(models.Model):
	fileExt=models.CharField(max_length=10)
	def __str__ (self):
		return '%s' %(self.fileExt)

class Platform(models.Model):
	platformName=models.CharField(max_length=20)

	def __str__ (self):
 		return '%s' %(self.platformName)

class AvailableService(MPTTModel):
	availServiceDescription=models.CharField(max_length=100)
	parent=TreeForeignKey('self',models.SET_NULL,null=True,blank=True) 

	def __str__(self):
		return self.availServiceDescription
	
	class Meta: 
		ordering = ["tree_id","lft"]
		verbose_name = ("AvailableService")
		verbose_name_plural = ("AvailableServices")

class Service(models.Model):
	serviceName=models.CharField(max_length=50)
	servicePosition=models.CharField(max_length=50)
	serviceCenter=models.ForeignKey(Center)
	serviceArea=models.CharField(max_length=50)
	serviceExtension=models.CharField(max_length=5)
	serviceEmail=models.EmailField(max_length=45)
	serviceSeqCenter=models.CharField(max_length=50)
	servicePlatform=models.ForeignKey(Platform)
	serviceRunSpecs=models.CharField(max_length=10)
	serviceFileExt=models.ForeignKey(FileExt)
	serviceAvailableService=TreeManyToManyField(AvailableService,verbose_name=("AvailableServices"))
	serviceFile=models.FileField(upload_to='documents/')
	serviceStatus=models.CharField(max_length=10)
	serviceNotes=models.TextField(max_length=500)
	
	def __str__ (self):
		return '%s' %(self.serviceName)

class Resolution(models.Model):
	resolutionServiceID=models.ForeignKey(Service)
	resolutionNumber=models.IntegerField
	resolutionServiceSRV=models.CharField(max_length=10)
	resolutionDate=models.DateField(auto_now_add=True)


class Delivery(models.Model):
	deliveryResolutionID=models.ForeignKey(Resolution)
	deliveryNumber=models.IntegerField()
	deliveryEstimatedDate=models.DateField()
	deliveryDate=models.DateField(auto_now_add=True)
	deliveryNotes=models.TextField()

