import os
from django.db import models
from django import forms
from django.utils.encoding import python_2_unicode_compatible
from django.utils.translation import ugettext_lazy as _
from mptt.models import MPTTModel
from mptt.fields import TreeForeignKey, TreeManyToManyField
from django.utils.timezone import now as timezone_now
from wetlab.models import RunProcess, Projects
from utils.models import Profile,Center
##jlgarcia 08/02/2018
from django.contrib.auth.models import User

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
	##jlgarcia 07/02/2018. Changed serviceUsername's foreign key to 'User'
	# workaround to store 'Service' instances in the database 
	#serviceUsername=models.ForeignKey(Profile ,on_delete=models.CASCADE)
	##jlgarcia 09/02/2018. Refactoring:'serviceUsername' --> 'serviceUserid'
	# which is its real nature
	#serviceUsername=models.ForeignKey(User ,on_delete=models.CASCADE, null=True)
	#User-requesting service#User-requesting service
	serviceUserId=models.ForeignKey(User ,on_delete=models.CASCADE, null=True)
	serviceSeqCenter=models.CharField(_("Sequencing center"),max_length=50,blank=False,null=True)
	##jlgarcia 01/02/2018.  serviceRunID not used in forms.py/serviceRequestForm()
	#serviceRunID=models.ForeignKey(RunProcess,null=True ,on_delete=models.CASCADE)

	##jlgarcia 02/02/2018. Addition of member 'serviceProjectNames' to support 
	# implementation of drop down menu to choose a project name of a list of projects
	# belonging to the logged-in user in the service request form
	# serviceProjectNames=models.ForeignKey(Projects,on_delete=models.CASCADE, verbose_name=_("User's project names"),blank=False,null=True)
	##jlgarcia 07/02/2018. Dropdown menu needs "ManyToManyField" support. 
	# "null has no effect" (according to docs.djangoproject.com, Model field reference)
	# and execution of makemigrations 
	#serviceProjectNames=models.ForeignKey(Projects,on_delete=models.CASCADE, verbose_name=_("User's project names"),blank=False,null=True)
	serviceProjectNames=models.ManyToManyField(Projects,verbose_name=_("User's projects"),blank=False)
	servicePlatform=models.ForeignKey(Platform ,on_delete=models.CASCADE , verbose_name=_("Sequencing platform"),blank=False,null=True)
	serviceRunSpecs=models.CharField(_("Run specifications"),max_length=10,blank=False,null=True)
	serviceFileExt=models.ForeignKey(FileExt ,on_delete=models.CASCADE ,verbose_name=_("File extension"),blank=False,null=True)
	serviceAvailableService=TreeManyToManyField(AvailableService,verbose_name=_("AvailableServices"))
	serviceFile=models.FileField(_("Service description file"),upload_to=service_files_upload)
	serviceStatus=models.CharField(_("Service status"),max_length=10,choices=STATUS_CHOICES)
	serviceNotes=models.TextField(_("Service Notes"),max_length=500)
	
	def __str__ (self):
		return '%s' %(self.serviceName)

class Resolution(models.Model):
	resolutionServiceID=models.ForeignKey(Service ,on_delete=models.CASCADE)
	resolutionNumber=models.IntegerField(_("Number of resolutions"))
	resolutionServiceSRV=models.CharField(_("Service identifier"),max_length=10)
	resolutionDate=models.DateField(_("Resolution date"),auto_now_add=True)


class Delivery(models.Model):
	deliveryResolutionID=models.ForeignKey(Resolution ,on_delete=models.CASCADE )
	deliveryNumber=models.IntegerField(_("Number of deliveries"))
	deliveryEstimatedDate=models.DateField(_("Delivery estimated date"))
	deliveryDate=models.DateField(_("Delivery date"),auto_now_add=True)
	deliveryNotes=models.TextField(_("Delivery notes"))

