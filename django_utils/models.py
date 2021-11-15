from django.db import models
import os
from django.db import models
from django import forms
from django.utils.translation import gettext_lazy as _
from django.contrib.auth.models import User,AbstractUser


class Center(models.Model):
  	centerName = models.CharField(_("Center"),max_length=50)
  	centerAbbr = models.CharField(_("Acronym"),max_length=25)

  	def __str__ (self):
  		return '%s' %(self.centerName)

  	def get_center_abbr(self):
  	  	return '%s' %(self.centerAbbr)

class ClassificationArea(models.Model):
    classificationAreaName = models.CharField(max_length = 80)
    classificationAreaDescription = models.CharField(max_length = 255, null=True, blank = True)

    def __str__ (self):
        return '%s' %(self.classificationAreaName)

    def get_classification_area_name(self):
        return '%s' %(self.classificationAreaName)

class Profile(models.Model):
    profileUserID = models.OneToOneField(
            User,
            on_delete=models.CASCADE)
    profileClassificationArea = models.ForeignKey(
            ClassificationArea,
            on_delete= models.CASCADE, null=True , blank=True)
    profileCenter = models.ForeignKey(
            Center,
            on_delete=models.CASCADE, verbose_name=_("Center"))
    profilePosition = models.CharField(_("Position"),max_length=50)
    profileArea = models.CharField(_("Area / Unit"),max_length=50)
    profileExtension = models.CharField(_("Phone extension"),max_length=5)

    def __str__ (self):
      	return '%s' %(self.pk)

    def get_clasification_area(self):
        if self.profileClassificationArea != None:
            return '%s' %(self.profileClassificationArea.get_classification_area_name())
        return 'Not available'

    def get_profile_center_abbr(self):
        if self.profileCenter != None:
            return '%s' %(self.profileCenter.get_center_abbr())
        return 'Not defined'
