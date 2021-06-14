from django.db import models
import os
from django.db import models
from django import forms
from django.utils.translation import gettext_lazy as _
from django.contrib.auth.models import User,AbstractUser


class Center(models.Model):
  	centerName=models.CharField(_("Center"),max_length=50)
  	centerAbbr=models.CharField(_("Acronym"),max_length=25)

  	def __str__ (self):
  		return '%s' %(self.centerName)

class Profile(models.Model):
     profileUserID = models.OneToOneField(User, on_delete=models.CASCADE)
     profilePosition=models.CharField(_("Position"),max_length=50)
     profileCenter=models.ForeignKey(Center, on_delete=models.CASCADE, verbose_name=_("Center"))
     profileArea=models.CharField(_("Area / Unit"),max_length=50)
     profileExtension=models.CharField(_("Phone extension"),max_length=5)

     def __str__ (self):
      	return '%s' %(self.pk)
