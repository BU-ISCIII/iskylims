from django.db import models
import os
from django.db import models
from django import forms
from django.utils.encoding import python_2_unicode_compatible
from django.utils.translation import ugettext_lazy as _
from django.contrib.auth.models import User,AbstractUser
from django.db.models.signals import post_save
from django.dispatch import receiver

# Create your models here.
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

# @receiver(post_save, sender=User)
# def create_user_profile(sender, instance, created, **kwargs):
# 	if created:
# 		Profile.objects.create(user=instance)
#
# @receiver(post_save, sender=User)
# def save_user_profile(sender, instance, **kwargs):
# 	instance.profile.save()
