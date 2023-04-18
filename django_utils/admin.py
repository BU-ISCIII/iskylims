from django.contrib import admin
from django.contrib.auth.admin import UserAdmin as BaseUserAdmin
from django.contrib.auth.models import User

from django_utils.models import *


# Register your models here.
class CenterAdmin(admin.ModelAdmin):
    list_display = ("center_name", "center_abbr")


class ProfileInline(admin.StackedInline):
    model = Profile
    can_delete = False
    verbose_name_plural = "Profiles"


# Define a new User admin
class UserAdmin(BaseUserAdmin):
    inlines = (ProfileInline,)


class ProfileAdmin(admin.ModelAdmin):
    list_display = [
        "profile_user_id",
        "profile_position",
        "profile_center",
        "profile_area",
        "profile_extension",
    ]


class ClassificationAreaAdmin(admin.ModelAdmin):
    list_display = ["classification_area_name", "classification_area_description"]


# Re-register UserAdmin
admin.site.unregister(User)
admin.site.register(User, UserAdmin)
admin.site.register(Center, CenterAdmin)
admin.site.register(Profile, ProfileAdmin)
admin.site.register(ClassificationArea, ClassificationAreaAdmin)
