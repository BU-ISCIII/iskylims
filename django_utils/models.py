from django.db import models

from django.utils.translation import gettext_lazy as _
from django.contrib.auth.models import User


class Center(models.Model):
    center_name = models.CharField(_("Center"), max_length=50)
    center_abbr = models.CharField(_("Acronym"), max_length=25)

    class Meta:
        db_table = "utils_center"

    def __str__(self):
        return "%s" % (self.center_abbr)

    def get_user_center_abbr(self):
        return "%s" % (self.center_abbr)

    def get_center_name(self):
        return "%s" % (self.center_name)


class ClassificationArea(models.Model):
    classification_area_name = models.CharField(max_length=80)
    classification_area_description = models.CharField(
        max_length=255, null=True, blank=True
    )

    class Meta:
        db_table = "utils_classification_area"

    def __str__(self):
        return "%s" % (self.classification_area_name)

    def get_classification_area_name(self):
        return "%s" % (self.classification_area_name)


class Profile(models.Model):
    profile_user_id = models.OneToOneField(
        User, related_name="profile", on_delete=models.CASCADE
    )
    profile_classification_area = models.ForeignKey(
        ClassificationArea, on_delete=models.CASCADE, null=True, blank=True
    )
    profile_center = models.ForeignKey(
        Center, on_delete=models.CASCADE, verbose_name=_("Center")
    )
    profile_position = models.CharField(_("Position"), max_length=50)
    profile_area = models.CharField(_("Area / Unit"), max_length=50)
    profile_extension = models.CharField(_("Phone extension"), max_length=5)

    class Meta:
        db_table = "utils_profile"

    def __str__(self):
        return "%s" % (self.pk)

    def get_clasification_area(self):
        if self.profile_classification_area is not None:
            return "%s" % (
                self.profile_classification_area.get_classification_area_name()
            )
        return "Not available"

    def get_user_center_abbr(self):
        if self.profile_center is not None:
            return "%s" % (self.profile_center.get_user_center_abbr())
        return "Not defined"
