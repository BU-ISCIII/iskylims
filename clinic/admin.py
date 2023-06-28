from django.contrib import admin

import clinic.models


# Register your models here.
class ClinicSampleRequestAdmin(admin.ModelAdmin):
    list_display = (
        "sample_core",
        "clinic_sample_state",
        "sample_request_user",
        "doctor_id",
        "confirmation_code",
        "priority",
        "service_date",
    )


class ClinicSampleStateAdmin(admin.ModelAdmin):
    list_display = ("clinic_state",)


class ServiceUnitsAdmin(admin.ModelAdmin):
    list_display = ("service_unit_name",)


class PatientDataAdmin(admin.ModelAdmin):
    list_display = (
        "patien_core",
        "address",
        "phone",
        "birthday",
        "smoker",
        "notification_preference",
    )


class DoctorAdmin(admin.ModelAdmin):
    list_display = ("doctor_name",)


class ConfigSettingAdmin(admin.ModelAdmin):
    list_display = ["configuration_name", "configuration_value"]


admin.site.register(clinic.models.ClinicSampleRequest, ClinicSampleRequestAdmin)
admin.site.register(clinic.models.ClinicSampleState, ClinicSampleStateAdmin)
admin.site.register(clinic.models.Doctor, DoctorAdmin)
admin.site.register(clinic.models.PatientData, PatientDataAdmin)
admin.site.register(clinic.models.ServiceUnits, ServiceUnitsAdmin)
admin.site.register(clinic.models.ConfigSetting, ConfigSettingAdmin)
