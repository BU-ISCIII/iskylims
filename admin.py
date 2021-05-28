from django.contrib import admin
from iSkyLIMS_clinic.models import *

# Register your models here.
class ClinicSampleRequestAdmin (admin.ModelAdmin):
    list_display= ('sampleCore', 'clinicSampleState', 'sampleRequestUser', 'doctor_id', 'confirmationCode', 'priority', 'serviceDate')

class ClinicSampleStateAdmin (admin.ModelAdmin):
    list_display = ('clinicState',)

class ServiceUnitsAdmin(admin.ModelAdmin):
    list_display = ('serviceUnitName',)

class PatientDataAdmin(admin.ModelAdmin):
    list_display = ('patienCore', 'address', 'phone', 'birthday', 'smoker', 'notificationPreference')

class DoctorAdmin(admin.ModelAdmin):
    list_display = ('doctorName',)

class ConfigSettingAdmin(admin.ModelAdmin):
    list_display = ['configurationName', 'configurationValue']

admin.site.register(ClinicSampleRequest, ClinicSampleRequestAdmin)
admin.site.register(ClinicSampleState,ClinicSampleStateAdmin)
admin.site.register(Doctor,DoctorAdmin)


admin.site.register(PatientData,PatientDataAdmin)
admin.site.register(ServiceUnits,ServiceUnitsAdmin)
admin.site.register(ConfigSetting,ConfigSettingAdmin)
