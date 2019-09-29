from django.contrib import admin
from iSkyLIMS_clinic.models import *

# Register your models here.
class ClinicSampleRequestAdmin (admin.ModelAdmin):
    list_display= ('sampleCore', 'clinicSampleState' ,'patient_id', 'doctor_id', 'confirmationCode', 'priority', 'coments')

class ClinicSampleStateAdmin (admin.ModelAdmin):
    list_display = ('clinicState',)

class ServiceUnitsAdmin(admin.ModelAdmin):
    list_display = ('serviceUnitName',)

class PatientAdmin(admin.ModelAdmin):
    list_display = ('patientName', 'numberOfHistory')

class DoctorAdmin(admin.ModelAdmin):
    list_display = ('doctorName',)

class SupiciousHistoryAdmin(admin.ModelAdmin):
    list_display = ('clinicSample_id', 'patient_id', 'description')

class ResultParametersAdmin(admin.ModelAdmin):
    list_display = ('sampleType_id', 'parameterName', 'parameterDescription', 'parameterOrder','parameterUsed', 'parameterMaxValue','parameterMinValue')



admin.site.register(ClinicSampleRequest, ClinicSampleRequestAdmin)
admin.site.register(ClinicSampleState,ClinicSampleStateAdmin)
admin.site.register(Doctor,DoctorAdmin)
admin.site.register(SupiciousHistory,SupiciousHistoryAdmin)
admin.site.register(ResultParameters,ResultParametersAdmin)
admin.site.register(Patient,PatientAdmin)
admin.site.register(ServiceUnits,ServiceUnitsAdmin)
