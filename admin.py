from django.contrib import admin
from iSkyLIMS_clinic.models import *

# Register your models here.
class ClinicSampleRequestAdmin (admin.ModelAdmin):
    list_display= ('sampleCore', 'clinicSampleState', 'sampleRequestUser','patient_id', 'doctor_id', 'protocol_id', 'confirmationCode', 'priority', 'comments')

class ClinicSampleStateAdmin (admin.ModelAdmin):
    list_display = ('clinicState',)

class ServiceUnitsAdmin(admin.ModelAdmin):
    list_display = ('serviceUnitName',)

class PatientAdmin(admin.ModelAdmin):
    list_display = ('patientName', 'numberOfHistory')

class DoctorAdmin(admin.ModelAdmin):
    list_display = ('doctorName',)

class SuspicionHistoryAdmin(admin.ModelAdmin):
    list_display = ('clinicSample_id', 'patient_id', 'description')



admin.site.register(ClinicSampleRequest, ClinicSampleRequestAdmin)
admin.site.register(ClinicSampleState,ClinicSampleStateAdmin)
admin.site.register(Doctor,DoctorAdmin)
admin.site.register(SuspicionHistory,SuspicionHistoryAdmin)

admin.site.register(Patient,PatientAdmin)
admin.site.register(ServiceUnits,ServiceUnitsAdmin)
