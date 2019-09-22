from django.contrib import admin
from iSkyLIMS_clinic.models import *

# Register your models here.
class ClinicSampleRequestAdmin (admin.ModelAdmin):
    list_display= ('sampleCore', 'clinicSampleState' ,'patient_id', 'doctor_id', 'confirmationCode', 'priority', 'coments')

class ClinicSampleStateAdmin (admin.ModelAdmin):
    list_display = ('clinicState',)

admin.site.register(ClinicSampleRequest, ClinicSampleRequestAdmin)
admin.site.register(ClinicSampleState,ClinicSampleStateAdmin)
