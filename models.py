from django.db import models
from iSkyLIMS_core.models import Samples
import datetime


class Patient (models.Model):
    patientName = models.CharField(max_length = 80)
    numberOfHistory = models.CharField(max_length = 80)

    def __str_ (self):
        return '%s' %(self.patientName)

class ServiceUnits (models.Model) :
    serviceUnitName = models.CharField(max_length = 80)

    def __str_ (self):
        return '%s' %(self.serviceUnitName)

class Doctor (models.Model):
    doctorName = models.CharField(max_length = 80)

    def __str_ (self):
        return '%s' %(self.doctorName)

class ClinicSampleState (models.Model):
    clinicState = models.CharField(max_length = 20)

    def __str__ (self):
        return '%s' %(self.clinicState)


class ClinicSampleRequest (models.Model):
    sampleCore = models.ForeignKey(
                Samples,
                on_delete = models.CASCADE)
    patient_id = models.ForeignKey(
                Patient,
                on_delete = models.CASCADE, null = True)
    doctor_id = models.ForeignKey(
                Doctor,
                on_delete = models.CASCADE, null = True)
    clinicSampleState = models.ForeignKey(
                ClinicSampleState,
                on_delete = models.CASCADE)
    orderInEntry = models.CharField(max_length = 8, null = True)
    confirmationCode = models.CharField(max_length = 80)
    priority = models.CharField(max_length = 10)
    coments = models.CharField(max_length = 255)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str_ (self):
        return '%s' %(self.sampleCore.get_name())

    def get_id (self):
        return '%s' %(self.id)

    def get_sample_name(self):
        return '%s' %(self.sampleCore.get_sample_name())

class SupiciousHistory (models.Model):
    clinicSample_id = models.ForeignKey(
            ClinicSampleRequest,
            on_delete = models.CASCADE, null = True)
    patient_id = models.ForeignKey(
            Patient,
            on_delete = models.CASCADE, null = True)
    description = models.CharField(max_length = 255)

    def __str_ (self):
        return '%s' %(self.description)
