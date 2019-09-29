from django.db import models
from iSkyLIMS_core.models import Samples, SampleType
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

    def get_name (self):
        return '%s' %(self.serviceUnitName)

class Doctor (models.Model):
    serviceUnit_id = models.ForeignKey(
                ServiceUnits,
                on_delete = models.CASCADE, null = True, blank = True)
    doctorName = models.CharField(max_length = 80)

    def __str_ (self):
        return '%s' %(self.doctorName)

    def get_name (self):
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
    serviceUnit_id = models.ForeignKey(
                ServiceUnits,
                on_delete = models.CASCADE, null = True, blank = True)
    orderInEntry = models.CharField(max_length = 8, null = True)
    confirmationCode = models.CharField(max_length = 80)
    priority = models.CharField(max_length = 10)
    coments = models.CharField(max_length = 255)
    serviceDate = models.DateTimeField(auto_now_add=False)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str_ (self):
        return '%s' %(self.sampleCore.get_name())

    def get_id (self):
        return '%s' %(self.id)

    def get_sample_name(self):
        return '%s' %(self.sampleCore.get_sample_name())

    def update(self, patient_data):
        self.clinicSampleState = ClinicSampleState.objects.get(clinicState__exact = 'Patient update')
        self.orderInEntry = patient_data['orderInEntry']
        self.confirmationCode = patient_data['confirmationCode']
        self.priority = patient_data['priority']
        self.coments = patient_data['coments']
        self.serviceDate = patient_data['serviceDate']
        self.doctor_id = patient_data['doctor_id']
        self.patient_id = patient_data['patient_id']
        self.serviceUnit_id = patient_data['serviceUnit_id']
        self.save()
        return self

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

class SampleResults (models.Model):
    sampleRequest_id = models.ForeignKey(
            ClinicSampleRequest,
            on_delete = models.CASCADE, null = True)
    sampleType_id =models.ForeignKey(
            SampleType,
            on_delete = models.CASCADE, null = True)


class ResultParametersManager(models.Manager) :
    def create_protocol_parameter (self, prot_param_data):
        new_prot_parameter = self.create(protocol_id =prot_param_data['protocol_id'],parameterName = prot_param_data['Parameter name'],
                    parameterDescription = prot_param_data['Description'], parameterOrder = prot_param_data['Order'],
                    parameterUsed = prot_param_data['Used'], parameterMaxValue = prot_param_data['Max Value'],
                    parameterMinValue = prot_param_data['Min Value'] )
        return new_prot_parameter

class ResultParameters (models.Model):
    sampleType_id = models.ForeignKey(
                SampleType,
                on_delete = models.CASCADE, null = True, blank = True)


    parameterName = models.CharField(max_length=255)
    parameterDescription = models.CharField(max_length= 400, null=True, blank=True)
    parameterOrder = models.IntegerField()
    parameterUsed = models.BooleanField()
    parameterMaxValue = models.CharField(max_length = 50, null = True, blank = True)
    parameterMinValue = models.CharField(max_length = 50, null = True, blank = True)

    def __str__ (self):
        return "%s" %(self.parameterName)

    def get_parameter_name (self):
        return "%s"  %(self.parameterName)

    def get_all_parameter_info(self):
        param_info = []
        param_info.append(self.parameterName)
        param_info.append(self.parameterOrder)
        param_info.append(self.parameterUsed)
        param_info.append(self.parameterMinValue)
        param_info.append(self.parameterMaxValue)
        param_info.append(self.parameterDescription)
        return param_info

    objects = ResultParametersManager()

class ResultParameterValueManager (models.Manager):
    def create_result_parameter_value (self, parameter_value):
        new_result_parameter_data = self.create(resultParameter_id = parameter_value['resultParameter_id'],
                result_id = parameter_value['result_id'],
                parameterValue = parameter_value['parameterValue'])
        return new_result_parameter_data

class ResultParameterValue (models.Model):
    resultParameter_id = models.ForeignKey(
                    ResultParameters,
                    on_delete= models.CASCADE)
    result_id = models.ForeignKey(
                SampleResults,
                on_delete= models.CASCADE)
    parameterValue = models.CharField(max_length=255)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__ (self):
        return '%s' %(self.parameterValue)

    def get_param_value(self):
        return '%s' %(self.parameterValue)

    objects = ResultParameterValueManager()
