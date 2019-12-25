from django.db import models
from iSkyLIMS_core.models import Samples, SampleType, Protocols, ProtocolParameters, PatientCore

from django.contrib.auth.models import User
import datetime


class PatientHistoryManager(models.Manager):
    def create_patient_history ( self, pat_hist_data):
        new_suspicion_history = self.create( patiendData_id = pat_hist_data['patiendData_id'],
                    entryDate = pat_hist_data['entryDate'], description = pat_hist_data['description'] )




class PatientData (models.Model):
    patienCore =  models.OneToOneField(
            PatientCore,
            on_delete=models.CASCADE,
            primary_key=True, )
    address = models.CharField(max_length = 255, null = True, blank = True)
    phone = models.CharField(max_length = 20, null = True, blank = True)
    email = models.CharField(max_length=50, null = True, blank = True)
    sex = models.CharField(max_length = 20, null = True, blank = True)
    birthday = models.DateTimeField(auto_now_add=False, null = True, blank = True)
    smoker =  models.CharField(max_length = 20, null = True, blank = True)
    notificationPreference =  models.CharField(max_length = 20, null = True, blank = True)

    def __str__ (self):
        return '%s' %(self.patienCore)

    def get_patient_code(self):
        return '%s' %(self.patienCore_id.get_patient_code())

    def get_patient_name(self):
        return '%s' %(self.patienCore_id.get_name())

    def get_patient_full_data(self):
        patient_data = []
        patient_data.append(self.address)
        patient_data.append(self.phone)
        patient_data.append(self.email)
        patient_data.append(self.sex)
        patient_data.append(self.birthday.strftime("%d , %B , %Y"))
        patient_data.append(self.smoker)
        patient_data.append(self.notificationPreference)

        return patient_data


class PatientHistory (models.Model):
    patiendData_id = models.ForeignKey (
            PatientData,
            on_delete = models.CASCADE, null = True, blank = True)
    entryDate = models.DateTimeField(auto_now_add=False, null = True, blank = True)
    description = models.CharField(max_length = 255)

    def __str__ (self):
        return '%s' %(self.description)

    def get_history_text(self):
        return '%s' %(self.description)

    objects = PatientHistoryManager()


class ServiceUnits (models.Model) :
    serviceUnitName = models.CharField(max_length = 80)

    def __str__ (self):
        return '%s' %(self.serviceUnitName)

    def get_name (self):
        return '%s' %(self.serviceUnitName)

class Doctor (models.Model):
    serviceUnit_id = models.ForeignKey(
                ServiceUnits,
                on_delete = models.CASCADE, null = True, blank = True)
    doctorName = models.CharField(max_length = 80)

    def __str__ (self):
        return '%s' %(self.doctorName)

    def get_name (self):
        return '%s' %(self.doctorName)

class ClinicSampleState (models.Model):
    clinicState = models.CharField(max_length = 20)

    def __str__ (self):
        return '%s' %(self.clinicState)

    def get_state(self):
        return '%s' %(self.clinicState)


class ClinicSampleRequest (models.Model):
    sampleCore = models.ForeignKey(
                Samples,
                on_delete = models.CASCADE)
    patientCore = models.ForeignKey(
                PatientCore,
                on_delete = models.CASCADE, null = True, blank = True)
    doctor_id = models.ForeignKey(
                Doctor,
                on_delete = models.CASCADE, null = True)
    clinicSampleState = models.ForeignKey(
                ClinicSampleState,
                on_delete = models.CASCADE)
    serviceUnit_id = models.ForeignKey(
                ServiceUnits,
                on_delete = models.CASCADE, null = True, blank = True)
    sampleRequestUser = models.ForeignKey(
                User,
                on_delete=models.CASCADE, null = True, blank = True)
    orderInEntry = models.CharField(max_length = 8, null = True)
    confirmationCode = models.CharField(max_length = 80, null = True, blank = True)
    priority = models.IntegerField(null = True, blank =True)

    #priority = models.CharField(max_length = 10, null = True, blank = True)
    comments = models.CharField(max_length = 255, null = True, blank = True)
    serviceDate = models.DateTimeField(auto_now_add=False, null = True, blank = True)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__ (self):
        return '%s' %(self.sampleCore)

    def get_id (self):
        return '%s' %(self.id)

    def get_comments(self):
        return '%s'  %(self.comments)

    def get_sample_name(self):
        return '%s' %(self.sampleCore.get_sample_name())

    def get_core_sample_obj(self):
        return self.sampleCore

    def get_patient_information(self):

        patient_info = []
        patient_info.append(self.sampleCore.get_sample_patient_name())
        patient_info.append(self.sampleCore.get_sample_patient_surname())
        patient_info.append(self.sampleCore.get_sample_patient_code())

        return patient_info

    def get_protocol(self):
        if self.protocol_id != None:
            return '%s' %(self.protocol_id.get_name())
        else:
            return 'Not Defined'

    def get_protocol_obj(self):
        return self.protocol_id

    def get_requested_by_information(self):
        requested_by = []
        if self.serviceUnit_id != None:
            requested_by.append(self.serviceUnit_id.get_name())
        else:
            requested_by.append('Not Defined')
        if self.doctor_id != None:
            requested_by.append(self.doctor_id.get_name())
        else:
            requested_by.append('Not Defined')

        return requested_by

    def get_sample_core_info(self):
        s_core_info = []
        s_core_info.append(self.orderInEntry)
        s_core_info.append(self.confirmationCode)
        s_core_info.append(self.priority)
        s_core_info.append(self.sampleCore.get_sample_origin())
        s_core_info.append(self.sampleCore.get_sample_type())
        s_core_info.append(self.sampleCore.get_species())
        s_core_info.append(self.sampleCore.get_extraction_date())
        s_core_info.append(self.sampleCore.get_register_user())

        return s_core_info

    def get_info_for_defined_state(self):
        s_core_info = []
        s_core_info.append(self.sampleCore.get_sample_name())
        s_core_info.append(self.sampleCore.get_sample_origin())
        s_core_info.append(self.sampleCore.get_sample_type())
        s_core_info.append(self.sampleCore.get_species())
        s_core_info.append(self.sampleCore.get_extraction_date())
        s_core_info.append(self.pk)
        return s_core_info

    def get_info_for_patient_sequencing_state(self):
        s_core_info = []
        s_core_info.append(self.sampleCore.get_sample_name())
        s_core_info.append(self.orderInEntry)
        s_core_info.append(self.confirmationCode)
        s_core_info.append(self.priority)
        s_core_info.append(self.patient_id.get_history_number())
        s_core_info.append(self.doctor_id.get_name())
        s_core_info.append(self.serviceUnit_id.get_name())
        return s_core_info


    def get_sample_info_for_list (self):
        s_core_info = []
        s_core_info.append(self.pk)
        s_core_info.append(self.sampleCore.get_sample_name())
        try:
            s_core_info.append(self.patient_id.get_history_number())
        except:
            s_core_info.append('History Number not available')
        try:
            s_core_info.append(self.doctor_id.get_name())
        except:
            s_core_info.append('Doctor name not available')
        s_core_info.append(self.priority)
        try:
            s_core_info.append(self.serviceDate.strftime('%Y-%m-%d'))
        except:
            s_core_info.append('Not available')
        s_core_info.append(self.clinicSampleState.get_state())
        return s_core_info

    def get_sample_info_for_protocol (self):
        s_core_info = []
        s_core_info.append(self.sampleCore.get_sample_name())
        s_core_info.append(self.patient_id.get_history_number())
        s_core_info.append(self.priority)
        return s_core_info

    def get_sample_core_state (self):
        return '%s' %(self.sampleCore.get_sample_state())

    def set_protocol (self, protocol_name):
        self.protocol_id = Protocols.objects.get(name__exact = protocol_name)
        self.save()
        return self

    def set_state(self, new_state):
        self.clinicSampleState = ClinicSampleState.objects.get(clinicState__exact = new_state)
        self.save()
        return self

    def update(self, patient_data):
        self.clinicSampleState = ClinicSampleState.objects.get(clinicState__exact = 'Patient update')
        self.orderInEntry = patient_data['orderInEntry']
        self.confirmationCode = patient_data['confirmationCode']
        self.priority = patient_data['priority']
        self.comments = patient_data['comments']
        self.serviceDate = patient_data['serviceDate']
        self.doctor_id = patient_data['doctor_id']
        self.patient_id = patient_data['patient_id']
        self.serviceUnit_id = patient_data['serviceUnit_id']
        self.save()
        return self


'''
class SampleResults (models.Model):
    sampleRequest_id = models.ForeignKey(
            ClinicSampleRequest,
            on_delete = models.CASCADE, null = True)
    sampleType_id =models.ForeignKey(
            SampleType,
            on_delete = models.CASCADE, null = True)



class ResultParametersManager(models.Manager) :
    def create_result_parameter (self, prot_param_data):
        new_result_parameter = self.create(parameterName = prot_param_data['Parameter name'],
                    parameterDescription = prot_param_data['Description'], parameterOrder = prot_param_data['Order'],
                    parameterUsed = prot_param_data['Used'], parameterMaxValue = prot_param_data['Max Value'],
                    parameterMinValue = prot_param_data['Min Value'] )
        return new_result_parameter

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
'''
class ResultParameterValueManager (models.Manager):
    def create_result_parameter_value (self, parameter_value):
        new_result_parameter_data = self.create(parameter_id = parameter_value['parameter_id'],
                clinicSample_id = parameter_value['clinicSample_id'],
                parameterValue = parameter_value['parameterValue'])
        return new_result_parameter_data

class ResultParameterValue (models.Model):
    '''
    resultParameter_id = models.ForeignKey(
                    ResultParameters,
                    on_delete= models.CASCADE)

    result_id = models.ForeignKey(
                SampleResults,
                on_delete= models.CASCADE)
    '''
    clinicSample_id = models.ForeignKey(
            ClinicSampleRequest,
            on_delete = models.CASCADE, null = True)
    parameter_id = models.ForeignKey(
            ProtocolParameters,
            on_delete= models.CASCADE, null = True)
    parameterValue = models.CharField(max_length=255)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__ (self):
        return '%s' %(self.parameterValue)

    def get_parameter_value(self):
        return '%s' %(self.parameterValue)

    objects = ResultParameterValueManager()
