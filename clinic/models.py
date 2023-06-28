# Generic imports
from django.contrib.auth.models import User
from django.db import models

# Local imports
import core.models


class PatientDataManager(models.Manager):
    def create_patient_opt_data(self, p_opt_data):
        new_p_opt_data = self.create(
            patien_core=p_opt_data["patienCore"],
            address=p_opt_data["address"],
            phone=p_opt_data["phone"],
            email=p_opt_data["email"],
            birthday=p_opt_data["birthday"],
            smoker=p_opt_data["smoker"],
            notification_preference=p_opt_data["notificationPreference"],
            comments=p_opt_data["comments"],
        )
        return new_p_opt_data


class PatientData(models.Model):
    patien_core = models.OneToOneField(
        core.models.PatientCore,
        on_delete=models.CASCADE,
        primary_key=True,
    )
    address = models.CharField(max_length=255, null=True, blank=True)
    phone = models.CharField(max_length=20, null=True, blank=True)
    email = models.CharField(max_length=50, null=True, blank=True)
    birthday = models.DateTimeField(auto_now_add=False, null=True, blank=True)
    smoker = models.CharField(max_length=20, null=True, blank=True)
    notification_preference = models.CharField(max_length=20, null=True, blank=True)
    comments = models.CharField(max_length=255, null=True, blank=True)

    def __str__(self):
        return "%s" % (self.patien_core)

    def get_patient_code(self):
        return "%s" % (self.patienCore_id.get_patient_code())

    def get_patient_name(self):
        return "%s" % (self.patienCore_id.get_name())

    def get_patient_full_data(self):
        patient_data = []
        patient_data.append(self.address)
        patient_data.append(self.phone)
        patient_data.append(self.email)
        patient_data.append(self.birthday.strftime("%d , %B , %Y"))
        patient_data.append(self.smoker)
        patient_data.append(self.notification_preference)
        patient_data.append(self.comments)
        return patient_data

    objects = PatientDataManager()


class PatientHistoryManager(models.Manager):
    def create_patient_history(self, pat_hist_data):
        self.create(
            patient_core=pat_hist_data["patientCore"],
            entry_date=pat_hist_data["entryDate"],
            description=pat_hist_data["description"],
        )


class PatientHistory(models.Model):
    patient_core = models.ForeignKey(
        core.models.PatientCore, on_delete=models.CASCADE, null=True, blank=True
    )
    entry_date = models.DateTimeField(auto_now_add=False, null=True, blank=True)
    description = models.CharField(max_length=255)

    def __str__(self):
        return "%s" % (self.description)

    def get_history_text(self):
        return "%s" % (self.description)

    objects = PatientHistoryManager()


class ServiceUnits(models.Model):
    service_unit_name = models.CharField(max_length=80)

    def __str__(self):
        return "%s" % (self.service_unit_name)

    def get_name(self):
        return "%s" % (self.service_unit_name)


class Doctor(models.Model):
    service_unit_id = models.ForeignKey(
        ServiceUnits, on_delete=models.CASCADE, null=True, blank=True
    )
    doctor_name = models.CharField(max_length=80)

    def __str__(self):
        return "%s" % (self.doctor_name)

    def get_name(self):
        return "%s" % (self.doctor_name)


class ClinicSampleState(models.Model):
    clinic_state = models.CharField(max_length=20)

    def __str__(self):
        return "%s" % (self.clinic_state)

    def get_state(self):
        return "%s" % (self.clinic_state)


class ClinicSampleRequestManager(models.Manager):
    def create_clinic_sample(self, c_sample_data):
        new_clinic_sample = self.create(
            sample_core=c_sample_data["sampleCore"],
            patient_core=c_sample_data["patientCore"],
            clinic_sample_state=ClinicSampleState.objects.get(
                clinic_state__exact=c_sample_data["state"]
            ),
            sample_request_user=c_sample_data["user"],
        )
        return new_clinic_sample


class ClinicSampleRequest(models.Model):
    sample_core = models.ForeignKey(core.models.Samples, on_delete=models.CASCADE)
    patient_core = models.ForeignKey(
        core.models.PatientCore, on_delete=models.CASCADE, null=True, blank=True
    )
    doctor_id = models.ForeignKey(
        Doctor, on_delete=models.CASCADE, null=True, blank=True
    )
    clinic_sample_state = models.ForeignKey(ClinicSampleState, on_delete=models.CASCADE)
    service_unit_id = models.ForeignKey(
        ServiceUnits, on_delete=models.CASCADE, null=True, blank=True
    )
    sample_request_user = models.ForeignKey(
        User, on_delete=models.CASCADE, null=True, blank=True
    )
    entry_order = models.CharField(max_length=8, null=True, blank=True)
    confirmation_code = models.CharField(max_length=80, null=True, blank=True)
    priority = models.IntegerField(null=True, blank=True)
    comments = models.CharField(max_length=255, null=True, blank=True)
    service_date = models.DateTimeField(auto_now_add=False, null=True, blank=True)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return "%s" % (self.sample_core)

    def get_id(self):
        return "%s" % (self.id)

    def get_comments(self):
        return "%s" % (self.comments)

    def get_sample_name(self):
        return "%s" % (self.sample_core.get_sample_name())

    def get_core_sample_obj(self):
        return self.sample_core

    def get_patient_information(self):
        patient_info = []
        patient_info.append(self.sample_core.get_sample_patient_name())
        patient_info.append(self.sample_core.get_sample_patient_surname())
        patient_info.append(self.sample_core.get_sample_patient_code())

        return patient_info

    def get_protocol(self):
        if self.protocol_id is not None:
            return "%s" % (self.protocol_id.get_name())
        else:
            return "Not Defined"

    def get_protocol_obj(self):
        return self.protocol_id

    def get_requested_by_information(self):
        requested_by = []
        if self.service_unit_id is not None:
            requested_by.append(self.service_unit_id.get_name())
        else:
            requested_by.append("Not Defined")
        if self.doctor_id is not None:
            requested_by.append(self.doctor_id.get_name())
        else:
            requested_by.append("Not Defined")

        return requested_by

    def get_sample_core_info(self):
        s_core_info = []
        s_core_info.append(self.sample_core.get_sample_type())
        s_core_info.append(self.sample_core.get_species())
        s_core_info.append(self.sample_core.get_entry_date())
        s_core_info.append(self.sample_core.get_register_user())

        return s_core_info

    def get_info_for_defined_state(self):
        s_core_info = []
        s_core_info.append(self.sample_core.get_sample_name())
        s_core_info.append(self.sample_core.get_sample_type())
        s_core_info.append(self.sample_core.get_species())
        s_core_info.append(self.sample_core.get_entry_date())
        s_core_info.append(self.pk)
        return s_core_info

    def get_info_for_patient_sequencing_state(self):
        s_core_info = []
        s_core_info.append(self.sample_core.get_sample_name())
        s_core_info.append(self.orderInEntry)
        s_core_info.append(self.confirmation_code)
        s_core_info.append(self.priority)
        s_core_info.append(self.doctor_id.get_name())
        s_core_info.append(self.service_unit_id.get_name())
        return s_core_info

    def get_sample_info_for_list(self):
        s_core_info = []
        s_core_info.append(self.pk)
        s_core_info.append(self.sample_core.get_sample_name())
        try:
            s_core_info.append(self.doctor_id.get_name())
        except Exception:
            s_core_info.append("Doctor name not available")
        s_core_info.append(self.priority)
        try:
            s_core_info.append(self.service_date.strftime("%Y-%m-%d"))
        except Exception:
            s_core_info.append("Not available")
        s_core_info.append(self.clinic_sample_state.get_state())
        return s_core_info

    def get_sample_info_for_protocol(self):
        s_core_info = []
        s_core_info.append(self.sample_core.get_sample_name())
        s_core_info.append(self.priority)
        return s_core_info

    def get_sample_core_state(self):
        return "%s" % (self.sample_core.get_sample_state())

    def set_protocol(self, protocol_name):
        self.protocol_id = core.models.Protocols.objects.get(name__exact=protocol_name)
        self.save()
        return self

    def set_state(self, new_state):
        self.clinic_sample_state = ClinicSampleState.objects.get(
            clinicState__exact=new_state
        )
        self.save()
        return self

    def update(self, patient_data):
        self.clinic_sample_state = ClinicSampleState.objects.get(
            clinicState__exact="Patient update"
        )
        self.orderInEntry = patient_data["orderInEntry"]
        self.confirmation_code = patient_data["confirmationCode"]
        self.priority = patient_data["priority"]
        self.comments = patient_data["comments"]
        self.service_date = patient_data["serviceDate"]
        self.doctor_id = patient_data["doctor_id"]
        self.patient_id = patient_data["patient_id"]
        self.service_unit_id = patient_data["serviceUnit_id"]
        self.save()
        return self

    objects = ClinicSampleRequestManager()


class ResultParameterValueManager(models.Manager):
    def create_result_parameter_value(self, parameter_value):
        new_result_parameter_data = self.create(
            parameter_id=parameter_value["parameter_id"],
            clinicSample_id=parameter_value["clinicSample_id"],
            parameterValue=parameter_value["parameterValue"],
        )
        return new_result_parameter_data


class ResultParameterValue(models.Model):
    clinicSample_id = models.ForeignKey(
        ClinicSampleRequest, on_delete=models.CASCADE, null=True
    )
    parameter_id = models.ForeignKey(
        core.models.ProtocolParameters, on_delete=models.CASCADE, null=True
    )
    parameterValue = models.CharField(max_length=255)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return "%s" % (self.parameterValue)

    def get_parameter_value(self):
        return "%s" % (self.parameterValue)

    objects = ResultParameterValueManager()


class FamilyRelatives(models.Model):
    relationship = models.CharField(max_length=255)

    def __str__(self):
        return "%s" % (self.relationship)

    def get_relation(self):
        return "%s" % (self.relationship)


class Family(models.Model):
    patien_core_id = models.ForeignKey(
        core.models.PatientCore, on_delete=models.CASCADE, null=True, blank=True
    )
    family_relatives_id = models.ForeignKey(
        FamilyRelatives, on_delete=models.CASCADE, null=True, blank=True
    )
    family_id = models.CharField(max_length=50)
    relative_1 = models.CharField(max_length=50)
    relative_2 = models.CharField(max_length=50)
    affected = models.BooleanField()
    family_comments = models.CharField(max_length=255, null=True, blank=True)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return "%s" % (self.family_id)


class ConfigSetting(models.Model):
    configuration_name = models.CharField(max_length=80)
    configuration_value = models.CharField(max_length=255, null=True, blank=True)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return "%s" % (self.configuration_name)

    def get_configuration_value(self):
        return "%s" % (self.configuration_value)
