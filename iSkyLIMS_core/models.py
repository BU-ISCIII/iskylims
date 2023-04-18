import datetime

from django.contrib.auth.models import User
from django.db import models


class StateInCountryManager(models.Manager):
    def create_new_state(self, data):
        new_state = self.create(state_name=data["state"], apps_name=data["apps_name"])
        return new_state


class StateInCountry(models.Model):
    state_name = models.CharField(max_length=80)
    apps_name = models.CharField(max_length=40, null=True)

    def __str__(self):
        return "%s" % (self.state_name)

    def get_state_name(self):
        return "%s" % (self.state_name)

    def get_state_id(self):
        return "%s" % (self.pk)

    objects = StateInCountryManager()


class Contact(models.Model):
    contact_name = models.CharField(max_length=80)
    contact_mail = models.CharField(max_length=40, null=True)

    def __str__(self):
        return "%s" % (self.contact_name)

    def get_contact_name(self):
        return "%s" % (self.contact_name)

    def get_contact_email(self):
        return "%s" % (self.contact_mail)


class CityManager(models.Manager):
    def create_new_city(self, data):
        if StateInCountry.objects.filter(pk__exact=data["state"]).exists():
            state_obj = StateInCountry.objects.filter(pk__exact=data["state"]).last()
        else:
            state_obj = None
        new_city = self.create(
            belongs_to_state=state_obj,
            city_name=data["cityName"],
            geo_loc_latitude=data["latitude"],
            geo_loc_longitude=data["longitude"],
            apps_name=data["apps_name"],
        )
        return new_city


class City(models.Model):
    belongs_to_state = models.ForeignKey(
        StateInCountry, on_delete=models.CASCADE, null=True, blank=True
    )
    city_name = models.CharField(max_length=80)
    geo_loc_latitude = models.CharField(max_length=80)
    geo_loc_longitude = models.CharField(max_length=80)
    apps_name = models.CharField(max_length=40, null=True)

    def __str__(self):
        return "%s" % (self.city_name)

    def get_city_name(self):
        return "%s" % (self.city_name)

    def get_city_id(self):
        return "%s" % (self.pk)

    def get_coordenates(self):
        return {"latitude": self.geo_loc_latitude, "longitude": self.geo_loc_longitude}

    def get_state(self):
        if self.belongs_to_state is None:
            return ""
        return "%s" % (self.belongs_to_state.get_state_name())

    objects = CityManager()


class LabRequestManager(models.Manager):
    def create_lab_request(self, data):
        city_obj = City.objects.filter(pk__exact=data["city"]).last()
        new_lab_request = self.create(
            lab_name=data["labName"],
            lab_name_coding=data["labNameCoding"],
            lab_unit=data["labUnit"],
            lab_contact_name=data["labContactName"],
            lab_phone=data["labPhone"],
            lab_email=data["labEmail"],
            address=data["address"],
            apps_name=data["apps_name"],
            lab_city=city_obj,
        )
        return new_lab_request


class LabRequest(models.Model):
    lab_city = models.ForeignKey(City, on_delete=models.CASCADE, null=True, blank=True)
    lab_name = models.CharField(max_length=80)
    lab_name_coding = models.CharField(max_length=50)
    lab_unit = models.CharField(max_length=50)
    lab_contact_name = models.CharField(max_length=50)
    lab_phone = models.CharField(max_length=20)
    lab_email = models.CharField(max_length=70)
    address = models.CharField(max_length=255)

    apps_name = models.CharField(max_length=40, null=True)

    def __str__(self):
        return "%s" % (self.lab_name)

    def get_name(self):
        return "%s" % (self.lab_name)

    def get_id(self):
        return "%s" % (self.pk)

    def get_lab_request_code(self):
        return "%s" % (self.lab_name_coding)

    def get_all_data(self):
        data = []
        data.append(self.lab_name)
        data.append(self.lab_name_coding)
        data.append(self.lab_unit)
        data.append(self.lab_contact_name)
        data.append(self.lab_phone)
        data.append(self.lab_email)
        data.append(self.address)
        return data

    def get_fields_and_data(self):
        data = {}
        if self.lab_city is None:
            data["state"] = ""
            data["city"] = ""
            data["latitude"] = ""
            data["longitude"] = ""
        else:
            data["state"] = self.lab_city.get_state()
            data["city"] = self.lab_city.get_city_name()
            data.update(self.lab_city.get_coordenates())
        data["Lab Name"] = self.lab_name
        data["Lab Coding"] = self.lab_name_coding
        return data

    def get_region(self):
        if self.lab_city is None:
            return ""
        else:
            return "%s" % (self.lab_city.get_state())

    objects = LabRequestManager()


class MoleculeTypeManager(models.Manager):
    def create_molecule_type(self, data):
        new_molecule_type = self.create(
            molecule_type=data["moleculeType"], apps_name=data["apps_name"]
        )
        return new_molecule_type


class MoleculeType(models.Model):
    molecule_type = models.CharField(max_length=30)
    apps_name = models.CharField(max_length=40, null=True)

    def __str__(self):
        return "%s" % (self.molecule_type)

    def get_name(self):
        return "%s" % (self.molecule_type)

    def get_id(self):
        return "%s" % (self.pk)

    objects = MoleculeTypeManager()


class ProtocolTypeManager(models.Manager):
    def create_protocol_type(self, data):
        if data["molecule"] is not None:
            molecule_obj = MoleculeType.objects.filter(
                molecule_type__iexact=data["molecule"],
                apps_name__exact=data["apps_name"],
            ).last()
        else:
            molecule_obj = None
        new_protocol_type = self.create(
            molecule=molecule_obj,
            protocol_type=data["protocol_type"],
            apps_name=data["apps_name"],
        )
        return new_protocol_type


class ProtocolType(models.Model):
    molecule = models.ForeignKey(
        MoleculeType, on_delete=models.CASCADE, null=True, blank=True
    )
    protocol_type = models.CharField(max_length=40)
    apps_name = models.CharField(max_length=40)

    def __str__(self):
        return "%s" % (self.protocol_type)

    def get_molecule_type(self):
        if self.molecule is None:
            return ""
        else:
            return "%s" % (self.molecule.get_name())

    def get_name(self):
        return "%s" % (self.protocol_type)

    def get_id(self):
        return "%s" % (self.pk)

    objects = ProtocolTypeManager()


class Protocols(models.Model):
    type = models.ForeignKey(ProtocolType, on_delete=models.CASCADE)
    name = models.CharField(max_length=40)
    description = models.CharField(max_length=160, null=True, blank=True)

    def __str__(self):
        return "%s" % (self.name)

    def get_name(self):
        return "%s" % (self.name)

    def get_type(self):
        return "%s" % (self.type.get_name())

    def get_protocol_id(self):
        return "%s" % (self.pk)


class ProtocolParametersManager(models.Manager):
    def create_protocol_parameter(self, prot_param_data):
        new_prot_parameter = self.create(
            protocol_id=prot_param_data["protocol_id"],
            parameter_name=prot_param_data["Parameter name"],
            parameter_description=prot_param_data["Description"],
            parameter_order=prot_param_data["Order"],
            parameter_used=prot_param_data["Used"],
            parameter_max_value=prot_param_data["Max Value"],
            parameter_min_value=prot_param_data["Min Value"],
            parameter_option_values=prot_param_data["Option Values"],
            parameter_type=prot_param_data["Parameter Type"],
        )
        return new_prot_parameter


class ProtocolParameters(models.Model):
    protocol_id = models.ForeignKey(Protocols, on_delete=models.CASCADE)
    parameter_name = models.CharField(max_length=255)
    parameter_description = models.CharField(max_length=400, null=True, blank=True)
    parameter_order = models.IntegerField()
    parameter_used = models.BooleanField()
    parameter_type = models.CharField(max_length=20, default="string")
    parameter_option_values = models.CharField(max_length=400, null=True, blank=True)
    parameter_max_value = models.CharField(max_length=50, null=True, blank=True)
    parameter_min_value = models.CharField(max_length=50, null=True, blank=True)

    def __str__(self):
        return "%s" % (self.parameter_name)

    def get_parameter_name(self):
        return "%s" % (self.parameter_name)

    def get_parameter_protocol_id(self):
        return "%s" % (self.pk)

    def get_parameter_option_values(self):
        return "%s" % (self.parameter_option_values)

    def get_parameter_type(self):
        return "%s" % (self.parameter_type)

    def get_all_parameter_info(self):
        param_info = []
        param_info.append(self.parameter_name)
        param_info.append(self.parameter_order)
        param_info.append(self.parameter_used)
        param_info.append(self.parameter_type)
        param_info.append(self.parameter_option_values)
        param_info.append(self.parameter_min_value)
        param_info.append(self.parameter_max_value)
        param_info.append(self.parameter_description)
        return param_info

    def get_protocol_fields_for_javascript(self):
        if self.parameter_used:
            used = "true"
        else:
            used = "false"
        if self.parameter_option_values is None:
            parameter_option_values = ""
        else:
            parameter_option_values = self.parameter_option_values
        field_data = []
        field_data.append(self.parameter_name)

        field_data.append(self.parameter_order)
        field_data.append(used)
        field_data.append(self.parameter_type)
        field_data.append(parameter_option_values)
        field_data.append(self.parameter_description)
        return field_data

    def update_protocol_fields(self, prot_param_data):
        self.parameter_name = prot_param_data["Parameter name"]
        self.parameter_description = prot_param_data["Description"]
        self.parameter_order = prot_param_data["Order"]
        self.parameter_used = prot_param_data["Used"]
        self.parameter_option_values = prot_param_data["Option Values"]
        self.parameter_type = prot_param_data["Parameter Type"]
        self.save()

    objects = ProtocolParametersManager()


class StatesForSample(models.Model):
    sample_state_name = models.CharField(max_length=50)

    def __str__(self):
        return "%s" % (self.sample_state_name)

    def get_sample_state(self):
        return "%s" % (self.sample_state_name)

    def get_id(self):
        return "%s" % (self.pk)


class StatesForMolecule(models.Model):
    molecule_state_name = models.CharField(max_length=50)

    def __str__(self):
        return "%s" % (self.molecule_state_name)

    def get_molecule_state(self):
        return "%s" % (self.molecule_state_name)


class SampleTypeManager(models.Manager):
    def create_sample_type(self, sample_type_data):
        new_sample_type = self.create(
            sample_type=sample_type_data["sampleType"],
            apps_name=sample_type_data["apps_name"],
            optional_fields=sample_type_data["optional_fields"],
        )
        return new_sample_type


class SampleType(models.Model):
    sample_type = models.CharField(max_length=50)
    apps_name = models.CharField(max_length=50)
    optional_fields = models.CharField(max_length=50, null=True, blank=True)
    generated_at = models.DateTimeField(auto_now_add=True, null=True)

    def __str__(self):
        return "%s" % (self.sample_type)

    def get_sample_type_id(self):
        return "%s" % (self.pk)

    def get_name(self):
        return "%s" % (self.sample_type)

    def get_optional_values(self):
        if self.optional_fields == "" or self.optional_fields is None:
            return []
        else:
            return list(map(int, self.optional_fields.split(",")))

    objects = SampleTypeManager()


class SpeciesManager(models.Manager):
    def create_new_specie(self, data):
        new_specie = self.create(species_name=data["name"], apps_name=data["apps_name"])
        return new_specie


class Species(models.Model):
    species_name = models.CharField(max_length=50)
    ref_genome_name = models.CharField(max_length=255, null=True, blank=True)
    ref_genome_size = models.CharField(max_length=100, null=True, blank=True)
    ref_genome_id = models.CharField(max_length=255, null=True, blank=True)
    apps_name = models.CharField(max_length=50, null=True)
    generated_at = models.DateTimeField(auto_now_add=True, null=True)

    def __str__(self):
        return "%s" % (self.species_name)

    def get_name(self):
        return "%s" % (self.species_name)

    def get_id(self):
        return "%s" % (self.pk)

    objects = SpeciesManager()


class SequencingPlatform(models.Model):
    platform_name = models.CharField(max_length=30)
    company_name = models.CharField(max_length=30)
    sequencing_technology = models.CharField(max_length=30)

    def __str__(self):
        return "%s" % (self.platform_name)

    def get_platform_name(self):
        return "%s" % (self.platform_name)

    def get_platform_id(self):
        return "%s" % (self.pk)

    def get_company_name(self):
        return "%s" % (self.company_name)


class CommercialKitsManager(models.Manager):
    def create_commercial_kit(self, kit_data):
        if "platform" in kit_data:
            platform_obj = SequencingPlatform.objects.get(
                pk__exact=kit_data["platform"]
            )
        else:
            platform_obj = None
        new_commercial_kit = self.create(
            name=kit_data["name"],
            provider=kit_data["provider"],
            cat_number=kit_data["cat_number"],
            description=kit_data["description"],
            platform_kits=platform_obj,
        )
        return new_commercial_kit


class CommercialKits(models.Model):
    protocol_kits = models.ManyToManyField(Protocols, blank=True)

    platform_kits = models.ForeignKey(
        SequencingPlatform, on_delete=models.CASCADE, null=True, blank=True
    )

    """ Add if finally we use another type of for reagentsKits
    reagentsKits = models.ForeignKey(
                    SequencingPlatform, related_name = 'reagentsKits',
                    on_delete= models.CASCADE, null = True, blank = True)
    """
    name = models.CharField(max_length=150)
    provider = models.CharField(max_length=30)
    cat_number = models.CharField(max_length=40, null=True, blank=True)
    description = models.CharField(max_length=255, null=True, blank=True)
    generated_at = models.DateTimeField(auto_now_add=True, null=True)

    def __str__(self):
        return "%s" % (self.name)

    def get_name(self):
        return "%s" % (self.name)

    def platform_kit_obj(self):
        return self.platform_kits

    def get_platform_name(self):
        if self.platform_kits is not None:
            return "%s" % (self.platform_kits.get_platform_name())
        else:
            return ""

    def get_protocol_objs(self):
        return self.protocol_kits.all()

    def get_provider_kit_name(self):
        return "%s" % (self.provider)

    def get_cat_number(self):
        return "%s" % (self.cat_number)

    def get_commercial_platform_basic_data(self):
        kit_basic_data = []
        kit_basic_data.append(self.name)
        kit_basic_data.append(self.provider)
        kit_basic_data.append(self.platform_kits.get_platform_name())
        return kit_basic_data

    def get_commercial_protocol_basic_data(self):
        kit_basic_data = []
        kit_basic_data.append(self.name)
        kit_basic_data.append(self.provider)
        protocols = []
        protocol_objs = self.protocol_kits.all()
        for protocol_obj in protocol_objs:
            protocols.append(protocol_obj.get_name())
        kit_basic_data.append(protocols)

        return kit_basic_data

    objects = CommercialKitsManager()


class UserLotCommercialKitsManager(models.Manager):
    def create_user_lot_commercial_kit(self, kit_data):
        expiration_date = datetime.datetime.strptime(
            kit_data["expirationDate"], "%Y-%m-%d"
        ).date()

        new_user_lot_commercial_kit = self.create(
            user=kit_data["user"],
            based_commercial=kit_data["basedCommercial"],
            chip_lot=kit_data["chipLot"],
            expiration_date=expiration_date,
        )
        return new_user_lot_commercial_kit


class UserLotCommercialKits(models.Model):
    user = models.ForeignKey(User, on_delete=models.CASCADE, null=True)
    based_commercial = models.ForeignKey(
        CommercialKits, on_delete=models.CASCADE, null=True
    )

    uses_number = models.IntegerField(null=True, default=0)
    chip_lot = models.CharField(max_length=50)
    latest_used_date = models.DateTimeField(null=True, blank=True)
    expiration_date = models.DateField(auto_now_add=False)
    run_out = models.BooleanField(default=False)
    generated_at = models.DateTimeField(auto_now_add=True, null=True)

    def __str__(self):
        return "%s" % (self.chip_lot)

    def get_basic_data(self):
        lot_data = []
        lot_data.append(self.based_commercial.get_name())
        lot_data.append(self.chip_lot)
        lot_data.append(self.expiration_date.strftime("%d %B %Y"))
        return lot_data

    def get_commercial_kit(self):
        return "%s" % (self.based_commercial.get_name())

    def get_commercial_obj(self):
        return self.based_commercial

    def get_lot_number(self):
        return "%s" % (self.chip_lot)

    def get_number_of_uses(self):
        return "%s" % (self.uses_number)

    def get_protocol_for_kit(self):
        return "%s" % (self.based_commercial.get_protocol())

    def get_protocol_obj_for_kit(self):
        return self.based_commercial.get_protocol_obj()

    def get_expiration_date(self):
        exp_date = self.expiration_date
        return "%s" % (exp_date.strftime("%d %B, %Y"))

    def get_user_lot_kit_id(self):
        return "%s" % (self.pk)

    def set_increase_use(self):
        self.latest_used_date = datetime.datetime.now()
        self.uses_number += 1
        self.save()
        return self

    def set_latest_use(self, date):
        self.latest_used_date = date
        self.save()
        return self

    def set_run_out(self):
        self.run_out = True
        self.save()
        return self

    objects = UserLotCommercialKitsManager()


class PatientProjectsManager(models.Manager):
    def create_project(self, project_data):
        new_project = self.create(
            project_mame=project_data["projectName"],
            project_manager=project_data["projectManager"],
            project_contact=project_data["projectContact"],
            project_description=project_data["projectDescription"],
            apps_name=project_data["apps_name"],
        )
        return new_project


class PatientProjects(models.Model):
    project_name = models.CharField(max_length=50)
    project_manager = models.CharField(max_length=50, null=True, blank=True)
    project_contact = models.CharField(max_length=50, null=True, blank=True)
    project_description = models.CharField(max_length=255, null=True, blank=True)
    apps_name = models.CharField(max_length=40)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return "%s" % (self.project_name)

    def get_project_id(self):
        return "%s" % (self.pk)

    def get_project_name(self):
        return "%s" % (self.project_name)

    def get_patient_project_data(self):
        p_data = []
        p_data.append(self.project_name)
        p_data.append(self.project_manager)
        p_data.append(self.project_contact)
        p_data.append(self.project_description)
        p_data.append(self.pk)
        return p_data

    objects = PatientProjectsManager()


class PatientProjectsFieldsManager(models.Manager):
    def create_project_fields(self, project_field_data):
        new_project_field = self.create(
            patient_projects_id=project_field_data["project_id"],
            project_field_name=project_field_data["Field name"],
            project_field_description=project_field_data["Description"],
            project_field_order=project_field_data["Order"],
            project_field_used=project_field_data["Used"],
        )
        return new_project_field


class PatientProjectsFields(models.Model):
    patient_projects_id = models.ForeignKey(
        PatientProjects, on_delete=models.CASCADE, null=True, blank=True
    )
    project_field_name = models.CharField(max_length=50)
    project_field_description = models.CharField(max_length=400, null=True, blank=True)
    project_field_order = models.IntegerField()
    project_field_used = models.BooleanField()

    def __str__(self):
        return "%s" % (self.project_field_name)

    def get_field_id(self):
        return "%s" % (self.id)

    def get_field_name(self):
        return "%s" % (self.project_field_name)

    def get_description(self):
        return "%s" % (self.project_field_description)

    def get_all_fields_info(self):
        if self.project_field_used:
            used = "Yes"
        else:
            used = "No"
        field_data = []
        field_data.append(self.project_field_name)

        field_data.append(self.project_field_order)
        field_data.append(used)
        field_data.append(self.project_field_description)
        return field_data

    objects = PatientProjectsFieldsManager()


class PatientSex(models.Model):
    sex = models.CharField(max_length=16)

    def __str__(self):
        return "%s" % (self.sex)

    def get_patient_sex(self):
        return "%s" % (self.sex)


class PatientCoreManager(models.Manager):
    def create_patient(self, p_data):
        try:
            sex = PatientSex.objects.get(sex__exact=p_data["patientSex"])
        except ValueError:
            sex = None
        new_patient = self.create(
            patient_name=p_data["patientName"],
            patient_surname=p_data["patientSurname"],
            patient_code=p_data["patientCode"],
            patient_sex=sex,
        )
        return new_patient


class PatientCore(models.Model):
    patient_projects = models.ManyToManyField(PatientProjects, blank=True)
    patient_name = models.CharField(max_length=255, null=True)
    patient_surname = models.CharField(max_length=255, null=True)
    patient_code = models.CharField(max_length=255, null=True)
    patient_sex = models.ForeignKey(
        PatientSex, on_delete=models.CASCADE, null=True, blank=True
    )

    def __str__(self):
        return "%s" % (self.patient_code)

    def get_patient_id(self):
        return "%s" % (self.id)

    def get_patient_name(self):
        return "%s" % (self.patient_name)

    def get_patient_surname(self):
        return "%s" % (self.patient_surname)

    def get_patient_code(self):
        return "%s" % (self.patient_code)

    def get_patient_sex(self):
        return "%s" % (self.patient_sex)

    objects = PatientCoreManager()


class PatientProjectFieldValueManager(models.Manager):
    def create_project_field_value(self, field_value):
        new_field_data = self.create(
            project_field_id=field_value["projectField_id"],
            patient_core_id=field_value["patientCore_id"],
            project_field_value=field_value["projectFieldValue"],
        )
        return new_field_data


class PatientProjectFieldValue(models.Model):
    patient_core_id = models.ForeignKey(
        PatientCore, on_delete=models.CASCADE, null=True, blank=True
    )
    project_field_id = models.ForeignKey(
        PatientProjectsFields, on_delete=models.CASCADE, null=True
    )
    project_field_value = models.CharField(max_length=255)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return "%s" % (self.project_field_value)

    def get_field_value(self):
        return "%s" % (self.project_field_value)

    objects = PatientProjectFieldValueManager()


class SampleProjectsManager(models.Manager):
    def create_sample_project(self, s_project_data):
        new_sample_project = self.create(
            sample_project_name=s_project_data["sampleProjectName"],
            sample_project_manager=s_project_data["sampleProjectManager"],
            sample_project_contact=s_project_data["sampleProjectContact"],
            sample_project_description=s_project_data["sampleProjectDescription"],
            apps_name=s_project_data["apps_name"],
        )
        return new_sample_project


class SampleProjects(models.Model):
    sample_project_name = models.CharField(max_length=255)
    sample_project_manager = models.CharField(max_length=50, null=True, blank=True)
    sample_project_contact = models.CharField(max_length=250, null=True, blank=True)
    sample_project_description = models.CharField(max_length=255, null=True, blank=True)
    apps_name = models.CharField(max_length=255)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return "%s" % (self.sample_project_name)

    def get_id(self):
        return "%s" % (self.pk)

    def get_sample_project_name(self):
        return "%s" % (self.sample_project_name)

    def get_info_to_display(self):
        s_project_info = []
        s_project_info.append(self.pk)
        s_project_info.append(self.sample_project_name)
        s_project_info.append(self.sample_project_manager)
        return s_project_info

    def get_full_info_to_display(self):
        s_project_info = []
        s_project_info.append(self.sample_project_name)
        s_project_info.append(self.sample_project_manager)
        s_project_info.append(self.sampleProjectContact)
        s_project_info.append(self.sample_project_description)
        return s_project_info

    objects = SampleProjectsManager()


class SampleProjectFieldClassificationManager(models.Manager):
    def create_sample_project_field_classification(self, data):
        if "classification_display" not in data:
            data["classification_display"] = data["classification_name"]
        new_s_p_classification = self.create(
            sample_projects_id=data["sample_project_id"],
            classification_name=data["classification_name"],
            classification_display=data["classification_display"],
        )
        return new_s_p_classification


class SampleProjectFieldClassification(models.Model):
    sample_projects_id = models.ForeignKey(
        SampleProjects, on_delete=models.CASCADE, null=True, blank=True
    )
    classification_name = models.CharField(max_length=80)
    classification_display = models.CharField(max_length=100)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return "%s" % (self.classification_name)

    def get_classification_name(self):
        return "%s" % (self.classification_name)

    def get_classification_display(self):
        return "%s" % (self.classification_display)

    objects = SampleProjectFieldClassificationManager()


class SampleProjectsFieldsManager(models.Manager):
    def create_sample_project_fields(self, project_field_data):
        new_project_field = self.create(
            sample_projects_id=project_field_data["sample_project_id"],
            Sample_project_field_classification_id=project_field_data[
                "SampleProjectFieldClassificationID"
            ],
            sample_project_field_name=project_field_data["Field name"],
            sample_project_field_description=project_field_data["Description"],
            sample_project_field_order=project_field_data["Order"],
            sample_project_field_used=project_field_data["Used"],
            sample_project_field_type=project_field_data["Field type"],
            sample_project_searchable=project_field_data["Searchable"],
            sample_project_option_list=project_field_data["Option Values"],
        )
        return new_project_field


class SampleProjectsFields(models.Model):
    sample_projects_id = models.ForeignKey(
        SampleProjects, on_delete=models.CASCADE, null=True, blank=True
    )
    sample_project_field_classification_id = models.ForeignKey(
        SampleProjectFieldClassification,
        on_delete=models.CASCADE,
        null=True,
        blank=True,
    )
    sample_project_field_name = models.CharField(max_length=80)
    sample_project_field_description = models.CharField(
        max_length=400, null=True, blank=True
    )
    sample_project_field_order = models.IntegerField()
    sample_project_field_used = models.BooleanField()
    sample_project_field_type = models.CharField(max_length=20)
    sample_project_option_list = models.CharField(max_length=255, null=True, blank=True)
    sample_project_searchable = models.BooleanField(default=False)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return "%s" % (self.sample_project_field_name)

    def get_field_id(self):
        return "%s" % (self.id)

    def get_field_options_list(self):
        if self.sample_project_field_type == "Options List":
            if SamplesProjectsTableOptions.objects.filter(
                sample_project_field=self
            ).exists():
                s_p_opt_objs = SamplesProjectsTableOptions.objects.filter(
                    sample_project_field=self
                )
                data = []
                for s_p_opt_obj in s_p_opt_objs:
                    data.append(s_p_opt_obj.get_option_value())
                return data
            return ""
        else:
            return ""

    def get_field_name(self):
        return "%s" % (self.sample_project_field_name)

    def get_field_type(self):
        return "%s" % (self.sample_project_field_type)

    def get_description(self):
        return "%s" % (self.sample_project_field_description)

    def get_classification_name(self):
        if self.sample_project_field_classification_id is not None:
            return self.sample_project_field_classification_id.get_classification_display()

    def get_sample_project_fields_name(self):
        if self.sample_project_field_used:
            used = "Yes"
        else:
            used = "No"
        if self.sample_project_field_classification_id is not None:
            classification = (
                self.sample_project_field_classification_id.get_classification_display()
            )
        else:
            classification = ""
        field_data = []
        field_data.append(self.sample_project_field_name)

        field_data.append(self.sample_project_field_order)
        field_data.append(used)
        field_data.append(self.sample_project_searchable)
        field_data.append(self.sample_project_field_type)
        field_data.append(",".join(self.get_field_options_list()))
        field_data.append(self.sample_project_field_description)
        field_data.append(classification)

        return field_data

    def get_sample_project_fields_for_javascript(self):
        if self.sample_project_field_used:
            used = "true"
        else:
            used = "false"
        if self.sample_project_searchable:
            searchable = "true"
        else:
            searchable = "false"
        field_data = []
        field_data.append(self.sample_project_field_name)

        field_data.append(self.sample_project_field_order)
        field_data.append(used)
        field_data.append(searchable)
        field_data.append(self.sample_project_field_type)
        field_data.append(",".join(self.get_field_options_list()))
        field_data.append(self.sample_project_field_description)
        return field_data

    def update_sample_project_fields(self, project_field_data):
        self.sample_project_field_name = project_field_data["Field name"]
        self.sample_project_field_description = project_field_data["Description"]
        self.sample_project_field_order = project_field_data["Order"]
        self.sample_project_field_used = project_field_data["Used"]
        self.sample_project_field_type = project_field_data["Field type"]
        self.sample_project_searchable = project_field_data["Searchable"]
        self.save()
        return self

    objects = SampleProjectsFieldsManager()


class SamplesProjectsTableOptionsManager(models.Manager):
    def create_new_s_proj_table_opt(self, data):
        new_s_proj_table_opt = self.create(
            sample_project_field=data["s_proj_obj"], option_value=data["opt_value"]
        )
        return new_s_proj_table_opt


class SamplesProjectsTableOptions(models.Model):
    sample_project_field = models.ForeignKey(
        SampleProjectsFields,
        on_delete=models.CASCADE,
        null=True,
        blank=True,
        related_name="opt_value_prop",
    )
    option_value = models.CharField(max_length=120)

    def __str__(self):
        return "%s" % (self.option_value)

    def get_option_value(self):
        return "%s" % (self.option_value)

    def get_option_and_pk(self):
        return [self.pk, self.option_value]

    objects = SamplesProjectsTableOptionsManager()


class SamplesManager(models.Manager):
    def create_sample(self, sample_data):
        if sample_data["labRequest"] != "":
            sample_data["labRequest"] = LabRequest.objects.get(
                lab_name_coding__exact=sample_data["labRequest"]
            )
        else:
            sample_data["labRequest"] = None
        if sample_data["species"] != "":
            sample_data["species"] = Species.objects.get(
                species_name__exact=sample_data["species"]
            )
        else:
            sample_data["species"] = None
        if "sampleEntryDate" in sample_data:
            if not isinstance(sample_data["sampleEntryDate"], datetime.date):
                try:
                    sample_entry_date = datetime.datetime.strptime(
                        sample_data["sampleEntryDate"], "%Y-%m-%d %H:%M:%S"
                    )
                except ValueError:
                    sample_entry_date = None
            else:
                sample_entry_date = sample_data["sampleEntryDate"]
        else:
            sample_entry_date = None

        if "collectionSampleDate" in sample_data:
            if not isinstance(sample_data["collectionSampleDate"], datetime.date):
                try:
                    collection_sample_date = datetime.datetime.strptime(
                        sample_data["collectionSampleDate"], "%Y-%m-%d %H:%M:%S"
                    )
                except ValueError:
                    collection_sample_date = None
            else:
                collection_sample_date = sample_data["collectionSampleDate"]
        else:
            collection_sample_date = None

        if "completedDate" in sample_data:
            if not isinstance(sample_data["completedDate"], datetime.date):
                try:
                    completed_date = datetime.datetime.strptime(
                        sample_data["completedDate"], "%Y-%m-%d %H:%M:%S"
                    )
                except ValueError:
                    completed_date = None
            else:
                completed_date = sample_data["completedDate"]
        else:
            completed_date = None

        new_sample = self.create(
            sample_state=StatesForSample.objects.get(
                sample_state_name__exact=sample_data["sampleState"]
            ),
            patient_core=sample_data["patient"],
            lab_request=sample_data["labRequest"],
            sample_project=sample_data["sampleProject"],
            sample_type=SampleType.objects.get(
                sample_type__exact=sample_data["sampleType"],
                apps_name__exact=sample_data["app_name"],
            ),
            sample_user=User.objects.get(username__exact=sample_data["user"]),
            sample_code_id=sample_data["sample_id"],
            sample_name=sample_data["sampleName"],
            unique_sample_id=sample_data["new_unique_value"],
            species=sample_data["species"],
            sample_location=sample_data["sampleLocation"],
            only_recorded=sample_data["onlyRecorded"],
            sample_entry_date=sample_entry_date,
            collection_sample_date=collection_sample_date,
            completed_date=completed_date,
        )

        return new_sample


class Samples(models.Model):
    sample_state = models.ForeignKey(
        StatesForSample,
        on_delete=models.CASCADE,
        verbose_name="Sample state",
        null=True,
    )

    patient_core = models.ForeignKey(
        PatientCore, on_delete=models.CASCADE, null=True, blank=True
    )
    lab_request = models.ForeignKey(
        LabRequest,
        on_delete=models.CASCADE,
        verbose_name="Laboratory",
        null=True,
        blank=True,
    )

    sample_type = models.ForeignKey(
        SampleType, on_delete=models.CASCADE, verbose_name="Sample type", null=True
    )

    sample_user = models.ForeignKey(
        User, on_delete=models.CASCADE, verbose_name="Username", null=True
    )

    species = models.ForeignKey(
        Species, on_delete=models.CASCADE, null=True, verbose_name="Species", blank=True
    )

    sample_project = models.ForeignKey(
        SampleProjects,
        on_delete=models.CASCADE,
        verbose_name="Sample Project",
        null=True,
        blank=True,
    )

    sample_name = models.CharField(max_length=255, null=True, verbose_name="Sample Name")
    sample_location = models.CharField(
        max_length=255, null=True, verbose_name="Sample location", blank=True
    )
    sample_entry_date = models.DateTimeField(
        auto_now_add=False, null=True, verbose_name="Sample defined date", blank=True
    )
    collection_sample_date = models.DateTimeField(
        auto_now_add=False, null=True, verbose_name="Sample collection date", blank=True
    )
    unique_sample_id = models.CharField(
        max_length=8, verbose_name="Unique sample id", null=True
    )
    sample_code_id = models.CharField(
        max_length=60, null=True, verbose_name="Sample code id"
    )
    reused_number = models.IntegerField(
        default=0, verbose_name="Number of type reused"
    )
    sequencing_date = models.DateTimeField(
        auto_now_add=False, null=True, blank=True, verbose_name="Sequencing date"
    )
    completed_date = models.DateTimeField(
        auto_now_add=False, null=True, blank=True, verbose_name="Completion date"
    )
    generated_at = models.DateTimeField(auto_now_add=True, verbose_name="Generated at")
    only_recorded = models.BooleanField(
        default=False, null=True, blank=True, verbose_name="Only recorded?"
    )

    def __str__(self):
        return "%s" % (self.sample_name)

    def get_sample_definition_information(self):
        sample_info = []
        if self.sample_entry_date is not None:
            recorded_date = self.sample_entry_date.strftime("%d , %B , %Y")
        else:
            recorded_date = ""
        sample_info.append(self.unique_sample_id)
        sample_info.append(self.sample_code_id)
        sample_info.append(self.sample_name)
        sample_info.append(recorded_date)
        sample_info.append(self.sample_type.get_name())
        return sample_info

    def get_info_in_defined_state(self):
        sample_info = []
        if self.sample_entry_date is not None:
            sample_entry_date = self.sample_entry_date.strftime("%d , %B , %Y")
        else:
            sample_entry_date = ""
        sample_info.append(sample_entry_date)
        sample_info.append(self.sample_code_id)
        sample_info.append(self.sample_name)
        sample_info.append(str(self.pk))
        return sample_info

    def get_info_for_searching(self):
        sample_info = []
        if self.sample_entry_date is not None:
            sample_entry_date = self.sample_entry_date.strftime("%d , %B , %Y")
        else:
            sample_entry_date = ""
        sample_info.append(str(self.pk))
        sample_info.append(self.sample_name)
        sample_info.append(self.sample_state.get_sample_state())
        sample_info.append(sample_entry_date)
        sample_info.append(self.sample_code_id)
        sample_info.append(self.sample_type.get_name())
        try:
            sample_info.append(self.species.get_name())
        except KeyError:
            sample_info.append("Not defined")
        return sample_info

    def get_info_for_display(self):
        if self.collection_sample_date:
            collection_sample_date = self.collection_sample_date.strftime("%d , %B , %Y")
        else:
            collection_sample_date = ""
        if self.sample_entry_date:
            sample_entry_date = self.sample_entry_date.strftime("%d , %B , %Y")
        else:
            sample_entry_date = ""
        sample_info = []
        sample_info.append(self.sample_name)
        sample_info.append(self.sample_code_id)
        sample_info.append(self.sample_state.get_sample_state())
        sample_info.append(self.generated_at.strftime("%d , %B , %Y"))
        sample_info.append(collection_sample_date)
        sample_info.append(sample_entry_date)
        sample_info.append(self.sample_type.get_name())
        sample_info.append(self.species.get_name())
        sample_info.append(self.reused_number)
        sample_info.append(self.sample_user.username)
        return sample_info

    def get_info_for_patient(self):
        sample_info = []
        if self.sample_entry_date:
            sample_entry_date = self.sample_entry_date.strftime("%d , %B , %Y")
        else:
            sample_entry_date = ""
        sample_info.append(str(self.pk))
        sample_info.append(self.sample_name)
        sample_info.append(self.lab_request.get_name())
        sample_info.append(sample_entry_date)
        sample_info.append(self.sample_type.get_name())
        sample_info.append(self.sample_state.get_sample_state())
        return sample_info

    def get_entry_date(self):
        if self.sample_entry_date:
            sample_entry_date = self.sample_entry_date.strftime("%d , %B , %Y")
        else:
            sample_entry_date = ""
        return "%s" % sample_entry_date

    def get_lab_request(self):
        return "%s" % (self.lab_request.get_name())

    def get_sample_code(self):
        return "%s" % (self.sample_code_id)

    def get_sample_id(self):
        return "%s" % (self.pk)

    def get_sample_name(self):
        return "%s" % (self.sample_name)

    def get_sample_patient_code(self):
        return "%s" % (self.patient_core.get_patient_code())

    def get_sample_patient_name(self):
        return "%s" % (self.patient_core.get_patient_name())

    def get_sample_patient_surname(self):
        return "%s" % (self.patient_core.get_patient_surname())

    def get_sample_patient_obj(self):
        return self.patient_core

    def get_sample_project(self):
        if self.sample_project is None:
            return "None"
        else:
            return "%s" % (self.sample_project.get_sample_project_name())

    def get_region(self):
        if self.lab_request is None:
            return ""
        else:
            return "%s" % (self.lab_request.get_region())

    def get_sample_state(self):
        return "%s" % (self.sample_state)

    def get_sample_type(self):
        return "%s" % (self.sample_type.get_name())

    def get_species(self):
        return "%s" % (self.species.get_name())

    def get_register_user(self):
        if self.sample_user is None:
            return "Not available"
        else:
            return "%s" % (self.sample_user)

    def get_registered_sample(self):
        recorded_date = self.generated_at.strftime("%B %d, %Y")
        return "%s" % recorded_date

    def get_sample_project_obj(self):
        return self.sample_project

    def is_only_recorded(self):
        return self.only_recorded

    def get_unique_sample_id(self):
        return "%s" % (self.unique_sample_id)

    def set_state(self, state_value):
        self.sample_state = StatesForSample.objects.get(
            sample_state_name__exact=state_value
        )
        if state_value == "Sequencing":
            self.sequencing_date = datetime.datetime.today()
        self.save()

    def set_increase_reuse(self):
        self.reused_number += 1
        self.save()

    objects = SamplesManager()


class SampleProjectsFieldsValueManager(models.Manager):
    def create_project_field_value(self, field_value):
        new_field_data = self.create(
            sample_id=field_value["sample_id"],
            sample_project_field_id=field_value["sampleProjecttField_id"],
            sample_project_field_value=field_value["sampleProjectFieldValue"],
        )
        return new_field_data


class SampleProjectsFieldsValue(models.Model):
    sample_id = models.ForeignKey(
        Samples,
        related_name="project_values",
        on_delete=models.CASCADE,
        null=True,
        blank=True,
        verbose_name="Sample Name",
    )
    sample_project_field_id = models.ForeignKey(
        SampleProjectsFields, on_delete=models.CASCADE, null=True
    )
    sample_project_field_value = models.CharField(max_length=255, null=True, blank=True)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return "%s" % (self.sample_project_field_value)

    def get_field_value(self):
        return "%s" % (self.sample_project_field_value)

    objects = SampleProjectsFieldsValueManager()


class MoleculeUsedForManager(models.Manager):
    def create_molecule_use_for(self, molecule_use_data):
        new_molecule_use = self.create(
            used_for=molecule_use_data["usedFor"],
            apps_name=molecule_use_data["apps_name"],
            massive_use=molecule_use_data["massiveUse"],
        )
        return new_molecule_use


class MoleculeUsedFor(models.Model):
    used_for = models.CharField(max_length=50)
    apps_name = models.CharField(max_length=50)
    massive_use = models.BooleanField(default=False)

    def __str__(self):
        return "%s" % (self.usedFor)

    def get_molecule_use_name(self):
        return "%s" % (self.usedFor)

    def get_massive(self):
        return "%s" % (self.massive_use)

    objects = MoleculeUsedForManager()


class MoleculePreparationManager(models.Manager):
    def create_molecule(self, molecule_data):
        molecule_used_obj = MoleculeType.objects.filter(
            molecule_type__exact=molecule_data["moleculeType"]
        ).last()

        protocol_type_obj = ProtocolType.objects.filter(
            molecule=molecule_used_obj, apps_name__exact=molecule_data["app_name"]
        ).last()
        protocol_used_obj = Protocols.objects.filter(
            name__exact=molecule_data["protocolUsed"], type__exact=protocol_type_obj
        ).last()
        new_molecule = self.create(
            protocol_used=protocol_used_obj,
            sample=molecule_data["sample"],
            molecule_type=molecule_used_obj,
            state=StatesForMolecule.objects.get(molecule_state_name__exact="Defined"),
            molecule_code_id=molecule_data["moleculeCodeId"],
            molecule_extraction_date=molecule_data["moleculeExtractionDate"],
            extraction_type=molecule_data["extractionType"],
            molecule_user=User.objects.get(username__exact=molecule_data["user"]),
        )

        return new_molecule


class MoleculePreparation(models.Model):
    protocol_used = models.ForeignKey(Protocols, on_delete=models.CASCADE)
    sample = models.ForeignKey(Samples, on_delete=models.CASCADE)
    molecule_type = models.ForeignKey(MoleculeType, on_delete=models.CASCADE)
    state = models.ForeignKey(StatesForMolecule, on_delete=models.CASCADE, null=True)
    molecule_user = models.ForeignKey(
        User, on_delete=models.CASCADE, null=True, blank=True
    )

    user_lot_kit_id = models.ForeignKey(
        UserLotCommercialKits, on_delete=models.CASCADE, null=True, blank=True
    )

    molecule_used_for = models.ForeignKey(
        MoleculeUsedFor, on_delete=models.CASCADE, null=True, blank=True
    )

    molecule_code_id = models.CharField(max_length=255)
    extraction_type = models.CharField(max_length=50)
    molecule_extraction_date = models.DateTimeField(auto_now_add=False, null=True)
    reused_number = models.IntegerField(default=0)
    used_for_massive_sequencing = models.BooleanField(null=True, blank=True)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return "%s" % (self.molecule_code_id)

    def get_info_for_display(self):
        extraction_date = self.molecule_extraction_date.strftime("%d, %B, %Y")
        molecule_info = []
        molecule_info.append(self.molecule_code_id)
        molecule_info.append(self.state.get_molecule_state())
        molecule_info.append(extraction_date)
        molecule_info.append(self.extraction_type)
        molecule_info.append(self.molecule_type.get_name())
        if self.molecule_used_for is None:
            molecule_info.append("Not defined yet")
        else:
            molecule_info.append(self.molecule_used_for.get_molecule_use_name())
        molecule_info.append(self.protocol_used.get_name())
        molecule_info.append(self.reused_number)
        return molecule_info

    def get_extraction_date(self):
        return "%s" % (self.molecule_extraction_date.strftime("%B %d, %Y"))

    def get_molecule_id(self):
        return "%s" % (self.pk)

    def get_molecule_code_id(self):
        return "%s" % (self.molecule_code_id)

    def get_molecule_information(self):
        data = []
        data.append(self.molecule_code_id)
        data.append(self.molecule_extraction_date.strftime("%B %d, %Y"))
        data.append(self.protocol_used.get_name())
        data.append(str(self.pk))
        return data

    def get_sample_name(self):
        return "%s" % (self.sample.get_sample_name())

    def get_sample_obj(self):
        return self.sample

    def get_protocol(self):
        return "%s" % (self.protocol_used.get_name())

    def get_protocol_obj(self):
        return self.protocol_used

    def get_state(self):
        return "%s" % (self.state)

    def get_used_for_massive(self):
        return self.used_for_massive_sequencing

    def get_user_lot_kit_obj(self):
        return self.user_lot_kit_id

    def set_molecule_use(self, use_for_molecule, app_name):
        self.molecule_used_for = MoleculeUsedFor.objects.get(
            used_for__exact=use_for_molecule, apps_name__exact=app_name
        )
        self.save()
        self.used_for_massive_sequencing = self.molecule_used_for.get_massive()
        self.save()
        return self

    def set_state(self, state_value):
        self.state = StatesForMolecule.objects.get(molecule_state_name__exact=state_value)
        self.save()

    def set_increase_reuse(self):
        self.reused_number += 1
        self.save()

    def set_user_lot_kit(self, lot_kit_name):
        self.user_lot_kit_id = UserLotCommercialKits.objects.get(
            chip_lot__exact=lot_kit_name
        )
        self.save()

    def set_user_lot_kit_obj(self, lot_kit_obj):
        self.user_lot_kit_id = lot_kit_obj
        self.save()

    objects = MoleculePreparationManager()


class MoleculeParameterValueManager(models.Manager):
    def create_molecule_parameter_value(self, parameter_value):
        new_molecule_parameter_data = self.create(
            molecule_parameter_id=parameter_value["moleculeParameter_id"],
            molecule_id=parameter_value["molecule_id"],
            parameter_value=parameter_value["parameterValue"],
        )
        return new_molecule_parameter_data


class MoleculeParameterValue(models.Model):
    molecule_parameter_id = models.ForeignKey(
        ProtocolParameters, on_delete=models.CASCADE
    )
    molecule_id = models.ForeignKey(MoleculePreparation, on_delete=models.CASCADE)
    parameter_value = models.CharField(max_length=255)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return "%s" % (self.parameter_value)

    def get_param_value(self):
        return "%s" % (self.parameter_value)

    objects = MoleculeParameterValueManager()


class SequencingConfigurationManager(models.Manager):
    def create_new_configuration(self, data):
        try:
            platform_obj = SequencingPlatform.objects.get(pk__exact=data["platformID"])
        except SequencingPlatform.DoesNotExist:
            platform_obj = None
        new_sequencer_configuration = self.create(
            platform_id=platform_obj, configuration_name=data["configurationName"]
        )
        return new_sequencer_configuration


class SequencingConfiguration(models.Model):
    platform_id = models.ForeignKey(
        SequencingPlatform, on_delete=models.CASCADE, null=True, blank=True
    )
    configuration_name = models.CharField(max_length=255)

    def __str__(self):
        return "%s" % (self.configuration_name)

    def get_configuration_name(self):
        return "%s" % (self.configuration_name)

    def get_platform_name(self):
        return "%s" % (self.platform_id.get_platform_name())

    def get_platform_obj(self):
        return self.platform_id

    objects = SequencingConfigurationManager()


class SequencerInLabManager(models.Manager):
    def create_sequencer_in_lab(self, sequencer_value):
        try:
            platform_obj = SequencingPlatform.objects.get(
                pk__exact=sequencer_value["platformID"]
            )
        except models.SequencingPlatform.DoesNotExist:
            platform_obj = None
        new_sequencer = self.create(
            platform_id=platform_obj,
            sequencer_name=sequencer_value["sequencerName"],
            sequencer_description=sequencer_value["sequencerDescription"],
            sequencer_location=sequencer_value["sequencerLocation"],
            sequencer_serial_number=sequencer_value["sequencerSerialNumber"],
            sequencer_state="In Use",
            sequencer_operation_start=sequencer_value["sequencerOperationStart"],
            sequencer_number_lanes=sequencer_value["sequencerNumberLanes"],
        )

        return new_sequencer


class SequencerInLab(models.Model):
    platform_id = models.ForeignKey(
        SequencingPlatform, on_delete=models.CASCADE, null=True, blank=True
    )
    sequencer_name = models.CharField(max_length=255)
    sequencer_description = models.CharField(max_length=255, null=True, blank=True)
    sequencer_location = models.CharField(max_length=255, null=True, blank=True)
    sequencer_serial_number = models.CharField(max_length=255, null=True, blank=True)
    sequencer_state = models.CharField(max_length=50, null=True, blank=True)
    sequencer_operation_start = models.DateField(
        auto_now_add=False, null=True, blank=True
    )
    sequencer_operation_end = models.DateField(auto_now_add=False, null=True, blank=True)
    sequencer_number_lanes = models.CharField(max_length=5, null=True, blank=True)

    def __str__(self):
        return "%s" % (self.sequencer_name)

    def get_sequencer_name(self):
        return "%s" % (self.sequencer_name)

    def get_number_of_lanes(self):
        return "%s" % (self.sequencer_number_lanes)

    def get_sequencer_id(self):
        return "%s" % (self.pk)

    def get_sequencing_platform_name(self):
        if self.platform_id is None:
            return "Not Defined"
        else:
            return "%s" % (self.platform_id.get_platform_name())

    def get_sequencing_platform_id(self):
        if self.platform_id is None:
            return ""
        else:
            return "%s" % (self.platform_id.get_platform_id())

    def get_all_sequencer_data(self):
        data = []
        if self.platform_id is not None:
            platform_name = self.platformID.get_platform_name()
        else:
            platform_name = "Not Defined"
        if self.sequencer_operation_start is not None:
            op_start = self.sequencer_operation_start.strftime("%B %d, %Y")
        else:
            op_start = "Not available date"
        if self.sequencer_operation_end is not None:
            op_end = self.sequencer_operation_end.strftime("%B %d, %Y")
        else:
            op_end = "Not available date"
        data.append(platform_name)
        data.append(self.sequencer_name)
        data.append(self.sequencer_description)
        data.append(self.sequencer_location)
        data.append(self.sequencer_serial_number)
        data.append(self.sequencer_state)
        data.append(op_start)
        data.append(op_end)
        data.append(self.sequencer_number_lanes)
        return data

    objects = SequencerInLabManager()


class OntologyMap(models.Model):
    label = models.CharField(max_length=255)
    ontology = models.CharField(max_length=50, null=True, blank=True)

    def __str__(self):
        return "%s" % (self.label)

    def get_label(self):
        return "%s" % (self.label)

    def get_ontology(self):
        return "%s" % (self.ontology)

    def get_label_and_ontology(self):
        return {self.label: self.ontology}
