from django.db import models
from django.contrib.auth.models import User
import datetime


class StateInCountryManager(models.Manager):
    def create_new_state(self, data):
        new_state = self.create(stateName=data["state"], apps_name=data["apps_name"])
        return new_state


class StateInCountry(models.Model):
    stateName = models.CharField(max_length=80)
    apps_name = models.CharField(max_length=40, null=True)

    def __str__(self):
        return "%s" % (self.stateName)

    def get_state_name(self):
        return "%s" % (self.stateName)

    def get_state_id(self):
        return "%s" % (self.pk)

    objects = StateInCountryManager()


class CityManager(models.Manager):
    def create_new_city(self, data):
        if StateInCountry.objects.filter(pk__exact=data["state"]).exists():
            state_obj = StateInCountry.objects.filter(pk__exact=data["state"]).last()
        else:
            state_obj = None
        new_city = self.create(
            belongsToState=state_obj,
            cityName=data["cityName"],
            geoLocLatitude=data["latitude"],
            geoLocLongitude=data["longitude"],
            apps_name=data["apps_name"],
        )
        return new_city


class City(models.Model):
    belongsToState = models.ForeignKey(
        StateInCountry, on_delete=models.CASCADE, null=True, blank=True
    )
    cityName = models.CharField(max_length=80)
    geoLocLatitude = models.CharField(max_length=80)
    geoLocLongitude = models.CharField(max_length=80)
    apps_name = models.CharField(max_length=40, null=True)

    def __str__(self):
        return "%s" % (self.cityName)

    def get_city_name(self):
        return "%s" % (self.cityName)

    def get_city_id(self):
        return "%s" % (self.pk)

    def get_coordenates(self):
        return {"latitude": self.geoLocLatitude, "longitude": self.geoLocLongitude}

    def get_state(self):
        if self.belongsToState is None:
            return ""
        return "%s" % (self.belongsToState.get_state_name())

    objects = CityManager()


class LabRequestManager(models.Manager):
    def create_lab_request(self, data):
        city_obj = City.objects.filter(pk__exact=data["city"]).last()
        new_lab_request = self.create(
            labName=data["labName"],
            labNameCoding=data["labNameCoding"],
            labUnit=data["labUnit"],
            labContactName=data["labContactName"],
            labPhone=data["labPhone"],
            labEmail=data["labEmail"],
            address=data["address"],
            apps_name=data["apps_name"],
            labCity=city_obj,
        )
        return new_lab_request


class LabRequest(models.Model):
    labCity = models.ForeignKey(City, on_delete=models.CASCADE, null=True, blank=True)
    labName = models.CharField(max_length=80)
    labNameCoding = models.CharField(max_length=50)
    labUnit = models.CharField(max_length=50)
    labContactName = models.CharField(max_length=50)
    labPhone = models.CharField(max_length=20)
    labEmail = models.CharField(max_length=70)
    address = models.CharField(max_length=255)

    apps_name = models.CharField(max_length=40, null=True)

    def __str__(self):
        return "%s" % (self.labName)

    def get_name(self):
        return "%s" % (self.labName)

    def get_id(self):
        return "%s" % (self.pk)

    def get_lab_request_code(self):
        return "%s" % (self.labNameCoding)

    def get_all_data(self):
        data = []
        data.append(self.labName)
        data.append(self.labNameCoding)
        data.append(self.labUnit)
        data.append(self.labContactName)
        data.append(self.labPhone)
        data.append(self.labEmail)
        data.append(self.address)
        return data

    def get_fields_and_data(self):
        data = {}
        if self.labCity is None:
            data["state"] = ""
            data["city"] = ""
            data["latitude"] = ""
            data["longitude"] = ""
        else:
            data["state"] = self.labCity.get_state()
            data["city"] = self.labCity.get_city_name()
            data.update(self.labCity.get_coordenates())
        data["Lab Name"] = self.labName
        data["Lab Coding"] = self.labNameCoding
        return data

    def get_region(self):
        if self.labCity is None:
            return ""
        else:
            return "%s" % (self.labCity.get_state())

    objects = LabRequestManager()


class MoleculeTypeManager(models.Manager):
    def create_molecule_type(self, data):
        new_molecule_type = self.create(
            moleculeType=data["moleculeType"], apps_name=data["apps_name"]
        )
        return new_molecule_type


class MoleculeType(models.Model):
    moleculeType = models.CharField(max_length=30)
    apps_name = models.CharField(max_length=40, null=True)

    def __str__(self):
        return "%s" % (self.moleculeType)

    def get_name(self):
        return "%s" % (self.moleculeType)

    def get_id(self):
        return "%s" % (self.pk)

    objects = MoleculeTypeManager()


class ProtocolTypeManager(models.Manager):
    def create_protocol_type(self, data):
        if data["molecule"] is not None:
            molecule_obj = MoleculeType.objects.filter(
                moleculeType__iexact=data["molecule"],
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
            parameterName=prot_param_data["Parameter name"],
            parameterDescription=prot_param_data["Description"],
            parameterOrder=prot_param_data["Order"],
            parameterUsed=prot_param_data["Used"],
            parameterMaxValue=prot_param_data["Max Value"],
            parameterMinValue=prot_param_data["Min Value"],
            parameterOptionValues=prot_param_data["Option Values"],
            parameterType=prot_param_data["Parameter Type"],
        )
        return new_prot_parameter


class ProtocolParameters(models.Model):
    protocol_id = models.ForeignKey(Protocols, on_delete=models.CASCADE)
    parameterName = models.CharField(max_length=255)
    parameterDescription = models.CharField(max_length=400, null=True, blank=True)
    parameterOrder = models.IntegerField()
    parameterUsed = models.BooleanField()
    parameterType = models.CharField(max_length=20, default="string")
    parameterOptionValues = models.CharField(max_length=400, null=True, blank=True)
    parameterMaxValue = models.CharField(max_length=50, null=True, blank=True)
    parameterMinValue = models.CharField(max_length=50, null=True, blank=True)

    def __str__(self):
        return "%s" % (self.parameterName)

    def get_parameter_name(self):
        return "%s" % (self.parameterName)

    def get_parameter_protocol_id(self):
        return "%s" % (self.pk)

    def get_parameter_option_values(self):
        return "%s" % (self.parameterOptionValues)

    def get_parameter_type(self):
        return "%s" % (self.parameterType)

    def get_all_parameter_info(self):
        param_info = []
        param_info.append(self.parameterName)
        param_info.append(self.parameterOrder)
        param_info.append(self.parameterUsed)
        param_info.append(self.parameterType)
        param_info.append(self.parameterOptionValues)
        param_info.append(self.parameterMinValue)
        param_info.append(self.parameterMaxValue)
        param_info.append(self.parameterDescription)
        return param_info

    def get_protocol_fields_for_javascript(self):
        if self.parameterUsed:
            used = "true"
        else:
            used = "false"
        if self.parameterOptionValues is None:
            parameterOptionValues = ""
        else:
            parameterOptionValues = self.parameterOptionValues
        field_data = []
        field_data.append(self.parameterName)

        field_data.append(self.parameterOrder)
        field_data.append(used)
        field_data.append(self.parameterType)
        field_data.append(parameterOptionValues)
        field_data.append(self.parameterDescription)
        return field_data

    def update_protocol_fields(self, prot_param_data):
        self.parameterName = prot_param_data["Parameter name"]
        self.parameterDescription = prot_param_data["Description"]
        self.parameterOrder = prot_param_data["Order"]
        self.parameterUsed = prot_param_data["Used"]
        self.parameterOptionValues = prot_param_data["Option Values"]
        self.parameterType = prot_param_data["Parameter Type"]
        self.save()

    objects = ProtocolParametersManager()


class StatesForSample(models.Model):
    sampleStateName = models.CharField(max_length=50)

    def __str__(self):
        return "%s" % (self.sampleStateName)

    def get_sample_state(self):
        return "%s" % (self.sampleStateName)

    def get_id(self):
        return "%s" % (self.pk)


class StatesForMolecule(models.Model):
    moleculeStateName = models.CharField(max_length=50)

    def __str__(self):
        return "%s" % (self.moleculeStateName)

    def get_molecule_state(self):
        return "%s" % (self.moleculeStateName)


class SampleTypeManager(models.Manager):
    def create_sample_type(self, sample_type_data):
        new_sample_type = self.create(
            sampleType=sample_type_data["sampleType"],
            apps_name=sample_type_data["apps_name"],
            optional_fields=sample_type_data["optional_fields"],
        )
        return new_sample_type


class SampleType(models.Model):
    sampleType = models.CharField(max_length=50)
    apps_name = models.CharField(max_length=50)
    optional_fields = models.CharField(max_length=50, null=True, blank=True)
    generatedat = models.DateTimeField(auto_now_add=True, null=True)

    def __str__(self):
        return "%s" % (self.sampleType)

    def get_sample_type_id(self):
        return "%s" % (self.pk)

    def get_name(self):
        return "%s" % (self.sampleType)

    def get_optional_values(self):
        if self.optional_fields == "" or self.optional_fields is None:
            return []
        else:
            return list(map(int, self.optional_fields.split(",")))

    objects = SampleTypeManager()


class SpeciesManager(models.Manager):
    def create_new_specie(self, data):
        new_specie = self.create(speciesName=data["name"], apps_name=data["apps_name"])
        return new_specie


class Species(models.Model):
    speciesName = models.CharField(max_length=50)
    refGenomeName = models.CharField(max_length=255, null=True, blank=True)
    refGenomeSize = models.CharField(max_length=100, null=True, blank=True)
    refGenomeID = models.CharField(max_length=255, null=True, blank=True)
    apps_name = models.CharField(max_length=50, null=True)
    generatedat = models.DateTimeField(auto_now_add=True, null=True)

    def __str__(self):
        return "%s" % (self.speciesName)

    def get_name(self):
        return "%s" % (self.speciesName)

    def get_id(self):
        return "%s" % (self.pk)

    objects = SpeciesManager()


class SequencingPlatform(models.Model):
    platformName = models.CharField(max_length=30)
    companyName = models.CharField(max_length=30)
    sequencingTecnology = models.CharField(max_length=30)

    def __str__(self):
        return "%s" % (self.platformName)

    def get_platform_name(self):
        return "%s" % (self.platformName)

    def get_platform_id(self):
        return "%s" % (self.pk)

    def get_company_name(self):
        return "%s" % (self.companyName)


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
            platformKits=platform_obj,
        )
        return new_commercial_kit


class CommercialKits(models.Model):

    protocolKits = models.ManyToManyField(Protocols, blank=True)

    platformKits = models.ForeignKey(
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
    generatedat = models.DateTimeField(auto_now_add=True, null=True)

    def __str__(self):
        return "%s" % (self.name)

    def get_name(self):
        return "%s" % (self.name)

    def platform_kit_obj(self):
        return self.platformKits

    def get_platform_name(self):
        if self.platformKits is not None:
            return "%s" % (self.platformKits.get_platform_name())
        else:
            return ""

    def get_protocol_objs(self):
        return self.protocolKits.all()

    def get_provider_kit_name(self):
        return "%s" % (self.provider)

    def get_cat_number(self):
        return "%s" % (self.cat_number)

    def get_commercial_platform_basic_data(self):
        kit_basic_data = []
        kit_basic_data.append(self.name)
        kit_basic_data.append(self.provider)
        kit_basic_data.append(self.platformKits.get_platform_name())
        return kit_basic_data

    def get_commercial_protocol_basic_data(self):
        kit_basic_data = []
        kit_basic_data.append(self.name)
        kit_basic_data.append(self.provider)
        protocols = []
        protocol_objs = self.protocolKits.all()
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
            basedCommercial=kit_data["basedCommercial"],
            chipLot=kit_data["chipLot"],
            expirationDate=expiration_date,
        )
        return new_user_lot_commercial_kit


class UserLotCommercialKits(models.Model):
    user = models.ForeignKey(User, on_delete=models.CASCADE, null=True)
    basedCommercial = models.ForeignKey(
        CommercialKits, on_delete=models.CASCADE, null=True
    )

    numberOfuses = models.IntegerField(null=True, default=0)
    chipLot = models.CharField(max_length=50)
    latestUsedDate = models.DateTimeField(null=True, blank=True)
    expirationDate = models.DateField(auto_now_add=False)
    runOut = models.BooleanField(default=False)
    generatedat = models.DateTimeField(auto_now_add=True, null=True)

    def __str__(self):
        return "%s" % (self.chipLot)

    def get_basic_data(self):
        lot_data = []
        lot_data.append(self.basedCommercial.get_name())
        lot_data.append(self.chipLot)
        lot_data.append(self.expirationDate.strftime("%d %B %Y"))
        return lot_data

    def get_commercial_kit(self):
        return "%s" % (self.basedCommercial.get_name())

    def get_commercial_obj(self):
        return self.basedCommercial

    def get_lot_number(self):
        return "%s" % (self.chipLot)

    def get_number_of_uses(self):
        return "%s" % (self.numberOfuses)

    def get_protocol_for_kit(self):
        return "%s" % (self.basedCommercial.get_protocol())

    def get_protocol_obj_for_kit(self):
        return self.basedCommercial.get_protocol_obj()

    def get_expiration_date(self):
        exp_date = self.expirationDate
        return "%s" % (exp_date.strftime("%d %B, %Y"))

    def get_user_lot_kit_id(self):
        return "%s" % (self.pk)

    def set_increase_use(self):
        self.latestUsedDate = datetime.datetime.now()
        self.numberOfuses += 1
        self.save()
        return self

    def set_latest_use(self, date):
        self.latestUsedDate = date
        self.save()
        return self

    def set_run_out(self):
        self.runOut = True
        self.save()
        return self

    objects = UserLotCommercialKitsManager()


class PatientProjectsManager(models.Manager):
    def create_project(self, project_data):
        new_project = self.create(
            projectName=project_data["projectName"],
            projectManager=project_data["projectManager"],
            projectContact=project_data["projectContact"],
            projectDescription=project_data["projectDescription"],
            apps_name=project_data["apps_name"],
        )
        return new_project


class PatientProjects(models.Model):

    projectName = models.CharField(max_length=50)
    projectManager = models.CharField(max_length=50, null=True, blank=True)
    projectContact = models.CharField(max_length=50, null=True, blank=True)
    projectDescription = models.CharField(max_length=255, null=True, blank=True)
    apps_name = models.CharField(max_length=40)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return "%s" % (self.projectName)

    def get_project_id(self):
        return "%s" % (self.pk)

    def get_project_name(self):
        return "%s" % (self.projectName)

    def get_patient_project_data(self):
        p_data = []
        p_data.append(self.projectName)
        p_data.append(self.projectManager)
        p_data.append(self.projectContact)
        p_data.append(self.projectDescription)
        p_data.append(self.pk)
        return p_data

    objects = PatientProjectsManager()


class PatientProjectsFieldsManager(models.Manager):
    def create_project_fields(self, project_field_data):
        new_project_field = self.create(
            patientProjects_id=project_field_data["project_id"],
            projectFieldName=project_field_data["Field name"],
            projectFieldDescription=project_field_data["Description"],
            projectFieldOrder=project_field_data["Order"],
            projectFieldUsed=project_field_data["Used"],
        )
        return new_project_field


class PatientProjectsFields(models.Model):
    patientProjects_id = models.ForeignKey(
        PatientProjects, on_delete=models.CASCADE, null=True, blank=True
    )
    projectFieldName = models.CharField(max_length=50)
    projectFieldDescription = models.CharField(max_length=400, null=True, blank=True)
    projectFieldOrder = models.IntegerField()
    projectFieldUsed = models.BooleanField()

    def __str__(self):
        return "%s" % (self.projectFieldName)

    def get_field_id(self):
        return "%s" % (self.id)

    def get_field_name(self):
        return "%s" % (self.projectFieldName)

    def get_description(self):
        return "%s" % (self.projectFieldDescription)

    def get_all_fields_info(self):
        if self.projectFieldUsed:
            used = "Yes"
        else:
            used = "No"
        field_data = []
        field_data.append(self.projectFieldName)

        field_data.append(self.projectFieldOrder)
        field_data.append(used)
        field_data.append(self.projectFieldDescription)
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
            patientName=p_data["patientName"],
            patientSurname=p_data["patientSurname"],
            patientCode=p_data["patientCode"],
            patientSex=sex,
        )
        return new_patient


class PatientCore(models.Model):
    patientProjects = models.ManyToManyField(PatientProjects, blank=True)
    patientName = models.CharField(max_length=255, null=True)
    patientSurname = models.CharField(max_length=255, null=True)
    patientCode = models.CharField(max_length=255, null=True)
    patientSex = models.ForeignKey(
        PatientSex, on_delete=models.CASCADE, null=True, blank=True
    )

    def __str__(self):
        return "%s" % (self.patientCode)

    def get_patient_id(self):
        return "%s" % (self.id)

    def get_patient_name(self):
        return "%s" % (self.patientName)

    def get_patient_surname(self):
        return "%s" % (self.patientSurname)

    def get_patient_code(self):
        return "%s" % (self.patientCode)

    def get_patient_sex(self):
        return "%s" % (self.patientSex)

    objects = PatientCoreManager()


class PatientProjectFieldValueManager(models.Manager):
    def create_project_field_value(self, field_value):
        new_field_data = self.create(
            projectField_id=field_value["projectField_id"],
            patientCore_id=field_value["patientCore_id"],
            projectFieldValue=field_value["projectFieldValue"],
        )
        return new_field_data


class PatientProjectFieldValue(models.Model):
    patientCore_id = models.ForeignKey(
        PatientCore, on_delete=models.CASCADE, null=True, blank=True
    )
    projectField_id = models.ForeignKey(
        PatientProjectsFields, on_delete=models.CASCADE, null=True
    )
    projectFieldValue = models.CharField(max_length=255)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return "%s" % (self.projectFieldValue)

    def get_field_value(self):
        return "%s" % (self.projectFieldValue)

    objects = PatientProjectFieldValueManager()


class SampleProjectsManager(models.Manager):
    def create_sample_project(self, s_project_data):
        new_sample_project = self.create(
            sampleProjectName=s_project_data["sampleProjectName"],
            sampleProjectManager=s_project_data["sampleProjectManager"],
            sampleProjectContact=s_project_data["sampleProjectContact"],
            sampleProjectDescription=s_project_data["sampleProjectDescription"],
            apps_name=s_project_data["apps_name"],
        )
        return new_sample_project


class SampleProjects(models.Model):
    sampleProjectName = models.CharField(max_length=255)
    sampleProjectManager = models.CharField(max_length=50, null=True, blank=True)
    sampleProjectContact = models.CharField(max_length=250, null=True, blank=True)
    sampleProjectDescription = models.CharField(max_length=255, null=True, blank=True)
    apps_name = models.CharField(max_length=255)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return "%s" % (self.sampleProjectName)

    def get_id(self):
        return "%s" % (self.pk)

    def get_sample_project_name(self):
        return "%s" % (self.sampleProjectName)

    def get_info_to_display(self):
        s_project_info = []
        s_project_info.append(self.pk)
        s_project_info.append(self.sampleProjectName)
        s_project_info.append(self.sampleProjectManager)
        return s_project_info

    def get_full_info_to_display(self):
        s_project_info = []
        s_project_info.append(self.sampleProjectName)
        s_project_info.append(self.sampleProjectManager)
        s_project_info.append(self.sampleProjectContact)
        s_project_info.append(self.sampleProjectDescription)
        return s_project_info

    objects = SampleProjectsManager()


class SampleProjectsFieldsManager(models.Manager):
    def create_sample_project_fields(self, project_field_data):
        new_project_field = self.create(
            sampleProjects_id=project_field_data["sample_project_id"],
            sampleProjectFieldName=project_field_data["Field name"],
            sampleProjectFieldDescription=project_field_data["Description"],
            sampleProjectFieldOrder=project_field_data["Order"],
            sampleProjectFieldUsed=project_field_data["Used"],
            sampleProjectFieldType=project_field_data["Field type"],
            sampleProjectSearchable=project_field_data["Searchable"],
            sampleProjectOptionList=project_field_data["Option Values"],
        )
        return new_project_field


class SampleProjectsFields(models.Model):
    sampleProjects_id = models.ForeignKey(
        SampleProjects, on_delete=models.CASCADE, null=True, blank=True
    )
    # sampleProjectsOptionValues = models.ManyToManyField(SamplesProjectsOptionValues, blank = True)

    sampleProjectFieldName = models.CharField(max_length=80)
    sampleProjectFieldDescription = models.CharField(
        max_length=400, null=True, blank=True
    )
    sampleProjectFieldOrder = models.IntegerField()
    sampleProjectFieldUsed = models.BooleanField()
    sampleProjectFieldType = models.CharField(max_length=20)
    sampleProjectOptionList = models.CharField(max_length=255, null=True, blank=True)
    sampleProjectSearchable = models.BooleanField(default=False)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return "%s" % (self.sampleProjectFieldName)

    def get_field_id(self):
        return "%s" % (self.id)

    def get_field_options_list(self):
        if self.sampleProjectFieldType == "Options List":
            if SamplesProjectsTableOptions.objects.filter(
                sampleProjectField=self
            ).exists():
                s_p_opt_objs = SamplesProjectsTableOptions.objects.filter(
                    sampleProjectField=self
                )
                data = []
                for s_p_opt_obj in s_p_opt_objs:
                    data.append(s_p_opt_obj.get_option_value())
                return data
            return ""
        else:
            return ""

    def get_field_name(self):
        return "%s" % (self.sampleProjectFieldName)

    def get_field_type(self):
        return "%s" % (self.sampleProjectFieldType)

    def get_description(self):
        return "%s" % (self.sampleProjectFieldDescription)

    def get_sample_project_fields_name(self):
        if self.sampleProjectFieldUsed:
            used = "Yes"
        else:
            used = "No"

        field_data = []
        field_data.append(self.sampleProjectFieldName)

        field_data.append(self.sampleProjectFieldOrder)
        field_data.append(used)
        field_data.append(self.sampleProjectSearchable)
        field_data.append(self.sampleProjectFieldType)
        field_data.append(",".join(self.get_field_options_list()))
        field_data.append(self.sampleProjectFieldDescription)
        return field_data

    def get_sample_project_fields_for_javascript(self):
        if self.sampleProjectFieldUsed:
            used = "true"
        else:
            used = "false"
        if self.sampleProjectSearchable:
            searchable = "true"
        else:
            searchable = "false"
        field_data = []
        field_data.append(self.sampleProjectFieldName)

        field_data.append(self.sampleProjectFieldOrder)
        field_data.append(used)
        field_data.append(searchable)
        field_data.append(self.sampleProjectFieldType)
        field_data.append(",".join(self.get_field_options_list()))
        field_data.append(self.sampleProjectFieldDescription)
        return field_data

    def update_sample_project_fields(self, project_field_data):
        # self.sampleProjects_id =project_field_data['sample_project_id']
        self.sampleProjectFieldName = project_field_data["Field name"]
        self.sampleProjectFieldDescription = project_field_data["Description"]
        self.sampleProjectFieldOrder = project_field_data["Order"]
        self.sampleProjectFieldUsed = project_field_data["Used"]
        self.sampleProjectFieldType = project_field_data["Field type"]
        self.sampleProjectSearchable = project_field_data["Searchable"]
        # self.sampleProjectOptionList = project_field_data["Option Values"]
        self.save()
        return self

    objects = SampleProjectsFieldsManager()


class SamplesProjectsTableOptionsManager(models.Manager):
    def create_new_s_proj_table_opt(self, data):
        new_s_proj_table_opt = self.create(
            sampleProjectField=data["s_proj_obj"], optionValue=data["opt_value"]
        )
        return new_s_proj_table_opt


class SamplesProjectsTableOptions(models.Model):
    sampleProjectField = models.ForeignKey(
        SampleProjectsFields,
        on_delete=models.CASCADE,
        null=True,
        blank=True,
        related_name="opt_value_prop",
    )
    optionValue = models.CharField(max_length=120)

    def __str__(self):
        return "%s" % (self.optionValue)

    def get_option_value(self):
        return "%s" % (self.optionValue)

    def get_option_and_pk(self):
        return [self.pk, self.optionValue]

    objects = SamplesProjectsTableOptionsManager()


class SamplesManager(models.Manager):
    def create_sample(self, sample_data):
        if sample_data["labRequest"] != "":
            sample_data["labRequest"] = LabRequest.objects.get(
                labNameCoding__exact=sample_data["labRequest"]
            )
        else:
            sample_data["labRequest"] = None
        if sample_data["species"] != "":
            sample_data["species"] = Species.objects.get(
                speciesName__exact=sample_data["species"]
            )
        else:
            sample_data["species"] = None
        if "completedDate" in sample_data:
            completedDate = sample_data["completedDate"]
        else:
            completedDate = None
        new_sample = self.create(
            sampleState=StatesForSample.objects.get(
                sampleStateName__exact=sample_data["sampleState"]
            ),
            patientCore=sample_data["patient"],
            labRequest=sample_data["labRequest"],
            sampleProject=sample_data["sampleProject"],
            sampleType=SampleType.objects.get(
                sampleType__exact=sample_data["sampleType"],
                apps_name__exact=sample_data["app_name"],
            ),
            sampleUser=User.objects.get(username__exact=sample_data["user"]),
            sampleCodeID=sample_data["sample_id"],
            sampleName=sample_data["sampleName"],
            uniqueSampleID=sample_data["new_unique_value"],
            species=sample_data["species"],
            sampleLocation=sample_data["sampleLocation"],
            onlyRecorded=sample_data["onlyRecorded"],
            sampleEntryDate=datetime.datetime.strptime(
                sample_data["sampleEntryDate"], "%Y-%m-%d %H:%M:%S"
            ),
            collectionSampleDate=datetime.datetime.strptime(
                sample_data["sampleEntryDate"], "%Y-%m-%d %H:%M:%S"
            ),
            completedDate=completedDate,
        )

        return new_sample


class Samples(models.Model):
    sampleState = models.ForeignKey(
        StatesForSample, on_delete=models.CASCADE, null=True
    )

    patientCore = models.ForeignKey(
        PatientCore, on_delete=models.CASCADE, null=True, blank=True
    )
    labRequest = models.ForeignKey(
        LabRequest, on_delete=models.CASCADE, null=True, blank=True
    )

    sampleType = models.ForeignKey(SampleType, on_delete=models.CASCADE, null=True)

    sampleUser = models.ForeignKey(User, on_delete=models.CASCADE, null=True)

    species = models.ForeignKey(
        Species, on_delete=models.CASCADE, null=True, blank=True
    )

    sampleProject = models.ForeignKey(
        SampleProjects, on_delete=models.CASCADE, null=True, blank=True
    )

    sampleName = models.CharField(max_length=255, null=True)
    sampleLocation = models.CharField(max_length=255, null=True, blank=True)
    sampleEntryDate = models.DateTimeField(auto_now_add=False, null=True)
    collectionSampleDate = models.DateTimeField(auto_now_add=False, null=True)
    uniqueSampleID = models.CharField(max_length=8, null=True)
    sampleCodeID = models.CharField(max_length=60, null=True)
    numberOfReused = models.IntegerField(default=0)
    sequencingDate = models.DateTimeField(auto_now_add=False, null=True, blank=True)
    completedDate = models.DateTimeField(auto_now_add=False, null=True, blank=True)
    generated_at = models.DateTimeField(auto_now_add=True)
    onlyRecorded = models.BooleanField(default=False, null=True, blank=True)

    def __str__(self):
        return "%s" % (self.sampleName)

    def get_sample_definition_information(self):
        sample_info = []
        recordeddate = self.sampleEntryDate.strftime("%d , %B , %Y")
        sample_info.append(self.uniqueSampleID)
        sample_info.append(self.sampleCodeID)
        sample_info.append(self.sampleName)
        sample_info.append(recordeddate)
        sample_info.append(self.sampleType.get_name())
        return sample_info

    def get_info_in_defined_state(self):
        sample_info = []
        sample_info.append(self.sampleEntryDate.strftime("%d , %B , %Y"))
        sample_info.append(self.sampleCodeID)
        sample_info.append(self.sampleName)
        sample_info.append(str(self.pk))
        return sample_info

    def get_info_for_searching(self):
        sample_info = []
        sample_info.append(str(self.pk))
        sample_info.append(self.sampleName)
        sample_info.append(self.sampleState.get_sample_state())
        sample_info.append(self.sampleEntryDate.strftime("%d , %B , %Y"))
        sample_info.append(self.sampleCodeID)
        sample_info.append(self.sampleType.get_name())
        try:
            sample_info.append(self.species.get_name())
        except KeyError:
            sample_info.append("Not defined")
        return sample_info

    def get_info_for_display(self):
        if self.collectionSampleDate:
            collectionSampleDate = self.collectionSampleDate.strftime("%d , %B , %Y")
        else:
            collectionSampleDate = ""
        if self.sampleEntryDate:
            sampleEntryDate = self.sampleEntryDate.strftime("%d , %B , %Y")
        else:
            sampleEntryDate = ""
        sample_info = []
        sample_info.append(self.sampleName)
        sample_info.append(self.sampleCodeID)
        sample_info.append(self.sampleState.get_sample_state())
        sample_info.append(self.generated_at.strftime("%d , %B , %Y"))
        sample_info.append(collectionSampleDate)
        sample_info.append(sampleEntryDate)
        sample_info.append(self.sampleType.get_name())
        sample_info.append(self.species.get_name())
        sample_info.append(self.numberOfReused)
        sample_info.append(self.sampleUser.username)
        return sample_info

    def get_info_for_patient(self):
        sample_info = []
        sample_info.append(str(self.pk))
        sample_info.append(self.sampleName)
        sample_info.append(self.labRequest.get_name())
        sample_info.append(self.sampleEntryDate.strftime("%d , %B , %Y"))
        sample_info.append(self.sampleType.get_name())
        sample_info.append(self.sampleState.get_sample_state())
        return sample_info

    def get_extraction_date(self):
        recordeddate = self.sampleEntryDate.strftime("%d , %B , %Y")
        return "%s" % (recordeddate)

    def get_lab_request(self):
        return "%s" % (self.labRequest.get_name())

    def get_sample_code(self):
        return "%s" % (self.sampleCodeID)

    def get_sample_id(self):
        return "%s" % (self.pk)

    def get_sample_name(self):
        return "%s" % (self.sampleName)

    def get_sample_patient_code(self):
        return "%s" % (self.patientCore.get_patient_code())

    def get_sample_patient_name(self):
        return "%s" % (self.patientCore.get_patient_name())

    def get_sample_patient_surname(self):
        return "%s" % (self.patientCore.get_patient_surname())

    def get_sample_patient_obj(self):
        return self.patientCore

    def get_sample_project(self):
        if self.sampleProject is None:
            return "None"
        else:
            return "%s" % (self.sampleProject.get_sample_project_name())

    def get_region(self):
        if self.labRequest is None:
            return ""
        else:
            return "%s" % (self.labRequest.get_region())

    def get_sample_state(self):
        return "%s" % (self.sampleState)

    def get_sample_type(self):
        return "%s" % (self.sampleType.get_name())

    def get_species(self):
        return "%s" % (self.species.get_name())

    def get_register_user(self):
        if self.sampleUser is None:
            return "Not available"
        else:
            return "%s" % (self.sampleUser)

    def get_registered_sample(self):
        recordeddate = self.generated_at.strftime("%B %d, %Y")
        return "%s" % (recordeddate)

    def get_sample_project_obj(self):
        return self.sampleProject

    def is_only_recorded(self):
        return self.onlyRecorded

    def get_unique_sample_id(self):
        return "%s" % (self.uniqueSampleID)

    def set_state(self, state_value):
        self.sampleState = StatesForSample.objects.get(
            sampleStateName__exact=state_value
        )
        if state_value == "Sequencing":
            self.sequencingDate = datetime.datetime.today()
        self.save()

    def set_increase_reuse(self):
        self.numberOfReused += 1
        self.save()

    objects = SamplesManager()


class SampleProjectsFieldsValueManager(models.Manager):
    def create_project_field_value(self, field_value):
        new_field_data = self.create(
            sample_id=field_value["sample_id"],
            sampleProjecttField_id=field_value["sampleProjecttField_id"],
            sampleProjectFieldValue=field_value["sampleProjectFieldValue"],
        )
        return new_field_data


class SampleProjectsFieldsValue(models.Model):
    sample_id = models.ForeignKey(
        Samples, on_delete=models.CASCADE, null=True, blank=True
    )
    sampleProjecttField_id = models.ForeignKey(
        SampleProjectsFields, on_delete=models.CASCADE, null=True
    )
    sampleProjectFieldValue = models.CharField(max_length=255, null=True, blank=True)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return "%s" % (self.sampleProjectFieldValue)

    def get_field_value(self):
        return "%s" % (self.sampleProjectFieldValue)

    objects = SampleProjectsFieldsValueManager()


class MoleculeUsedForManager(models.Manager):
    def create_molecule_use_for(self, molecule_use_data):
        new_molecule_use = self.create(
            usedFor=molecule_use_data["usedFor"],
            apps_name=molecule_use_data["apps_name"],
            massiveUse=molecule_use_data["massiveUse"],
        )
        return new_molecule_use


class MoleculeUsedFor(models.Model):
    usedFor = models.CharField(max_length=50)
    apps_name = models.CharField(max_length=50)
    massiveUse = models.BooleanField(default=False)

    def __str__(self):
        return "%s" % (self.usedFor)

    def get_molecule_use_name(self):
        return "%s" % (self.usedFor)

    def get_massive(self):
        return "%s" % (self.massiveUse)

    objects = MoleculeUsedForManager()


class MoleculePreparationManager(models.Manager):
    def create_molecule(self, molecule_data):

        molecule_used_obj = MoleculeType.objects.filter(
            moleculeType__exact=molecule_data["moleculeType"]
        ).last()

        protocol_type_obj = ProtocolType.objects.filter(
            molecule=molecule_used_obj, apps_name__exact=molecule_data["app_name"]
        ).last()
        protocol_used_obj = Protocols.objects.filter(
            name__exact=molecule_data["protocolUsed"], type__exact=protocol_type_obj
        ).last()
        new_molecule = self.create(
            protocolUsed=protocol_used_obj,
            sample=molecule_data["sample"],
            moleculeType=molecule_used_obj,
            state=StatesForMolecule.objects.get(moleculeStateName__exact="Defined"),
            moleculeCodeId=molecule_data["moleculeCodeId"],
            moleculeExtractionDate=molecule_data["moleculeExtractionDate"],
            extractionType=molecule_data["extractionType"],
            moleculeUser=User.objects.get(username__exact=molecule_data["user"]),
        )

        return new_molecule


class MoleculePreparation(models.Model):
    protocolUsed = models.ForeignKey(Protocols, on_delete=models.CASCADE)
    sample = models.ForeignKey(Samples, on_delete=models.CASCADE)
    moleculeType = models.ForeignKey(MoleculeType, on_delete=models.CASCADE)
    state = models.ForeignKey(StatesForMolecule, on_delete=models.CASCADE, null=True)
    moleculeUser = models.ForeignKey(
        User, on_delete=models.CASCADE, null=True, blank=True
    )

    userLotKit_id = models.ForeignKey(
        UserLotCommercialKits, on_delete=models.CASCADE, null=True, blank=True
    )

    moleculeUsedFor = models.ForeignKey(
        MoleculeUsedFor, on_delete=models.CASCADE, null=True, blank=True
    )

    moleculeCodeId = models.CharField(max_length=255)
    extractionType = models.CharField(max_length=50)
    moleculeExtractionDate = models.DateTimeField(auto_now_add=False, null=True)
    numberOfReused = models.IntegerField(default=0)
    usedForMassiveSequencing = models.BooleanField(null=True, blank=True)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return "%s" % (self.moleculeCodeId)

    def get_info_for_display(self):
        extraction_date = self.moleculeExtractionDate.strftime("%d, %B, %Y")
        molecule_info = []
        molecule_info.append(self.moleculeCodeId)
        molecule_info.append(self.state.get_molecule_state())
        molecule_info.append(extraction_date)
        molecule_info.append(self.extractionType)
        molecule_info.append(self.moleculeType.get_name())
        if self.moleculeUsedFor is None:
            molecule_info.append("Not defined yet")
        else:
            molecule_info.append(self.moleculeUsedFor.get_molecule_use_name())
        molecule_info.append(self.protocolUsed.get_name())
        molecule_info.append(self.numberOfReused)
        return molecule_info

    def get_extraction_date(self):
        return "%s" % (self.moleculeExtractionDate.strftime("%B %d, %Y"))

    def get_molecule_id(self):
        return "%s" % (self.pk)

    def get_molecule_code_id(self):
        return "%s" % (self.moleculeCodeId)

    def get_molecule_information(self):
        data = []
        data.append(self.moleculeCodeId)
        data.append(self.moleculeExtractionDate.strftime("%B %d, %Y"))
        data.append(self.protocolUsed.get_name())
        data.append(str(self.pk))
        return data

    def get_sample_name(self):
        return "%s" % (self.sample.get_sample_name())

    def get_sample_obj(self):
        return self.sample

    def get_protocol(self):
        return "%s" % (self.protocolUsed.get_name())

    def get_protocol_obj(self):
        return self.protocolUsed

    def get_state(self):
        return "%s" % (self.state)

    def get_used_for_massive(self):
        return self.usedForMassiveSequencing

    def get_user_lot_kit_obj(self):
        return self.userLotKit_id

    def set_molecule_use(self, use_for_molecule, app_name):
        self.moleculeUsedFor = MoleculeUsedFor.objects.get(
            usedFor__exact=use_for_molecule, apps_name__exact=app_name
        )
        self.save()
        self.usedForMassiveSequencing = self.moleculeUsedFor.get_massive()
        self.save()
        return self

    def set_state(self, state_value):
        self.state = StatesForMolecule.objects.get(moleculeStateName__exact=state_value)
        self.save()

    def set_increase_reuse(self):
        self.numberOfReused += 1
        self.save()

    def set_user_lot_kit(self, lot_kit_name):
        self.userLotKit_id = UserLotCommercialKits.objects.get(
            chipLot__exact=lot_kit_name
        )
        self.save()

    def set_user_lot_kit_obj(self, lot_kit_obj):
        self.userLotKit_id = lot_kit_obj
        self.save()

    objects = MoleculePreparationManager()


class MoleculeParameterValueManager(models.Manager):
    def create_molecule_parameter_value(self, parameter_value):
        new_molecule_parameter_data = self.create(
            moleculeParameter_id=parameter_value["moleculeParameter_id"],
            molecule_id=parameter_value["molecule_id"],
            parameterValue=parameter_value["parameterValue"],
        )
        return new_molecule_parameter_data


class MoleculeParameterValue(models.Model):
    moleculeParameter_id = models.ForeignKey(
        ProtocolParameters, on_delete=models.CASCADE
    )
    molecule_id = models.ForeignKey(MoleculePreparation, on_delete=models.CASCADE)
    parameterValue = models.CharField(max_length=255)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return "%s" % (self.parameterValue)

    def get_param_value(self):
        return "%s" % (self.parameterValue)

    objects = MoleculeParameterValueManager()


class SequencingConfigurationManager(models.Manager):
    def create_new_configuration(self, data):
        try:
            platform_obj = SequencingPlatform.objects.get(pk__exact=data["platformID"])
        except models.SequencingPlatform.DoesNotExist:
            platform_obj = None
        new_sequencer_configuration = self.create(
            platformID=platform_obj, configurationName=data["configurationName"]
        )
        return new_sequencer_configuration


class SequencingConfiguration(models.Model):
    platformID = models.ForeignKey(
        SequencingPlatform, on_delete=models.CASCADE, null=True, blank=True
    )
    configurationName = models.CharField(max_length=255)

    def __str__(self):
        return "%s" % (self.configurationName)

    def get_configuration_name(self):
        return "%s" % (self.configurationName)

    def get_platform_name(self):
        return "%s" % (self.platformID.get_platform_name())

    def get_platform_obj(self):
        return self.platformID

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
            platformID=platform_obj,
            sequencerName=sequencer_value["sequencerName"],
            sequencerDescription=sequencer_value["sequencerDescription"],
            sequencerLocation=sequencer_value["sequencerLocation"],
            sequencerSerialNumber=sequencer_value["sequencerSerialNumber"],
            sequencerState="In Use",
            sequencerOperationStart=sequencer_value["sequencerOperationStart"],
            sequencerNumberLanes=sequencer_value["sequencerNumberLanes"],
        )

        return new_sequencer


class SequencerInLab(models.Model):
    platformID = models.ForeignKey(
        SequencingPlatform, on_delete=models.CASCADE, null=True, blank=True
    )
    sequencerName = models.CharField(max_length=255)
    sequencerDescription = models.CharField(max_length=255, null=True, blank=True)
    sequencerLocation = models.CharField(max_length=255, null=True, blank=True)
    sequencerSerialNumber = models.CharField(max_length=255, null=True, blank=True)
    sequencerState = models.CharField(max_length=50, null=True, blank=True)
    sequencerOperationStart = models.DateField(
        auto_now_add=False, null=True, blank=True
    )
    sequencerOperationEnd = models.DateField(auto_now_add=False, null=True, blank=True)
    sequencerNumberLanes = models.CharField(max_length=5, null=True, blank=True)

    def __str__(self):
        return "%s" % (self.sequencerName)

    def get_sequencer_name(self):
        return "%s" % (self.sequencerName)

    def get_number_of_lanes(self):
        return "%s" % (self.sequencerNumberLanes)

    def get_sequencer_id(self):
        return "%s" % (self.pk)

    def get_sequencing_platform_name(self):
        if self.platformID is None:
            return "Not Defined"
        else:
            return "%s" % (self.platformID.get_platform_name())

    def get_sequencing_platform_id(self):
        if self.platformID is None:
            return ""
        else:
            return "%s" % (self.platformID.get_platform_id())

    def get_all_sequencer_data(self):

        data = []
        if self.platformID is not None:
            platform_name = self.platformID.get_platform_name()
        else:
            platform_name = "Not Defined"
        if self.sequencerOperationStart is not None:
            op_start = self.sequencerOperationStart.strftime("%B %d, %Y")
        else:
            op_start = "Not available date"
        if self.sequencerOperationEnd is not None:
            op_end = self.sequencerOperationEnd.strftime("%B %d, %Y")
        else:
            op_end = "Not available date"
        data.append(platform_name)
        data.append(self.sequencerName)
        data.append(self.sequencerDescription)
        data.append(self.sequencerLocation)
        data.append(self.sequencerSerialNumber)
        data.append(self.sequencerState)
        data.append(op_start)
        data.append(op_end)
        data.append(self.sequencerNumberLanes)
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
