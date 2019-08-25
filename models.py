from django.db import models
from django.contrib.auth.models import User
import datetime





class Laboratory (models.Model):
    labName = models.CharField(max_length=50)
    labCoding = models.CharField(max_length=10)
    labLocation = models.CharField(max_length=255)

    def __str__ (self):
        return '%s' %(self.labName)

    def get_name(self):
        return '%s' %(self.labName)

    def get_lab_code (self):
        return '%s' %(self.labCoding)


class MoleculeType (models.Model):
    moleculeType = models.CharField(max_length = 30)

    def __str__ (self):
        return '%s' %(self.moleculeType)

    def get_name (self):
        return '%s' %(self.moleculeType)


class ProtocolType (models.Model):
    molecule = models.ForeignKey(
                    MoleculeType,
                    on_delete= models.CASCADE, null = True, blank = True)
    protocol_type = models.CharField(max_length = 40)

    def __str__ (self) :
        return '%s' %(self.protocol_type)

    def get_name (self):
        return '%s' %(self.protocol_type)

class Protocols (models.Model):
    type =  models.ForeignKey(
                    ProtocolType,
                    on_delete= models.CASCADE)
    name = models.CharField(max_length = 40)
    description = models.CharField(max_length = 160, null = True, blank = True)

    def __str__ (self) :
        return '%s' %(self.name)

    def get_name (self):
        return '%s' %(self.name)

class ProtocolParameters (models.Model):
    protocol_id = models.ForeignKey(
                    Protocols,
                    on_delete= models.CASCADE)
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

class StatesForSample (models.Model):
    sampleStateName = models.CharField(max_length=50)

    def __str__ (self):
        return '%s' %(self.sampleStateName)


class StatesForMolecule (models.Model):
    moleculeStateName = models.CharField(max_length=50)

    def __str__ (self):
        return '%s' %(self.moleculeStateName)


class SampleType (models.Model):
    sampleType = models.CharField(max_length=50)


    def __str__ (self):
        return '%s' %(self.sampleType)

    def get_name(self):
        return '%s' %(self.sampleType)


class Species (models.Model):
    spicesName = models.CharField(max_length=50)
    refGenomeName = models.CharField(max_length=255)
    refGenomeSize = models.CharField(max_length=100)
    refGenomeID = models.CharField(max_length=255)

    def __str__ (self):
        return '%s' %(self.spicesName)

    def get_name(self):
        return '%s' %(self.spicesName)

class SamplesManager (models.Manager):

    def create_sample (self, sample_data):
        new_sample = self.create(sampleState = StatesForSample.objects.get(sampleStateName__exact = 'Defined'),
                            laboratory = Laboratory.objects.get(labName__exact = sample_data['laboratory']),
                            patientCodeName = sample_data['patientCodeName'],
                            sampleType = SampleType.objects.get(sampleType__exact = sample_data['sampleType']) ,
                            sampleUser = User.objects.get(username__exact = sample_data['user']),
                            sampleCodeID = sample_data['sample_id'] , sampleName =  sample_data['sampleName'],
                            uniqueSampleID = sample_data['new_unique_value'],
                            labSampleName = sample_data['labSampleName'],
                            species = Species.objects.get(spicesName__exact = sample_data['species']),
                            sampleExtractionDate = datetime.datetime.strptime(sample_data['extractionDate'],'%Y-%m-%d %H:%M:%S'))
        return new_sample

class Samples (models.Model):
    sampleState = models.ForeignKey(
                StatesForSample,
                on_delete = models.CASCADE, null = True)
    laboratory = models.ForeignKey(
                Laboratory,
                on_delete = models.CASCADE, null = True)
    sampleType = models.ForeignKey(
                SampleType,
                on_delete = models.CASCADE, null = True)


    sampleUser = models.ForeignKey(
                User,
                on_delete=models.CASCADE, null = True)

    species = models.ForeignKey(
                Species,
                on_delete=models.CASCADE, null = True)
    sampleName = models.CharField(max_length=255, null = True)
    labSampleName = models.CharField(max_length=255, null = True, blank = True)
    sampleExtractionDate = models.DateTimeField(auto_now_add = False, null =True)
    uniqueSampleID = models.CharField(max_length=8, null = True)
    patientCodeName = models.CharField(max_length=255, null = True)
    sampleCodeID = models.CharField(max_length=60, null = True)
    #reused = models.BooleanField(null = True, blank = True)
    numberOfReused = models.IntegerField(default=0)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__ (self):
        return '%s' %(self.sampleName)

    def get_sample_definition_information (self):
        recordeddate=self.sampleExtractionDate.strftime("%d , %B , %Y")
        return '%s;%s;%s' %(recordeddate, self.sampleCodeID, self.sampleType.sampleType)

    def get_info_in_defined_state(self):
        sample_info = []
        sample_info.append(self.sampleExtractionDate.strftime("%d , %B , %Y"))
        sample_info.append(self.sampleCodeID)
        sample_info.append(self.sampleType.get_name())
        sample_info.append(str(self.pk))
        return sample_info

    def get_extraction_date (self):
        recordeddate=self.sampleExtractionDate.strftime("%B %d, %Y")
        return '%s' %(recordeddate)

    def get_sample_code (self):
        return '%s' %(self.sampleCodeID)

    def get_sample_id(self):
        return '%s' %(self.id)

    def get_sample_type (self):
        return '%s' %(self.sampleType.get_name())


    def get_register_user(self):
        if self.sampleUser is None:
            return 'Not avilable'
        else:
            return '%s' %(self.sampleUser)

    def get_registered_sample (self):
        recordeddate=self.generated_at.strftime("%B %d, %Y")
        return '%s' %(recordeddate)

    def set_state (self, state_value):
        self.sampleState = StatesForSample.objects.get(sampleStateName__exact = state_value)
        self.save()

    objects = SamplesManager()

class MoleculePreparationManager (models.Manager):

    def  create_molecule (self, molecule_data) :
        protocol_used = molecule_data['protocolUsed']
        molecule_used = MoleculeType.objects.get(moleculeType__exact = molecule_data['moleculeUsed'])
        new_molecule = self.create( protocolUsed = protocol_used,
        sample =  molecule_data['sample'],
        moleculeUsed =  molecule_used,
        state = StatesForMolecule.objects.get(moleculeStateName__exact = 'Defined'),
        moleculeCodeId =  molecule_data['moleculeCodeId'],
        moleculeExtractionDate = molecule_data['moleculeExtractionDate'],
        extractionType =  molecule_data['extractionType'],
        numberOfReused = molecule_data['numberOfReused'])

        return new_molecule

class MoleculePreparation (models.Model):
    protocolUsed =  models.ForeignKey(
                Protocols,
                on_delete = models.CASCADE)
    sample = models.ForeignKey(
                Samples,
                on_delete = models.CASCADE)
    moleculeUsed = models.ForeignKey(
                MoleculeType,
                on_delete = models.CASCADE)
    state =  models.ForeignKey(
                StatesForMolecule,
                on_delete = models.CASCADE, null = True)
    moleculeCodeId = models.CharField(max_length=255)
    extractionType = models.CharField(max_length=50)
    moleculeExtractionDate = models.DateTimeField(auto_now_add = False, null =True)
    numberOfReused = models.IntegerField(default=0)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__ (self):
        return "%s" %(self.moleculeCodeId)

    def get_id (self):
        return "%s" %(self.pk)



    def get_molecule_code_id (self):
        return '%s' %(self.moleculeCodeId)

    def get_molecule_information(self):
        data =[]
        data.append(self.moleculeCodeId)
        data.append(self.moleculeExtractionDate.strftime("%B %d, %Y"))
        data.append(self.protocolUsed.get_name())
        data.append(str(self.pk))
        return data

    def get_sample_obj(self):
        return self.sample

    def get_protocol (self):
        return '%s' %(self.protocolUsed.get_name())

    def set_state (self, state_value):
        self.state = StatesForMolecule.objects.get(moleculeStateName__exact = state_value)
        self.save()
    objects = MoleculePreparationManager()


class MoleculeParameterValueManager (models.Manager):
    def create_molecule_parameter_value (self, molecule_parameter_value):
        new_molecule_parameter_data = self.create(moleculeParameter_id = molecule_parameter_value['moleculeParameter_id'],
                molecule_id = molecule_parameter_value['molecule_id'],
                parameterValue = molecule_parameter_value['parameterValue'])
        return new_molecule_parameter_data

class MoleculeParameterValue (models.Model):
    moleculeParameter_id = models.ForeignKey(
                    ProtocolParameters,
                    on_delete= models.CASCADE)
    molecule_id = models.ForeignKey(
                MoleculePreparation,
                on_delete= models.CASCADE)
    parameterValue = models.CharField(max_length=255)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__ (self):
        return '%s' %(self.parameterValue)

    def get_parameter_information (self):
        return '%s' %(self.parameterValue)

    objects = MoleculeParameterValueManager()
