from django.db import models
from django.contrib.auth.models import User
import datetime
from iSkyLIMS_core.core_config import COLLECTION_INDEX_KITS_DIRECTORY




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

    def get_molecule_type (self):
        return '%s' %(self.molecule.get_name())

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

    def get_type (self):
        return '%s' %(self.type.get_name())


class ProtocolParametersManager(models.Manager) :
    def create_protocol_parameter (self, prot_param_data):
        new_prot_parameter = self.create(protocol_id =prot_param_data['protocol_id'],parameterName = prot_param_data['Parameter name'],
                    parameterDescription = prot_param_data['Description'], parameterOrder = prot_param_data['Order'],
                    parameterUsed = prot_param_data['Used'], parameterMaxValue = prot_param_data['Max Value'],
                    parameterMinValue = prot_param_data['Min Value'] )
        return new_prot_parameter

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

    def get_all_parameter_info(self):
        param_info = []
        param_info.append(self.parameterName)
        param_info.append(self.parameterOrder)
        param_info.append(self.parameterUsed)
        param_info.append(self.parameterMinValue)
        param_info.append(self.parameterMaxValue)
        param_info.append(self.parameterDescription)
        return param_info

    objects = ProtocolParametersManager()

class StatesForSample (models.Model):
    sampleStateName = models.CharField(max_length=50)

    def __str__ (self):
        return '%s' %(self.sampleStateName)

    def get_sample_state (self):
        return '%s' %(self.sampleStateName)


class StatesForMolecule (models.Model):
    moleculeStateName = models.CharField(max_length=50)

    def __str__ (self):
        return '%s' %(self.moleculeStateName)

    def get_molecule_state (self):
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



class CollectionIndexKit (models.Model):
    collectionIndexName = models.CharField(max_length=125)
    version = models.CharField(max_length=80,null=True)
    plateExtension = models.CharField(max_length=125, null=True)
    adapter1 = models.CharField(max_length=125,null=True)
    adapter2 = models.CharField(max_length=125, null=True)
    collectionLibraryFile =  models.FileField(upload_to = COLLECTION_INDEX_KITS_DIRECTORY )
    generatedat = models.DateTimeField(auto_now_add=True, null=True)

    def __srt__ (self):
        return '%s'(self.collectionLibraryName)

    def get_collection_index_information (self):
        if self.adapter2 =='':
            adapter2 = 'Not used on this collection'
        else:
            adapter2 = self.adapter2
        collection_info = []
        collection_info.append(self.indexLibraryName)
        collection_info.append(self.version)
        collection_info.append(self.plateExtension)
        collection_info.append(self.adapter1)
        collection_info.append(self.adapter2)
        collection_info.append(self.collectionLibraryFile)
        return collection_info



class CollectionIndexValues (models.Model):
    collectionIndexKitid = models.ForeignKey(
        CollectionIndexKit,
        on_delete=models.CASCADE)
    defaultWell = models.CharField(max_length=10,null=True)
    index_7 = models.CharField(max_length=25,null=True)
    i_7_seq = models.CharField(max_length=25,null=True)
    index_5 = models.CharField(max_length=25,null=True)
    i_5_seq = models.CharField(max_length=25,null=True)


    def get_index_value_information (self):
        collection_info = []
        collection_info.append(self.defaultWell)
        collection_info.append(self.index_7)
        collection_info.append(self.i_7_seq)
        collection_info.append(self.index_5)
        collection_info.append(self.i_5_seq)

        return collection_info


class ComercialKits (models.Model):
    protocol_id = models.ForeignKey(
                    ProtocolType,
                    on_delete= models.CASCADE, null = True)
    molecule_id = models.ForeignKey(
                    MoleculeType,
                    on_delete= models.CASCADE, null = True, blank = True)
    name = models.CharField(max_length =60)
    provider = models.CharField(max_length =30)
    maximumUses = models.IntegerField(null = True, default = 0)

    cat_number = models.CharField(max_length = 40, null = True, blank = True)
    generatedat = models.DateTimeField(auto_now_add=True, null=True)

    def __str__ (self):
        return '%s' %(self.name)

    def get_name (self):
        return '%s' %(self.name)

    def get_protocol (self):
        return '%s' %(self.protocol_id.get_name())


class UserComercialKits (models.Model):
    user = models.ForeignKey(
                User,
                on_delete=models.CASCADE, null = True)
    basedComercial = models.ForeignKey(
                    ComercialKits,
                    on_delete= models.CASCADE, null = True)
    nickName =  models.CharField(max_length = 50, null = True, blank = True)
    numberOfuses = models.IntegerField(null = True, default = 0)
    chipLot = models.CharField(max_length = 50)
    latestUsedDate = models.DateTimeField(null = True, blank = True)
    expirationDate = models.DateField(auto_now_add=False)
    generatedat = models.DateTimeField(auto_now_add=True, null=True)

    def __str__ (self):
        return '%s' %(self.nickName)

    def get_comercial_kit(self):
        return '%s' %(self.basedComercial.get_name())

    def get_nick_name (self):
        return '%s' %(self.nickName)

    def get_protocol_for_kit (self):
        return '%s' %(self.basedComercial.get_protocol())




class SamplesManager (models.Manager):

    def create_sample (self, sample_data):
        ## Set to null in case that there is not Laboratory defined
        if sample_data['Laboratory'] != '':
            sample_data['Laboratory'] = Laboratory.objects.get(labName__exact = sample_data['Laboratory'])
        else:
            sample_data['Laboratory'] = None
        new_sample = self.create(sampleState = StatesForSample.objects.get(sampleStateName__exact = 'Defined'),
                            laboratory = sample_data['Laboratory'],
                            sampleType = SampleType.objects.get(sampleType__exact = sample_data['Type of Sample']) ,
                            sampleUser = User.objects.get(username__exact = sample_data['user']),
                            sampleCodeID = sample_data['sample_id'] , sampleName =  sample_data['Sample Name'],
                            uniqueSampleID = sample_data['new_unique_value'],
                            species = Species.objects.get(spicesName__exact = sample_data['Species']),
                            sampleExtractionDate = datetime.datetime.strptime(sample_data['Sample Extraction date'],'%Y-%m-%d %H:%M:%S'))

        return new_sample

class Samples (models.Model):
    sampleState = models.ForeignKey(
                StatesForSample,
                on_delete = models.CASCADE, null = True)
    laboratory = models.ForeignKey(
                Laboratory,
                on_delete = models.CASCADE, null = True, blank = True)
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
    #labSampleName = models.CharField(max_length=255, null = True, blank = True)
    sampleExtractionDate = models.DateTimeField(auto_now_add = False, null =True)
    uniqueSampleID = models.CharField(max_length=8, null = True)
    #patientCodeName = models.CharField(max_length=255, null = True)
    sampleCodeID = models.CharField(max_length=60, null = True)
    #reused = models.BooleanField(null = True, blank = True)
    numberOfReused = models.IntegerField(default=0)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__ (self):
        return '%s' %(self.sampleName)

    def get_sample_definition_information (self):
        sample_info = []
        recordeddate=self.sampleExtractionDate.strftime("%d , %B , %Y")
        sample_info.append(self.uniqueSampleID)
        sample_info.append(self.sampleCodeID)
        sample_info.append(self.sampleName)
        sample_info.append(recordeddate)
        sample_info.append(self.sampleType)
        return sample_info

    def get_info_in_defined_state(self):
        sample_info = []
        sample_info.append(self.sampleExtractionDate.strftime("%d , %B , %Y"))
        sample_info.append(self.sampleCodeID)
        sample_info.append(self.sampleName)
        sample_info.append(str(self.pk))
        return sample_info

    def get_info_for_searching (self):
        sample_info = []
        sample_info.append(str(self.pk))
        sample_info.append(self.sampleName)
        sample_info.append(self.sampleState.get_sample_state())
        sample_info.append(self.sampleExtractionDate.strftime("%d , %B , %Y"))
        sample_info.append(self.sampleCodeID)
        sample_info.append(self.sampleType.get_name())
        sample_info.append(self.species.get_name())
        return sample_info

    def get_info_for_display (self):
        sample_info = []
        sample_info.append(self.sampleName)
        sample_info.append(self.sampleCodeID)
        sample_info.append(self.sampleState.get_sample_state())
        sample_info.append(self.generated_at.strftime("%d , %B , %Y"))
        sample_info.append(self.sampleType.get_name())
        sample_info.append(self.species.get_name())
        sample_info.append(self.numberOfReused)
        sample_info.append(self.sampleUser.username)
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

    def get_sample_name (self):
        return '%s' %(self.sampleName)

    def get_register_user(self):
        if self.sampleUser is None:
            return 'Not available'
        else:
            return '%s' %(self.sampleUser)

    def get_registered_sample (self):
        recordeddate=self.generated_at.strftime("%B %d, %Y")
        return '%s' %(recordeddate)

    def set_state (self, state_value):
        self.sampleState = StatesForSample.objects.get(sampleStateName__exact = state_value)
        self.save()

    def set_increase_reuse(self):
        self.numberOfReused += 1
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
    '''
    usedKit =  models.ForeignKey(
                MoleculeUserComercialKits,
                on_delete = models.CASCADE, null = True)
    '''
    moleculeCodeId = models.CharField(max_length=255)
    extractionType = models.CharField(max_length=50)
    moleculeExtractionDate = models.DateTimeField(auto_now_add = False, null =True)
    numberOfReused = models.IntegerField(default=0)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__ (self):
        return "%s" %(self.moleculeCodeId)

    def get_id (self):
        return "%s" %(self.pk)

    def get_info_for_display(self):
        molecule_info =[]
        molecule_info.append(self.moleculeCodeId)
        molecule_info.append(self.state.get_molecule_state())
        molecule_info.append(self.moleculeExtractionDate)
        molecule_info.append(self.extractionType)
        molecule_info.append(self.moleculeUsed.get_name())
        molecule_info.append(self.protocolUsed.get_name())
        molecule_info.append(self.numberOfReused)
        return molecule_info


    def get_extraction_date (self):
        return "%s" %(self.moleculeExtractionDate.strftime("%B %d, %Y"))


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

    def get_protocol_obj(self):
        return self.protocolUsed

    def get_state (self):
        return '%s' %(self.state)

    def set_state (self, state_value):
        self.state = StatesForMolecule.objects.get(moleculeStateName__exact = state_value)
        self.save()

    def set_increase_reuse(self):
        self.numberOfReused += 1
        self.save()

    objects = MoleculePreparationManager()


class MoleculeParameterValueManager (models.Manager):
    def create_molecule_parameter_value (self, parameter_value):
        new_molecule_parameter_data = self.create(moleculeParameter_id = parameter_value['moleculeParameter_id'],
                molecule_id = parameter_value['molecule_id'],
                parameterValue = parameter_value['parameterValue'])
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

    def get_param_value(self):
        return '%s' %(self.parameterValue)

    objects = MoleculeParameterValueManager()
