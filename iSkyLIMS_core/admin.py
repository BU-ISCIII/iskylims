from django.contrib import admin
from iSkyLIMS_core.models import *

class StateInCountryAdmin(admin.ModelAdmin):
    list_display = ["stateName", 'apps_name']

class CityAdmin(admin.ModelAdmin):
    list_display = ["cityName", "apps_name"]


class CommercialKitsAdmin (admin.ModelAdmin):
    list_display =( 'name','provider' , 'platformKits' , 'cat_number')

class LabRequestAdmin (admin.ModelAdmin):
    list_display = ['labName' , 'labNameCoding' ,'labContactName', 'labPhone', 'labEmail' ]

class MoleculeTypeAdmin( admin.ModelAdmin):
    list_display = ('moleculeType',)

class MoleculePreparationAdmin(admin.ModelAdmin):
    list_display = ('moleculeCodeId', 'state','sample', 'moleculeType', 'extractionType', 'protocolUsed',
                    'moleculeExtractionDate', 'moleculeUsedFor','numberOfReused')
    list_filter = ('generated_at',)
    search_fields = ('sample__startswith',)

class MoleculeUsedForAdmin(admin.ModelAdmin):
    list_display = ['usedFor' ,'apps_name', 'massiveUse']

class PatientCoreAdmin(admin.ModelAdmin):
    list_display = ('patientCode', 'patientName', 'patientSurname', 'patientSex')
    search_fields = ("patientCode__icontains", )

class PatientSexAdmin(admin.ModelAdmin):
    list_display = ('sex',)

class PatientProjectsAdmin(admin.ModelAdmin):
    list_display = ('projectName', 'projectDescription')

class ProtocolsAdmin(admin.ModelAdmin):
    list_display = ('name', 'type', 'description')

class ProtocolTypeAdmin( admin.ModelAdmin):
    list_display = ('protocol_type', 'molecule', 'apps_name')


class ProtocolParametersAdmin (admin.ModelAdmin):
    list_display = ('protocol_id', 'parameterName', 'parameterOrder', 'parameterUsed', 'parameterMinValue', 'parameterMaxValue', 'parameterDescription')
    list_filter = ('protocol_id',)
    search_fields = ('parameterName__startswith',)

class MoleculeParameterValueAdmin (admin.ModelAdmin):
    list_display = ('moleculeParameter_id', 'molecule_id','parameterValue')

class SamplesAdmin(admin.ModelAdmin):
    list_display = ('sampleCodeID', 'sampleName', 'sampleState', 'labRequest', 'sampleType', 'sampleUser', 'species','sampleProject',
                    'sampleEntryDate', 'uniqueSampleID',  'numberOfReused','sequencingDate' )
    list_filter = ('generated_at',)
    search_fields = ('sampleName__icontains',)

class SampleTypeAdmin(admin.ModelAdmin):
    list_display = ('sampleType','apps_name' ,'optional_fields')
'''
class SampleProjectBelongsAdmin(admin.ModelAdmin):
    list_display = ('projectName', 'projectManager', 'projectDescription', 'contactEmail', 'contactPhone', 'contactComments')
'''
class SampleProjectsAdmin(admin.ModelAdmin):
    list_display = ['sampleProjectName', 'sampleProjectManager', 'sampleProjectContact', 'sampleProjectDescription', 'apps_name']


class SampleProjectsFieldsAdmin(admin.ModelAdmin):
    list_display  = ['sampleProjects_id', 'sampleProjectFieldName', 'sampleProjectFieldType', 'sampleProjectOptionList', 'sampleProjectFieldOrder', 'sampleProjectFieldUsed']

class SampleProjectsFieldsValueAdmin(admin.ModelAdmin):
    list_display = ['sample_id', 'sampleProjecttField_id', 'sampleProjectFieldValue']
    list_filter = ['sample_id']
    search_fields = ('sampleProjecttField_id__sampleProjectFieldName', )

class SequencingPlatformAdmin(admin.ModelAdmin):
    list_display = ('platformName', 'companyName', 'sequencingTecnology')

class SequencingConfigurationAdmin(admin.ModelAdmin):
    list_display = ['platformID', 'configurationName']

class SequencerInLabAdmin(admin.ModelAdmin):
    list_display = ('sequencerName','platformID', 'sequencerDescription', 'sequencerLocation', 'sequencerSerialNumber', 'sequencerState', 'sequencerOperationStart', 'sequencerOperationEnd','sequencerNumberLanes')

class SpeciesAdmin (admin.ModelAdmin):
    list_display= ('speciesName', 'refGenomeName', 'refGenomeSize' , 'refGenomeID' )

class StatesForMoleculeAdmin (admin.ModelAdmin):
    list_display = ('moleculeStateName',)

class StatesForSampleAdmin (admin.ModelAdmin):
    list_display = ('sampleStateName',)

class UserLotCommercialKitsAdmin (admin.ModelAdmin):
    list_display = ('user', 'basedCommercial', 'runOut', 'numberOfuses', 'chipLot', 'latestUsedDate','expirationDate')


admin.site.register(StateInCountry, StateInCountryAdmin)
admin.site.register(City, CityAdmin)
admin.site.register(CommercialKits, CommercialKitsAdmin)
admin.site.register(LabRequest, LabRequestAdmin)

admin.site.register(MoleculeType,MoleculeTypeAdmin)
admin.site.register(MoleculeUsedFor,MoleculeUsedForAdmin)
admin.site.register(MoleculePreparation,MoleculePreparationAdmin)
admin.site.register(PatientCore,PatientCoreAdmin)
admin.site.register(PatientSex,PatientSexAdmin)
admin.site.register(PatientProjects,PatientProjectsAdmin)
admin.site.register(ProtocolType,ProtocolTypeAdmin)
admin.site.register(Protocols,ProtocolsAdmin)
admin.site.register(ProtocolParameters,ProtocolParametersAdmin)
admin.site.register(MoleculeParameterValue, MoleculeParameterValueAdmin)
admin.site.register(SampleType, SampleTypeAdmin)
admin.site.register(Samples, SamplesAdmin)
admin.site.register(SampleProjects,SampleProjectsAdmin)
admin.site.register(SampleProjectsFields,SampleProjectsFieldsAdmin)
admin.site.register(SampleProjectsFieldsValue,SampleProjectsFieldsValueAdmin)

admin.site.register(SequencingPlatform, SequencingPlatformAdmin)
admin.site.register(SequencingConfiguration, SequencingConfigurationAdmin)
admin.site.register(SequencerInLab, SequencerInLabAdmin)

admin.site.register(Species, SpeciesAdmin)
# admin.site.register(SampleProjectBelongs, SampleProjectBelongsAdmin)
admin.site.register(StatesForMolecule, StatesForMoleculeAdmin)
admin.site.register(StatesForSample, StatesForSampleAdmin)
admin.site.register(UserLotCommercialKits, UserLotCommercialKitsAdmin)
