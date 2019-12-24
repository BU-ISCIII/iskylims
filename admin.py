from django.contrib import admin
from iSkyLIMS_core.models import *

class CommercialKitsAdmin (admin.ModelAdmin):
    list_display =( 'name','provider', 'protocol_id' ,'maximumUses','cat_number')

class SamplesOriginAdmin (admin.ModelAdmin):
    list_display = ['originName' , 'originNameCoding' , 'location' ]

class MoleculeTypeAdmin( admin.ModelAdmin):
    list_display = ('moleculeType',)

class MoleculePreparationAdmin(admin.ModelAdmin):
    list_display = ('moleculeCodeId', 'state','sample', 'moleculeUsed', 'extractionType', 'protocolUsed',
                    'moleculeExtractionDate', 'numberOfReused')

class ProtocolsAdmin(admin.ModelAdmin):
    list_display = ('type', 'name',  'description')

class ProtocolTypeAdmin( admin.ModelAdmin):
    list_display = ('protocol_type', 'molecule', 'apps_name')


class ProtocolParametersAdmin (admin.ModelAdmin):
    list_display = ('protocol_id', 'parameterName', 'parameterOrder', 'parameterUsed', 'parameterMinValue', 'parameterMaxValue', 'parameterDescription')

class MoleculeParameterValueAdmin (admin.ModelAdmin):
    list_display = ('moleculeParameter_id', 'molecule_id','parameterValue')

class SamplesAdmin(admin.ModelAdmin):
    list_display = ('sampleCodeID', 'sampleState', 'samplesOrigin', 'sampleType', 'sampleUser', 'species', 'sampleName',
                    'sampleEntryDate', 'uniqueSampleID',  'numberOfReused','sequencingDate' )

class SampleTypeAdmin(admin.ModelAdmin):
    list_display = ('sampleType',)

class SampleProjectBelongsAdmin(admin.ModelAdmin):
    list_display = ('projectName', 'projectManager', 'projectDescription', 'contactEmail', 'contactPhone', 'contactComments')


class SpeciesAdmin (admin.ModelAdmin):
    list_display= ('speciesName', 'refGenomeName', 'refGenomeSize' , 'refGenomeID' )

class StatesForMoleculeAdmin (admin.ModelAdmin):
    list_display = ('moleculeStateName',)

class StatesForSampleAdmin (admin.ModelAdmin):
    list_display = ('sampleStateName',)

class UserLotCommercialKitsAdmin (admin.ModelAdmin):
    list_display = ('user', 'basedCommercial', 'nickName', 'numberOfuses', 'chipLot', 'latestUsedDate','expirationDate')


admin.site.register(CommercialKits, CommercialKitsAdmin)
admin.site.register(SamplesOrigin, SamplesOriginAdmin)

admin.site.register(MoleculeType,MoleculeTypeAdmin)
admin.site.register(MoleculePreparation,MoleculePreparationAdmin)
admin.site.register(ProtocolType,ProtocolTypeAdmin)
admin.site.register(Protocols,ProtocolsAdmin)
admin.site.register(ProtocolParameters,ProtocolParametersAdmin)
admin.site.register(MoleculeParameterValue, MoleculeParameterValueAdmin)
admin.site.register(SampleType, SampleTypeAdmin)
admin.site.register(Samples, SamplesAdmin)
admin.site.register(Species, SpeciesAdmin)
admin.site.register(SampleProjectBelongs, SampleProjectBelongsAdmin)
admin.site.register(StatesForMolecule, StatesForMoleculeAdmin)
admin.site.register(StatesForSample, StatesForSampleAdmin)
admin.site.register(UserLotCommercialKits, UserLotCommercialKitsAdmin)
