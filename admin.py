from django.contrib import admin
from iSkyLIMS_wetlab.models import *




class AppAdmin(admin.ModelAdmin):
    list_display = ('run_name','description','csv_file','name')
    def file_link(self, obj):
        if obj.file:
            return "<a href='%s' download>Download</a>" % (obj.file.url,)
        else:
            return "No attachment"
    file_link.allow_tags = True
    file_link.short_description = 'File Download'


class RunErrorsAdmin (admin.ModelAdmin):
    list_display = ('errorCode', 'errorText')


class RunStatesAdmin (admin.ModelAdmin):
    list_display = ('runStateName',)
'''
class StatesForSampleAdmin (admin.ModelAdmin):
    list_display = ('sampleStateName',)
'''
class RunningParametersAdmin(admin.ModelAdmin):
    list_display = ('runName_id','RunID','ExperimentName','RunStartDate')
    #import pdb; pdb.set_trace()

class RunProcessAdmin (admin.ModelAdmin):
    list_display = ('runName','sequencerModel','sampleSheet','generatedat','run_date','runState','runError', 'state','index_library','samples','centerRequestedBy','useSpaceImgMb','useSpaceFastaMb','useSpaceOtherMb')

class ProjectsAdmin (admin.ModelAdmin):
    list_display= ('runprocess_id','user_id','LibraryKit_id','projectName','libraryKit','baseSpaceFile','generatedat','project_run_date')
    #list_display= ('runprocess_id','projectName','procState','libraryKit','baseSpaceFile')

class BaseSpaceLibraryNameAdmin(admin.ModelAdmin):
    list_display = ('libraryName','generatedat')

class IndexLibraryKitAdmin(admin.ModelAdmin):
    list_display = ('indexLibraryName', 'version', 'plateExtension', 'adapter1', 'adapter2', 'indexLibraryFile','generatedat')


class SamplesInProjectAdmin (admin.ModelAdmin):
    list_display = ('project_id','sampleName','barcodeName','pfClusters','percentInProject','yieldMb','qualityQ30')


class StatesForLibraryPreparationAdmin (admin.ModelAdmin):
    list_display = ('libPrepState',)


### New for library preparation
'''
class ProtocolInLabAdmin (admin.ModelAdmin):
    list_display = ('protocolName',)

class LaboratoryAdmin (admin.ModelAdmin):
    list_display = ['labName' , 'labCoding' , 'labLocation' ]
'''
#class LibraryProtocolParametersAdmin (admin.ModelAdmin):
#    list_display = ('parameterName', 'parameterDescription', 'parameterOrder','parameterUsed','parameterMaxValue', 'parameterMinValue')
'''


class ReagentsCommercialKitsAdmin (admin.ModelAdmin):
    list_display = ('name', 'protocol_id', 'provider', 'reagentLibraryName','chipLot','usedDate', 'expirationDate')








class ProtocolParametersAdmin (admin.ModelAdmin):
    list_display = ('protocol_id', 'parameterName', 'parameterOrder', 'parameterUsed', 'parameterMinValue', 'parameterMaxValue', 'parameterDescription')

'''
admin.site.register(RunningParameters , RunningParametersAdmin)
admin.site.register(RunProcess , RunProcessAdmin)
admin.site.register(BaseSpaceLibraryName,BaseSpaceLibraryNameAdmin)
admin.site.register(IndexLibraryKit,IndexLibraryKitAdmin)
admin.site.register(Projects, ProjectsAdmin)
admin.site.register(RunErrors, RunErrorsAdmin)
admin.site.register(RunStates, RunStatesAdmin)
admin.site.register(StatesForLibraryPreparation, StatesForLibraryPreparationAdmin)





'''
admin.site.register(Laboratory, LaboratoryAdmin)
#admin.site.register(LibraryProtocolParameters, LibraryProtocolParametersAdmin)
admin.site.register(NucleotidesComercialKits,NucleotidesComercialKitsAdmin)
admin.site.register(ReagentsCommercialKits,ReagentsCommercialKitsAdmin)


admin.site.register(NAProtocolParameters,NAProtocolParametersAdmin)

admin.site.register(NAProtParamData, NAProtParamDataAdmin)
admin.site.register(SampleType, SampleTypeAdmin)
admin.site.register(Species, SpeciesAdmin)
admin.site.register(StatesForSample, StatesForSampleAdmin)



### new
admin.site.register(MoleculeType,MoleculeTypeAdmin)
admin.site.register(ProtocolType,ProtocolTypeAdmin)
admin.site.register(Protocols,ProtocolsAdmin)
admin.site.register(ProtocolParameters,ProtocolParametersAdmin)
'''
