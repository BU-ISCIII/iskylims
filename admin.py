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

class LibraryPreparationAdmin (admin.ModelAdmin):
    list_display = ('registerUser', 'molecule_id', 'sample_id', 'protocol_id', 'libPrepState', 
            'user_sample_sheet', 'libPrepCodeID', 'userSampleID', 'projectInSampleSheet','samplePlate',
            'sampleWell', 'i7IndexID', 'i7Index', 'i5IndexID', 'i5Index', 'singlePairedEnd', 'lengthRead',
            'numberOfReused')

class LibParameterValueAdmin (admin.ModelAdmin):
    list_display = ('parameter_id', 'library_id', 'parameterValue')


class LibraryPoolAdmin(admin.ModelAdmin):
    list_display = ('registerUser', 'poolState', 'poolName', 'poolCodeID')

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
    list_display = ('runName','sequencerModel','sampleSheet','generatedat','run_date','runError', 'state','index_library','samples','centerRequestedBy','useSpaceImgMb','useSpaceFastaMb','useSpaceOtherMb')

class ProjectsAdmin (admin.ModelAdmin):
    list_display= ('runprocess_id','user_id','LibraryKit_id','projectName','libraryKit','baseSpaceFile','generatedat','project_run_date')
    #list_display= ('runprocess_id','projectName','procState','libraryKit','baseSpaceFile')
'''
class BaseSpaceLibraryNameAdmin(admin.ModelAdmin):
    list_display = ('libraryName','generatedat')
'''
class CollectionIndexKitAdmin(admin.ModelAdmin):
    list_display = ('collectionIndexName', 'version', 'plateExtension', 'adapter1', 'adapter2', 'collectionIndexFile','generatedat')


class SamplesInProjectAdmin (admin.ModelAdmin):
    list_display = ('project_id','sampleName','barcodeName','pfClusters','percentInProject','yieldMb','qualityQ30')


class StatesForLibraryPreparationAdmin (admin.ModelAdmin):
    list_display = ('libPrepState',)

class StatesForPoolAdmin (admin.ModelAdmin):
    list_display = ('poolState',)



admin.site.register(LibraryPreparation, LibraryPreparationAdmin)
admin.site.register(LibParameterValue, LibParameterValueAdmin)

admin.site.register(LibraryPool, LibraryPoolAdmin)

admin.site.register(RunningParameters , RunningParametersAdmin)
admin.site.register(RunProcess , RunProcessAdmin)
#admin.site.register(BaseSpaceLibraryName,BaseSpaceLibraryNameAdmin)
admin.site.register(CollectionIndexKit,CollectionIndexKitAdmin)
admin.site.register(Projects, ProjectsAdmin)
admin.site.register(RunErrors, RunErrorsAdmin)
admin.site.register(RunStates, RunStatesAdmin)
admin.site.register(StatesForLibraryPreparation, StatesForLibraryPreparationAdmin)
admin.site.register(StatesForPool, StatesForPoolAdmin)
