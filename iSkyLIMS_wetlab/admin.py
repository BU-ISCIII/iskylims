from django.contrib import admin
from iSkyLIMS_wetlab.models import *


class LibraryPreparationAdmin (admin.ModelAdmin):
    list_display = ( 'libPrepCodeID','registerUser', 'molecule_id', 'sample_id', 'protocol_id', 'libPrepState',
            'user_sample_sheet', 'userSampleID', 'projectInSampleSheet','samplePlate',
            'sampleWell', 'i7IndexID', 'i7Index', 'i5IndexID', 'i5Index', 'numberOfReused')
    search_fields = ("libPrepCodeID__icontains", )


class LibParameterValueAdmin (admin.ModelAdmin):
    list_display = ('parameter_id', 'library_id', 'parameterValue')


class LibraryPoolAdmin(admin.ModelAdmin):
    list_display = ('registerUser', 'poolState', 'poolName', 'platform' ,'poolCodeID', 'runProcess_id')

class libPreparationUserSampleSheetAdmin(admin.ModelAdmin):
    list_display = ['registerUser', 'collectionIndexKit_id', 'sequencingConfiguration', 'sampleSheet', 'application','instrument', 'adapter1', 'adapter2', 'assay','reads']

class AdditionaKitsLibraryPreparationAdmin(admin.ModelAdmin):
    list_display = ['kitName', 'protocol_id','commercialKit_id']

class AdditionalUserLotKitAdmin(admin.ModelAdmin):
    list_display = ['lib_prep_id', 'additionalLotKits' , 'userLotKit_id']

class RunErrorsAdmin (admin.ModelAdmin):
    list_display = ('errorCode', 'errorText')


class RunStatesAdmin (admin.ModelAdmin):
    list_display = ('runStateName',)

class RunningParametersAdmin(admin.ModelAdmin):
    list_display = ('runName_id','RunID','ExperimentName','RunStartDate')


class RunProcessAdmin (admin.ModelAdmin):
    list_display = ('runName','state','usedSequencer','sampleSheet','run_date','runError', 'index_library','samples','centerRequestedBy','useSpaceImgMb','useSpaceFastaMb','useSpaceOtherMb')
    list_filter= ('generatedat',)
    search_fields = ("runName__icontains", )


class ProjectsAdmin (admin.ModelAdmin):
    list_display= ('projectName','user_id','LibraryKit_id','libraryKit','baseSpaceFile','generatedat','project_run_date')
    list_filter= ('generatedat',)
    search_fields = ("projectName__icontains", )

class CollectionIndexKitAdmin(admin.ModelAdmin):
    list_display = ('collectionIndexName', 'version', 'plateExtension', 'adapter1', 'adapter2', 'collectionIndexFile','generatedat')


class SamplesInProjectAdmin (admin.ModelAdmin):
    list_display = ['sampleName','project_id' ,'runProcess_id', 'user_id','barcodeName','pfClusters','percentInProject','yieldMb','qualityQ30']
    list_filter= ('generated_at',)
    search_fields = ("sampleName__icontains", )

class StatesForLibraryPreparationAdmin (admin.ModelAdmin):
    list_display = ('libPrepState',)

class StatesForPoolAdmin (admin.ModelAdmin):
    list_display = ('poolState',)

class StatsRunSummaryAdmin (admin.ModelAdmin):
    list_display = ('runprocess_id', 'level', 'yieldTotal', 'projectedTotalYield', 'aligned', 'errorRate', 'intensityCycle', 'biggerQ30', 'stats_summary_run_date')


class StatsRunReadAdmin(admin.ModelAdmin):
    list_display = ('runprocess_id', 'read', 'lane', 'tiles', 'density', 'cluster_PF', 'phas_prephas', 'reads', 'reads_PF', 'q30', 'yields', 'cyclesErrRated', 'aligned',
                    'errorRate', 'errorRate35',  'errorRate50', 'errorRate75', 'errorRate100' ,'intensityCycle' ,'stats_read_run_date')

class RawDemuxStatsAdmin(admin.ModelAdmin):
    list_display = ('runprocess_id', 'project_id', 'defaultAll' ,'rawYield', 'rawYieldQ30', 'rawQuality', 'PF_Yield', 'PF_YieldQ30' ,'PF_QualityScore')

class RawTopUnknowBarcodesdmin(admin.ModelAdmin):
    list_display = ('runprocess_id', 'lane_number' , 'top_number', 'count', 'sequence')


class StatsLaneSummaryAdmin(admin.ModelAdmin):
    list_display = ('runprocess_id', 'project_id', 'defaultAll' ,'lane', 'pfCluster', 'percentLane', 'perfectBarcode', 'oneMismatch' ,'yieldMb', 'biggerQ30', 'meanQuality' )

class StatsFlSummaryAdmin(admin.ModelAdmin):
    list_display = ('runprocess_id', 'project_id', 'defaultAll' , 'flowRawCluster', 'flowPfCluster', 'flowYieldMb', 'sampleNumber')



class GraphicsStatsAdmin(admin.ModelAdmin):
    list_display = ('runprocess_id', 'folderRunGraphic', 'cluserCountGraph', 'flowCellGraph', 'intensityByCycleGraph', 'heatMapGraph', 'histogramGraph', 'sampleQcGraph')

class SambaConnectionDataAdmin(admin.ModelAdmin):
    list_display = ['SAMBA_IP_SERVER','SAMBA_HOST_NAME', 'SAMBA_PORT_SERVER', 'SAMBA_USER_ID', 'SAMBA_USER_PASSWORD']


class ConfigSettingAdmin(admin.ModelAdmin):
    list_display = ['configurationName', 'configurationValue']

admin.site.register(LibraryPreparation, LibraryPreparationAdmin)
admin.site.register(LibParameterValue, LibParameterValueAdmin)

admin.site.register(LibraryPool, LibraryPoolAdmin)
admin.site.register(libPreparationUserSampleSheet, libPreparationUserSampleSheetAdmin)

admin.site.register(AdditionaKitsLibraryPreparation, AdditionaKitsLibraryPreparationAdmin)
admin.site.register(AdditionalUserLotKit, AdditionalUserLotKitAdmin)

admin.site.register(RunningParameters , RunningParametersAdmin)
admin.site.register(RunProcess , RunProcessAdmin)
#admin.site.register(BaseSpaceLibraryName,BaseSpaceLibraryNameAdmin)
admin.site.register(CollectionIndexKit,CollectionIndexKitAdmin)
admin.site.register(Projects, ProjectsAdmin)
admin.site.register(RunErrors, RunErrorsAdmin)
admin.site.register(RunStates, RunStatesAdmin)
admin.site.register(StatesForLibraryPreparation, StatesForLibraryPreparationAdmin)
admin.site.register(StatesForPool, StatesForPoolAdmin)
admin.site.register(SamplesInProject,SamplesInProjectAdmin)

admin.site.register(StatsRunSummary, StatsRunSummaryAdmin)
admin.site.register(StatsRunRead, StatsRunReadAdmin)
admin.site.register(RawDemuxStats, RawDemuxStatsAdmin)
admin.site.register(RawTopUnknowBarcodes, RawTopUnknowBarcodesdmin)
admin.site.register(StatsLaneSummary, StatsLaneSummaryAdmin)
admin.site.register(StatsFlSummary, StatsFlSummaryAdmin)
admin.site.register(GraphicsStats, GraphicsStatsAdmin)

admin.site.register(SambaConnectionData,SambaConnectionDataAdmin)
admin.site.register(ConfigSetting, ConfigSettingAdmin)
