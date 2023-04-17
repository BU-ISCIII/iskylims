from django.contrib import admin
from iSkyLIMS_wetlab.models import *


class LibraryPreparationAdmin(admin.ModelAdmin):
    list_display = (
        "lib_prep_code_id",
        "register_user",
        "molecule_id",
        "sample_id",
        "protocol_id",
        "lib_prep_state",
        "user_sample_sheet",
        "user_sample_id",
        "project_in_samplesheet",
        "sample_plate",
        "sample_well",
        "i7_index_id",
        "i7_index",
        "i5_index_id",
        "i5_index",
        "reused_number",
    )
    search_fields = ("lib_prep_code_id__icontains",)


class LibParameterValueAdmin(admin.ModelAdmin):
    list_display = ("parameter_id", "library_id", "parameter_value")


class LibraryPoolAdmin(admin.ModelAdmin):
    list_display = (
        "register_user",
        "pool_state",
        "pool_name",
        "platform",
        "pool_code_id",
        "run_process_id",
    )


class libPreparationUserSampleSheetAdmin(admin.ModelAdmin):
    list_display = [
        "register_user",
        "collection_index_kit_id",
        "sequencing_configuration",
        "sample_sheet",
        "application",
        "instrument",
        "adapter_1",
        "adapter_2",
        "assay",
        "reads",
    ]


class AdditionaKitsLibraryPreparationAdmin(admin.ModelAdmin):
    list_display = ["kit_name", "protocol_id", "commercial_kit_id"]


class AdditionalUserLotKitAdmin(admin.ModelAdmin):
    list_display = ["lib_prep_id", "additional_lot_kits", "user_lot_kit_id"]


class RunErrorsAdmin(admin.ModelAdmin):
    list_display = ("error_code", "error_text")


class RunStatesAdmin(admin.ModelAdmin):
    list_display = ("run_state_name",)


class RunningParametersAdmin(admin.ModelAdmin):
    list_display = ("run_name_id", "run_id", "experiment_name", "run_start_date")


class RunProcessAdmin(admin.ModelAdmin):
    list_display = (
        "run_name",
        "state",
        "used_sequencer",
        "sample_sheet",
        "run_date",
        "run_error",
        "index_library",
        "samples",
        "center_requested_by",
        "use_space_img_mb",
        "use_space_fasta_mb",
        "use_space_other_mb",
    )
    list_filter = ("generated_at",)
    search_fields = ("run_name__icontains",)


class ProjectsAdmin(admin.ModelAdmin):
    list_display = (
        "project_name",
        "user_id",
        "library_kit_id",
        "library_kit",
        "base_space_file",
        "generated_at",
        "project_run_date",
    )
    list_filter = ("generated_at",)
    search_fields = ("project_name__icontains",)


class CollectionIndexKitAdmin(admin.ModelAdmin):
    list_display = (
        "collection_index_name",
        "version",
        "plate_extension",
        "adapter_1",
        "adapter_2",
        "collection_index_file",
        "generated_at",
    )


class SamplesInProjectAdmin(admin.ModelAdmin):
    list_display = [
        "sample_name",
        "project_id",
        "run_process_id",
        "user_id",
        "barcode_name",
        "pf_clusters",
        "percent_in_project",
        "yield_mb",
        "quality_q30",
    ]
    list_filter = ("generated_at",)
    search_fields = ("sample_name__icontains",)


class StatesForLibraryPreparationAdmin(admin.ModelAdmin):
    list_display = ("lib_prep_state",)


class StatesForPoolAdmin(admin.ModelAdmin):
    list_display = ("pool_state",)


class StatsRunSummaryAdmin(admin.ModelAdmin):
    list_display = (
        "runprocess_id",
        "level",
        "yield_total",
        "projected_total_yield",
        "aligned",
        "error_rate",
        "intensity_cycle",
        "bigger_q30",
        "stats_summary_run_date",
    )


class StatsRunReadAdmin(admin.ModelAdmin):
    list_display = (
        "runprocess_id",
        "read",
        "lane",
        "tiles",
        "density",
        "cluster_pf",
        "phas_prephas",
        "reads",
        "reads_pf",
        "q30",
        "yields",
        "cycles_err_rated",
        "aligned",
        "error_rate",
        "error_rate_35",
        "error_rate_50",
        "error_rate_75",
        "error_rate_100",
        "intensity_cycle",
        "stats_read_run_date",
    )


class RawDemuxStatsAdmin(admin.ModelAdmin):
    list_display = (
        "runprocess_id",
        "project_id",
        "default_all",
        "raw_yield",
        "raw_yield_q30",
        "raw_quality",
        "pf_yield",
        "pf_yield_q30",
        "pf_quality_score",
    )


class RawTopUnknowBarcodesdmin(admin.ModelAdmin):
    list_display = ("runprocess_id", "lane_number", "top_number", "count", "sequence")


class StatsLaneSummaryAdmin(admin.ModelAdmin):
    list_display = (
        "runprocess_id",
        "project_id",
        "default_all",
        "lane",
        "pf_cluster",
        "percent_lane",
        "perfect_barcode",
        "one_mismatch",
        "yield_mb",
        "bigger_q30",
        "mean_quality",
    )


class StatsFlSummaryAdmin(admin.ModelAdmin):
    list_display = (
        "runprocess_id",
        "project_id",
        "default_all",
        "flow_raw_cluster",
        "flow_pf_cluster",
        "flow_yield_mb",
        "sample_number",
    )


class GraphicsStatsAdmin(admin.ModelAdmin):
    list_display = (
        "runprocess_id",
        "folder_run_graphic",
        "cluser_count_graph",
        "flow_cell_graph",
        "intensity_by_cycle_graph",
        "heatmap_graph",
        "histogram_graph",
        "sample_qc_graph",
    )


class SambaConnectionDataAdmin(admin.ModelAdmin):
    list_display = [
        "ip_server",
        "host_name",
        "port_server",
        "user_id",
        "user_password",
    ]


class ConfigSettingAdmin(admin.ModelAdmin):
    list_display = ["configuration_name", "configuration_value"]


class RunConfigurationTestAdmin(admin.ModelAdmin):
    list_display = ["run_test_name", "run_test_folder"]


admin.site.register(LibPrepare, LibraryPreparationAdmin)
admin.site.register(LibParameterValue, LibParameterValueAdmin)
admin.site.register(LibraryPool, LibraryPoolAdmin)
admin.site.register(libUserSampleSheet, libPreparationUserSampleSheetAdmin)
admin.site.register(AdditionaKitsLibPrepare, AdditionaKitsLibraryPreparationAdmin)
admin.site.register(AdditionalUserLotKit, AdditionalUserLotKitAdmin)
admin.site.register(RunningParameters, RunningParametersAdmin)
admin.site.register(RunProcess, RunProcessAdmin)
admin.site.register(CollectionIndexKit, CollectionIndexKitAdmin)
admin.site.register(Projects, ProjectsAdmin)
admin.site.register(RunErrors, RunErrorsAdmin)
admin.site.register(RunStates, RunStatesAdmin)
admin.site.register(LibPrepareStates, StatesForLibraryPreparationAdmin)
admin.site.register(PoolStates, StatesForPoolAdmin)
admin.site.register(SamplesInProject, SamplesInProjectAdmin)
admin.site.register(StatsRunSummary, StatsRunSummaryAdmin)
admin.site.register(StatsRunRead, StatsRunReadAdmin)
admin.site.register(RawDemuxStats, RawDemuxStatsAdmin)
admin.site.register(RawTopUnknowBarcodes, RawTopUnknowBarcodesdmin)
admin.site.register(StatsLaneSummary, StatsLaneSummaryAdmin)
admin.site.register(StatsFlSummary, StatsFlSummaryAdmin)
admin.site.register(GraphicsStats, GraphicsStatsAdmin)
admin.site.register(SambaConnectionData, SambaConnectionDataAdmin)
admin.site.register(ConfigSetting, ConfigSettingAdmin)
admin.site.register(RunConfigurationTest, RunConfigurationTestAdmin)
