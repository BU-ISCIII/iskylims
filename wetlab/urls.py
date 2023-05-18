# Generic imports
from django.conf import settings
from django.conf.urls.static import static
from django.urls import path

# Local imports
import wetlab.views

urlpatterns = [
    path("", wetlab.views.index, name="index"),
    path(
        "AddBasespaceLibrary/",
        wetlab.views.add_basespace_library,
        name="add_basespace_library",
    ),
    path("addCommercialKit", wetlab.views.add_commercial_kit, name="add_commercial_kit"),
    path(
        "addCollectionIndexKit",
        wetlab.views.add_collection_index_kit,
        name="add_collection_index_kit",
    ),
    path(
        "addUserLotCommercialKit",
        wetlab.views.add_user_lot_commercial_kit,
        name="add_user_lot_commercial_kit",
    ),
    path("AnnualReport/", wetlab.views.annual_report, name="annual_report"),
    path(
        "change_project_libKit=<int:project_id>",
        wetlab.views.change_project_libKit,
        name="change_project_libKit",
    ),
    path("ChangeRunName=<int:run_id>", wetlab.views.change_run_name, name="change_run_name"),
    path("compareSamples", wetlab.views.compare_samples, name="compare_samples"),
    path("createNewRun/", wetlab.views.create_new_run, name="create_new_run"),
    path("createNextSeqRun/", wetlab.views.create_nextseq_run, name="create_nextseq_run"),
    path("createPool/", wetlab.views.create_pool, name="create_pool"),
    path("createProtocol", wetlab.views.create_protocol, name="create_protocol"),
    path("configurationEmail", wetlab.views.configuration_email, name="configuration_email"),
    path("configurationSamba", wetlab.views.configuration_samba, name="configuration_samba"),
    path("configurationTest/", wetlab.views.configuration_test, name="configuration_test"),
    path(
        "defineAdditionalKits=<int:protocol_id>",
        wetlab.views.define_additional_kits,
        name="define_additional_kits",
    ),
    path("defineMoleculeUses", wetlab.views.define_molecule_uses, name="define_molecule_uses"),
    path(
        "defineProtocolParameters=<int:protocol_id>",
        wetlab.views.define_protocol_parameters,
        name="define_protocol_parameters",
    ),
    path(
        "defineSampleProjects",
        wetlab.views.define_sample_projects,
        name="define_sample_projects",
    ),
    path(
        "defineTypeOfSamples",
        wetlab.views.define_type_of_samples,
        name="define_type_of_samples",
    ),
    path(
        "defineSampleProjectFields=<int:sample_project_id>",
        wetlab.views.define_sample_projects_fields,
        name="define_sample_projects_fields",
    ),
    path(
        "DisplayCollectionIndex=<int:collection_index_id>/",
        wetlab.views.display_collection_index,
        name="display_collection_index",
    ),
    path("displaySample=<int:sample_id>/", wetlab.views.display_sample, name="display_sample"),
    path(
        "displaySampleInRun=<int:sample_run_id>/",
        wetlab.views.display_sample_in_run,
        name="display_sample_in_run",
    ),
    path(
        "displayProject=<int:project_id>/",
        wetlab.views.display_project,
        name="display_project",
    ),
    path(
        "displayProtocol=<int:protocol_id>",
        wetlab.views.display_protocol,
        name="display_protocol",
    ),
    path("displayRun=<int:run_id>/", wetlab.views.display_run, name="display_run"),
    path("displaySample=<int:sample_id>/", wetlab.views.display_sample, name="display_sample"),
    path(
        "displaySampleInRun=<int:sample_run_project_id>/",
        wetlab.views.display_sample_in_run,
        name="display_sample_in_run",
    ),
    path(
        "displaySampleProject=<int:sample_project_id>/",
        wetlab.views.display_sample_project,
        name="display_sample_project",
    ),
    path(
        "displayTypeOfSample=<int:sample_type_id>/",
        wetlab.views.display_type_of_sample,
        name="display_type_of_sample",
    ),
    path(
        "displayUserLotKit=<int:user_kit_id>/",
        wetlab.views.display_user_lot_kit,
        name="display_user_lot_kit",
    ),
    path(
        "skipCancelSituation", wetlab.views.skip_cancel_situation, name="skip_cancel_situation"
    ),
    path(
        "handlingLibraryPreparations",
        wetlab.views.handling_library_preparations,
        name="handling_library_preparations",
    ),
    path("handlingMolecules", wetlab.views.handling_molecules, name="handling_molecules"),
    path("initialSettings", wetlab.views.initial_settings, name="initial_settings"),
    path(
        "lastRunBySequencer/", wetlab.views.last_run_by_sequencer, name="last_run_by_sequencer"
    ),
    path("incompletedRuns", wetlab.views.incompleted_runs, name="incompleted_runs"),
    path(
        "modifyAdditionalKits=<int:protocol_id>",
        wetlab.views.modify_additional_kits,
        name="modify_additional_kits",
    ),
    path(
        "modifyProtocolFields=<int:protocol_id>/",
        wetlab.views.modify_protocol_fields,
        name="modify_protocol_fields",
    ),
    path(
        "modifySampleProjectFields=<int:sample_project_id>/",
        wetlab.views.modify_sample_project_fields,
        name="modify_sample_project_fields",
    ),
    path("MonthlyReport/", wetlab.views.monthly_report, name="montly_report"),
    path(
        "pendingSamplePreparations",
        wetlab.views.pending_sample_preparations,
        name="pending_sample_preparations",
    ),
    path("pendingToUpdate/", wetlab.views.pending_to_update, name="pending_to_update"),
    path("QuarterReport/", wetlab.views.quarter_report, name="quarter_report"),
    path("recordSamples", wetlab.views.record_samples, name="record_samples"),
    path(
        "repeatLibraryPreparation",
        wetlab.views.repeat_library_preparation,
        name="repeat_library_preparation",
    ),
    path(
        "repeatMoleculeExtraction",
        wetlab.views.repeat_molecule_extraction,
        name="repeat_molecule_extraction",
    ),
    path("repeatPool", wetlab.views.repeat_pool, name="repeat_pool"),
    path("retryErrorRun", wetlab.views.retry_error_run, name="retry_error_run"),
    path(
        "searchCollectionIndexLibrary",
        wetlab.views.search_collection_index_library,
        name="search_collection_index_library",
    ),
    path("searchProject/", wetlab.views.search_project, name="search_project"),
    path("searchRun/", wetlab.views.search_run, name="search_run"),
    path("searchSample/", wetlab.views.search_sample, name="search_sample"),
    path("searchUserLotKit", wetlab.views.search_user_lot_kit, name="search_user_lot_kit"),
    path("sequencerInventory", wetlab.views.sequencer_inventory, name="sequencer_inventory"),
    path("setMoleculeValues", wetlab.views.set_molecule_values, name="set_molecule_values"),
    # path('setLibraryValues', wetlab.views.set_library_values, name = 'set_library_values'),
    path(
        "sequencerConfiguration",
        wetlab.views.sequencer_configuration,
        name="sequencer_configuration",
    ),
    path("StatsExperiment/", wetlab.views.stats_experiment, name="stats_experiment"),
    path("StatsLibrary/", wetlab.views.stats_per_library, name="stats_per_library"),
    path("StatsPerSequencer/", wetlab.views.stats_per_sequencer, name="stats_per_sequencer"),
    path(
        "StatsPerResearcher/", wetlab.views.stats_per_researcher, name="stats_per_researcher"
    ),
    path("StatsPerTime/", wetlab.views.stats_per_time, name="stats_per_time"),
    path(
        "userCommercialKitInventory/",
        wetlab.views.user_commercial_kit_inventory,
        name="user_commercial_kit_inventory",
    ),
]
# + static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
