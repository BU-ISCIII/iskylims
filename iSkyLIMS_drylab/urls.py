# Generic imports
from django.urls import path
from django.conf import settings
from django.shortcuts import redirect
from django.conf.urls.static import static
from django.views.generic import ListView, DetailView

# Local imports
import iSkyLIMS_drylab.views

urlpatterns = [
    path("", iSkyLIMS_drylab.views.index, name="index"),
    path("addResolution", iSkyLIMS_drylab.views.add_resolution, name="add_resolution"),
    path(
        "addInProgress", iSkyLIMS_drylab.views.add_in_progress, name="add_in_progress"
    ),
    path("addDelivery", iSkyLIMS_drylab.views.add_delivery, name="add_delivery"),
    path(
        "addSamplesInService",
        iSkyLIMS_drylab.views.add_samples_in_service,
        name="add_samples_in_service",
    ),
    path(
        "counselingRequest",
        iSkyLIMS_drylab.views.counseling_request,
        name="counseling_service",
    ),
    path(
        "configurationEmail",
        iSkyLIMS_drylab.views.configuration_email,
        name="configuration_email",
    ),
    path(
        "configurationTest",
        iSkyLIMS_drylab.views.configuration_test,
        name="configuration_test",
    ),
    path(
        "deleteSamplesInService",
        iSkyLIMS_drylab.views.delete_samples_in_service,
        name="deleteSamplesInService",
    ),
    path(
        "detailPipeline=<int:pipeline_id>",
        iSkyLIMS_drylab.views.detail_pipeline,
        name="detail_pipeline",
    ),
    path(
        "definePipelineService",
        iSkyLIMS_drylab.views.define_pipeline_service,
        name="define_pipeline_service",
    ),
    path(
        "display_service=<int:service_id>/",
        iSkyLIMS_drylab.views.display_service,
        name="display_service",
    ),
    path(
        "infrastructureRequest",
        iSkyLIMS_drylab.views.infrastructure_request,
        name="infrastructure_service",
    ),
    path(
        "managePipelines",
        iSkyLIMS_drylab.views.manage_pipelines,
        name="manage_pipelines",
    ),
    path("openSessions", iSkyLIMS_drylab.views.open_sessions, name="open_sessions"),
    path(
        "pendingServices",
        iSkyLIMS_drylab.views.pending_services,
        name="peding_services",
    ),
    path(
        "requestSequencingService",
        iSkyLIMS_drylab.views.request_sequencing_service,
        name="request_sequencing_service",
    ),
    path("searchService", iSkyLIMS_drylab.views.search_service, name="search_service"),
    # path('service_request_<str:serviceRequestType>',views.service_request, name='service_request'),
    path("serviceInWaitingInfo", iSkyLIMS_drylab.views.service_in_waiting_info, name="service_in_waiting_info"),
    path("statsByUser", iSkyLIMS_drylab.views.stats_by_user, name="stats_by_user"),
    path(
        "statsByServicesRequest",
        iSkyLIMS_drylab.views.stats_by_services_request,
        name="stats_by_services_request",
    ),
    path("userLogin", iSkyLIMS_drylab.views.user_login, name="user_login"),
    # path('multipleFiles', views.multiple_files, name='multiple_files'),
    path(
        "uploadServiceFileDelete=<int:file_id>",
        iSkyLIMS_drylab.views.upload_service_file_delete,
        name="upload_service_file_delete",
    ),
]

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
