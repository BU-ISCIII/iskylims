from django.urls import path
from django.conf import settings
from django.shortcuts import redirect

from . import views
from django.conf.urls.static import static
from django.views.generic import ListView, DetailView

urlpatterns = [
    path("", views.index, name="index"),
    path("addResolution",views.add_resolution, name="add_resolution"),
    path("addInProgress", views.add_in_progress, name="add_in_progress"),
    path("addDelivery", views.add_delivery, name= "add_delivery"),
    path("addSamplesInService", views.add_samples_in_service, name = "add_samples_in_service"),
    path("counselingRequest",views.counseling_request, name="counseling_service"),
    path("configurationEmail", views.configuration_email, name="configuration_email"),
    path("configurationTest",views.configuration_test, name="configuration_test"),
    path("deleteSamplesInService", views.delete_samples_in_service, name= "deleteSamplesInService"),
    path("detailPipeline=<int:pipeline_id>", views.detail_pipeline, name = "detail_pipeline"),
    path("definePipelineService", views.define_pipeline_service, name = "define_pipeline_service"),
    path("displayService=<int:service_id>/",views.display_service, name= "display_service"),
    path("infrastructureRequest",views.infrastructure_request, name="infrastructure_service"),
    path("managePipelines", views.manage_pipelines, name = "manage_pipelines"),
    path("openSessions", views.open_sessions, name="open_sessions"),
    path("pendingServices", views.pending_services, name ="peding_services"),
    path("requestSequencingService", views.request_sequencing_service, name = "request_sequencing_service"),
    path("searchService", views.search_service, name="search_service"),
    #path("service_request_<str:serviceRequestType>",views.service_request, name="service_request"),
    path("serviceInWaitingInfo", views.service_in_waiting_info, name="service_in_waiting_info"),
    path("statsByUser",views.stats_by_user, name = "stats_by_user"),
    path("statsByServicesRequest",views.stats_by_services_request, name = "stats_by_services_request"),
    path("userLogin", views.user_login, name = "user_login"),

    #path("multipleFiles", views.multiple_files, name="multiple_files"),
    path("uploadServiceFileDelete=<int:file_id>", views.upload_service_file_delete, name="upload_service_file_delete"),

]

if settings.DEBUG:
     urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
