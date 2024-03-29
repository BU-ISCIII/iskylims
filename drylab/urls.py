# Generic imports
from django.conf import settings
from django.conf.urls.static import static
from django.urls import path

# Local imports
import drylab.views

urlpatterns = [
    path("", drylab.views.index, name="index"),
    path("add-resolution", drylab.views.add_resolution, name="add_resolution"),
    path("add-in-progress", drylab.views.add_in_progress, name="add_in_progress"),
    path("add-delivery", drylab.views.add_delivery, name="add_delivery"),
    path(
        "add-samples-service",
        drylab.views.add_samples_service,
        name="add_samples_service_in_service",
    ),
    path(
        "counseling-request",
        drylab.views.counseling_request,
        name="counseling_service",
    ),
    path(
        "configuration-email",
        drylab.views.configuration_email,
        name="configuration_email",
    ),
    path(
        "delete-samples-service",
        drylab.views.delete_samples_service,
        name="delete_samples_service",
    ),
    path(
        "detail-pipeline=<int:pipeline_id>",
        drylab.views.detail_pipeline,
        name="detail_pipeline",
    ),
    path(
        "define-pipeline",
        drylab.views.define_pipeline,
        name="define_pipeline_service",
    ),
    path(
        "display-service=<int:service_id>/",
        drylab.views.display_service,
        name="display_service",
    ),
    path(
        "infrastructure-request",
        drylab.views.infrastructure_request,
        name="infrastructure_service",
    ),
    path(
        "manage-pipelines",
        drylab.views.manage_pipelines,
        name="manage_pipelines",
    ),
    path(
        "pending-services",
        drylab.views.pending_services,
        name="peding_services",
    ),
    path(
        "sequencing-request",
        drylab.views.request_seq_service,
        name="request_seq_service",
    ),
    path("search-service", drylab.views.search_service, name="search_service"),
    path(
        "add-on-hold",
        drylab.views.add_on_hold,
        name="add_on_hold",
    ),
    path("stats-by-user", drylab.views.stats_by_user, name="stats_by_user"),
    path(
        "stats-services-time",
        drylab.views.stats_by_services_request,
        name="stats_by_services_request",
    ),
    path(
        "upload-file-delete=<int:file_id>",
        drylab.views.upload_service_file_delete,
        name="upload_service_file_delete",
    ),
]

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
