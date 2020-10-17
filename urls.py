from django.urls import path
from django.conf import settings
from django.shortcuts import redirect

from . import views
from django.conf.urls.static import static
from django.views.generic import ListView, DetailView

urlpatterns = [
    path('', views.index, name='index'),
    path('addResolution',views.add_resolution, name='add_resolution'),
    path('addinProgress=<int:resolution_id>/', views.add_in_progress, name='add_in_progress'),
    path('addDelivery=<int:resolution_id>/', views.add_delivery, name= 'add_delivery'),
    path('counseling_request',views.counseling_request, name='counseling_service'),
    path('configurationEmail', views.configuration_email, name='configuration_email'),
    path('configurationTest',views.configuration_test, name='configuration_test'),
    path('detailPipeline=<int:pipeline_id>', views.detail_pipeline, name = 'detail_pipeline'),
    path('definePipelineService', views.define_pipeline_service, name = 'define_pipeline_service'),
    path('display_service=<int:service_id>/',views.display_service, name= 'display_service'),
    path('infrastructure_request',views.infrastructure_request, name='infrastructure_service'),
    path('managePipelines', views.manage_pipelines, name = 'manage_pipelines'),
    path('openSessions', views.open_sessions, name='open_sessions'),
    path('pendingServices', views.pending_services, name ='peding_services'),
    path('searchService', views.search_service, name='search_service'),
    path('service_request_<str:serviceRequestType>',views.service_request, name='service_request'),
    path('statsByDateUser',views.stats_by_date_user, name = 'stats_by_date_user'),
    path('statsByServicesRequest',views.stats_by_services_request, name = 'stats_by_services_request'),
    path('statsBySamplesProcessed',views.stats_by_samples_processed, name = 'stats_by_samples_processed'),
    path('statsTimeDelivery', views.stats_time_delivery, name = 'stats_time_delivery'),
    path('userLogin', views.user_login, name = 'user_login'),
    path('newRequestService', views.new_request_service, name = 'new_request_service'),

]

if settings.DEBUG:
     urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
