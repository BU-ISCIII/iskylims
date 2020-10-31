from django.urls import path
from iSkyLIMS_drylab.api import views
#from . import views
from rest_framework.urlpatterns import format_suffix_patterns
app_name = 'iSkyLIMS_drylab_api'

urlpatterns = [
#     path('jobslist',views.jobs_list, name='jobs_list'),
#     path('<int:state>/jobsliststate',views.jobs_list_state, name='jobs_list_state'),

    #path('',views.jobs_list, name='jobs_list'),
    #path('list/<int:state>',views.jobs_list_state, name='jobs_list_state'),
    #path('update/<str:service>',views.api_update_job, name='api_update_job'),
    #path('updatestate/<str:service>',views.api_update_state, name='api_update_state'),
    #path('drylab/api/pipeline/<int:pipeline>',views.get_pipeline, name='get_pipeline'),
    #path('samples/<str:service>', views.get_samplesinproject, name='get_samplesinproject'),


    path('services/', views.service_list,  name='service_list'),
    path('resolution', views.resolution_data, name='resolution_data'),
    path('samplesInService', views.samples_in_service, name = 'samples_in_service'),
    path('updateResolution',views.update_resolution_to_in_progress, name = 'update_resolution_to_in_progress'),


    #path('runfolder/<int:project>', views.get_runfolder,  name='get_runfolder'),
]
urlpatterns = format_suffix_patterns(urlpatterns)
