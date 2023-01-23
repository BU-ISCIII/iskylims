from django.urls import path
from iSkyLIMS_drylab.api import views

app_name = 'iSkyLIMS_drylab_api'

urlpatterns = [
    path('services', views.service_list, name='service_list'),
    path('serviceFullData', views.service_full_data, name='service_full_data'),
    path('resolution', views.resolution_data, name='resolution_data'),
    path('samplesInService', views.samples_in_service, name='samples_in_service'),
    path('updateResolution', views.update_resolution, name='update_resolution'),
    path('createDelivery', views.create_delivery, name='create_delivery'),
    ]
# urlpatterns = format_suffix_patterns(urlpatterns)
