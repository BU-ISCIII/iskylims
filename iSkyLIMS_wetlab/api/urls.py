from django.urls import path
from iSkyLIMS_wetlab.api import views
from rest_framework.urlpatterns import format_suffix_patterns

app_name = 'iSkyLIMS_wetlab_api'


urlpatterns = [
    path("sampleList/", views.sample_list, name="service_list"),
    path("createSampleData", views.create_sample_data, name="create_sample_data")
    ]
urlpatterns = format_suffix_patterns(urlpatterns)
