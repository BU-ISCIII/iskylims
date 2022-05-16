from django.urls import path
from iSkyLIMS_wetlab.api import views

app_name = 'iSkyLIMS_wetlab_api'


urlpatterns = [
    path("laboratoryData", views.get_lab_information_contact, name="get_lab_information_contact"),
    path("createSampleData", views.create_sample_data, name="create_sample_data"),
    path("samplefields", views.sample_fields, name="sample_fields"),
    path("sampleFields/", views.sample_fields, name="sample_fields"),
    path("sampleProjectFields", views.sample_project_fields, name="sample_project_fields"),
    path("updateLab", views.update_lab, name="update_lab"),
    ]
# urlpatterns = format_suffix_patterns(urlpatterns)
