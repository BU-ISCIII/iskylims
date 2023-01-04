from django.urls import path
from iSkyLIMS_wetlab.api import views

app_name = "iSkyLIMS_wetlab_api"


urlpatterns = [
    path(
        "fetchRunInformation", views.fetch_run_information, name="fetch_run_information"
    ),
    path(
        "laboratoryData",
        views.get_lab_information_contact,
        name="get_lab_information_contact",
    ),
    path(
        "fetchSamplesOnParameter",
        views.fetch_samples_on_parameter,
        name="fetch_samples_on_parameter",
    ),
    path(
        "fetchSampleInformation",
        views.fetch_sample_information,
        name="fetch_sample_information",
    ),
    path("createSampleData", views.create_sample_data, name="create_sample_data"),
    # path("samplefields", views.sample_fields, name="sample_fields"),
    path("sampleFields/", views.sample_fields, name="sample_fields"),
    path(
        "sampleProjectFields", views.sample_project_fields, name="sample_project_fields"
    ),
    path(
        "statisticsInformation",
        views.statistic_information,
        name="statistic_information",
    ),
    path(
        "summarizeDataInformation",
        views.summarize_data_information,
        name="summarize_data_information",
    ),
    path("updateLab", views.update_lab, name="update_lab"),
]
# urlpatterns = format_suffix_patterns(urlpatterns)
