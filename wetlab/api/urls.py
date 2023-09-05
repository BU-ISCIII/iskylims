from django.urls import path
from wetlab.api import views

app_name = "wetlab_api"


urlpatterns = [
    path("run-info", views.fetch_run_information, name="fetch_run_information"),
    path(
        "lab-data",
        views.get_lab_information_contact,
        name="get_lab_information_contact",
    ),
    path(
        "sample-info",
        views.fetch_sample_information,
        name="fetch_sample_information",
    ),
    path("create-sample", views.create_sample_data, name="create_sample_data"),
    path("sample-fields", views.sample_fields, name="sample_fields"),
    path("projects-fields", views.sample_project_fields, name="sample_project_fields"),
    path(
        "stats-info",
        views.statistic_information,
        name="statistic_information",
    ),
    path(
        "summarized-info",
        views.summarize_data_information,
        name="summarize_data_information",
    ),
    path("update-lab", views.update_lab, name="update_lab"),
]
