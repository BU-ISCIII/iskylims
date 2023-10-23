from django.urls import path

from drylab.api import views

app_name = "drylab_api"

urlpatterns = [
    path("services", views.service_list, name="service_list"),
    path("service-data", views.service_full_data, name="service_full_data"),
    path("resolution", views.resolution_data, name="resolution_data"),
    path("samples_in_service", views.samples_in_service, name="samples_in_service"),
    path("update-state", views.update_state, name="update_state"),
    path("create-delivery", views.create_delivery, name="create_delivery"),
]
