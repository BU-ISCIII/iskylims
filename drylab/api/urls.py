from django.urls import path

from drylab.api import views

app_name = "drylab_api"

urlpatterns = [
    path("services", views.service_list, name="service_list"),
    path("serviceFullData", views.service_full_data, name="service_full_data"),
    path("resolution", views.resolution_data, name="resolution_data"),
    path("samples_in_service", views.samples_in_service, name="samples_in_service"),
    path("updateState", views.update_state, name="update_state"),
    path("createDelivery", views.create_delivery, name="create_delivery"),
]
# urlpatterns = format_suffix_patterns(urlpatterns)
