"""URLs for the standard deployment health endpoint."""

from django.urls import path

from .views import health_check


urlpatterns = [
    path("", health_check, name="deployment-health"),
]
