from django.urls import path
from django.conf import settings

from . import views
from django.conf.urls.static import static
from django.views.generic import ListView, DetailView



urlpatterns = [
    path('',views.index, name = 'index'),
    path('defineNewSamples', views.define_new_samples, name = 'define_new_samples'),
    path('definePatientInformation', views.define_patient_information, name = 'define_patient_information'),

]

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
