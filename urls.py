from django.urls import path
from django.conf import settings

from . import views
from django.conf.urls.static import static
from django.views.generic import ListView, DetailView



urlpatterns = [
    path('',views.index, name = 'index'),
    path('defineNewSamples', views.define_new_samples, name = 'define_new_samples'),
    path('definePatientInformation', views.define_patient_information, name = 'define_patient_information'),
    path('defineResultProtocol', views.define_result_protocol, name = 'define_result_protocol'),
    path('defineResultProtocolParameters=<int:result_protocol_id>', views.define_result_protocol_parameters, name = 'define_result_protocol_parameters'),
    path('displayResultProtocol=<int:result_protocol_id>', views.display_result_protocol, name = 'display_result_protocol'),
    path('displaySampleInfo=<int:sample_c_id>', views.display_sample_info, name = 'display_sample_info'),
    path('searchSample', views.search_sample, name = 'search_sample'),

]

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
