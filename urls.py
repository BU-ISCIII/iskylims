from django.urls import path
from django.conf import settings

from . import views
from django.conf.urls.static import static
from django.views.generic import ListView, DetailView



urlpatterns = [
    path('',views.index, name = 'index'),
    path('addResultData', views.add_result_data, name = 'add_result_data'),
    path('addCommercialKit', views.add_commercial_kit, name='add_commercial_kit'),
    path('addUserLotCommercialKit', views.add_user_lot_commercial_kit, name ='add_user_lot_commercial_kit'),
    path('createProtocol', views.create_protocol, name = 'create_protocol'),
    path('defineProtocolParameters=<int:protocol_id>', views.define_protocol_parameters, name = 'define_protocol_parameters'),
    path('defineNewSamples', views.define_new_samples, name = 'define_new_samples'),
    path('definePatientInformation', views.define_patient_information, name = 'define_patient_information'),
    path('displayProtocol=<int:protocol_id>', views.display_protocol, name = 'display_protocol'),
    path('defineResultProtocol', views.define_result_protocol, name = 'define_result_protocol'),
    path('defineResultProtocolParameters=<int:result_protocol_id>', views.define_result_protocol_parameters, name = 'define_result_protocol_parameters'),
    path('displayResultProtocol=<int:result_protocol_id>', views.display_result_protocol, name = 'display_result_protocol'),
    path('displayPatientInformation=<int:patient_id>/', views.display_patient_information, name = 'display_patient_information'),
    path('displaySampleInfo=<int:sample_c_id>', views.display_sample_info, name = 'display_sample_info'),
    path('pendingToUpdate', views.pending_to_update, name = 'pending_to_update'),
    path('searchSample', views.search_sample, name = 'search_sample'),
    path('searchPatient', views.search_patient, name = 'search_patient'),
    path('setMoleculeValues', views.set_molecule_values, name = 'set_molecule_values'),

]

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
