from django.urls import path
from django.conf import settings

from . import views
from django.conf.urls.static import static
from django.views.generic import ListView, DetailView



urlpatterns = [
    path('',views.index, name = 'index'),
    #path('addResultData', views.add_result_data, name = 'add_result_data'),
    path('addCommercialKit', views.add_commercial_kit, name='add_commercial_kit'),
    path('addUserLotCommercialKit', views.add_user_lot_commercial_kit, name ='add_user_lot_commercial_kit'),
    path('assignProject', views.assign_project, name='assign_project'),
    path('createNewPatientProject', views.create_new_patient_project, name= 'create_new_patient_project'),
    path('createProtocol', views.create_protocol, name = 'create_protocol'),
    path('createSampleProjects', views.create_sample_projects, name = 'create_sample_projects'),
    path('defineExtractionMolecules', views.define_extraction_molecules, name = 'define_extraction_molecules'),
    path('defineProjectFields=<int:project_id>', views.define_project_fields, name = 'define_project_fields'),
    path('defineNewPatient', views.define_new_patient, name = 'define_new_patient'),
    path('defineNewPatientHistory' ,  views.define_new_patient_history, name = 'define_new_patient_history'),
    path('defineNewSamples', views.define_new_samples, name = 'define_new_samples'),
    path('defineProtocolParameters=<int:protocol_id>', views.define_protocol_parameters, name = 'define_protocol_parameters'),
    path('defineSampleProjectFields=<int:sample_project_id>', views.define_sample_projects_fields, name = 'define_sample_projects_fields'),
    path('displayPatientProject=<int:project_id>', views.display_patient_project, name = 'display_patient_project'),
    path('displayProtocol=<int:protocol_id>', views.display_protocol, name = 'display_protocol'),
    #path('defineResultProtocol', views.define_result_protocol, name = 'define_result_protocol'),
    #path('defineResultProtocolParameters=<int:result_protocol_id>', views.define_result_protocol_parameters, name = 'define_result_protocol_parameters'),
    #path('displayResultProtocol=<int:result_protocol_id>', views.display_result_protocol, name = 'display_result_protocol'),
    path('displayPatientInformation=<int:patient_id>/', views.display_patient_information, name = 'display_patient_information'),
    path('displaySampleClinicInfo=<int:sample_c_id>', views.display_sample_clinic_info, name = 'display_sample_clinic_info'),
    path('displaySampleProject=<int:sample_project_id>/', views.display_sample_project, name = 'display_sample_project'),
    path('pendingToUpdate', views.pending_to_update, name = 'pending_to_update'),
    path('searchSample', views.search_sample, name = 'search_sample'),
    path('searchPatient', views.search_patient, name = 'search_patient'),
    path('setMoleculeValues', views.set_molecule_values, name = 'set_molecule_values'),

]

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
