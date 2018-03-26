##Adaptation of code to django 2.0:
# "from django.conf.urls import include , url" replaced by "from django.urls import path" 
# "url(r'^" replaced by "path('"

from django.urls import path
from django.conf import settings

from . import views
from django.conf.urls.static import static
from django.views.generic import ListView, DetailView

urlpatterns = [
     path('', views.index, name='index'),
     path('service_request_<str:serviceRequestType>',views.service_request, name='service_request'),
     path('counseling_request',views.counseling_request, name='counseling_service'),
     path('infrastructure_request',views.infrastructure_request, name='infrastructure_service'),
     path('pendingServices', views.pending_services, name ='peding_services'),
     path('display_service=<int:service_id>/',views.display_service, name= 'display_service'),
] 

if settings.DEBUG:
     urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)


