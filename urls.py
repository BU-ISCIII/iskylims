from django.conf.urls import include , url                                                        
from django.conf import settings                                                                  
from django.contrib.auth.views import LoginView
from . import views                                                                               
from django.conf.urls.static import static                                                        
from django.views.generic import ListView, DetailView                                             
                                                                                                   
urlpatterns = [                                                                                   
    url(r'^$', LoginView.as_view(template_name='iSkyLIMS_home/index.html'), name="index"),
   	#url(r'^infrastructure_request',views.infrastructure_request, name='infrastructure_service'),  
]                                                                                                 
                                                                                                   
if settings.DEBUG:                                                                                
     urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)                  
