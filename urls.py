from django.conf.urls import include , url                                                        
from django.conf import settings                                                                  
                                                                                                   
from . import views                                                                               
from django.conf.urls.static import static                                                        
from django.views.generic import ListView, DetailView                                             
#from polls.models import Document                                                                
                                                                                                   
                                                                                                   
urlpatterns = [                                                                                   
 	url(r'^user_creation',views.user_creation, name='user_creation'),                             
  	url(r'^user_edit',views.user_edit, name='user_edit'),                                         
 ]                                                                                                 
 #+ static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)                                
                                                                                                   
if settings.DEBUG:                                                                                
     urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)                  
