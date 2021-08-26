## Adaptation of code to django 2.0:
# "from django.conf.urls import include , url" replaced by "from django.urls import path"
# "url(r'^" replaced by "path('"

from django.urls import path
from django.conf import settings                                                                  
                                                                                                   
from . import views                                                                               
from django.conf.urls.static import static                                                        
from django.views.generic import ListView, DetailView                                             
                                                                                                   
urlpatterns = [                                                                                   
 	path('user_creation',views.user_creation, name='user_creation'),                             
  	path('user_edit',views.user_edit, name='user_edit'),                                         
 ]                                                                                                 
 #+ static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)                                
                                                                                                   
if settings.DEBUG:                                                                                
     urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)                  
