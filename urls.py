from django.conf.urls import include , url
from django.conf import settings

from . import views
from django.conf.urls.static import static
from django.views.generic import ListView, DetailView
#from polls.models import Document


urlpatterns = [
    # ex: /polls/
    url(r'^$', views.index, name='index'),
    url(r'^service_request',views.service_request, name='request_service')
] 
#+ static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)


