from django.conf import settings
from django.conf.urls.static import static
from django.urls import path

from . import views

urlpatterns = [
    path("", views.index, name="index"),
    path("addNewContacts", views.add_new_contacts, name="add_new_contacts"),
    path("contact", views.contact, name="contact"),
    path("thanks", views.thanks, name="thanks"),
]

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
