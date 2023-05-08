from django.conf import settings
from django.conf.urls.static import static
from django.contrib.auth.views import LoginView
from django.urls import path

from . import views

urlpatterns = [
    path("", views.index, name="index"),
    path("addNewContacts", views.add_new_contacts, name="add_new_contacts"),
    # path('', LoginView.as_view(template_name='iSkyLIMS_home/index.html'), name="index"),
    path("contact", views.contact, name="contact"),
    path("thanks", views.thanks, name="thanks"),
    path(
        "about-us",
        LoginView.as_view(template_name="iSkyLIMS_home/about_us.html"),
        name="about-us",
    ),
]

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
