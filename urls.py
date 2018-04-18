# Adaptation of code to django 2.0:
# "from django.conf.urls import include , url" replaced by "import" for "path"
# "url(r'^" replaced by "path('"

from django.urls import path
from django.conf import settings
from django.contrib.auth.views import LoginView
from . import views
from django.conf.urls.static import static
from django.views.generic import ListView, DetailView

urlpatterns = [
    path('', LoginView.as_view(template_name='iSkyLIMS_home/index.html'), name="index"),
<<<<<<< HEAD
    path('contact_email',views.contact_email, name='contact_email'),
    path('thanks', views.thanks, name='thanks'),
=======
    path('about-us', LoginView.as_view(template_name='iSkyLIMS_home/about_us.html'), name="about-us"),
>>>>>>> 274ad31b320993f385bcb3d73def3e0b8eb55d6d
]

if settings.DEBUG:
     urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)

