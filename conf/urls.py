from django.urls import include, path
from django.contrib import admin
from django.contrib.auth.views import LoginView

urlpatterns = [
    path('',include('iSkyLIMS_core.urls')),
    path('background', LoginView.as_view(template_name='iSkyLIMS_core/background.html'), name="background"),
    path('bu-isciii', LoginView.as_view(template_name='iSkyLIMS_core/bu_isciii.html'), name="about-us"),
    path('contact', LoginView.as_view(template_name='iSkyLIMS_core/contact.html'), name="contact"),
    path('admin/', admin.site.urls),
    path('wetlab/', include('iSkyLIMS_wetlab.urls')),
    path('drylab/',include('iSkyLIMS_drylab.urls')),
    path('utils/',include('django_utils.urls')),
    path('clinic/',include ('iSkyLIMS_clinic.urls')),
    path('accounts/', include('django.contrib.auth.urls')),
    
    # REST FRAMEWORK URLS
    path('drylab/api/', include('iSkyLIMS_drylab.api.urls')),
]
