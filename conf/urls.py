from django.urls import include, path
from django.contrib import admin
from django.contrib.auth.views import LoginView
from drf_yasg.views import get_schema_view
from drf_yasg import openapi

from drf_yasg.generators import OpenAPISchemaGenerator

class BothHttpAndHttpsSchemaGenerator(OpenAPISchemaGenerator):
    def get_schema(self, request=None, public=False):
        schema = super().get_schema(request, public)
        schema.schemes = ["http", "https"]
        return schema
        
        
schema_view = get_schema_view(openapi.Info(
    title="iSkyLIMS API",
    default_version='v0.0.1',
    description="iSkyLIMS API",
    ),
    generator_class=BothHttpAndHttpsSchemaGenerator, # Here
    public=True,
    
)


urlpatterns = [
    path('',include('iSkyLIMS_core.urls')),
    path('background', LoginView.as_view(template_name='iSkyLIMS_core/background.html'), name="background"),
    path('bu-isciii', LoginView.as_view(template_name='iSkyLIMS_core/bu_isciii.html'), name="about-us"),
    path('contact', LoginView.as_view(template_name='iSkyLIMS_core/contact.html'), name="contact"),
    path('admin/', admin.site.urls),
    path('wetlab/', include('iSkyLIMS_wetlab.urls')),
    path('drylab/',include('drylab.urls')),
    path('utils/',include('django_utils.urls')),
    path('clinic/',include ('iSkyLIMS_clinic.urls')),
    path('accounts/', include('django.contrib.auth.urls')),
    
    # REST FRAMEWORK URLS
    path('drylab/api/', include('drylab.api.urls')),
    path('wetlab/api/', include('iSkyLIMS_wetlab.api.urls')),
    path('swagger/', schema_view.with_ui('swagger', cache_timeout=0))
]
