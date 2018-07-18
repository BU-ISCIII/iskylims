"""iSkyLIMS URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.11/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.conf.urls import url, include
    2. Add a URL to urlpatterns:  url(r'^blog/', include('blog.urls'))
"""

"""iSkyLIMS URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/2.0/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
# Adapting to Django 2.0 (aligned to comment above)
# "from django.conf.urls import include, url" -->  " from django.urls import include, path"


from django.urls import include, path
from django.contrib import admin
from django.contrib.auth.views import LoginView

urlpatterns = [
    path('',include('iSkyLIMS_home.urls')),
    path('background', LoginView.as_view(template_name='iSkyLIMS_home/background.html'), name="background"),
    path('bu-isciii', LoginView.as_view(template_name='iSkyLIMS_home/bu_isciii.html'), name="about-us"),
    path('faqs', LoginView.as_view(template_name='iSkyLIMS_home/faqs.html'), name="faqs"),
    path('contact', LoginView.as_view(template_name='iSkyLIMS_home/contact.html'), name="contact"),
    path('admin/', admin.site.urls),
    path('wetlab/', include('iSkyLIMS_wetlab.urls')),
    path('drylab/',include('iSkyLIMS_drylab.urls')),
    path('utils/',include('django_utils.urls')),
    path('accounts/', include('django.contrib.auth.urls')),
    path('tutorials/', include('zinnia.urls')),
]
