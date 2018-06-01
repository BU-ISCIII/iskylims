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
    path('',include('iSkyLIMS-home.urls')),
    path('background', LoginView.as_view(template_name='iSkyLIMS-home/background.html'), name="background"),
    path('bu-isciii', LoginView.as_view(template_name='iSkyLIMS-home/bu_isciii.html'), name="about-us"),
    path('faqs', LoginView.as_view(template_name='iSkyLIMS-home/faqs.html'), name="faqs"),
    path('contact', LoginView.as_view(template_name='iSkyLIMS-home/contact.html'), name="contact"),
    path('admin/', admin.site.urls),
    path('iSkyLIMS-wetlab/', include('iSkyLIMS-wetlab.urls')),
    path('iSkyLIMS-drylab/',include('iSkyLIMS-drylab.urls')),
    path('django-utils/',include('django-utils.urls')),
    path('accounts/', include('django.contrib.auth.urls')),
    path('tutorials/', include('zinnia.urls')),
]
