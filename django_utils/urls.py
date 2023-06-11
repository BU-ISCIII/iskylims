from django.conf import settings
from django.conf.urls.static import static
from django.urls import path

import django_utils.views

urlpatterns = [
    path("user_creation", django_utils.views.user_creation, name="user_creation"),
    path("user_edit", django_utils.views.user_edit, name="user_edit"),
]
# + static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
