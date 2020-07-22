from django.contrib import admin
from django.urls import path
from rest_framework.routers import DefaultRouter
from api.views import PipelinesViewSet

router = DefaultRouter()
router.register('api', PipelinesViewSet)

urlpatterns =  router.urls
