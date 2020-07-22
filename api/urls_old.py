from django.contrib import admin
from django.urls import path
from rest_framework.routers import DefaultRouter
from iSkyLIMS_drylab.api.views import PipelineExternalDataJobsViewSet

router = DefaultRouter()
router.register('apidrylab', PipelineExternalDataJobsViewSet)

urlpatterns =  router.urls
