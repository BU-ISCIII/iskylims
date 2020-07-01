from django.shortcuts import render
from rest_framework import viewsets
from iSkyLIMS_drylab.models import Pipelines
from .serializers import PipelinesSerializer
# Create your views here.
class PipelinesViewSet(viewsets.ModelViewSet):
    serializer_class = PipelinesSerializer
    queryset = Pipelines.objects.all()
