
from rest_framework import serializers
from iSkyLIMS_drylab.models import PipelineExternalDataJobs
from iSkyLIMS_drylab.models import ParameterPipeline
from iSkyLIMS_drylab.models import Service
from iSkyLIMS_wetlab.models import SamplesInProject


class ParameterPipelineSerializer (serializers.ModelSerializer):
    class Meta:
         model = ParameterPipeline
         fields = '__all__'


class ServiceSerializer (serializers.ModelSerializer):
     class Meta:
         model = Service
         fields = '__all__'

class PipelineExternalDataJobsSerializer (serializers.ModelSerializer):
     class Meta:
         model = PipelineExternalDataJobs
         fields = '__all__'

class PipelineExternalDataJobsBSerializer (serializers.ModelSerializer):
    class Meta:
         model = PipelineExternalDataJobs
         fields =  ['serviceRequestNumber','jobState']
