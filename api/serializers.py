from rest_framework import serializers
from iSkyLIMS_drylab.models import PipelineExternalDataJobs

class PipelineExternalDataJobsSerializer (serializers.ModelSerializer):
    class Meta:
         model = PipelineExternalDataJobs
         fields = '__all__'

class PipelineExternalDataJobsBSerializer (serializers.ModelSerializer):
    class Meta:
         model = PipelineExternalDataJobs
         fields =  ['serviceRequestNumber','jobState']
