from rest_framework import serializers
from iSkyLIMS_drylab.models import Pipelines

class PipelinesSerializer (serializers.ModelSerializer):
    class Meta:
         model = Pipelines
         fields = '__all__'
