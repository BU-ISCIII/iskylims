
from rest_framework import serializers
from django.contrib.auth.models import User
#from iSkyLIMS_drylab.models import ParameterPipeline
from iSkyLIMS_drylab.models import Service , Resolution, RequestedSamplesInServices, AvailableService
#from iSkyLIMS_wetlab.models import SamplesInProject

'''
class ParameterPipelineSerializer (serializers.ModelSerializer):
    class Meta:
         model = ParameterPipeline
         fields = '__all__'
'''

class UserIDSerializer(serializers.ModelSerializer):
    class Meta:
        model = User
        fields = ['username','first_name','last_name', 'email']

class AvailableServicesSerializer(serializers.ModelSerializer):
    class Meta:
        model = AvailableService
        fields = ['availServiceDescription' , 'serviceId']

class ServiceSerializer (serializers.ModelSerializer):
    serviceFileExt = serializers.StringRelatedField(many=False)
    serviceUserId = UserIDSerializer(many=False)
    serviceAvailableService = serializers.StringRelatedField(many=True)

    class Meta:
        model = Service
        #fields = '__all__'

        fields = ['pk', 'serviceRequestNumber','serviceStatus', 'serviceUserId','serviceCreatedOnDate',
            'serviceSeqCenter', 'serviceAvailableService', 'serviceFileExt' , 'serviceNotes']



class UpdateResolutionSerializer(serializers.ModelSerializer):
    resolutionState = serializers.StringRelatedField(many = False)
    class Meta:
        model = Resolution

        fields = ['resolutionNumber','resolutionState']

    def update (self,state_obj):
        self.resolutionState = state_obj
        self.save()
        return self

'''
class PipelineExternalDataJobsSerializer (serializers.ModelSerializer):
     class Meta:
         model = PipelineExternalDataJobs
         fields = '__all__'

class PipelineExternalDataJobsBSerializer (serializers.ModelSerializer):
    class Meta:
         model = PipelineExternalDataJobs
         fields =  ['serviceRequestNumber','jobState']

'''
class ResolutionSerializer(serializers.ModelSerializer):
    #resolutionServiceID = serializers.StringRelatedField(many = False)
    resolutionPipelines = serializers.StringRelatedField(many = True)
    #resolutionServiceID = ServiceSerializer(many = False)
    availableServices = AvailableServicesSerializer(many=True)
    class Meta:
        model = Resolution
        fields = ['pk', 'resolutionNumber', 'resolutionFullNumber', 'resolutionServiceID' , 'resolutionDate', 'resolutionEstimatedDate' ,
            'resolutionOnQueuedDate' , 'resolutionOnInProgressDate' , 'resolutionDeliveryDate' , 'resolutionNotes', 'resolutionPipelines','availableServices']

class RequestedSamplesInServicesSerializer(serializers.ModelSerializer):

    class Meta:
        model  = RequestedSamplesInServices
        fields = ['runName', 'projectName', 'sampleName' , 'samplePath']
