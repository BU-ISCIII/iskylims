from rest_framework import serializers
from django.contrib.auth.models import User
from django_utils.models import Center, Profile

# from iSkyLIMS_drylab.models import ParameterPipeline
from iSkyLIMS_drylab.models import (
    Service,
    Resolution,
    RequestedSamplesInServices,
    Delivery
)

# from iSkyLIMS_wetlab.models import SamplesInProject

"""
class ParameterPipelineSerializer (serializers.ModelSerializer):
    class Meta:
         model = ParameterPipeline
         fields = '__all__'
"""


class CreateDeliveryPostSerializer(serializers.ModelSerializer):
    class Meta:
        model = Delivery
        fields = ["deliveryResolutionID", "deliveryDate",
                  "deliveryNotes", "executionStartDate"]


class UserIDSerializer(serializers.ModelSerializer):

    class Meta:
        model = User
        fields = ["username", "first_name", "last_name", "email"]


"""
class AvailableServicesSerializer(serializers.ModelSerializer):
    #serviceId = CustomAvailableServiceField(many=True, read_only=True)
    class Meta:
        model = AvailableService
        fields = ['availServiceDescription', 'parent', 'level', 'pk']
        #fields = ['availServiceDescription', 'serviceId']
        #attributes :availServiceDescription, :serviceId

    def to_representation(self, instance):
        data = super().to_representation(instance)
        import pdb; pdb.set_trace()
        if not data['serviceId']:
            removed = data.pop("servideId", None)
            import pdb; pdb.set_trace()
            # data['serviceId'] = ""
        return data
"""


class ServiceSerializer(serializers.ModelSerializer):
    serviceFileExt = serializers.StringRelatedField(many=False)
    serviceUserId = UserIDSerializer(many=False)
    serviceAvailableService = serializers.StringRelatedField(many=True)

    class Meta:
        model = Service
        fields = [
            "pk",
            "serviceRequestNumber",
            "serviceStatus",
            "serviceUserId",
            "serviceCreatedOnDate",
            "serviceSeqCenter",
            "serviceAvailableService",
            "serviceFileExt",
            "serviceNotes",
        ]


class UpdateResolutionSerializer(serializers.ModelSerializer):
    resolutionState = serializers.StringRelatedField(many=False)

    class Meta:
        model = Resolution

        fields = ["resolutionNumber", "resolutionState"]

    def update(self, state_obj):
        self.resolutionState = state_obj
        self.save()
        return self


class CustomAvailableServiceField(serializers.RelatedField):
    def to_representation(self, service):
        data = {"availServiceDescription": service.availServiceDescription}
        if service.serviceId:
            data["serviceId"] = service.serviceId
        else:
            data["serviceId"] = None
        return data


class ResolutionSerializer(serializers.ModelSerializer):
    # resolutionServiceID = serializers.StringRelatedField(many = False)
    resolutionPipelines = serializers.StringRelatedField(many=True)
    # esolutionServiceID = ServiceSerializer(many = False)
    # availableServices = AvailableServicesSerializer(many=True)
    availableServices = CustomAvailableServiceField(many=True, read_only=True)

    class Meta:
        model = Resolution
        fields = [
            "pk",
            "resolutionNumber",
            "resolutionFullNumber",
            "resolutionServiceID",
            "resolutionDate",
            "resolutionEstimatedDate",
            "resolutionOnQueuedDate",
            "resolutionOnInProgressDate",
            "resolutionDeliveryDate",
            "resolutionNotes",
            "resolutionPipelines",
            "availableServices",
        ]


class RequestedSamplesInServicesSerializer(serializers.ModelSerializer):
    class Meta:
        model = RequestedSamplesInServices
        fields = ["runName", "projectName", "sampleName", "samplePath"]
