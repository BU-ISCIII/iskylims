from rest_framework import serializers
from django.contrib.auth.models import User
from django_utils.models import Profile
from iSkyLIMS_drylab.models import (
    Service,
    Resolution,
    RequestedSamplesInServices,
    Delivery,
    Pipelines
    )


class CreateDeliveryPostSerializer(serializers.ModelSerializer):
    class Meta:
        model = Delivery
        fields = [
            "deliveryResolutionID",
            "pipelinesInDelivery",
            "deliveryDate",
            "executionStartDate",
            "executionEndDate",
            "permanentUsedSpace",
            "temporaryUsedSpace"
            ]


class UpdateResolutionStateSerializer(serializers.ModelSerializer):
    #resolutionState = serializers.StringRelatedField(many=False)

    class Meta:
        model = Resolution

        fields = ["resolutionNumber", "resolutionState"]

    def update(self, instance, validated_data):
        instance.resolutionNumber = validated_data["resolutionNumber"]
        instance.resolutionState = validated_data["resolutionState"]
        instance.save()
        return instance

class UpdateServiceStateSerializer(serializers.ModelSerializer):

    class Meta:
        model = Service

        fields = ["serviceRequestNumber", "serviceStatus"]

class ProfileUserSerializer(serializers.ModelSerializer):
    profileClassificationArea= serializers.StringRelatedField()
    profileCenter = serializers.StringRelatedField()

    class Meta:
        model = Profile
        fields = ["profileClassificationArea", "profileCenter"]

class UserIDSerializer(serializers.ModelSerializer):
    profile = ProfileUserSerializer(many=False)

    class Meta:
        model = User
        fields = ["username", "first_name", "last_name", "email", "profile"]

class CustomAvailableServiceField(serializers.RelatedField):
    def to_representation(self, service):
        data = {"availServiceDescription": service.availServiceDescription}
        if service.serviceId:
            data["serviceId"] = service.serviceId
        else:
            data["serviceId"] = None
        return data

class PipelinesSerializer(serializers.ModelSerializer):

    class Meta:
        model=Pipelines
        fields = [
            "pipelineName",
            "PipelineVersion"
        ]

class DeliverySerializer(serializers.ModelSerializer):
    deliveryResolutionID = serializers.StringRelatedField(many=False)
    pipelinesInDelivery = PipelinesSerializer(many=True)

    class Meta:
        model= Delivery
        fields = [
            "deliveryResolutionID",
            "pipelinesInDelivery",
            "deliveryDate",
            "executionStartDate",
            "executionEndDate",
            "permanentUsedSpace",
            "temporaryUsedSpace"
        ]

class ResolutionSerializer(serializers.ModelSerializer):
    resolutionPipelines = serializers.StringRelatedField(many=True)
    availableServices = CustomAvailableServiceField(many=True, read_only=True)
    resolutionServiceID = serializers.StringRelatedField(many=False)
    delivery = DeliverySerializer(many=True)

    class Meta:
        model = Resolution
        fields = [
            "resolutionNumber",
            "resolutionFullNumber",
            "resolutionDate",
            "resolutionEstimatedDate",
            "resolutionServiceID",
            "resolutionOnQueuedDate",
            "resolutionOnInProgressDate",
            "resolutionDeliveryDate",
            "resolutionNotes",
            "resolutionPipelines",
            "availableServices",
            "delivery"
            ]

class RequestedSamplesInServicesSerializer(serializers.ModelSerializer):
    class Meta:
        model = RequestedSamplesInServices
        fields = ["runName", "projectName", "sampleName", "samplePath"]


class ServiceListSerializer(serializers.ModelSerializer):

    class Meta:
        model = Service
        fields = [
            "serviceRequestNumber",
            "serviceStatus",
            "serviceCreatedOnDate",
            "serviceOnDeliveredDate",
            ]

class ServiceSerializer(serializers.ModelSerializer):
    serviceFileExt = serializers.StringRelatedField(many=False)
    serviceUserId = UserIDSerializer(many=False)
    serviceAvailableService = serializers.StringRelatedField(many=True)
    resolutions = ResolutionSerializer(source="filtered_resolutions", many=True)
    samples = RequestedSamplesInServicesSerializer(many=True)

    class Meta:
        model = Service
        fields = [
            "serviceRequestNumber",
            "serviceStatus",
            "serviceUserId",
            "serviceCreatedOnDate",
            "serviceOnDeliveredDate",
            "serviceSeqCenter",
            "serviceAvailableService",
            "serviceFileExt",
            "serviceNotes",
            "resolutions",
            "samples"
            ]
