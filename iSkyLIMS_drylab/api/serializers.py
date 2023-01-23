from rest_framework import serializers
from django.contrib.auth.models import User
from iSkyLIMS_drylab.models import (
    Service,
    Resolution,
    RequestedSamplesInServices,
    Delivery
    )


class CreateDeliveryPostSerializer(serializers.ModelSerializer):
    class Meta:
        model = Delivery
        fields = ["deliveryResolutionID", "deliveryDate",
                  "deliveryNotes", "executionStartDate"]


class UserIDSerializer(serializers.ModelSerializer):

    class Meta:
        model = User
        fields = ["username", "first_name", "last_name", "email"]

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

    class Meta:
        model = Service
        fields = [
            "pk",
            "serviceRequestNumber",
            "serviceStatus",
            "serviceUserId",
            "serviceCreatedOnDate",
            "serviceOnDeliveredDate",
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
