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
            "delivery_resolutionID",
            "pipelines_in_delivery",
            "delivery_date",
            "executionStartDate",
            "execution_end_date",
            "permanent_used_space",
            "temporary_used_space",
            "deliveryNotes"
            ]


class UpdateResolutionStateSerializer(serializers.ModelSerializer):
    #resolution_state = serializers.StringRelatedField(many=False)

    class Meta:
        model = Resolution

        fields = ["resolutionNumber", "resolutionOnInProgressDate", "resolutionDeliveryDate","resolution_state"]

    def update(self, instance, validated_data):
        instance.resolutionNumber = validated_data["resolutionNumber"]
        instance.resolution_state = validated_data["resolution_state"]
        if "resolutionOnInProgressDate" in validated_data:
            instance.resolutionOnInProgressDate = validated_data["resolutionOnInProgressDate"]
        if "resolutionDeliveryDate" in validated_data:
            instance.resolutionDeliveryDate = validated_data["resolutionDeliveryDate"]
        instance.save()
        return instance

class UpdateServiceStateSerializer(serializers.ModelSerializer):

    class Meta:
        model = Service

        fields = ["service_request_number", "service_delivered_date", "serviceStatus"]

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
        data = {"avail_service_description": service.avail_service_description}
        if service.serviceId:
            data["serviceId"] = service.serviceId
        else:
            data["serviceId"] = None
        return data

class PipelinesSerializer(serializers.ModelSerializer):

    class Meta:
        model=Pipelines
        fields = [
            "pipeline_name",
            "pipeline_version"
        ]

class DeliverySerializer(serializers.ModelSerializer):
    delivery_resolutionID = serializers.StringRelatedField(many=False)
    pipelines_in_delivery = PipelinesSerializer(many=True)

    class Meta:
        model= Delivery
        fields = [
            "delivery_resolutionID",
            "pipelines_in_delivery",
            "delivery_date",
            "executionStartDate",
            "execution_end_date",
            "permanent_used_space",
            "temporary_used_space",
            "deliveryNotes"
    ]

class ResolutionSerializer(serializers.ModelSerializer):
    resolution_state=serializers.StringRelatedField(many=False)
    resolution_pipelines = serializers.StringRelatedField(many=True)
    available_services = CustomAvailableServiceField(many=True, read_only=True)
    resolution_serviceID = serializers.StringRelatedField(many=False)
    delivery = DeliverySerializer(many=True)

    class Meta:
        model = Resolution
        fields = [
            "resolutionNumber",
            "resolution_full_number",
            "resolution_state",
            "resolutionDate",
            "resolutionEstimatedDate",
            "resolution_serviceID",
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
        fields = ["run_name", "project_name", "sample_name", "sample_path"]


class ServiceListSerializer(serializers.ModelSerializer):

    class Meta:
        model = Service
        fields = [
            "service_request_number",
            "serviceStatus",
            "service_created_date",
            "service_delivered_date",
            ]

class ServiceSerializer(serializers.ModelSerializer):
    serviceFileExt = serializers.StringRelatedField(many=False)
    serviceUserId = UserIDSerializer(many=False)
    service_available_service = serializers.StringRelatedField(many=True)
    resolutions = ResolutionSerializer(source="filtered_resolutions", many=True)
    samples = RequestedSamplesInServicesSerializer(many=True)

    class Meta:
        model = Service
        fields = [
            "service_request_number",
            "serviceStatus",
            "serviceUserId",
            "service_created_date",
            "service_delivered_date",
            "service_seq_center",
            "service_available_service",
            "serviceFileExt",
            "service_notes",
            "resolutions",
            "samples"
            ]
