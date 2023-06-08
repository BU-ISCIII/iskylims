from collections import OrderedDict

from django.contrib.auth.models import User
from rest_framework import serializers

from django_utils.models import Profile
from drylab.models import (Delivery, Pipelines, RequestedSamplesInServices,
                           Resolution, Service)


class CreateDeliveryPostSerializer(serializers.ModelSerializer):
    class Meta:
        model = Delivery
        fields = [
            "delivery_resolution_id",
            "pipelines_in_delivery",
            "execution_start_date",
            "execution_end_date",
            "permanent_used_space",
            "temporary_used_space",
            "delivery_notes",
        ]


class UpdateResolutionStateSerializer(serializers.ModelSerializer):
    class Meta:
        model = Resolution

        fields = [
            "resolution_number",
            "resolution_in_progress_date",
            "resolution_delivery_date",
            "resolution_state",
        ]

    def update(self, instance, validated_data):
        instance.resolution_number = validated_data["resolution_number"]
        instance.resolution_state = validated_data["resolution_state"]
        if "resolution_inprogress_date" in validated_data:
            instance.resolution_in_progress_date = validated_data[
                "resolution_inprogress_date"
            ]
        if "resolution_delivery_date" in validated_data:
            instance.resolution_delivery_date = validated_data["resolution_delivery_date"]
        instance.save()
        return instance


class UpdateServiceStateSerializer(serializers.ModelSerializer):
    class Meta:
        model = Service

        fields = ["service_request_number", "service_delivered_date", "service_state"]


class ProfileUserSerializer(serializers.ModelSerializer):
    profile_classification_area = serializers.StringRelatedField()
    profile_center = serializers.StringRelatedField()

    class Meta:
        model = Profile
        fields = ["profile_classification_area", "profile_center"]


class UserIDSerializer(serializers.ModelSerializer):
    profile = ProfileUserSerializer(many=False)

    class Meta:
        model = User
        fields = ["username", "first_name", "last_name", "email", "profile"]


class CustomAvailableServiceField(serializers.RelatedField):
    def to_representation(self, service):
        data = {"avail_service_description": service.avail_service_description}
        if service.service_id:
            data["service_id"] = service.service_id
        else:
            data["service_id"] = None
        return data


class PipelinesSerializer(serializers.ModelSerializer):
    class Meta:
        model = Pipelines
        fields = ["pipeline_name", "pipeline_version"]


class DeliverySerializer(serializers.ModelSerializer):
    delivery_resolution_id = serializers.StringRelatedField(many=False)
    pipelines_in_delivery = PipelinesSerializer(many=True)

    class Meta:
        model = Delivery
        fields = [
            "delivery_resolution_id",
            "pipelines_in_delivery",
            "execution_start_date",
            "execution_end_date",
            "permanent_used_space",
            "temporary_used_space",
            "delivery_notes",
        ]


class ResolutionSerializer(serializers.ModelSerializer):
    resolution_state = serializers.StringRelatedField(many=False)
    resolution_pipelines = serializers.StringRelatedField(many=True)
    available_services = CustomAvailableServiceField(many=True, read_only=True)
    resolution_service_id = serializers.StringRelatedField(many=False)
    delivery = DeliverySerializer(many=True)

    class Meta:
        model = Resolution
        fields = [
            "resolution_number",
            "resolution_full_number",
            "resolution_state",
            "resolution_date",
            "resolution_estimated_date",
            "resolution_service_id",
            "resolution_queued_date",
            "resolution_in_progress_date",
            "resolution_delivery_date",
            "resolution_notes",
            "resolution_pipelines",
            "available_services",
            "delivery",
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
            "service_state",
            "service_created_date",
            "service_delivered_date",
        ]


class ServiceSerializer(serializers.ModelSerializer):
    service_user_id = UserIDSerializer(many=False)
    service_available_service = serializers.StringRelatedField(many=True)
    resolutions = ResolutionSerializer(source="filtered_resolutions", many=True)
    samples = RequestedSamplesInServicesSerializer(many=True)
    sample_name = serializers.CharField(source="sample_id.sampleName")

    def to_representation(self, instance):
        output_label = self.context["output_label"]
        data = super(ServiceSerializer, self).to_representation(instance)
        data_update = OrderedDict()
        for key in self.fields:
            if output_label:
                # output both label and value
                label_value = {}
                label_value["label"] = self.fields[key].label
                label_value["value"] = data[key]
                data_update[key] = label_value
            else:
                # normal key output using field name
                data_update[key] = data[key]
        return data_update

    class Meta:
        model = Service
        fields = [
            "service_request_number",
            "service_state",
            "service_user_id",
            "service_created_date",
            "service_delivered_date",
            "service_center",
            "service_available_service",
            "service_notes",
            "resolutions",
            "samples",
        ]
