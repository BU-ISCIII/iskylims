from rest_framework import serializers

from iSkyLIMS_core.models import (
    LabRequest,
    Samples,
    SampleProjectsFieldsValue,
    SampleProjectsFields
)


class CreateSampleSerializer(serializers.ModelSerializer):
    class Meta:
        model = Samples
        fields = [
            "patientCore",
            "sampleName",
            "labRequest",
            "sampleType",
            "species",
            "sampleProject",
            "sampleEntryDate",
            "collectionSampleDate",
            "sampleLocation",
            "onlyRecorded",
            "sampleCodeID",
            "uniqueSampleID",
            "sampleUser",
            "sampleState",
            "completedDate"
        ]


class CreateProjectDataSerializer(serializers.ModelSerializer):
    class Meta:
        model = SampleProjectsFieldsValue
        fields = ["sample_id", "sampleProjecttField_id", "sampleProjectFieldValue"]


class LabRequestSerializer(serializers.ModelSerializer):
    class Meta:
        model = LabRequest
        # fields = "__all__"
        exclude = ["id", "labCity", "apps_name"]

class SampleProjectFieldSerializer(serializers.ModelSerializer):
    class Meta:
        model = SampleProjectsFields
        fields = ["sampleProjectFieldName"]


class SampleFields(object):
    def __init__(self, sample_fields):
        self.sample_fields = sample_fields


class SampleFieldsSerializer(serializers.Serializer):
    sample_fields = serializers.CharField(max_length=800)
