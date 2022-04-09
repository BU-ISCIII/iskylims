from rest_framework import serializers

from iSkyLIMS_core.models import Samples, SampleProjectsFieldsValue


class CreateSampleSerializer(serializers.ModelSerializer):
    class Meta:
        model = Samples
        fields = [
            "sampleState",
            "patientCore",
            "labRequest",
            "sampleType",
            "sampleUser",
            "sampleCodeID",
            "uniqueSampleID",
            "species",
            "sampleLocation",
            "sampleEntryDate",
            "uniqueSampleID",
            "sampleCodeID",
        ]


class CreateProjectDataSerializer(serializers.ModelSerializer):
    class Meta:
        model = SampleProjectsFieldsValue
        fields = ["sample_id", "sampleProjecttField_id", "sampleProjectFieldValue"]
