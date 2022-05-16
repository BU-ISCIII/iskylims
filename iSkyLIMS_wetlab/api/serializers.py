from rest_framework import serializers

from iSkyLIMS_core.models import (
    LabRequest,
    Samples,
    SampleProjectsFieldsValue,
    SampleProjectsFields,
    SamplesProjectsTableOptions,
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
            "completedDate",
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

    def update(self, data):
        self.labName = data["labName"]
        self.labNameCoding = data["labNameCoding"]
        self.labUnit = data["labUnit"]
        self.labContactName = data["labContactName"]
        self.labPhone = data["labPhone"]
        self.labEmail = data["labEmail"]
        self.address = data["address"]
        self.apps_name = data["apps_name"]
        self.save()
        return self


class SamplProjetPropertyOptionSerializer(serializers.ModelSerializer):
    class Meta:
        model = SamplesProjectsTableOptions
        fields = ("id", "optionValue")


class SampleProjectFieldSerializer(serializers.ModelSerializer):
    # sampleProjectOptionList = SamplProjetPropertyOptionSerializer(SamplesProjectsTableOptions.ForeignKey(SampleProjectsFields, related_name='optionValue'))
    sampleProjectOptionList = SamplProjetPropertyOptionSerializer(
        source="opt_value_prop", many=True
    )

    class Meta:
        model = SampleProjectsFields
        depth = 1
        fields = [
            "sampleProjectFieldName",
            "sampleProjectFieldUsed",
            "sampleProjectFieldType",
            "sampleProjectOptionList",
            "sampleProjectFieldDescription",
        ]


class SampleFields(object):
    def __init__(self, sample_fields):
        self.sample_fields = sample_fields


class SampleFieldsSerializer(serializers.Serializer):
    sample_fields = serializers.CharField(max_length=800)
