from rest_framework import serializers

from iSkyLIMS_core.models import (
    LabRequest,
    Samples,
    SampleType,
    SampleProjectsFieldsValue,
    SampleProjectsFields,
    SamplesProjectsTableOptions,
)

from iSkyLIMS_wetlab.models import SamplesInProject, Projects, RunProcess
from django.contrib.auth.models import User

# from rest_framework.parsers import JSONParser


class UserIDSerializer(serializers.ModelSerializer):
    class Meta:
        model = User
        fields = ["username"]


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


class CreateSampleTypeSerializer(serializers.ModelSerializer):
    class Meta:
        model = SampleType
        fields = "__all__"


class CreateProjectDataSerializer(serializers.ModelSerializer):
    class Meta:
        model = SampleProjectsFieldsValue
        fields = ["sample_id", "sampleProjecttField_id", "sampleProjectFieldValue"]


class LabRequestSerializer(serializers.ModelSerializer):
    class Meta:
        model = LabRequest
        # fields = ["labContactName", "labPhone", "labEmail"]
        fields = "__all__"

    def update(self, data):
        self.labContactName = data["lab_contact_name"]
        self.labPhone = data["lab_contact_telephone"]
        self.labEmail = data["lab_contact_email"]
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


class ProjectsSerializers(serializers.ModelSerializer):
    class Meta:
        model = Projects
        fields = ["projectName"]


class RunProcessSerializers(serializers.ModelSerializer):
    class Meta:
        model = RunProcess
        fields = ["runName"]


class SampleRunInfoSerializers(serializers.ModelSerializer):
    project_id = ProjectsSerializers(many=False)
    user_id = UserIDSerializer(many=False)
    runProcess_id = RunProcessSerializers(many=False)
    generated_at = serializers.DateTimeField(format="%Y-%m-%d")

    class Meta:
        model = SamplesInProject
        # fields = "__all__"
        exclude = ["sampleInCore"]
