from rest_framework import serializers
from collections import OrderedDict

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

class ProjectValuesSerializers(serializers.ModelSerializer):
    sampleProjecttField_id=serializers.StringRelatedField(many=False)

    def to_representation(self,instance):
        data = super(ProjectValuesSerializers, self).to_representation(instance)
        data_update = dict()
        # Convert dict according to db table fields, to id and values inside the table.
        # WARNING: This will have to be change if db fields change in models.py
        data_update[data["sampleProjecttField_id"]] = data["sampleProjectFieldValue"]

        return data_update

    class Meta:
        model = SampleProjectsFieldsValue
        fields = ["sampleProjecttField_id", "sampleProjectFieldValue"]

class SampleSerializer(serializers.ModelSerializer):
    labRequest = serializers.StringRelatedField(many=False, label="Laboratory")
    sampleProject = serializers.StringRelatedField(many=False, label="Sample Project")
    sampleEntryDate = serializers.DateTimeField(format="%Y-%m-%d", label="Recorded sample date")
    collectionSampleDate = serializers.DateTimeField(format="%Y-%m-%d", label="Collection sample date")
    project_values = ProjectValuesSerializers(many=True)

    def to_representation(self,instance):
        data = super(SampleSerializer, self).to_representation(instance)
        data_update = OrderedDict()
        field_values = dict()

        for key in self.fields:

            # Append all dictionaries into one.
            # origin: [{'id1':value1}{id1:value2}]
            # dest: {'id1':value1, 'id2:value2'}
            if key == "project_values":

                for item in data["project_values"]:
                    field_values.update(item)

                data_update[self.fields[key].label] = field_values
            else:
                # Change id to label for api rest output
                data_update[self.fields[key].label] = data[key]

        return data_update

    class Meta:
        model = Samples
        fields = [
            "sampleName",
            "labRequest",
            "sampleProject",
            "sampleEntryDate",
            "collectionSampleDate",
            "project_values",
        ]


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

