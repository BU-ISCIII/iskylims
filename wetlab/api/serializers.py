from django.contrib.auth.models import User
from rest_framework import serializers
from collections import OrderedDict

import core.models
import wetlab.models


class UserIDSerializer(serializers.ModelSerializer):
    class Meta:
        model = User
        fields = ["username"]


class CreateProjectDataSerializer(serializers.ModelSerializer):
    class Meta:
        model = core.models.SampleProjectsFieldsValue
        fields = ["sample_id", "sample_project_field_id", "sample_project_field_value"]


class SampleProjectOptionSerializer(serializers.ModelSerializer):
    class Meta:
        model = core.models.SamplesProjectsTableOptions
        fields = ("id", "option_value")


class SampleProjectFieldSerializer(serializers.ModelSerializer):
    sample_project_option_list = SampleProjectOptionSerializer(
        source="opt_value_prop", many=True
    )

    class Meta:
        model = core.models.SampleProjectsFields
        # depth = 1
        fields = [
            "sample_project_field_name",
            "sample_project_field_used",
            "sample_project_field_type",
            "sample_project_option_list",
            "sample_project_field_description",
        ]


class ProjectsSerializer(serializers.ModelSerializer):
    sample_project_fields = SampleProjectFieldSerializer(
        source="project_fields_options", many=True
    )

    class Meta:
        model = core.models.SampleProjects
        fields = [
            "sample_project_name",
            "sample_project_manager",
            "sample_project_contact",
            "sample_project_fields",
        ]


class ProjectValuesSerializers(serializers.ModelSerializer):
    sample_project_field_id = serializers.StringRelatedField(many=False)

    def to_representation(self, instance):
        data = super(ProjectValuesSerializers, self).to_representation(instance)
        data_update = dict()
        # Convert dict according to db table fields, to id and values inside the table.
        # WARNING: This will have to be change if db fields change in models.py
        data_update[data["sample_project_field_id"]] = data[
            "sample_project_field_value"
        ]

        return data_update

    class Meta:
        model = core.models.SampleProjectsFieldsValue
        fields = ["sample_project_field_id", "sample_project_field_value"]


class SampleProjectParameterSerializer(serializers.ModelSerializer):
    sample_name = serializers.CharField(source="sample_id.sample_name")

    def to_representation(self, instance):
        data = super(SampleProjectParameterSerializer, self).to_representation(instance)
        data_update = OrderedDict()
        for key in self.fields:
            # change parameter label name
            if key == "sample_project_field_value":
                data_update[self.context["parameter"]] = data[key]
            else:
                # Change id to label for api rest output
                data_update[self.fields[key].label] = data[key]

        return data_update

    class Meta:
        model = core.models.SampleProjectsFieldsValue
        fields = ["sample_name", "sample_project_field_value"]


class SampleSerializer(serializers.ModelSerializer):
    lab_request = serializers.StringRelatedField(many=False, label="Laboratory")
    sample_project = serializers.StringRelatedField(many=False, label="Sample Project")
    project_values = ProjectValuesSerializers(many=True)
    sample_type = serializers.StringRelatedField(many=False, label="Sample type")
    species = serializers.StringRelatedField(many=False, label="Species")

    def to_representation(self, instance):
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
        model = core.models.Samples
        fields = [
            "sample_name",
            "lab_request",
            "sample_type",
            "species",
            "sample_project",
            "sample_entry_date",
            "collection_sample_date",
            "project_values",
        ]


class CreateSampleSerializer(serializers.ModelSerializer):
    class Meta:
        model = core.models.Samples
        fields = [
            "patient_core",
            "sample_name",
            "lab_request",
            "sample_type",
            "species",
            "sample_project",
            "sample_entry_date",
            "collection_sample_date",
            "sample_location",
            "only_recorded",
            "sample_code_id",
            "unique_sample_id",
            "sample_user",
            "sample_state",
            "completed_date",
        ]


class SampleParameterSerializer(serializers.ModelSerializer):
    req_param = serializers.SerializerMethodField("parameter_value")

    def to_representation(self, instance):
        data = super(SampleParameterSerializer, self).to_representation(instance)
        data_update = OrderedDict()
        for key in self.fields:
            # change parameter label name
            if key == "req_param":
                data_update[self.context["parameter"]] = data[key]
            else:
                # Change id to label for api rest output
                data_update[self.fields[key].label] = data[key]

        return data_update

    def parameter_value(self, obj):
        param = self.context["parameter"]
        if param:
            req_parameter = "obj." + param
            value = eval(req_parameter)
            if "date" in param.lower():
                if value is not None:
                    return value.strftime("%Y-%m-%d")
                else:
                    return value
            else:
                return value.__str__()
        return False

    class Meta:
        model = core.models.Samples
        fields = ["sample_name", "req_param"]


class CreateSampleTypeSerializer(serializers.ModelSerializer):
    class Meta:
        model = core.models.SampleType
        fields = "__all__"


class LabRequestSerializer(serializers.ModelSerializer):
    # lab_request = serializers.StringRelatedField(many=False, label="Laboratory")
    lab_city = serializers.StringRelatedField(many=False)
    
    def to_representation(self, instance):
        data = super(LabRequestSerializer, self).to_representation(instance)
        data_update = dict()
        for key in self.fields:

            data_update[self.fields[key].label] = data[key]

        return data_update

    class Meta:
        model = core.models.LabRequest
        fields = [
            "lab_name",
            "lab_name_coding",
            "lab_unit",
            "lab_contact_name",
            "lab_phone",
            "lab_email",
            "address",
            "lab_city",
        ]

    

    def update(self, data):
        self.labContactName = data["lab_contact_name"]
        self.labPhone = data["lab_contact_telephone"]
        self.labEmail = data["lab_contact_email"]
        self.save()
        return self


class SampleFields(object):
    def __init__(self, sample_fields):
        self.sample_fields = sample_fields


class SampleFieldsSerializer(serializers.Serializer):
    sample_fields = serializers.CharField(max_length=800)


class RunProcessSerializers(serializers.ModelSerializer):
    class Meta:
        model = wetlab.models.RunProcess
        fields = ["run_name"]


class SampleRunInfoSerializers(serializers.ModelSerializer):
    project_id = serializers.StringRelatedField(many=False)
    user_id = UserIDSerializer(many=False)
    run_process_id = RunProcessSerializers(many=False)
    generated_at = serializers.DateTimeField(format="%Y-%m-%d")

    class Meta:
        model = wetlab.models.SamplesInProject
        # fields = "__all__"
        exclude = ["sample_in_core"]
