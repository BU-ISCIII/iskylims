from django.http import QueryDict
from drf_yasg import openapi
from drf_yasg.utils import swagger_auto_schema
from rest_framework import status
from rest_framework.authentication import BasicAuthentication, SessionAuthentication
from rest_framework.decorators import (
    api_view,
    authentication_classes,
    permission_classes,
)
from rest_framework.permissions import IsAuthenticated
from rest_framework.response import Response

import core.models
import wetlab.api.serializers
import wetlab.api.utils.lab
import wetlab.api.utils.sample
import wetlab.models

sample_project_fields = openapi.Parameter(
    "project",
    openapi.IN_QUERY,
    description="Project name to fetch the sample project fields defined. Example Relecov",
    type=openapi.TYPE_STRING,
)

laboratory = openapi.Parameter(
    "laboratory",
    openapi.IN_QUERY,
    description="Laboratory name form to fetch contact information. Example Harlem Hospital Center",
    type=openapi.TYPE_STRING,
)

sample_in_run = openapi.Parameter(
    "samples",
    openapi.IN_QUERY,
    description="Sample name list to fetch run information",
    type=openapi.TYPE_STRING,
)

sample_state = openapi.Parameter(
    "sample_state",
    openapi.IN_QUERY,
    description="Filter result to the sample which are in this state",
    type=openapi.TYPE_STRING,
)

start_date = openapi.Parameter(
    "start_date",
    openapi.IN_QUERY,
    description="Start date from starting collecting samples",
    type=openapi.TYPE_STRING,
)

end_date = openapi.Parameter(
    "end_date",
    openapi.IN_QUERY,
    description="Start date from starting collecting samples",
    type=openapi.TYPE_STRING,
)

sample_parameter = openapi.Parameter(
    "sample_parameter",
    openapi.IN_QUERY,
    description="Get the samples grouped by the parameter",
    type=openapi.TYPE_STRING,
)

region = openapi.Parameter(
    "region",
    openapi.IN_QUERY,
    description="Filter the samples for the selected region",
    type=openapi.TYPE_STRING,
)

sample_list = openapi.Parameter(
    "samples",
    openapi.IN_QUERY,
    description="List of samples to collect information",
    type=openapi.TYPE_STRING,
)

sample_information = openapi.Parameter(
    "sample",
    openapi.IN_QUERY,
    description="Fecthing information from sample",
    type=openapi.TYPE_STRING,
)

sample_parameter = openapi.Parameter(
    "parameter",
    openapi.IN_QUERY,
    description="Fecthing only parameter information from sample",
    type=openapi.TYPE_STRING,
)


sample_project_name = openapi.Parameter(
    "sample_project_name",
    openapi.IN_QUERY,
    description="Sample Project Name",
    type=openapi.TYPE_STRING,
)

project_fields = openapi.Parameter(
    "project_field",
    openapi.IN_QUERY,
    description="Project fields to make the query. Maximum number of fiels are 2",
    type=openapi.TYPE_STRING,
)

sample_project_field = openapi.Parameter(
    "sample_project_field",
    openapi.IN_QUERY,
    description="Name of the sample project Field. Requires project_name",
    type=openapi.TYPE_STRING,
)


@swagger_auto_schema(
    method="post",
    operation_description="Create a new sample Data in iSkyLIMS",
    request_body=openapi.Schema(
        type=openapi.TYPE_OBJECT,
        properties={
            "sample": openapi.Schema(
                type=openapi.TYPE_STRING, description="Sample name"
            ),
            "sample_state": openapi.Schema(
                type=openapi.TYPE_STRING, description="Sample state"
            ),
            "patient_core": openapi.Schema(
                type=openapi.TYPE_STRING, description="Code assigned to the patient"
            ),
            "lab_request": openapi.Schema(
                type=openapi.TYPE_STRING,
                description="Laboratory that request the sample",
            ),
            "sample_type": openapi.Schema(
                type=openapi.TYPE_STRING, description="Type of the sample"
            ),
            "species": openapi.Schema(
                type=openapi.TYPE_STRING,
                description="Specie that the sample belongs to",
            ),
            "sample_location": openapi.Schema(
                type=openapi.TYPE_STRING,
                description="Location where the sample is stored",
            ),
            "sample_entry_date": openapi.Schema(
                type=openapi.TYPE_STRING,
                description="Date when sample is received in the lab",
            ),
            "sample_collection_date": openapi.Schema(
                type=openapi.TYPE_STRING,
                description="Date when the sample is collected from the specimen",
            ),
            "only_recorded": openapi.Schema(
                type=openapi.TYPE_STRING,
                description="Select if sample is just recorded or if DNA/RNA manipulation will be done in the lab ",
            ),
            "project": openapi.Schema(
                type=openapi.TYPE_STRING, description="Project name"
            ),
        },
    ),
    responses={
        201: "Successful create information",
        400: "Bad Request",
        500: "Internal Server Error",
    },
)
@authentication_classes([SessionAuthentication, BasicAuthentication])
@api_view(["POST"])
@permission_classes([IsAuthenticated])
def create_sample_data(request):
    if request.method == "POST":
        data = request.data
        if isinstance(data, QueryDict):
            data = data.dict()
        if "sample_name" not in data or "sample_project" not in data:
            return Response(status=status.HTTP_400_BAD_REQUEST)
        if core.models.Samples.objects.filter(
            sample_name__iexact=data["sample_name"]
        ).exists():
            error = {"ERROR": "sample already defined"}
            return Response(error, status=status.HTTP_400_BAD_REQUEST)
        split_data = wetlab.api.utils.sample.split_sample_data(data)
        if not isinstance(split_data, dict):
            return Response(split_data, status=status.HTTP_400_BAD_REQUEST)
        apps_name = __package__.split(".")[0]
        inst_req_sample = wetlab.api.utils.sample.include_instances_in_sample(
            split_data["s_data"], split_data["lab_data"], apps_name
        )
        if not isinstance(inst_req_sample, dict):
            return Response(inst_req_sample, status=status.HTTP_400_BAD_REQUEST)
        split_data["s_data"] = inst_req_sample
        split_data["s_data"]["sample_user"] = request.user.pk
        # Adding coding for sample
        split_data["s_data"].update(
            wetlab.api.utils.sample.include_coding(
                request.user.username, split_data["s_data"]["sample_name"]
            )
        )
        sample_serializer = wetlab.api.serializers.CreateSampleSerializer(
            data=split_data["s_data"]
        )
        if not sample_serializer.is_valid():
            return Response(
                sample_serializer.errors, status=status.HTTP_400_BAD_REQUEST
            )
        new_sample_id = sample_serializer.save().get_sample_id()
        for d_field in split_data["p_data"]:
            d_field["sample_id"] = new_sample_id
            s_project_serializer = wetlab.api.serializers.CreateProjectDataSerializer(
                data=d_field
            )
            if not s_project_serializer.is_valid():
                return Response(
                    s_project_serializer.errors, status=status.HTTP_400_BAD_REQUEST
                )
            s_project_serializer.save()

        return Response("Successful upload information", status=status.HTTP_201_CREATED)
    return Response(status=status.HTTP_400_BAD_REQUEST)


@swagger_auto_schema(
    method="get",
    operation_description="Get the stored Run information available in iSkyLIMS for the list of samples",
    manual_parameters=[sample_in_run],
)
@api_view(["GET"])
def fetch_run_information(request):
    if "samples" in request.GET:
        samples = request.GET["samples"]
        # sample_run_info = get_run_info_for_sample(apps_name, samples)
        s_list = samples.strip().split(",")
        s_data = []
        for sample in s_list:
            sample = sample.strip()
            if wetlab.models.SamplesInProject.objects.filter(
                sample_name__iexact=sample
            ).exists():
                s_found_objs = wetlab.models.SamplesInProject.objects.filter(
                    sample_name__iexact=sample
                )
                for s_found_obj in s_found_objs:
                    s_data.append(
                        wetlab.api.serializers.SampleRunInfoSerializers(
                            s_found_obj, many=False
                        ).data
                    )
            else:
                s_data.append({"sample_name": sample, "Run data": "Not found"})
        return Response(s_data, status=status.HTTP_200_OK)
    return Response(status=status.HTTP_400_BAD_REQUEST)


@swagger_auto_schema(
    method="get",
    manual_parameters=[sample_information, sample_project_name, sample_parameter],
)
@api_view(["GET"])
def fetch_sample_information(request):
    sample_data = {}
    if "sample" in request.GET:
        sample = request.GET["sample"]
        if not core.modesl.Samples.objects.filter(sample_name__iexact=sample).exists():
            return Response(status=status.HTTP_204_NO_CONTENT)
        sample_obj = core.models.Samples.objects.filter(
            sample_name__iexact=sample
        ).last()
        sample_data = wetlab.api.serializers.SampleSerializer(
            sample_obj, many=False
        ).data
    else:
        if "sample_project_name" in request.GET:
            project_name = request.GET["sample_project_name"]
            if "parameter" not in request.GET:
                return Response(status=status.HTTP_400_BAD_REQUEST)
            param = request.GET["parameter"]
            if core.models.SampleProjects.objects.filter(
                sample_project_name__iexact=project_name
            ).exists():
                project_obj = core.models.SampleProjects.objects.filter(
                    sample_project_name__iexact=project_name
                ).last()

                if not core.models.SampleProjectsFields.objects.filter(
                    sample_projects_id=project_obj,
                    sample_project_field_name__iexact=param,
                ).exists():
                    return Response(status=status.HTTP_400_BAD_REQUEST)
                s_p_field_obj = core.models.SampleProjectsFields.objects.filter(
                    sample_projects_id=project_obj,
                    sample_project_field_name__iexact=param,
                ).last()
                sample_obj = core.models.SampleProjectsFieldsValue.objects.filter(
                    sample_project_field_id=s_p_field_obj
                )
                sample_data = wetlab.api.serializers.SampleProjectParameterSerializer(
                    sample_obj, many=True, context={"parameter": param}
                ).data
                return Response(sample_data, status=status.HTTP_200_OK)
            else:
                return Response(status=status.HTTP_400_BAD_REQUEST)
        if "parameter" in request.GET:
            param = request.GET["parameter"]
            sample_obj = core.models.Samples.objects.all()
            # check if parameter exists
            try:
                eval("sample_obj[0]." + param)
            except AttributeError:
                return Response(status=status.HTTP_400_BAD_REQUEST)
            sample_data = wetlab.api.serializers.SampleParameterSerializer(
                sample_obj, many=True, context={"parameter": param}
            ).data
        else:
            sample_obj = core.models.Samples.objects.prefetch_related("project_values")
            sample_data = wetlab.api.serializers.SampleSerializer(
                sample_obj, many=True
            ).data
    return Response(sample_data, status=status.HTTP_200_OK)


@swagger_auto_schema(
    method="get",
    operation_description="Send request to gen the fields that are required when storing a new sample using the API",
)
@api_view(["GET"])
def sample_fields(request):
    apps_name = __package__.split(".")[0]

    sample_fields = wetlab.api.utils.sample.get_sample_fields(apps_name)

    if "ERROR" in sample_fields:
        return Response(sample_fields, status=status.HTTP_202_ACCEPTED)
    else:
        return Response(sample_fields, status=status.HTTP_200_OK)


@swagger_auto_schema(
    method="get",
    operation_description="Use this request to get the field' s names that are required for the sample project",
    manual_parameters=[sample_project_fields],
)
@api_view(["GET"])
def sample_project_fields(request):
    if "project" in request.GET:
        project = request.GET["project"].strip()
        if core.models.SampleProjects.objects.filter(
            sample_project_name__iexact=project
        ).exists():
            s_project_obj = core.models.SampleProjects.objects.filter(
                sample_project_name__iexact=project
            ).last()
            s_project_field_objs = core.models.SampleProjectsFields.objects.filter(
                sample_projects_id=s_project_obj
            )
            s_project_serializer = wetlab.api.serializers.SampleProjectFieldSerializer(
                s_project_field_objs, many=True
            )
            return Response(s_project_serializer.data, status=status.HTTP_200_OK)
        else:
            return Response(status=status.HTTP_204_NO_CONTENT)
    return Response(status=status.HTTP_400_BAD_REQUEST)


@swagger_auto_schema(method="get", manual_parameters=[laboratory])
@api_view(["GET"])
def get_lab_information_contact(request):
    if "laboratory" in request.GET:
        lab_name = request.GET["laboratory"].strip()
        if core.models.LabRequest.objects.filter(lab_name__iexact=lab_name).exists():
            lab_req_obj = core.models.LabRequest.objects.filter(
                lab_name__iexact=lab_name
            ).last()
            lab_req_serializer = core.models.LabRequestSerializer(
                lab_req_obj, many=False
            )
            return Response(lab_req_serializer.data, status=status.HTTP_200_OK)
        else:
            return Response(status=status.HTTP_204_NO_CONTENT)
    return Response(status=status.HTTP_400_BAD_REQUEST)


@swagger_auto_schema(
    method="get",
    operation_description="",
    manual_parameters=[
        sample_list,
        sample_state,
        start_date,
        end_date,
        region,
        laboratory,
        sample_project_name,
        sample_project_field,
    ],
)
@api_view(["GET"])
def summarize_data_information(request):
    summarize_data = wetlab.api.utils.sample.summarize_samples(request.GET)
    if "ERROR" in summarize_data:
        return Response(status=status.HTTP_400_BAD_REQUEST)
    if len(summarize_data) == 0:
        return Response(status=status.HTTP_204_NO_CONTENT)
    return Response(summarize_data, status=status.HTTP_200_OK)


@swagger_auto_schema(
    method="get",
    operation_description="",
    manual_parameters=[
        sample_project_name,
        project_fields,
    ],
)
@api_view(["GET"])
def statistic_information(request):
    statistics_data = wetlab.api.utils.sample.collect_statistics_information(
        request.GET
    )
    if "ERROR" in statistics_data:
        return Response(status=status.HTTP_400_BAD_REQUEST)
    if len(statistics_data) == 0:
        return Response(status=status.HTTP_204_NO_CONTENT)
    return Response(statistics_data, status=status.HTTP_200_OK)


@swagger_auto_schema(
    method="put",
    operation_description="Update laboratory contact information",
    request_body=openapi.Schema(
        type=openapi.TYPE_OBJECT,
        properties={
            "lab_name": openapi.Schema(
                type=openapi.TYPE_STRING,
                description="Name of the Laboratory",
            ),
            "lab_contact_name": openapi.Schema(
                description="Name for laboratory contact",
                type=openapi.TYPE_STRING,
            ),
            "lab_contact_telephone": openapi.Schema(
                description="Phone number of contact",
                type=openapi.TYPE_STRING,
            ),
            "lab_contact_email": openapi.Schema(
                description="Contact email",
                type=openapi.TYPE_STRING,
            ),
        },
    ),
    responses={
        201: "Successful create information",
        204: "Laboratory not defined",
        400: "Bad Request",
        500: "Internal Server Error",
    },
)
@authentication_classes([SessionAuthentication, BasicAuthentication])
@api_view(["PUT"])
@permission_classes([IsAuthenticated])
def update_lab(request):
    if request.method == "PUT":
        data = request.data
        if isinstance(data, QueryDict):
            data = data.dict()
        if "lab_name" in data:
            lab_obj = wetlab.api.utils.lab.get_laboratory_instance(data["lab_name"])
            if lab_obj is None:
                return Response(status=status.HTTP_204_NO_CONTENT)
            wetlab.api.serializers.LabRequestSerializer.update(lab_obj, data)

            return Response(
                "Successful Update information", status=status.HTTP_201_CREATED
            )
    return Response(status=status.HTTP_400_BAD_REQUEST)
