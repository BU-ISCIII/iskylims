from rest_framework.authentication import SessionAuthentication, BasicAuthentication
from rest_framework.permissions import IsAuthenticated

from rest_framework.decorators import (
    authentication_classes,
    permission_classes,
    api_view,
)
from rest_framework import status
from rest_framework.response import Response
from django.http import QueryDict

from iSkyLIMS_core.models import SampleProjects, SampleProjectsFields, LabRequest
from .serializers import (
    CreateSampleSerializer,
    CreateProjectDataSerializer,
    SampleProjectFieldSerializer,
    # SampleFieldsSerializer,
    SampleFields,
    LabRequestSerializer,
)

from drf_yasg.utils import swagger_auto_schema
from drf_yasg import openapi

from .utils.sample_request_handling import (
    split_sample_data,
    include_instances_in_sample,
    include_codding,
    get_sample_fields,
)

sample_project_fields = openapi.Parameter(
    "project",
    openapi.IN_QUERY,
    description="Project name to fetch the sample project fields defined. Example Relecov",
    type=openapi.TYPE_STRING,
)

sample_fields = openapi.Parameter(
    "project",
    openapi.IN_QUERY,
    description="Fetch the sample fields",
    type=openapi.TYPE_STRING,
)

laboratory = openapi.Parameter(
    "laboratory",
    openapi.IN_QUERY,
    description="Laboratory name form to fetch contact information. Example Harlem Hospital Center",
    type=openapi.TYPE_STRING,
)


@swagger_auto_schema(
    method="post",
    request_body=openapi.Schema(
        type=openapi.TYPE_OBJECT,
        properties={
            "sample": openapi.Schema(
                type=openapi.TYPE_STRING, description="Sample name"
            ),
            "sampleState": openapi.Schema(
                type=openapi.TYPE_STRING, description="Sample state"
            ),
            "patientCore": openapi.Schema(
                type=openapi.TYPE_STRING, description="Code assigned to the patient"
            ),
            "labRequest": openapi.Schema(
                type=openapi.TYPE_STRING,
                description="Laboratory that request the sample",
            ),
            "sampleType": openapi.Schema(
                type=openapi.TYPE_STRING, description="Type of the sample"
            ),
            "species": openapi.Schema(
                type=openapi.TYPE_STRING,
                description="Specie that the sample belongs to",
            ),
            "sampleLocation": openapi.Schema(
                type=openapi.TYPE_STRING,
                description="Location where the sample is stored",
            ),
            "sampleEntryDate": openapi.Schema(
                type=openapi.TYPE_STRING,
                description="Date when sample is received in the lab",
            ),
            "sampleCollectionDate": openapi.Schema(
                type=openapi.TYPE_STRING,
                description="Date when the sample is collected from the specimen",
            ),
            "onlyRecorded": openapi.Schema(
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
        if "sampleName" not in data and "sampleProject" not in data:
            return Response(status=status.HTTP_400_BAD_REQUEST)
        split_data = split_sample_data(data)
        if not isinstance(split_data, dict):
            return Response(split_data, status=status.HTTP_400_BAD_REQUEST)
        inst_req_sample = include_instances_in_sample(split_data["s_data"])
        if not isinstance(inst_req_sample, dict):
            return Response(inst_req_sample, status=status.HTTP_400_BAD_REQUEST)
        split_data["s_data"] = inst_req_sample
        split_data["s_data"]["sampleUser"] = request.user.pk
        # Adding coding for sample
        split_data["s_data"].update(
            include_codding(request.user.username, split_data["s_data"]["sampleName"])
        )
        sample_serializer = CreateSampleSerializer(data=split_data["s_data"])
        if not sample_serializer.is_valid():
            return Response(
                sample_serializer.errors, status=status.HTTP_400_BAD_REQUEST
            )
        new_sample_id = sample_serializer.save().get_sample_id()
        for d_field in split_data["p_data"]:
            d_field["sample_id"] = new_sample_id
            s_project_serializer = CreateProjectDataSerializer(data=d_field)
            if not s_project_serializer.is_valid():
                return Response(
                    s_project_serializer.errors, status=status.HTTP_400_BAD_REQUEST
                )
            s_project_serializer.save()

        return Response("Successful upload information", status=status.HTTP_201_CREATED)


@swagger_auto_schema(method="get")
@api_view(["GET"])
def sample_fields(request):
    apps_name = __package__.split(".")[0]
    #sample_fields = SampleFieldsSerializer(SampleFields(get_sample_fields(apps_name)))

    sample_fields = get_sample_fields(apps_name)

    if "ERROR" in sample_fields:
        return Response(sample_fields, status=status.HTTP_202_ACCEPTED)
        # return Response(sample_fields.data, status=status.HTTP_204_NO_CONTENT)
    else:
        return Response(sample_fields, status=status.HTTP_200_OK)
    return Response(status=status.HTTP_400_BAD_REQUEST)


@swagger_auto_schema(method="get", manual_parameters=[sample_project_fields])
@api_view(["GET"])
def sample_project_fields(request):
    if "project" in request.GET:
        project = request.GET["project"].strip()
        if SampleProjects.objects.filter(sampleProjectName__iexact=project).exists():
            s_project_obj = SampleProjects.objects.filter(
                sampleProjectName__iexact=project
            ).last()
            s_project_field_objs = SampleProjectsFields.objects.filter(
                sampleProjects_id=s_project_obj
            )
            s_project_serializer = SampleProjectFieldSerializer(
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
        if LabRequest.objects.filter(labName__iexact=lab_name).exists():
            lab_req_obj = LabRequest.objects.filter(labName__iexact=lab_name).last()
            lab_req_serializer = LabRequestSerializer(lab_req_obj, many=False)
            return Response(lab_req_serializer.data, status=status.HTTP_200_OK)
        else:
            return Response(status=status.HTTP_204_NO_CONTENT)
    return Response(status=status.HTTP_400_BAD_REQUEST)


@api_view(["PUT"])
def update_lab(request):
    if request.method == "PUT":
        data = request.data
        if isinstance(data, QueryDict):
            data = data.dict()
        data["apps_name"] = "wetlab"
        serializer = LabRequestSerializer(data=data)
        if not serializer.is_valid():
            return Response(status=status.HTTP_400_BAD_REQUEST)
        lab_name = data["labName"].strip()
        if not LabRequest.objects.filter(labName__iexact=lab_name).exists():
            return Response(status=status.HTTP_400_BAD_REQUEST)
        lab_req_obj = LabRequest.objects.filter(labName__iexact=lab_name).last()
        LabRequestSerializer.update(lab_req_obj, data=data)

        return Response("Successful Update information", status=status.HTTP_201_CREATED)
    return Response(status=status.HTTP_400_BAD_REQUEST)
