from rest_framework.authentication import SessionAuthentication, BasicAuthentication
from rest_framework.permissions import IsAuthenticated

from rest_framework.decorators import authentication_classes, permission_classes
from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response
from django.http import QueryDict

from .serializers import CreateSampleSerializer, CreateProjectDataSerializer
from .utils.request_handling import (
    split_sample_data,
    include_instances_in_sample,
    include_codding,
)


@api_view(["GET"])
def sample_list(request):
    param_requests = request.GET.keys()
    for param_request in param_requests:
        if param_request not in ["date", "state", "onlyRecorded"]:
            return Response(status=status.HTTP_400_BAD_REQUEST)


@api_view(["POST"])
@authentication_classes([SessionAuthentication, BasicAuthentication])
@permission_classes([IsAuthenticated])
def create_sample_data(request):
    if request.method == "POST":
        data = request.data

        if isinstance(data, QueryDict):
            data = data.dict()
        if "sample" not in data and "project" not in data:
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
        import pdb; pdb.set_trace()

        if not sample_serializer.is_valid():
            return Response(
                sample_serializer.errors, status=status.HTTP_400_BAD_REQUEST
            )
        # new_sample_obj = sample_serializer.save()
        project_serializer = CreateProjectDataSerializer(data=split_data["p_data"])
        if not project_serializer.is_valid():
            return Response(
                project_serializer.errors, status=status.HTTP_400_BAD_REQUEST
            )
        p_data_obj = project_serializer.save()

        sample_obj = serializer.save()
