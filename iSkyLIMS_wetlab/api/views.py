from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response
from .serializers import *


@api_view(["GET"])
def sample_list(request):
    param_requests = request.GET.keys()
    for param_request in param_requests:
        if param_request not in ["date", "state", "onlyRecorded"]:
            return Response(status=status.HTTP_400_BAD_REQUEST)



@api_view(["POST"])
def create_sample_data(request):
    if request.method == "POST":
        data = request.data
        import pdb; pdb.set_trace()
        if isinstance(data, QueryDict):
            data = data.dict()
        if "sample" not in data and "project" not in data:
            return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)
        if "project" in data:

            serializer = CreateSampleDataInProject(data=data)
            if not serializer.is_valid():
                return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

            sample_obj = serializer.save()
