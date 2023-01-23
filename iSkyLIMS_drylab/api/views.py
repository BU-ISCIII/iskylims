from datetime import datetime
from rest_framework.authentication import SessionAuthentication, BasicAuthentication
from rest_framework.permissions import IsAuthenticated
from django.db.models import Prefetch

from rest_framework.decorators import (
    authentication_classes,
    permission_classes,
    api_view,
)
from rest_framework import status
from rest_framework.response import Response

from iSkyLIMS_drylab.models import (
    Service,
    Resolution,
    ResolutionStates,
    RequestedSamplesInServices,
)
from django_utils.models import Profile
from django.http import QueryDict

from drf_yasg.utils import swagger_auto_schema
from drf_yasg import openapi
from .serializers import (
    ServiceListSerializer,
    ServiceSerializer,
    ResolutionSerializer,
    RequestedSamplesInServicesSerializer,
    UpdateResolutionSerializer,
    CreateDeliveryPostSerializer,
)

from iSkyLIMS_drylab.utils.handling_resolutions import send_resolution_in_progress_email


def check_valid_date_format(date):
    try:
        datetime.strptime(date, "%Y-%m-%d")
        return True
    except ValueError:
        return False


"""
try:
   from iSkyLIMS_wetlab.utils.api.wetlab_api  import *
   wetlab_api_available = True
except:
   wetlab_api_available = False

def get_projectsid(service):
   serv = Service.objects.get(serviceRequestNumber__exact = service)
   service_id = serv.id
   project_list = RequestedProjectInServices.objects.filter(projectService = service_id)
   projects_id=[]
   for project in project_list:
       projects_id.append(str(project.get_requested_external_project_id()))
   return projects_id
"""
service_state_param = openapi.Parameter(
    "state",
    openapi.IN_QUERY,
    description="State parameter is optional. The allowed values are: [approved/rejected/queued/in_progress/delivered/archived/recorded]",
    enum=[
        "approved",
        "rejected",
        "queued",
        "in_progress",
        "delivered",
        "archived",
        "recorded",
    ],
    type=openapi.TYPE_STRING,
)
date_from_param = openapi.Parameter(
    "date_from",
    openapi.IN_QUERY,
    description="Date parameter is optional.It will limit the results from the date specified in the parameter. Example 2022-01-01",
    type=openapi.TYPE_STRING,
)

date_until_param = openapi.Parameter(
    "date_until",
    openapi.IN_QUERY,
    description="Date parameter is optional.It will limit the results up to the date specified in the parameter. Example 2022-01-01",
    type=openapi.TYPE_STRING,
)

resolution_state_param = openapi.Parameter(
    "state",
    openapi.IN_QUERY,
    description="State parameter is optional. The allowed values are: [Recorded/ In Progress/ Delivery/ Cancelled]",
    enum=["Recorded", "In Progress", "Delivery", "Cancelled"],
    type=openapi.TYPE_STRING,
)
resolution_number_param = openapi.Parameter(
    "resolution",
    openapi.IN_QUERY,
    description="Resolution parameter resolution is optional. Example SRVCNM123.1",
    type=openapi.TYPE_STRING,
)
service_name_param = openapi.Parameter(
    "service",
    openapi.IN_QUERY,
    description="Service parameter is mandatory. Example SRVCNM123",
    type=openapi.TYPE_STRING,
)
resolution_state_mand = openapi.Parameter(
    "status",
    openapi.IN_QUERY,
    description="Parameter status is mandatory and can take the following possible values:[approved/rejected/queued/in_progress/delivered/archived/recorded]",
    type=openapi.TYPE_STRING,
)


@swagger_auto_schema(
    method="get", manual_parameters=[service_state_param, date_from_param, date_until_param]
)
@api_view(["GET"])
def service_list(request):
    param_requests = request.GET.keys()
    for param_request in param_requests:
        if param_request not in ["date_from", "date_until", "state"]:
            return Response(status=status.HTTP_400_BAD_REQUEST)
    if "date_from" in request.GET:
        date_from = request.GET["date_from"].strip()
        if not check_valid_date_format(date_from):
            return Response(status=status.HTTP_400_BAD_REQUEST)
    if "date_until" in request.GET:
        date_until = request.GET["date_until"].strip()
        if not check_valid_date_format(date_until):
            return Response(status=status.HTTP_400_BAD_REQUEST)
    if "state" in request.GET:
        state = request.GET["state"].strip()
        if not Service.objects.filter(serviceStatus__exact=state).exists():
            return Response(status=status.HTTP_204_NO_CONTENT)

    service_objs = Service.objects.all()
    if "state" in request.GET:
        service_objs = service_objs.filter(serviceStatus__iexact=state).order_by(
            "serviceRequestNumber"
        )
    if "date_from" in request.GET and "date_until" in request.GET:
        service_objs = service_objs.filter(
            serviceOnDeliveredDate__range=(date_from, date_until)
        ).order_by("serviceRequestNumber")
        if len(service_objs) == 0:
            return Response(status=status.HTTP_204_NO_CONTENT)
    elif "date_from" in request.GET:
        date_until = datetime.today()
        service_objs = service_objs.filter(
            serviceOnDeliveredDate__range=(date_from, date_until)
        ).order_by("serviceRequestNumber")
        if len(service_objs) == 0:
            return Response(status=status.HTTP_204_NO_CONTENT)

    services_list_serializer = ServiceListSerializer(service_objs, many=True)

    return Response(services_list_serializer.data, status=status.HTTP_200_OK)


@swagger_auto_schema(
    method="get", manual_parameters=[resolution_state_param, resolution_number_param]
)
@api_view(["GET"])
def resolution_data(request):
    if "resolution" in request.GET:
        resolution = request.GET["resolution"].strip()
        if Resolution.objects.filter(resolutionNumber__exact=resolution).exists():
            resolution_obj = Resolution.objects.filter(
                resolutionNumber__exact=resolution
            ).last()
        else:
            return Response(status=status.HTTP_204_NO_CONTENT)
    elif "state" in request.GET:
        if Resolution.objects.filter(
            resolutionState__resolutionStateName__exact=request.GET["state"]
        ).exists():
            resolution_objs = Resolution.objects.filter(
                resolutionState__resolutionStateName__exact=request.GET["state"]
            )
        else:
            return Response(status=status.HTTP_204_NO_CONTENT)
    else:
        return Response(status=status.HTTP_400_BAD_REQUEST)

    resolution_serializer = ResolutionSerializer(resolution_obj, many=False)
    return Response(resolution_serializer.data, status=status.HTTP_200_OK)


@swagger_auto_schema(method="get", manual_parameters=[service_name_param])
@api_view(["GET"])
def samples_in_service(request):
    if "service" in request.GET:
        if RequestedSamplesInServices.objects.filter(
            samplesInService__serviceRequestNumber__iexact=request.GET["service"]
        ).exists():
            sample_objs = RequestedSamplesInServices.objects.filter(
                samplesInService__serviceRequestNumber__iexact=request.GET["service"]
            )
            sample_serializers = RequestedSamplesInServicesSerializer(
                sample_objs, many=True
            )
            return Response(sample_serializers.data, status=status.HTTP_200_OK)
        else:
            return Response(status=status.HTTP_204_NO_CONTENT)
    else:
        return Response(status=status.HTTP_400_BAD_REQUEST)


@swagger_auto_schema(method="get", manual_parameters=[service_name_param, resolution_number_param])
@api_view(["GET"])
def service_full_data(request):
    if "service" in request.GET:
        service = request.GET["service"].strip()
        service_obj = Service.objects.prefetch_related(
            Prefetch('resolutions', queryset=Resolution.objects.all(), to_attr='filtered_resolutions')
        )
    elif "resolution" in request.GET:
        resolution = request.GET["resolution"].strip()
        service = Resolution.objects.filter(resolutionNumber__iexact= resolution).last().resolutionServiceID
        service_obj = Service.objects.prefetch_related(
            Prefetch('resolutions', queryset=Resolution.objects.filter(resolutionNumber__iexact=resolution), to_attr='filtered_resolutions')
        )
    else:
        return Response(status=status.HTTP_400_BAD_REQUEST)

    if Service.objects.filter(serviceRequestNumber__iexact=service).exists():
        service_full_data = {}

        service_obj = service_obj.filter(
            serviceRequestNumber__iexact=service
        ).last()

        service_full_data = ServiceSerializer(
            service_obj, many=False
        ).data
        return Response(service_full_data, status=status.HTTP_200_OK)
    else:
        return Response(status=status.HTTP_204_NO_CONTENT)

@swagger_auto_schema(
    method="put", manual_parameters=[resolution_number_param, resolution_state_param]
)
@api_view(["PUT"])
def update_resolution(request):
    if ("resolution" in request.query_params) and ("state" in request.query_params):
        resolution = request.query_params["resolution"].strip()
        if Resolution.objects.filter(resolutionNumber__exact=resolution).exists():
            resolution_obj = Resolution.objects.get(resolutionNumber__exact=resolution)
            state = request.query_params["state"].strip()
            try:
                state_obj = ResolutionStates.objects.get(
                    resolutionStateName__iexact=state
                )
            except Exception:
                return Response(status=status.HTTP_400_BAD_REQUEST)

            UpdateResolutionSerializer.update(resolution_obj, state_obj)
            updated_resolution_serializer = UpdateResolutionSerializer(resolution_obj)
            service_obj = resolution_obj.get_service_obj()

            email_data = {}
            email_data["user_email"] = service_obj.get_user_email()
            email_data["user_name"] = service_obj.get_username()
            email_data["resolution_number"] = resolution_obj.get_resolution_number()
            send_resolution_in_progress_email(email_data)
            return Response(
                updated_resolution_serializer.data, status=status.HTTP_200_OK
            )

        else:
            return Response(status=status.HTTP_204_NO_CONTENT)
    else:
        return Response(status=status.HTTP_400_BAD_REQUEST)


@swagger_auto_schema(
    method="post",
    request_body=openapi.Schema(
        type=openapi.TYPE_OBJECT,
        properties={
            "delvery": openapi.Schema(
                type=openapi.TYPE_STRING, description="Delivery ID"
            )
        },
    ),
    responses={200: "Successful delivery creation", 400: "Bad Request"},
)
@authentication_classes([SessionAuthentication, BasicAuthentication])
@api_view(["POST"])
@permission_classes([IsAuthenticated])
def create(request):
    if request.method == "POST":
        data = request.data
        if isinstance(data, QueryDict):
            data = data.dict()
        if "delivery" in data:
            serializer = CreateDeliveryPostSerializer(data=data)
            if not serializer.is_valid():
                return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

            serializer.save()
            return Response("Successful delivery creation", status=status.HTTP_200_OK)
    return Response(status=status.HTTP_400_BAD_REQUEST)
