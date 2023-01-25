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
    Pipelines,
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
    UpdateResolutionStateSerializer,
    UpdateServiceStateSerializer,
    CreateDeliveryPostSerializer,
)

from iSkyLIMS_drylab.utils.handling_resolutions import (
    send_resolution_in_progress_email,
)

from iSkyLIMS_drylab.utils.handling_deliveries import (
    send_delivery_service_email,
)

def check_valid_date_format(date):
    try:
        datetime.strptime(date, "%Y-%m-%d")
        return True
    except ValueError:
        return False


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
def update_state(request):
    def all_resolutions_delivered(service,state):
        if Resolution.objects.filter(resolutionServiceID=service).exclude(resolutionState=state).exists():
            return False
        else:
            return True

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

            data_resolution = {
                "resolutionNumber" : resolution,
                "resolutionState" : state_obj.pk
            }

            resolution_serializer = UpdateResolutionStateSerializer(resolution_obj, data=data_resolution)
            resolution_serializer.is_valid(raise_exception=True)
            resolution_serializer.save()

            service_obj = resolution_obj.get_service_obj()

            email_data = {}
            email_data["user_email"] = service_obj.get_user_email()
            email_data["user_name"] = service_obj.get_username()
            email_data["resolution_number"] = resolution_obj.get_resolution_number()

            if state == "In Progress":
                 data_service = {
                     "serviceRequestNumber" : service_obj.serviceRequestNumber,
                     "serviceStatus" : "in_progress"
                 }
                 serializer_services = UpdateServiceStateSerializer(service_obj,data=data_service)
                 if serializer_services.is_valid(raise_exception=True):
                    serializer_services.save()
                 # Send email in progress
                 send_resolution_in_progress_email(email_data)

            elif state == "Delivery" and all_resolutions_delivered(service_obj,state_obj):
                 data_service = {
                     "serviceRequestNumber" : service_obj.serviceRequestNumber,
                     "serviceStatus" : "delivered"
                 }
                 serializer_services = UpdateServiceStateSerializer(service_obj,data=data_service)
                 if serializer_services.is_valid(raise_exception=True):
                    serializer_services.save()
                 # Send email
                 send_delivery_service_email(email_data)

            return Response(
                resolution_serializer.data, status=status.HTTP_200_OK
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
            "resolutionNumber": openapi.Schema(
                type=openapi.TYPE_STRING, description="resolutionNumber. pe. SRVCNM123.1"
            ),
            "pipelinesInDelivery": openapi.Schema(
                type=openapi.TYPE_ARRAY, items = openapi.Items(type="string"), description="resolutionNumber. pe. ['viralrecon']"
            ),
            "deliveryDate": openapi.Schema(
                type=openapi.TYPE_STRING, description="delivery date. pe. 2022-01-01"
            ),
            "executionStartDate": openapi.Schema(
                type=openapi.TYPE_STRING, description="execution start date. pe. 2022-01-01"
            ),
            "executionEndDate": openapi.Schema(
                type=openapi.TYPE_STRING, description="execution end date. pe. 2022-01-01"
            ),
            "permanentUsedSpace": openapi.Schema(
                type=openapi.TYPE_INTEGER, description="permanent used space in GB. pe. 134"
            ),
            "temporaryUsedSpace": openapi.Schema(
                type=openapi.TYPE_INTEGER, description="temporary used space in GB. pe. SRVCNM123.1"
            ),
        },
    ),
    responses={200: "Successful delivery creation", 400: "Bad Request"},
)
@authentication_classes([SessionAuthentication, BasicAuthentication])
@api_view(["POST"])
@permission_classes([IsAuthenticated])
def create_delivery(request):
    if request.method == "POST":
        data = request.data
        if isinstance(data, QueryDict):
            data = data.dict()

        resolution_pk = Resolution.objects.filter(resolutionNumber__exact=data["resolutionNumber"]).last().pk
        pipelines = [ Pipelines.objects.filter(pipelineName__exact=pip).last().pk for pip in data["pipelinesInDelivery"]]

        data.pop("resolutionNumber")
        data["deliveryResolutionID"] = resolution_pk
        data["pipelinesInDelivery"] = pipelines

        serializer = CreateDeliveryPostSerializer(data=data)

        if not serializer.is_valid():
            return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)
        serializer.save()

        return Response(serializer.data, status=status.HTTP_200_OK)

    return Response(status=status.HTTP_400_BAD_REQUEST)
