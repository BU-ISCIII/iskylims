from datetime import datetime

from django.db.models import Prefetch
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

from drylab.models import (
    Delivery,
    Pipelines,
    RequestedSamplesInServices,
    Resolution,
    ResolutionStates,
    Service,
)
from drylab.utils.deliveries import send_delivery_service_email
from drylab.utils.resolutions import send_resolution_in_progress_email

from .serializers import (
    CreateDeliveryPostSerializer,
    RequestedSamplesInServicesSerializer,
    ResolutionSerializer,
    ServiceListSerializer,
    ServiceSerializer,
    UpdateResolutionStateSerializer,
    UpdateServiceStateSerializer,
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
    description="""State parameter is optional.
                   The allowed values are: [approved/rejected/queued/in_progress/delivered/archived/recorded]""",
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
    description="""Date parameter is optional.It will limit the results from the date specified in the parameter.
                Example 2022-01-01""",
    type=openapi.TYPE_STRING,
)

date_until_param = openapi.Parameter(
    "date_until",
    openapi.IN_QUERY,
    description="""Date parameter is optional.It will limit the results up to the date specified in the parameter.
                Example 2022-01-01""",
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
    description="""Parameter status is mandatory.
                It can take the following possible values:
                [approved/rejected/queued/in_progress/delivered/archived/recorded]""",
    type=openapi.TYPE_STRING,
)


@swagger_auto_schema(
    method="get",
    manual_parameters=[service_state_param, date_from_param, date_until_param],
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
        if not Service.objects.filter(service_state__exact=state).exists():
            return Response(status=status.HTTP_204_NO_CONTENT)

    service_objs = Service.objects.all()
    if "state" in request.GET:
        service_objs = service_objs.filter(service_state__iexact=state).order_by(
            "service_request_number"
        )
    if "date_from" in request.GET and "date_until" in request.GET:
        service_objs = service_objs.filter(
            service_delivered_date__range=(date_from, date_until)
        ).order_by("service_request_number")
        if len(service_objs) == 0:
            return Response(status=status.HTTP_204_NO_CONTENT)
    elif "date_from" in request.GET:
        date_until = datetime.today()
        service_objs = service_objs.filter(
            service_delivered_date__range=(date_from, date_until)
        ).order_by("service_request_number")
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
        if Resolution.objects.filter(resolution_number__exact=resolution).exists():
            resolution_obj = Resolution.objects.filter(
                resolution_number__exact=resolution
            ).last()
        else:
            return Response(status=status.HTTP_204_NO_CONTENT)
    elif "state" in request.GET:
        if Resolution.objects.filter(
            resolution_state__resolution_stateName__exact=request.GET["state"]
        ).exists():
            resolution_obj = Resolution.objects.filter(
                resolution_state__resolution_stateName__exact=request.GET["state"]
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
            samples_in_service__service_request_number__iexact=request.GET["service"]
        ).exists():
            sample_objs = RequestedSamplesInServices.objects.filter(
                samples_in_service__service_request_number__iexact=request.GET[
                    "service"
                ]
            )
            sample_serializers = RequestedSamplesInServicesSerializer(
                sample_objs, many=True
            )
            return Response(sample_serializers.data, status=status.HTTP_200_OK)
        else:
            return Response(status=status.HTTP_204_NO_CONTENT)
    else:
        return Response(status=status.HTTP_400_BAD_REQUEST)


@swagger_auto_schema(
    method="get", manual_parameters=[service_name_param, resolution_number_param]
)
@api_view(["GET"])
def service_full_data(request):
    if "service" in request.GET:
        service = request.GET["service"].strip()
        service_obj = Service.objects.prefetch_related(
            Prefetch(
                "resolutions",
                queryset=Resolution.objects.all(),
                to_attr="filtered_resolutions",
            )
        )
    elif "resolution" in request.GET:
        resolution = request.GET["resolution"].strip()
        service = (
            Resolution.objects.filter(resolution_number__iexact=resolution)
            .last()
            .resolution_service_id
        )
        service_obj = Service.objects.prefetch_related(
            Prefetch(
                "resolutions",
                queryset=Resolution.objects.filter(resolution_number__iexact=resolution),
                to_attr="filtered_resolutions",
            )
        )
    else:
        return Response(status=status.HTTP_400_BAD_REQUEST)

    if Service.objects.filter(service_request_number__iexact=service).exists():
        service_full_data = {}

        service_obj = service_obj.filter(service_request_number__iexact=service).last()

        service_full_data = ServiceSerializer(service_obj, many=False).data
        return Response(service_full_data, status=status.HTTP_200_OK)
    else:
        return Response(status=status.HTTP_204_NO_CONTENT)


@swagger_auto_schema(
    method="put", manual_parameters=[resolution_number_param, resolution_state_param]
)
@api_view(["PUT"])
def update_state(request):
    def all_resolutions_delivered(service, state):
        if (
            Resolution.objects.filter(resolution_service_id=service)
            .exclude(resolution_state=state)
            .exists()
        ):
            return False
        else:
            return True

    if ("resolution" in request.query_params) and ("state" in request.query_params):
        resolution = request.query_params["resolution"].strip()
        if Resolution.objects.filter(resolution_number__exact=resolution).exists():
            resolution_obj = Resolution.objects.get(resolution_number__exact=resolution)
            state = request.query_params["state"].strip()
            try:
                state_obj = ResolutionStates.objects.get(
                    resolution_state_name__iexact=state
                )
            except Exception:
                return Response(status=status.HTTP_400_BAD_REQUEST)

            service_obj = resolution_obj.get_service_obj()

            email_data = {}
            email_data["user_email"] = service_obj.get_user_email()
            email_data["user_name"] = service_obj.get_username()
            email_data["resolution_number"] = resolution_obj.get_resolution_number()

            if state == "In Progress":
                data_resolution = {
                    "resolution_number": resolution,
                    "resolution_state": state_obj.pk,
                    "resolution_inprogress_date": datetime.today().strftime("%Y-%m-%d"),
                }
                data_service = {
                    "service_request_number": service_obj.service_request_number,
                    "service_state": "in_progress",
                }
                # Send email in progress
                send_resolution_in_progress_email(email_data)

            elif state == "Delivery":
                data_resolution = {
                    "resolution_number": resolution,
                    "resolution_state": state_obj.pk,
                    "resolution_delivery_date": datetime.today().strftime("%Y-%m-%d"),
                }
                data_service = {
                    "service_request_number": service_obj.service_request_number,
                    "service_state": "delivered",
                    "service_delivered_date": datetime.today().strftime("%Y-%m-%d"),
                }
                # Send email
                send_delivery_service_email(email_data)

            resolution_serializer = UpdateResolutionStateSerializer(
                resolution_obj, data=data_resolution
            )
            if resolution_serializer.is_valid(raise_exception=True):
                resolution_serializer.save()

            if (
                state == "Delivery"
                and all_resolutions_delivered(service_obj, state_obj)
            ) or state == "In Progress":
                serializer_services = UpdateServiceStateSerializer(
                    service_obj, data=data_service
                )
                if serializer_services.is_valid(raise_exception=True):
                    serializer_services.save()

            return Response(resolution_serializer.data, status=status.HTTP_200_OK)

        else:
            return Response(status=status.HTTP_204_NO_CONTENT)
    else:
        return Response(status=status.HTTP_400_BAD_REQUEST)


@swagger_auto_schema(
    method="post",
    request_body=openapi.Schema(
        type=openapi.TYPE_OBJECT,
        properties={
            "resolution_number": openapi.Schema(
                type=openapi.TYPE_STRING,
                description="resolution_number. pe. SRVCNM123.1",
            ),
            "pipelines_in_delivery": openapi.Schema(
                type=openapi.TYPE_ARRAY,
                items=openapi.Items(type="string"),
                description="pipelines in delivery. pe. ['viralrecon']",
            ),
            "delivery_date": openapi.Schema(
                type=openapi.TYPE_STRING, description="delivery date. pe. 2022-01-01"
            ),
            "execution_start_date": openapi.Schema(
                type=openapi.TYPE_STRING,
                description="execution start date. pe. 2022-01-01",
            ),
            "execution_end_date": openapi.Schema(
                type=openapi.TYPE_STRING,
                description="execution end date. pe. 2022-01-01",
            ),
            "permanent_used_space": openapi.Schema(
                type=openapi.TYPE_INTEGER,
                description="permanent used space in GB. pe. 134",
            ),
            "temporary_used_space": openapi.Schema(
                type=openapi.TYPE_INTEGER,
                description="temporary used space in GB. pe. SRVCNM123.1",
            ),
            "delivery_notes": openapi.Schema(
                type=openapi.TYPE_INTEGER,
                description="delivery notes with brief results description.",
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

        resolution_pk = (
            Resolution.objects.filter(resolution_number__exact=data["resolution_number"])
            .last()
            .pk
        )

        if "pipelines_in_delivery" in data:
            pipelines = [
                Pipelines.objects.filter(pipeline_name__exact=pip).last().pk
                for pip in data["pipelines_in_delivery"]
            ]
            data["pipelines_in_delivery"] = pipelines

        data.pop("resolution_number")
        data["delivery_resolution_id"] = resolution_pk

        if Delivery.objects.filter(
            delivery_resolution_id__exact=resolution_pk
        ).exists():
            delivery_obj = Delivery.objects.filter(
                delivery_resolution_id__exact=resolution_pk
            ).last()
            serializer = CreateDeliveryPostSerializer(delivery_obj, data=data)
        else:
            serializer = CreateDeliveryPostSerializer(data=data)

        if not serializer.is_valid():
            return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)
        serializer.save()

        return Response(serializer.data, status=status.HTTP_200_OK)

    return Response(status=status.HTTP_400_BAD_REQUEST)
