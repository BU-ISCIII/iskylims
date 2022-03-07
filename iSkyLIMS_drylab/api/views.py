from datetime import datetime

from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response
from iSkyLIMS_drylab.models import Service, Resolution, ResolutionStates
from django_utils.models import Profile,Center

from .serializers import *

from iSkyLIMS_drylab.utils.handling_resolutions import send_resolution_in_progress_email


def check_valid_date_format(date):
    try:
        datetime.strptime(date, "%Y")
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


@api_view(["GET"])
def service_list(request):
    param_requests = request.GET.keys()
    for param_request in param_requests:
        if param_request not in ["date", "state"]:
            return Response(status=status.HTTP_400_BAD_REQUEST)
    if "date" in request.GET:
        date = request.GET["date"].strip()
        if not check_valid_date_format(date):
            return Response(status=status.HTTP_400_BAD_REQUEST)
        else:
            end_date = date + "-12-31"
            date += "-01-01"
    if "state" in request.GET:
        state = request.GET["state"].strip()
        if not Service.objects.filter(serviceStatus__exact=state).exists():
            return Response(status=status.HTTP_204_NO_CONTENT)

    service_objs = Service.objects.all()
    if "state" in request.GET:
        service_objs = service_objs.filter(serviceStatus__iexact=state).order_by(
            "serviceRequestNumber"
        )
    if "date" in request.GET:
        service_objs = service_objs.filter(
            serviceOnDeliveredDate__range=(date, end_date)
        ).order_by("serviceRequestNumber")
    if len(service_objs) == 0:
        return Response(status=status.HTTP_204_NO_CONTENT)

    services_serializer = ServiceSerializer(service_objs, many=True)
    for item in range(len(services_serializer.data)):
        user_id = services_serializer.data[item]["serviceUserId"]["username"]

        profile_obj = Profile.objects.get(profileUserID__username__exact=user_id)
        services_serializer.data[item]["serviceUserId"][
            "Center"
        ] = profile_obj.get_profile_center_abbr()
        services_serializer.data[item]["serviceUserId"][
            "Area"
        ] = profile_obj.get_clasification_area()
        
        if Resolution.objects.filter(resolutionServiceID__pk__exact = services_serializer.data[item]['pk']).exists():
            services_serializer.data[item]["ServiceFolderName"] = Resolution.objects.filter(resolutionServiceID__pk__exact = services_serializer.data[item]['pk']).last().resolutionFullNumber
        else:
            services_serializer.data[item]["ServiceFolderName"] = None

    for item in range(len(services_serializer.data)):
        if Resolution.objects.filter(
            resolutionServiceID__pk__exact=services_serializer.data[item]["pk"]
        ).exists():
            resolution_objs = Resolution.objects.filter(
                resolutionServiceID__pk__exact=services_serializer.data[item]["pk"]
            )
            resolution_list = []
            for resolution_obj in resolution_objs:
                resolution_data = {}
                resolution_data["id"] = resolution_obj.get_resolution_id()
                resolution_data["number"] = resolution_obj.get_resolution_number()
                resolution_data["state"] = resolution_obj.get_resolution_state()
                resolution_list.append(resolution_data)
            services_serializer.data[item]["resolutions"] = resolution_list

    return Response(services_serializer.data, status=status.HTTP_200_OK)


@api_view(["GET"])
def resolution_data(request):
    if "resolution" in request.GET:
        resolution = request.GET["resolution"].strip()
        if Resolution.objects.filter(resolutionNumber__exact=resolution).exists():
            resolution_obj = Resolution.objects.filter(
                resolutionNumber__exact=rresolution
            ).last()
            resolution_serializer = ResolutionSerializer(resolution_obj, many=False)
            return Response(resolution_serializer.data, status=status.HTTP_200_OK)
        else:
            return Response(status=status.HTTP_204_NO_CONTENT)
    elif "state" in request.GET:
        if Resolution.objects.filter(
            resolutionState__resolutionStateName__exact=request.GET["state"]
        ).exists():
            resolution_objs = Resolution.objects.filter(
                resolutionState__resolutionStateName__exact=request.GET["state"]
            )
            resolution_serializer = ResolutionSerializer(resolution_objs, many=True)
            return Response(resolution_serializer.data, status=status.HTTP_200_OK)
        else:
            return Response(status=status.HTTP_204_NO_CONTENT)
    else:
        return Response(status=status.HTTP_400_BAD_REQUEST)


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


@api_view(["GET"])
def service_full_data(request):
    if "service" in request.GET:
        service = request.GET["service"].strip()
        if Service.objects.filter(serviceRequestNumber__iexact=service).exists():
            service_full_data = {}
            service_obj = Service.objects.filter(
                serviceRequestNumber__iexact=service
            ).last()
            service_full_data["Service"] = ServiceSerializer(
                service_obj, many=False
            ).data

            for item in range(len(services_serializer.data)):
                user_id = services_serializer.data[item]["serviceUserId"]["username"]
                profile_obj = Profile.objects.get(
                    profileUserID__username__exact=user_id
                )
                services_serializer.data[item]["serviceUserId"][
                    "Center"
                ] = profile_obj.get_profile_center_abbr()
                services_serializer.data[item]["serviceUserId"][
                    "Area"
                ] = profile_obj.get_clasification_area()

            if Resolution.objects.filter(resolutionServiceID=service_obj).exists():
                resolution_objs = Resolution.objects.filter(
                    resolutionServiceID=service_obj
                )
                service_full_data["Resolutions"] = ResolutionSerializer(
                    resolution_objs, many=True
                ).data
            if RequestedSamplesInServices.objects.filter(
                samplesInService__serviceRequestNumber__iexact=service
            ).exists():
                sample_objs = RequestedSamplesInServices.objects.filter(
                    samplesInService__serviceRequestNumber__iexact=service
                )
                service_full_data["Samples"] = RequestedSamplesInServicesSerializer(
                    sample_objs, many=True
                ).data

            return Response(service_full_data, status=status.HTTP_200_OK)
        else:
            return Response(status=status.HTTP_204_NO_CONTENT)
    else:
        return Response(status=status.HTTP_400_BAD_REQUEST)


@api_view(["GET"])
def resolution_full_data(request):
    if "resolution" in request.GET:
        resolution = request.GET["resolution"].strip()
        if Resolution.objects.filter(resolutionNumber__iexact=resolution).exists():
            resolution_full_data = {}
            resolution_obj = Resolution.objects.filter(
                resolutionNumber__iexact=resolution
            ).last()
            service_obj = resolution_obj.get_service_obj()
            resolution_full_data["Service"] = ServiceSerializer(
                service_obj, many=False
            ).data
            resolution_full_data["Resolutions"] = ResolutionSerializer(
                resolution_obj, many=False
            ).data
            if RequestedSamplesInServices.objects.filter(
                samplesInService=service_obj
            ).exists():
                sample_objs = RequestedSamplesInServices.objects.filter(
                    samplesInService=service_obj
                )
                resolution_full_data["Samples"] = RequestedSamplesInServicesSerializer(
                    sample_objs, many=True
                ).data

            return Response(resolution_full_data, status=status.HTTP_200_OK)
        else:
            return Response(status=status.HTTP_204_NO_CONTENT)
    else:
        return Response(status=status.HTTP_400_BAD_REQUEST)


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


@api_view(["POST"])
def create(request):
    if request.method == "POST":
        data = request.data
        if isinstance(data, QueryDict):
            data = data.dict()
        if "delivery" in data:
            serializer = CreateDeliveryPostSerializer(data=data)
            if not serializer.is_valid():
                return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

            robot_action_obj = serializer.save()
