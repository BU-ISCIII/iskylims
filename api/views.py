from django.shortcuts import render
from django.shortcuts import get_object_or_404
from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response
from iSkyLIMS_drylab import drylab_config
from iSkyLIMS_drylab.models import *

from .serializers import *

from iSkyLIMS_drylab.utils.handling_resolutions import send_resolution_in_progress_email

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

@api_view(['GET',])
def service_list(request ):

    if 'state' in request.GET:
        state = request.GET['state']
        if Service.objects.filter(serviceStatus__exact = state).exists():
            service_objs = Service.objects.filter(serviceStatus__exact = state).order_by('serviceRequestNumber')
        else:
            return Response(status = status.HTTP_204_NO_CONTENT)
    else:
        if Service.objects.all().exists():
            service_objs = Service.objects.all()
        else:
            return Response(status = status.HTTP_204_NO_CONTENT)

    services_serializer = ServiceSerializer(service_objs, many=True)

    for item in range(len(services_serializer.data)):
        if Resolution.objects.filter(resolutionServiceID__pk__exact = services_serializer.data[item]['pk']).exists():
            resolution_objs = Resolution.objects.filter(resolutionServiceID__pk__exact = services_serializer.data[item]['pk'])
            resolution_list = []
            for resolution_obj in resolution_objs:
                resolution_data = {}
                resolution_data['id'] = resolution_obj.get_resolution_id()
                resolution_data['number'] = resolution_obj.get_resolution_number()
                resolution_data['state'] = resolution_obj.get_resolution_state()
                resolution_list.append(resolution_data)
            services_serializer.data[item]['resolutions'] = resolution_list

    return Response(services_serializer.data, status = status.HTTP_200_OK)

@api_view(['GET'])
def resolution_data(request):
    if 'resolution' in request.GET:
        if Resolution.objects.filter(resolutionNumber__exact = request.GET['resolution']).exists():
            resolution_obj = Resolution.objects.filter(resolutionNumber__exact =request.GET['resolution']).last()
            resolution_serializer = ResolutionSerializer(resolution_obj, many=False)
            return Response(resolution_serializer.data, status = status.HTTP_200_OK)
        else:
            return Response(status = status.HTTP_204_NO_CONTENT)
    else:
        return Response(status = status.HTTP_400_BAD_REQUEST)



@api_view(['GET'])
def samples_in_service(request):
    if 'service' in request.GET:
        if RequestedSamplesInServices.objects.filter(samplesInService__serviceRequestNumber__iexact = request.GET['service']).exists():
            sample_objs = RequestedSamplesInServices.objects.filter(samplesInService__serviceRequestNumber__iexact = request.GET['service'])
            sample_serializers = RequestedSamplesInServicesSerializer(sample_objs, many = True)

            return Response(sample_serializers.data, status = status.HTTP_200_OK)
        else:
            return Response(status = status.HTTP_204_NO_CONTENT)
    else:
        return Response(status = status.HTTP_400_BAD_REQUEST)

@api_view(['PUT'])
def update_resolution_to_in_progress(request):
    if 'resolution' in request.data:
        if Resolution.objects.filter(resolutionNumber__exact = request.data['resolution']).exists():
            resolution_obj = Resolution.objects.get(resolutionNumber__exact =request.data['resolution'])
            state_obj = ResolutionStates.objects.get(resolutionStateName__exact = 'In Progress')

            updated_resolution = UpdateResolutionSerializer.update(resolution_obj,state_obj)
            updated_resolution_serializer = UpdateResolutionSerializer(resolution_obj)
            service_obj = resolution_obj.get_service_obj()
            if drylab_config.EMAIL_USER_CONFIGURED :
                email_data = {}
                email_data['user_email'] = service_obj.get_user_email()
                email_data['user_name'] = service_obj.get_username()
                email_data['resolution_number'] = resolution_obj.get_resolution_number()
                send_resolution_in_progress_email(email_data)

            return Response(updated_resolution_serializer.data, status = status.HTTP_200_OK)

        else:
            return Response(status = status.HTTP_204_NO_CONTENT)
    else:
        return Response(status = status.HTTP_400_BAD_REQUEST)








'''
@api_view(['GET',])
def jobs_list(request):

   jobs = PipelineExternalDataJobs.objects.all()
   serializer = PipelineExternalDataJobsSerializer(jobs,many=True)
   return Response(serializer.data)

# view to load all the records from the model whose jobState is state
@api_view(['GET',])
def jobs_list_state(request,state):
   try:
       jobs = PipelineExternalDataJobs.objects.filter(jobState=state)
   except PipelineExternalDataJobs.DoesNotExist:
       return Response(status=status.HTTP_404__FOUND)
   #jobsstate = get_object_or_404(jobs, jobState=state).last()
   serializer = PipelineExternalDataJobsSerializer(jobs,many=True)
   return Response(serializer.data)
'''

@api_view(['GET'],)
def get_pipeline(request,pipeline):
   try:
      pipeline =  ParameterPipeline.objects.filter(id =pipeline)
   except ParameterPipeline.DoesNotExist:
      return Response(status=status.HTTP_404_FOUND)
   serializer = ParameterPipelineSerializer(pipeline,many=True)
   return Response(serializer.data)





# view to update the field passed in request, from the service passed by argument
@api_view(['PATCH',])
def api_update_state(request, service):
   try:
      pipejob = PipelineExternalDataJobs.objects.get(serviceRequestNumber=service)
   except PipelineExternalDataJobs.DoesNotExist:
      return Response({'message': 'Pipe does not exist'},status = status.HTTP_404_NOT_FOUND)

   #pipejob_data = JSONParser().parse(request)
   pipejob_serializer = PipelineExternalDataJobsBSerializer(pipejob,data=request.data)
   if pipejob_serializer.is_valid():
       pipejob_serializer.save()
       return Response(pipejob_serializer.data)
   return Response(pipejob_serializer.errors, status = status.HTTP_400_BAD_REQUEST)
