from django.shortcuts import render
from django.shortcuts import get_object_or_404
from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response
from iSkyLIMS_drylab.models import PipelineExternalDataJobs
from iSkyLIMS_drylab.models import ParameterPipeline
from iSkyLIMS_wetlab.models import SamplesInProject
from iSkyLIMS_wetlab.models import Projects
from iSkyLIMS_drylab.models import Service
from iSkyLIMS_drylab.models import RequestedProjectInServices
from .serializers import PipelineExternalDataJobsSerializer
from .serializers import PipelineExternalDataJobsBSerializer
from .serializers import ParameterPipelineSerializer
from .serializers import ServiceSerializer

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
def service_list(request):
    service = Service.objects.all()
    serializer = ServiceSerializer(service, many=True)
    return Response(serializer.data)

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

@api_view(['GET'],)
def get_pipeline(request,pipeline):
   try:
      pipeline =  ParameterPipeline.objects.filter(id =pipeline)
   except ParameterPipeline.DoesNotExist:
      return Response(status=status.HTTP_404_FOUND)
   serializer = ParameterPipelineSerializer(pipeline,many=True)
   return Response(serializer.data)

@api_view(['GET'],)
def get_samplesinproject(request,service):
   #serv = Service.objects.get(serviceRequestNumber__exact = service)
   #service_id = serv.id
   #project_list = RequestedProjectInServices.objects.filter(projectService = service_id)
   #projects_id=[]
   #for project in project_list:
 #   projects_id.append(str(project.get_requested_external_project_id()))
   projects_id = get_projectsid(service)
   samples =  get_samples_projects(projects_id)

   #serializer = samples
   return Response(samples)

@api_view(['GET'],)
def get_runfolder(request,project):
    runfolder =  get_run_folder_from_user_project(project)
    return Response(runfolder)

#view to update the entire record
@api_view(['PUT',])
def api_update_job(request, service):
   try:
      pipejob = PipelineExternalDataJobs.objects.get(serviceRequestNumber=service)
   except PipelineExternalDataJobs.DoesNotExist:
      return Response({'message': 'Pipe does not exist'},status = status.HTTP_404_NOT_FOUND)

   #pipejob_data = JSONParser().parse(request)
   pipejob_serializer = PipelineExternalDataJobsSerializer(pipejob,data=request.data)
   if pipejob_serializer.is_valid():
       pipejob_serializer.save()
       return Response(pipejob_serializer.data)
   return Response(pipejob_serializer.errors, status = status.HTTP_400_BAD_REQUEST)

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
