from django.shortcuts import render
from django.shortcuts import get_object_or_404
from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response
from iSkyLIMS_drylab.models import PipelineExternalDataJobs
from .serializers import PipelineExternalDataJobsSerializer
from .serializers import PipelineExternalDataJobsBSerializer
# Create your views here.
#class PipelinesViewSet(viewsets.ModelViewSet):
#    serializer_class = PipelinesSerializer
#    queryset = Pipelines.objects.all()

@api_view(['GET',])
def jobs_list(request):

   jobs = PipelineExternalDataJobs.objects.all()
   serializer = PipelineExternalDataJobsSerializer(jobs,many=True)
   return Response(serializer.data)

@api_view(['GET',])
def jobs_list_state(request,state):
   try:
       jobs = PipelineExternalDataJobs.objects.filter(jobState=state)
   except PipelineExternalDataJobs.DoesNotExist:
       return Response(status=status.HTTP_404__FOUND)
   #jobsstate = get_object_or_404(jobs, jobState=state).last()
   serializer = PipelineExternalDataJobsSerializer(jobs,many=True)
   return Response(serializer.data)

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
