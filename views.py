# -*- coding: utf-8 -*-
from django.shortcuts import get_object_or_404, render, redirect
from .models import *
from .forms import *
from .graphics import *
import os

from django.contrib.auth.decorators import login_required
from django.contrib.auth.models import Group
from django.conf import settings
from utils.fusioncharts.fusioncharts import FusionCharts

####### Import libraries for static files
#from django.shortcuts import render_to_response
#from django.shortcuts import RequestContext

#pdb.set_trace()

@login_required
def index(request):
    return render(request, 'drylab/index.html')

@login_required
def service_request(request, serviceRequestType):
	if serviceRequestType == 'internal_sequencing':
		if request.method == "POST":
			##Create a form instance with POST data
			# more information in: 
			# https://docs.djangoproject.com/en/2.0/topics/forms/modelforms/
			form = ServiceRequestFormInternalSequencing(data=request.POST,files=request.FILES)
			#import pdb; pdb.set_trace()
			if form.is_valid():
				# Create but dont save for following modif		
				new_service = form.save(commit=False)
				# Modification of 'serviceStatus'
				new_service.serviceStatus = "recorded"
				# Local ("in these premises :-)") sequencing
				new_service.serviceSeqCenter = "GENOMIC_SEQ_UNIT"
				# Previous 'serviceUsername'refactored to 'serviceUserId' which describes better its real nature
				new_service.serviceUserId = User.objects.get(id=request.user.id)
				#import pdb; pdb.set_trace()
				# Save the new instance
				new_service.save()
				# Save the many-to-many data for the form
				# (required for cases like this one where form.save(commit=False)) 
				
				form.save_m2m()
				return render(request,'drylab/info_page.html',{'content':['Your service request has been successfully recorded.','You will be contacted shortly.']})
		else: #No POST
			form = ServiceRequestFormInternalSequencing()
			## Addition of serviceProjectName for
			# implementation of drop down menu to choose a project name of a list of projects
			# belonging to the logged-in user in the service request form
			form.fields['serviceProjectNames'].queryset = Projects.objects.filter(user_id__exact = request.user.id)
			form.fields['serviceAvailableService'].queryset = AvailableService.objects.filter(availServiceDescription__exact="Genomic data analysis").get_descendants(include_self=True) 
			return render(request, 'drylab/RequestForm.html' , { 'form' : form })
	

	if serviceRequestType == 'external_sequencing':
		if request.method == "POST":
			form = ServiceRequestFormExternalSequencing(data=request.POST,files=request.FILES)
			if form.is_valid():
				new_service = form.save(commit=False)
				new_service.serviceStatus = "recorded"
				new_service.serviceUserId = User.objects.get(id=request.user.id)
				new_service.save()
				form.save_m2m()
				return render(request,'drylab/info_page.html',{'content':['Your service request has been successfully recorded.','You will be contacted shortly.']})
		else: #No POST
			form = ServiceRequestFormExternalSequencing()
			form.fields['serviceAvailableService'].queryset = AvailableService.objects.filter(availServiceDescription__exact="Genomic data analysis").get_descendants(include_self=True) 
			return render(request, 'drylab/RequestForm.html' , { 'form' : form })


@login_required
def service_request_external_sequencing(request):
	if request.method == "POST":
		form = ServiceRequestFormExternalSequencing(data=request.POST,files=request.FILES)
		if form.is_valid():
			new_service = form.save(commit=False)
			new_service.serviceStatus = "recorded"
			new_service.serviceUserId = User.objects.get(id=request.user.id)
			new_service.save()
			form.save_m2m()
			return render(request,'utils/info_page.html',{'content':['Your service request has been successfully recorded.','You will be contacted shortly.']})
	else:
		form = ServiceRequestFormExternalSequencing()
	
	form.fields['serviceAvailableService'].queryset = AvailableService.objects.filter(availServiceDescription__exact="Genomic data analysis").get_descendants(include_self=True) 
	return render(request, 'drylab/RequestForm.html' , { 'form' : form })


@login_required
def counseling_request(request):
	if request.method == "POST":
		form = ServiceRequestForm_extended(data=request.POST,files=request.FILES)

		if form.is_valid():
			new_service = form.save(commit=False)
			new_service.serviceStatus = "recorded"
			new_service.save()
			form.save_m2m()
			return render(request,'drylab/info_page.html',{'content':['Your service request has been successfully recorded.','You will be contacted shortly.']})
		else:
			#import pdb; pdb.set_trace()
			return render(request,'drylab/error_page.html',{'content':['Your service request cannot be recorded.',
												'Check that all information is provided correctly.']})
	else:
		form = ServiceRequestForm_extended()

	form.fields['serviceAvailableService'].queryset = AvailableService.objects.filter(availServiceDescription__exact="Bioinformatics consulting and training").get_descendants(include_self=True)
	return render(request, 'drylab/RequestForm.html' , { 'form' : form })

@login_required
def infrastructure_request(request):
	if request.method == "POST":
		form = ServiceRequestForm_extended(data=request.POST or None,files=request.FILES)
		#import pdb; pdb.set_trace()
		if form.is_valid():
			new_service = form.save(commit=False)
			new_service.serviceStatus = "recorded"
			#import pdb; pdb.set_trace()
			new_service.save()
			form.save_m2m()
			return render(request,'drylab/info_page.html',{'content':['Your service request has been successfully recorded.','You will be contacted shortly.']})
		else:
			#import pdb; pdb.set_trace()
			return render(request,'drylab/error_page.html',{'content':['Your service request cannot be recorded.',
												'Check that all information is provided correctly.']})
	else:
		form = ServiceRequestForm_extended()

	form.fields['serviceAvailableService'].queryset = AvailableService.objects.filter(availServiceDescription__exact="User support").get_descendants(include_self=True)
	#pdb.set_trace()
	#form.helper[1].update_atrributes(hidden="true")
	return render(request, 'drylab/RequestForm.html' , { 'form' : form })


@login_required
def display_service (request, service_id):
	if request.user.is_authenticated:
		try:
			groups = Group.objects.get(name='Admin_iSkyLIMS')
			if groups not in request.user.groups.all():
				return render (request,'drylab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
		except:
			return render (request,'drylab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
	else:
		#redirect to login webpage
		return redirect ('/accounts/login')
	if Service.objects.filter(pk=service_id).exists():
		service= Service.objects.get(pk=service_id)
		display_service_details = {}
		display_service_details['service_name'] = 'PRUEBA--CAMBIAR'
		# get the list of projects
		projects_in_service = {}
		projects_class = service.serviceProjectNames.all()
		for project in projects_class:
			project_id = project.id
			projects_in_service[project_id]=project.get_project_name()
		display_service_details['projects'] = projects_in_service
		display_service_details['user_name'] = service.serviceUserId.username
		display_service_details['file'] = os.path.join(settings.MEDIA_URL,str(service.serviceFile))
		display_service_details['state'] = service.serviceStatus
		display_service_details['service_dates'] = service.get_service_dates()
		
		# get all services 
		req_services_class = service.serviceAvailableService.all()
		# get only the level 2 services
		#service.serviceAvailableService.filter(level = 2)
		req_service_names = []
		for req_service in req_services_class:
			req_service_names.append(req_service.availServiceDescription)
		display_service_details['requested_services'] = req_services_class
		
		
		import pdb; pdb.set_trace()
		return render (request,'drylab/display_service.html',{'display_service': display_service_details})
	else:
		return render (request,'drylab/error_page.html', {'content':['The service that you are trying to get does not exist ','Contact with your administrator .']})



@login_required
def pending_services (request):
	if request.user.is_authenticated:
		try:
			groups = Group.objects.get(name='Admin_iSkyLIMS')
			if groups not in request.user.groups.all():
				return render (request,'drylab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
		except:
			return render (request,'drylab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
	else:
		#redirect to login webpage
		return redirect ('/accounts/login')

	pending_services_details = {}
	recorded, approved, queued, in_progress = {}, {}, {}, {}
	if Service.objects.filter(serviceStatus__exact = 'recorded').exists():
		services_in_request = Service.objects.filter(serviceStatus__exact = 'recorded').order_by('-serviceCreatedOnDate')
		for services in services_in_request:
			recorded[services.id]= [services.get_service_information().split(';')]
		pending_services_details['recorded'] = recorded
	if Service.objects.filter(serviceStatus__exact = 'approved').exists():
		services_in_approved = Service.objects.filter(serviceStatus__exact = 'approved').order_by('-serviceCreatedOnDate')
		for services in services_in_approved:
			approved[services.id]= services.get_service_information().split(';')
		pending_services_details['approved'] = approved
	if Service.objects.filter(serviceStatus__exact = 'queued').exists():
		services_in_queued = Service.objects.filter(serviceStatus__exact = 'queued').order_by('-serviceCreatedOnDate')
		for services in services_in_queued:
			queued[services.id]= services.get_service_information().split(';')
		pending_services_details['queued'] = queued
	if Service.objects.filter(serviceStatus__exact = 'in_progress').exists():
		services_in_progress = Service.objects.filter(serviceStatus__exact = 'in_progress').order_by('-serviceCreatedOnDate')
		for services in services_in_progress:
			in_progress[services.id]= services.get_service_information().split(';')
		pending_services_details['in_progress'] = in_progress

	number_of_services = {}
	number_of_services ['RECORDED'] = len (recorded)
	number_of_services ['APPROVED'] = len (approved)
	number_of_services ['QUEUED'] = len (queued)
	number_of_services ['IN PROGRESS'] = len (in_progress)
	data_source = graphic_3D_pie('Number of Pending Services', '', '', '','fint',number_of_services)
	graphic_pending_services = FusionCharts("pie3d", "ex1" , "425", "350", "chart-1", "json", data_source)
	pending_services_details ['graphic_pending_services'] = graphic_pending_services.render()
	
	#import pdb ; pdb.set_trace()
	
	return render (request, 'drylab/pendingServices.html', {'pending_services': pending_services_details})
