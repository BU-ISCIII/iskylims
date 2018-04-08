# -*- coding: utf-8 -*-
from django.shortcuts import get_object_or_404, render, redirect
from .models import *
from .forms import *
from .graphics import *
import os, re

from django.contrib.auth.decorators import login_required
from django.contrib.auth.models import Group
from django.conf import settings
from utils.fusioncharts.fusioncharts import FusionCharts

import datetime

####### Import libraries for static files
#from django.shortcuts import render_to_response
#from django.shortcuts import RequestContext

#pdb.set_trace()

@login_required
def index(request):
	if Service.objects.all().exclude(serviceStatus = 'delivered').exists():
		ongoing_services = Service.objects.all().exclude(serviceStatus = 'delivered').order_by('serviceCreatedOnDate')
		service_list = []
		for service in ongoing_services:
			service_info = []
			service_info.append(service.serviceRequestNumber)
			service_info.append(service.serviceStatus)
			if service.serviceStatus == 'recorded':
				service_delivery_date = 'Not defined yet'
			else:
				if Resolution.objects.filter(resolutionServiceID = service).exists():
					service_delivery_date = Resolution.objects.get(resolutionServiceID = service).last()
				else:
					service_delivery_date = 'Not defined yet'
			service_info.append(service_delivery_date)
			service_list.append(service_info)
		#import pdb; pdb.set_trace()
		return render(request, 'drylab/index.html',{'service_list': service_list})
	else:
		return render(request, 'drylab/index.html')

def increment_service_number ( user_name):
	# check user center
	#import pdb; pdb.set_trace()
	user_center = Profile.objects.get(profileUserID = user_name).profileCenter.centerAbbr
	# get latest service used for user's center
	if Service.objects.filter(serviceRequestNumber__icontains = user_center).exists():
		number_request = Service.objects.filter(serviceRequestNumber__icontains = user_center).last().serviceRequestNumber
		service_center= re.match('(^SRV[A-Z]+)(\d+)',number_request)
		last_sequence_number = int(service_center.group(2))
		new_service_number = str(last_sequence_number +1).zfill(3)
		service_number = service_center.group(1)  + new_service_number
	else:
		service_number = 'SRV'+ user_center + '001'
	return service_number 

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
				new_service.serviceRequestNumber = increment_service_number(request.user)
				#import pdb; pdb.set_trace()
				# Save the new instance
				new_service.save()
				# Save the many-to-many data for the form
				# (required for cases like this one where form.save(commit=False)) 
				
				form.save_m2m()
				#import pdb; pdb.set_trace()
				return render(request,'drylab/info_page.html',{'content':['Your service request has been successfully recorded.',
								'The sequence number assigned for your request is: ', new_service.serviceRequestNumber, 
								'Keep this number safe for refering your request','You will be contacted shortly.']})
		else: #No POST
			form = ServiceRequestFormInternalSequencing()
			## Addition of serviceProjectName for
			# implementation of drop down menu to choose a project name of a list of projects
			# belonging to the logged-in user in the service request form
			form.fields['serviceProjectNames'].queryset = Projects.objects.filter(user_id__exact = request.user.id)
			form.fields['serviceAvailableService'].queryset = AvailableService.objects.filter(availServiceDescription__exact="Genomic data analysis").get_descendants(include_self=True) 
			return render(request, 'drylab/RequestForm.html' , { 'form' : form , 'request_internal': 'request_internal'})
	

	if serviceRequestType == 'external_sequencing':
		if request.method == "POST":
			form = ServiceRequestFormExternalSequencing(data=request.POST,files=request.FILES)
			if form.is_valid():
				new_service = form.save(commit=False)
				new_service.serviceStatus = "recorded"
				new_service.serviceUserId = User.objects.get(id=request.user.id)
				new_service.serviceRequestNumber = increment_service_number(request.user)
				new_service.save()
				form.save_m2m()
				return render(request,'drylab/info_page.html',{'content':['Your service request has been successfully recorded.',
								'The sequence number assigned for your request is: ', new_service.serviceRequestNumber, 
								'Keep this number safe for refering your request','You will be contacted shortly.']})
		else: #No POST
			form = ServiceRequestFormExternalSequencing()
			form.fields['serviceAvailableService'].queryset = AvailableService.objects.filter(availServiceDescription__exact="Genomic data analysis").get_descendants(include_self=True) 
			return render(request, 'drylab/RequestForm.html' , { 'form' : form ,  'request_external': 'request_external' })


@login_required
def service_request_external_sequencing(request):
	if request.method == "POST":
		form = ServiceRequestFormExternalSequencing(data=request.POST,files=request.FILES)
		if form.is_valid():
			new_service = form.save(commit=False)
			new_service.serviceStatus = "recorded"
			new_service.serviceUserId = User.objects.get(id=request.user.id)
			new_service.serviceRequestNumber = increment_service_number(request.user)
			new_service.save()
			form.save_m2m()
			return render(request,'utils/info_page.html',{'content':['Your service request has been successfully recorded.',
								'The sequence number assigned for your request is: ', new_service.serviceRequestNumber, 
								'Keep this number safe for refering your request','You will be contacted shortly.']})
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
			new_service.serviceUserId = User.objects.get(id=request.user.id)
			new_service.serviceRequestNumber = increment_service_number(request.user)
			new_service.save()
			form.save_m2m()
			return render(request,'drylab/info_page.html',{'content':['Your service request has been successfully recorded.',
								'The sequence number assigned for your request is: ', new_service.serviceRequestNumber, 
								'Keep this number safe for refering your request','You will be contacted shortly.']})
		else:
			#import pdb; pdb.set_trace()
			return render(request,'drylab/error_page.html',{'content':['Your service request cannot be recorded.',
												'Check that all information is provided correctly.']})
	else:
		form = ServiceRequestForm_extended()

	form.fields['serviceAvailableService'].queryset = AvailableService.objects.filter(availServiceDescription__exact="Bioinformatics consulting and training").get_descendants(include_self=True)
	return render(request, 'drylab/RequestForm.html' , { 'form' : form ,  'consulting_request': 'consulting_request'})

@login_required
def infrastructure_request(request):
	if request.method == "POST":
		form = ServiceRequestForm_extended(data=request.POST or None,files=request.FILES)
		#import pdb; pdb.set_trace()
		if form.is_valid():
			new_service = form.save(commit=False)
			new_service.serviceStatus = "recorded"
			new_service.serviceUserId = User.objects.get(id=request.user.id)
			new_service.serviceRequestNumber = increment_service_number(request.user)
			#import pdb; pdb.set_trace()
			new_service.save()
			form.save_m2m()
			return render(request,'drylab/info_page.html',{'content':['Your service request has been successfully recorded.',
								'The sequence number assigned for your request is: ', new_service.serviceRequestNumber, 
								'Keep this number safe for refering your request','You will be contacted shortly.']})
		else:
			#import pdb; pdb.set_trace()
			return render(request,'drylab/error_page.html',{'content':['Your service request cannot be recorded.',
												'Check that all information is provided correctly.']})
	else:
		form = ServiceRequestForm_extended()

	form.fields['serviceAvailableService'].queryset = AvailableService.objects.filter(availServiceDescription__exact="User support").get_descendants(include_self=True)
	#pdb.set_trace()
	#form.helper[1].update_atrributes(hidden="true")
	return render(request, 'drylab/RequestForm.html' , { 'form' : form , 'infrastructure_request': 'infrastructure_request'})

def get_service_information (service_id):
	service= Service.objects.get(pk=service_id)
	display_service_details = {}
	text_for_dates = ['Service Date Creation', 'Approval Service Date', 'Rejected Service Date']
	service_dates = []
	display_service_details['service_name'] = service.serviceRequestNumber
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
	display_service_details['service_notes'] = service.serviceNotes
	dates_for_services = service.get_service_dates()
	for i in range(len(dates_for_services)):
		service_dates.append([text_for_dates[i],dates_for_services[i]])
	display_service_details['service_dates'] = service_dates
	
	# get all services
	display_service_details['nodes']= service.serviceAvailableService.all()
	
	if service.serviceStatus == 'recorded':
		display_service_details['allowed_action_appr_reject'] = True
	elif service.serviceStatus == 'approved':
		display_service_details['allowed_action_queued'] = True
	elif service.serviceStatus == 'queued':
		display_service_details['allowed_action_inProgress'] = True
	else:
		display_service_details['not_allowed_actions'] = True

	return display_service_details

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
		
		
		if request.method == 'POST':
			#import pdb; pdb.set_trace()
			service= Service.objects.get(pk=service_id)
			if request.POST['action'] == 'action_appr_reject' and 'inlineRadioOptions' in request.POST :
				
				if request.POST['inlineRadioOptions'] == 'accepted':
					# update state and dates for accepted service
					service.serviceStatus = 'approved'
					service.serviceOnApprovedDate = datetime.date.today()
				else:
					service.serviceStatus = 'rejected'
					service.serviceOnRejectedDate = datetime.date.today()
				service.save()
			elif request.POST['action'] == 'action_queued' and 'inlineRadioOptions' in request.POST:  
				if request.POST['inlineRadioOptions'] == 'queued':
					# update state and dates for accepted service
					service.serviceStatus = 'queued'
					service.save()
		# displays the service information with the latest changes done using the forms
		display_service_details = get_service_information(service_id)
			
			
		#import pdb; pdb.set_trace()
		return render (request,'drylab/display_service.html',{'display_service': display_service_details})
	else:
		return render (request,'drylab/error_page.html', {'content':['The service that you are trying to get does not exist ','Contact with your administrator .']})


@login_required
def search_service (request):
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
	if request.method == 'POST' and request.POST['action'] == 'searchservice':
		
		
		service_number_request = request.POST['servicenumber']
		service_state = request.POST['servicestate']
		start_date=request.POST['startdate']
		end_date=request.POST['enddate']
		user_name = request.POST['username']
		if service_number_request == '' and service_state == '' and start_date == '' and end_date == '' and user_name =='':
			 return render( request,'drylab/searchService.html',{'services_state_list':STATUS_CHOICES})
		### check the right format of start and end date
		if start_date != '':
			try:
				datetime.datetime.strptime(start_date, '%Y-%m-%d')
			except:
				return render (request,'drylab/error_page.html', {'content':['The format for the "Start Date Search" Field is incorrect ',
																			'ADVICE:', 'Use the format  (DD-MM-YYYY)']})
		if end_date != '':
			try:
				datetime.datetime.strptime(end_date, '%Y-%m-%d')
			except:
				return render (request,'drylab/error_page.html', {'content':['The format for the "End Date Search" Field is incorrect ',
																			'ADVICE:', 'Use the format  (DD-MM-YYYY)']})
		if service_number_request == '' and service_state == '':
			services_found = Service.objects.all()

		if service_number_request != '':
			# check if the requested service in the form matches exactly with the existing service in DB
			if Service.objects.filter(serviceRequestNumber__exact = service_number_request).exists():
				#import pdb; pdb.set_trace()
				services_found = Service.objects.get(serviceRequestNumber__exact = service_number_request)
				redirect_page = '/drylab/display_service=' + str(services_found.id)
				return redirect (redirect_page)
			if Service.objects.filter(serviceRequestNumber__icontains = service_number_request).exists():
				services_found = Service.objects.filter(serviceRequestNumber__icontains = service_number_request)
			else:
				return render (request,'drylab/error_page.html', {'content':['No matches have been found for the service number ', service_number_request ]})

		if service_state != '':
			if service_number_request =='':
				services_found = Service.objects.all()
			if services_found.filter(serviceStatus__exact = service_state).exists():
				services_found = services_found.filter(serviceStatus__exact = service_state)
			else:
				return render (request,'drylab/error_page.html', {'content':['No matches have been found for the service number in state', service_state ]})
		if start_date !='' and end_date != '':
			if services_found.filter(serviceCreatedOnDate__range=(start_date, end_date)).exists():
				 services_found = services_found.filter(serviceCreatedOnDate__range=(start_date, end_date))
			else:
				return render (request,'drylab/error_page.html', {'content':['There are no services containing ', service_number,
														' created between ', start_date, 'and the ', end_date]})
		if start_date !='' and end_date == '':
			if services_found.filter(serviceCreatedOnDate__gte = start_date).exists():
				services_found = services_found.filter(serviceCreatedOnDate__gte = start_date)
			else:
				return render (request,'drylab/error_page.html', {'content':['There are no services containing ', service_number,
														' created before ', start_date]})
		if start_date =='' and end_date != '':
			if services_found.filter(serviceCreatedOnDate__lte = end_date).exists():
				services_found = services_found.filter(serviceCreatedOnDate__lte = end_date)
			else:
				return render (request,'drylab/error_page.html', {'content':['There are no services containing ', service_number,
														' finish before ', end_date]})
		
		#If only 1 service mathes the user conditions, then get the user information
		if len(services_found) == 1 :
			redirect_page = '/drylab/display_service=' + str(services_found[0].id)
			return redirect (redirect_page)
		else:
			display_multiple_services ={}
			s_list  = {}
			for service_item in services_found:
				service_id = service_item.id
				service_number = service_item.serviceRequestNumber
				service_status = service_item.serviceStatus
				service_center = service_item.serviceSeqCenter
				s_list [service_id]=[[service_number, service_status, service_center]]
			display_multiple_services['s_list'] = s_list
			return render (request,'drylab/searchService.html', {'display_multiple_services': display_multiple_services})
	#import pdb; pdb.set_trace()
	return render( request,'drylab/searchService.html',{'services_state_list':STATUS_CHOICES})


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
			approved[services.id]= [services.get_service_information().split(';')]
		pending_services_details['approved'] = approved
	if Service.objects.filter(serviceStatus__exact = 'queued').exists():
		services_in_queued = Service.objects.filter(serviceStatus__exact = 'queued').order_by('-serviceCreatedOnDate')
		for services in services_in_queued:
			queued[services.id]= [services.get_service_information().split(';')]
		pending_services_details['queued'] = queued
	if Service.objects.filter(serviceStatus__exact = 'in_progress').exists():
		services_in_progress = Service.objects.filter(serviceStatus__exact = 'in_progress').order_by('-serviceCreatedOnDate')
		for services in services_in_progress:
			in_progress[services.id]= [services.get_service_information().split(';')]
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

@login_required
def add_resolution (request, service_id):
#def add_resolution (request):
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
		
	if request.method == "POST":
		
		form = addResolutionService(data=request.POST)
		import pdb ; pdb.set_trace()
		if form.is_valid():
			
			new_resolution = form.save(commit=False)
			#service_id =  new_resolution.
			#new_resolution.serviceStatus = "queued"
			#new_service.serviceUserId = User.objects.get(id=request.user.id)
			#new_service.serviceRequestID = increment_service_number(request.user.id)
			#new_resolution.save()
			#form.save_m2m()
			return render(request,'utils/info_page.html',{'content':['Your resolution proposal has been successfully recorded.','You will be contacted shortly.']})
	else:
		if Service.objects.filter(pk=service_id).exists():
			service_id= Service.objects.get(pk=service_id)
			form = addResolutionService()
			#form.fields['resolutionServiceID'].queryset = service_id
			import pdb ; pdb.set_trace()

			return render(request, 'drylab/addResolution.html' , { 'form' : form ,'prueba':'pepe'})

def get_current_users():
	from django.contrib.sessions.models import Session
	from django.utils import timezone
	active_sessions = Session.objects.filter(expire_date__gte=timezone.now())
	user_id_list = []
	for session in active_sessions:
		data = session.get_decoded()
		user_id_list.append(data.get('_auth_user_id', None))
	# Query all logged in users based on id list
	return User.objects.filter(id__in=user_id_list)
    
'''
Define in settings.py 
SESSION_EXPIRE_AT_BROWSER_CLOSE = True

'''
@login_required
def open_sessions (request):
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
	user_connected = {}
	if get_current_users().exists():
		user_list_connected = get_current_users()
		user_data = []
		for user in user_list_connected:
			user_data.append([user.username, user.first_name, user.last_name, user.email])

		user_connected['user_data']= user_data
			
		user_connected['number_of_users'] = user_list_connected.count()
		import pdb ; pdb.set_trace()
	return render (request, 'drylab/openSessions.html', {'user_connected': user_connected })

@login_required
def user_login (request):
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
	user_data = []
	login_data = {}
	user_list = User.objects.all().order_by('-last_login')
	for user in user_list:
		user_data.append([user.username, user.first_name, user.last_name, user.email, user.last_login])
	login_data['user_data'] = user_data

	return render(request, 'drylab/userLogin.html', {'login_data': login_data})

