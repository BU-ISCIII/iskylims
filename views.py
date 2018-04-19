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
from django.core.mail import send_mail

import datetime
import statistics

####### Import libraries for static files
#from django.shortcuts import render_to_response
#from django.shortcuts import RequestContext

#pdb.set_trace()

@login_required
def index(request):
	if Service.objects.all().exclude(serviceStatus = 'delivered').exclude(serviceStatus = 'rejected').exists():
		ongoing_services = Service.objects.all().exclude(serviceStatus = 'delivered').exclude(serviceStatus = 'rejected').order_by('serviceCreatedOnDate')
		service_list = []
		for service in ongoing_services:
			service_info = []
			service_info.append(service.serviceRequestNumber)
			service_info.append(service.serviceStatus)
			if service.serviceStatus == 'recorded':
				service_delivery_date = 'Not defined yet'
			else:
				if Resolution.objects.filter(resolutionServiceID = service).exists():
					#import pdb; pdb.set_trace()
					service_delivery_date = Resolution.objects.filter(resolutionServiceID = service).last().resolutionEstimatedDate.strftime("%d %B, %Y")
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
				## Send mail to user and bioinfo admins
				subject = 'Service ' + new_service.serviceRequestNumber + " has been recorded"
				body_message = 'Dear ' + request.user.username + "\n Your service " + new_service.serviceRequestNumber + " has been recorded. You will recieved the resolution of the request as soon as possible.\n Kind regards \n BU-ISCIII \n bioinformatica@isciii.es"
				from_user = 'bioinformatica@isciii.es'
				to_user = [request.user.email,'bioinformatica@isciii.es']
				send_mail (subject, body_message, from_user, to_user)

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
				## Send email
				subject = 'Service ' + new_service.serviceRequestNumber + " has been recorded"
				body_message = 'Dear ' + request.user.username + "\n Your service " + new_service.serviceRequestNumber + " has been recorded. You will recieved the resolution of the request as soon as possible.\n Kind regards \n BU-ISCIII \n bioinformatica@isciii.es"
				from_user = 'bioinformatica@isciii.es'
				to_user = [request.user.email,'bioinformatica@isciii.es']
				send_mail (subject, body_message, from_user, to_user)

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
			## Send email
			subject = 'Service ' + new_service.serviceRequestNumber + " has been recorded"
			body_message = 'Dear ' + request.user.username + "\n Your service " + new_service.serviceRequestNumber + " has been recorded. You will recieved the resolution of the request as soon as possible.\n Kind regards \n BU-ISCIII \n bioinformatica@isciii.es"
			from_user = 'bioinformatica@isciii.es'
			to_user = [request.user.email,'bioinformatica@isciii.es']
			send_mail (subject, body_message, from_user, to_user)
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
			## Send email
			subject = 'Service ' + new_service.serviceRequestNumber + " has been recorded"
			body_message = 'Dear ' + request.user.username + "\n Your service " + new_service.serviceRequestNumber + " has been recorded. You will recieved the resolution of the request as soon as possible.\n Kind regards \n BU-ISCIII \n bioinformatica@isciii.es"
			from_user = 'bioinformatica@isciii.es'
			to_user = [request.user.email,'bioinformatica@isciii.es']
			send_mail (subject, body_message, from_user, to_user)
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
	if service.serviceStatus != 'approved'and service.serviceStatus != 'recorded':
		# get the proposal for the delivery date
		#import pdb; pdb.set_trace()
		resolution_folder = Resolution.objects.filter(resolutionServiceID = service).last().resolutionFullNumber
		display_service_details['resolution_folder'] = resolution_folder
		resolution_date = Resolution.objects.filter(resolutionServiceID = service).last().resolutionEstimatedDate
		display_service_details['estimated_delivery_date'] = resolution_date

	# get all services
	display_service_details['nodes']= service.serviceAvailableService.all()
	# adding actions fields
	if service.serviceStatus != 'rejected' or service.serviceStatus != 'archived':
		display_service_details['add_resolution_action'] = service_id
	if service.serviceStatus == 'queued':
		resolution_id = Resolution.objects.filter(resolutionServiceID = service).last().id
		display_service_details['add_in_progress_action'] = resolution_id
	if service.serviceStatus == 'in_progress':
		resolution_id = Resolution.objects.filter(resolutionServiceID = service).last().id
		display_service_details['add_delivery_action'] = resolution_id


	if Resolution.objects.filter(resolutionServiceID = service).exists():
		resolution_list = Resolution.objects.filter(resolutionServiceID = service)
		resolution_info =[]
		for resolution_item in resolution_list :
			resolution_info.append([resolution_item.get_resolution_information()])
		display_service_details['resolutions'] = resolution_info
	#import pdb; pdb.set_trace()
	if Resolution.objects.filter(resolutionServiceID = service).exists():
		resolution_id = Resolution.objects.filter(resolutionServiceID = service).last().id
		if Delivery.objects.filter(deliveryResolutionID = resolution_id).exists():
			delivery = Delivery.objects.get(deliveryResolutionID = resolution_id)
			display_service_details['delivery'] = [delivery.get_delivery_information()]

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

		form = AddResolutionService(data=request.POST)
		#import pdb ; pdb.set_trace()
		if form.is_valid():
			service_acepted_rejected = request.POST['radio_buttons']
			new_resolution = form.save(commit=False)
			#import pdb ; pdb.set_trace()
			service_reference = Service.objects.get(pk=service_id)
			if service_acepted_rejected == 'accepted':
				service_reference.serviceStatus = "approved"
				service_reference.serviceOnApprovedDate = datetime.date.today()
				number_list = []
				number_list.append(str(service_reference.serviceRequestNumber))
				number_list.append(str(datetime.date.today()).replace('-',''))
				number_list.append(new_resolution.resolutionFullNumber)
				number_list.append(str(service_reference.serviceUserId))
				number_list.append('S')
				new_resolution.resolutionFullNumber = '_'.join(number_list)
			else:
				service_reference.serviceStatus = "rejected"
				service_reference.serviceOnRejectedDate = datetime.date.today()
			service_reference.save()
			if Resolution.objects.filter(resolutionServiceID = service_reference).exists():
				resolution_count =  Resolution.objects.filter(resolutionServiceID = service_reference).count()
				resolution_number = str(service_reference.serviceRequestNumber) + '.' + str(resolution_count +1 )
			else:
				resolution_number = str(service_reference.serviceRequestNumber) + '.1'
			new_resolution.resolutionServiceID = service_reference
			new_resolution.resolutionNumber = resolution_number
			new_resolution.resolutionOnQueuedDate = datetime.date.today()
			if service_reference.serviceStatus == "approved":
				service_reference.serviceStatus = "queued"
				service_reference.save()
			new_resolution.save()
			form.save_m2m()
			## Send email
			service_user_mail = service_reference.serviceUserId.email
			subject = 'Service ' + service_reference.serviceRequestNumber + " has been updated"
			if service_acepted_rejected == "accepted":
			    body_message = 'Dear ' + service_reference.serviceUserId.username + "\n A new resolution has been added for your service: " + resolution_number + "\n. Your service has been "+ service_reference.serviceStatus + " and your delivery estimated date is " + new_resolution.resolutionEstimatedDate.strftime('%d %B %Y') + ".\n Your service is now queued and you will be notified when it is updated. \n Kind regards \n BU-ISCIII \n bioinformatica@isciii.es"
			else:
			    body_message= 'Dear ' + service_reference.serviceUserId.username + "\n A new resolution has been added for your service: " + resolution_number + "\n. Your service  has been "+ service_reference.serviceStatus + " because it does not fullfil our requirements or is not in our services portfolio. If you have any question please contact us. \n Kind regards \n BU-ISCIII \n bioinformatica@isciii.es"
			from_user = 'bioinformatica@isciii.es'
			to_user = [service_user_mail,'bioinformatica@isciii.es']
			send_mail (subject, body_message, from_user, to_user)
			#import pdb ; pdb.set_trace()
			return render(request,'drylab/info_page.html',{'content':['Your resolution proposal has been successfully recorded with Resolution Number.', resolution_number]})
	else:
		if Service.objects.filter(pk=service_id).exists():
			service_id= Service.objects.get(pk=service_id)
			service_number = service_id.serviceRequestNumber
			form = AddResolutionService()
			#import pdb ; pdb.set_trace()

			return render(request, 'drylab/addResolution.html' , { 'form' : form ,'prueba':'pepe'})
@login_required
def add_in_progress (request, resolution_id):
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

	if Resolution.objects.filter(pk = resolution_id).exists():
		resolution = Resolution.objects.get(pk = resolution_id)
		service_to_update = resolution.resolutionServiceID
		# update the service status and in_porgress date
		#import pdb ; pdb.set_trace()
		service_to_update.serviceStatus = 'in_progress'
		service_to_update.save()
		resolution.resolutionOnInProgressDate = datetime.date.today()
		resolution.save()
		service_user_mail = service_to_update.serviceUserId.email
		## Send email
		subject = 'Service ' + service_to_update.serviceRequestNumber + " has been updated"
		body_message = 'Dear ' + service_to_update.serviceUserId.username + "\n Your service with resolution id: " + resolution.resolutionNumber + " is now in progress." + "\n Kind regards \n BU-ISCIII \n bioinformatica@isciii.es"
		from_user = 'bioinformatica@isciii.es'
		to_user = [service_user_mail,'bioinformatica@isciii.es']
		send_mail (subject, body_message, from_user, to_user)
		return render (request,'drylab/info_page.html',{'content':['Your resolution  request ', resolution.resolutionNumber,
								'has been successfully upated to In Progress state']})
	else:
		#import pdb ; pdb.set_trace()
		return render (request,'drylab/error_page.html', {'content':['The resolution that you are trying to upadate does not exists ','Contact with your administrator .']})
	return

@login_required
def add_delivery (request , resolution_id):
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
	if request.method == 'POST' :
		form = AddDeliveryService(data=request.POST)
		if form.is_valid():
			#import pdb ; pdb.set_trace()
			resolution_id = Resolution.objects.get(pk = resolution_id)
			new_delivery = form.save(commit=False)
			new_delivery.deliveryDate = datetime.date.today()
			new_delivery.deliveryResolutionID = resolution_id
			new_delivery.save()
			form.save_m2m()
			# Update the status service to delivery
			service_id = resolution_id.resolutionServiceID
			service_id.serviceStatus = 'delivered'
			service_id.serviceOnDeliveredDate = datetime.date.today()
			service_id.save()
			service_user_mail = service_id.serviceUserId.email
			## Send email
			subject = 'Service ' + service_id.serviceRequestNumber + " has been updated"
			body_message = 'Dear ' + service_id.serviceUserId.username + "\n. Your service with resolution id: " + resolution_id.resolutionNumber + " is finished. A mail with instructions for downloading the results will be shortly sent to you." + "\n Kind regards \n BU-ISCIII \n bioinformatica@isciii.es"
			from_user = 'bioinformatica@isciii.es'
			to_user = [service_user_mail,'bioinformatica@isciii.es']
			send_mail (subject, body_message, from_user, to_user)
			return render(request,'drylab/info_page.html',{'content':['The service is now on Delivery status ']})
	else:
		if Resolution.objects.filter(pk = resolution_id).exists():
			#import pdb ; pdb.set_trace()
			form = AddDeliveryService()
			delivery_info = {}
			return render (request, 'drylab/addDelivery.html', {'form':form, 'delivery_info': delivery_info})
		else:
			#import pdb ; pdb.set_trace()
			return render (request,'drylab/error_page.html', {'content':['The resolution that you are trying to upadate does not exists ','Contact with your administrator .']})
	return

@login_required
def stats_by_date_user (request):
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
	if request.method == 'POST':
		form = ByDateUserStats(data=request.POST)
		if form.is_valid():
			# validate the input data in the form
			user_name = form['user_name'].data
			start_date = form['start_date'].data
			end_date = form['end_date'].data
			if User.objects.filter(username__icontains = user_name).exists():
				matched_names = User.objects.filter(username__icontains = user_name)
				if len(matched_names) > 1:
					name_list =[]
					for names in matched_names:
						name_list.append(names.username)
					name_string = '  ,  '.join(name_list)
					return render (request,'drylab/error_page.html', {'content':['Too many matches have been found for the user name field', user_name,
																				'ADVICE:', 'Please write down one of the following user name and repeate again the search',
																				name_string,]})
				else:
					user_name_id = matched_names[0].id
					user_name = matched_names[0].username
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
			#import pdb ; pdb.set_trace()
			services_user = Service.objects.filter(serviceUserId__exact = user_name_id).order_by('-serviceRequestNumber')
			if start_date != '' and end_date !='':
				if services_user.filter(serviceCreatedOnDate__range=(start_date,end_date)).exists():
					services_user = services_user.filter(serviceCreatedOnDate__range(start_date,end_date))
				else:
					return render (request,'drylab/error_page.html', {'content':['There are no services created by ', user_name , 'For the time of period of between:',
																start_date , 'and', end_date]})
			if start_date !='' and end_date == '':
				if services_user.filter(serviceCreatedOnDate__gte = end_date).exists():
					services_user = services_user.filter(serviceCreatedOnDate__lte = end_date)
				else:
					return render (request,'drylab/error_page.html', {'content':['There are no services created by ', user_name , 'Starting from ', start_date ]})
			if start_date =='' and end_date != '':
				if services_user.filter(serviceCreatedOnDate__lte = end_date).exists():
					services_user = services_user.filter(serviceCreatedOnDate__lte = end_date)
				else:
					return render (request,'drylab/error_page.html', {'content':['There are no services created by ', user_name , 'Finish before ', end_date ]})


			stats_info = {}
			service_by_user=[]
			for service_item in services_user:
				service_by_user.append(service_item.get_stats_information())

			#import pdb ; pdb.set_trace()
			stats_info ['user_name'] = user_name
			stats_info ['service_by_user'] = service_by_user

			# perform calculation time media delivery for user
			if services_user.filter(serviceStatus__exact = 'delivered').exists():
				delivery_services = services_user.filter(serviceStatus__exact = 'delivered')

				delivery_time_in_days = []
				for service_item in delivery_services :
					delivery_time_in_days.append(int (service_item.get_time_to_delivery()))

				stats_info['time_mean_for_user']=  format(statistics.mean (delivery_time_in_days), '.2f')

			else:
				# there are not delivery services for the user in the specified period of time
				pass
			# preparing graphic for status of the services
			number_of_services = {}
			if services_user.filter(serviceStatus__exact = 'recorded').exists():
				number_of_services ['RECORDED'] = len (services_user.filter(serviceStatus__exact = 'recorded'))
			else:
				number_of_services ['RECORDED'] = 0
			if services_user.filter(serviceStatus__exact = 'queued').exists():
				number_of_services ['QUEUED'] = len (services_user.filter(serviceStatus__exact = 'queued'))
			else:
				number_of_services ['QUEUED'] = 0
			if services_user.filter(serviceStatus__exact = 'in_progress').exists():
				number_of_services ['IN PROGRESS'] = len (services_user.filter(serviceStatus__exact = 'in_progress'))
			else:
				number_of_services ['IN PROGRESS'] = 0
			if services_user.filter(serviceStatus__exact = 'delivered').exists():
				number_of_services ['DELIVERED'] = len (services_user.filter(serviceStatus__exact = 'delivered'))
			else:
				number_of_services ['DELIVERED'] = 0

			data_source = graphic_3D_pie('Status of Requested Services of:', user_name, '', '','fint',number_of_services)
			graphic_by_user_date_services = FusionCharts("pie3d", "ex1" , "600", "350", "chart-1", "json", data_source)
			stats_info ['graphic_by_user_date_services'] = graphic_by_user_date_services.render()

			# getting statistics of the created services

			service_dict ={}
			for service_available in services_user :
				service_list = service_available.serviceAvailableService.filter(level=3)
				for service in service_list:
					service_name = service.availServiceDescription
					if service_name in service_dict:
						service_dict [service_name] += 1
					else:
						service_dict [service_name] = 1
			#creating the graphic for requested services
			data_source = column_graphic_dict('Requested Services by:', user_name, '', '','fint',service_dict)
			graphic_requested_services = FusionCharts("column3d", "ex2" , "600", "350", "chart-2", "json", data_source)
			stats_info ['graphic_requested_services'] = graphic_requested_services.render()

			# getting statistics for requested per time
			service_time_dict ={}
			for service_per_time in services_user :
				date_service = service_per_time.serviceCreatedOnDate.strftime("%m_%Y")
				if date_service in service_time_dict:
					service_time_dict[date_service] +=1
				else:
					service_time_dict[date_service] =1
			# sorting the dictionary to get
			#creating the graphic for monthly requested services
			service_time_tupla =[]
			for key , value in sorted(service_time_dict.items()):
				#import pdb ; pdb.set_trace()
				service_time_tupla.append([key,service_time_dict[key]])
			data_source = column_graphic_tupla('Requested Services by:', user_name, '', '','fint',service_time_tupla)
			graphic_date_requested_services = FusionCharts("column3d", "ex3" , "600", "350", "chart-3", "json", data_source)
			stats_info ['graphic_date_requested_services'] = graphic_date_requested_services.render()

			#import pdb ; pdb.set_trace()
			return render (request, 'drylab/statsByDateUser.html', {'stats_info':stats_info})
	else:
		form = ByDateUserStats()
		return render(request, 'drylab/statsByDateUser.html', {'form':form})


@login_required
def stats_by_services_request (request):
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
	if request.method == 'POST':
		form = ByServicesRequest(data=request.POST)
		if form.is_valid():
			# validate the input data in the form
			start_date = form['start_date'].data
			end_date = form['end_date'].data
			if start_date != '':
				try:
					start_date_format = datetime.datetime.strptime(start_date, '%Y-%m-%d')
				except:
					return render (request,'drylab/error_page.html', {'content':['The format for the "Start Date Search" Field is incorrect ',
																				'ADVICE:', 'Use the format  (DD-MM-YYYY)']})
			if end_date != '':
				try:
					end_date_format = datetime.datetime.strptime(end_date, '%Y-%m-%d')
				except:
					return render (request,'drylab/error_page.html', {'content':['The format for the "End Date Search" Field is incorrect ',
																				'ADVICE:', 'Use the format  (DD-MM-YYYY)']})

			if Service.objects.filter(serviceCreatedOnDate__range=(start_date,end_date)).exists():
				services_found = Service.objects.filter(serviceCreatedOnDate__range=(start_date,end_date)). order_by('-serviceCreatedOnDate')
				services_stats_info = {}
				#preparing stats for services request by users
				user_services ={}
				for service in services_found:
					user = service.serviceUserId.username
					if user in user_services :
						user_services[user] +=1
					else:
						user_services[user] = 1
				#import pdb ; pdb.set_trace()
				period_of_time_selected = str(' For the period between ' + start_date + ' and ' + end_date)
				#creating the graphic for requested services
				data_source = column_graphic_dict('Requested Services by:', period_of_time_selected , 'User names', 'Number of Services','fint',user_services)
				graphic_requested_services = FusionCharts("column3d", "ex1" , "525", "350", "chart-1", "json", data_source)
				services_stats_info ['graphic_requested_services'] = graphic_requested_services.render()
				#preparing stats for status of the services
				status_services ={}
				for service in services_found:
					#user_id = service.serviceUserId.id
					#import pdb ; pdb.set_trace()
					status = service.serviceStatus
					if status in status_services :
						status_services[status] +=1
					else:
						status_services[status] = 1
				#creating the graphic for status services
				data_source = graphic_3D_pie('Status of Requested Services', period_of_time_selected ,'', '','fint',status_services)
				graphic_status_requested_services = FusionCharts("pie3d", "ex2" , "525", "350", "chart-2", "json", data_source)
				services_stats_info ['graphic_status_requested_services'] = graphic_status_requested_services.render()

				#preparing stats for request by Area
				user_area_dict ={}
				for service in services_found:
					user_id = service.serviceUserId.id
					user_area = Profile.objects.get(profileUserID = user_id).profileArea

					if user_area in user_area_dict:
						user_area_dict[user_area] +=1
					else:
						user_area_dict[user_area] = 1
				#creating the graphic for areas
				data_source = column_graphic_dict('Services requested per Area', period_of_time_selected, 'Area ', 'Number of Services','fint',user_area_dict)
				graphic_area_services = FusionCharts("column3d", "ex3" , "600", "350", "chart-3", "json", data_source)
				services_stats_info ['graphic_area_services'] = graphic_area_services.render()

				#preparing stats for services request by Center
				user_center_dict ={}
				for service in services_found:
					user_id = service.serviceUserId.id
					user_center = Profile.objects.get(profileUserID = user_id).profileCenter.centerAbbr

					if user_center in user_center_dict:
						user_center_dict[user_center] +=1
					else:
						user_center_dict[user_center] = 1
				#creating the graphic for areas
				data_source = column_graphic_dict('Services requested per Center', period_of_time_selected, 'Center ', 'Number of Services','fint',user_center_dict)
				graphic_center_services = FusionCharts("column3d", "ex4" , "600", "350", "chart-4", "json", data_source)
				services_stats_info ['graphic_center_services'] = graphic_center_services.render()
				#import pdb ; pdb.set_trace()

				################################################
				## Preparing the statistics per period of time
				################################################
				# calculating the period to be displayed the graphic (per month o per year)
				delta_dates = (end_date_format - start_date_format).days
				if delta_dates > 366 :
					period_year_month = '%Y'
				else:
					period_year_month = '%m_%Y'

				## Preparing the statistics for Center on period of time
				user_services_period ={}
				center_period = {}
				time_values_dict = {}
				for service in services_found:
					user_id = service.serviceUserId.id
					date_service = service.serviceCreatedOnDate.strftime(period_year_month)
					user_center = Profile.objects.get(profileUserID = user_id).profileCenter.centerAbbr
					if not date_service in time_values_dict:
						time_values_dict[date_service] = 1
					if user_center in user_services_period:
						if date_service in user_services_period[user_center] :
							user_services_period[user_center][date_service] +=1
						else:
							user_services_period[user_center][date_service] = 1
					else:
						user_services_period[user_center]= {}
						user_services_period[user_center][date_service] = 1
				time_values =[]
				for key , values in sorted(time_values_dict.items()):
					time_values.append(key)
				# fill with zero for the centers that have no sevice during some period
				for center , value in user_services_period.items():
					for d_period in time_values:
						if not d_period in user_services_period[center]:
							user_services_period[center][d_period] = 0

				data_source = column_graphic_per_time ('Services requested by center ',period_of_time_selected,  'date', 'number of services', time_values , user_services_period)
				graphic_center_services_per_time = FusionCharts("mscolumn3d", "ex5" , "525", "350", "chart-5", "json", data_source)
				services_stats_info ['graphic_center_services_per_time'] = graphic_center_services_per_time.render()
				#import pdb ; pdb.set_trace()

				## Preparing the statistics for Area on period of time
				user_area_services_period ={}
				area_period = {}
				time_values_dict = {}
				for service in services_found:
					user_id = service.serviceUserId.id
					date_service = service.serviceCreatedOnDate.strftime(period_year_month)
					user_area =  Profile.objects.get(profileUserID = user_id).profileArea
					if not date_service in time_values_dict:
						time_values_dict[date_service] = 1
					if user_area in user_area_services_period:
						if date_service in user_area_services_period[user_area] :
							user_area_services_period[user_area][date_service] +=1
						else:
							user_area_services_period[user_area][date_service] = 1
					else:
						user_area_services_period[user_area]= {}
						user_area_services_period[user_area][date_service] = 1
				time_values =[]
				for key , values in sorted(time_values_dict.items()):
					time_values.append(key)
				# fill with zero for the centers that have no sevice during some period
				for area , value in user_area_services_period.items():
					for d_period in time_values:
						if not d_period in user_area_services_period[area]:
							user_area_services_period[area][d_period] = 0
				#import pdb ; pdb.set_trace()
				data_source = column_graphic_per_time ('Services requested by Area ',period_of_time_selected,  'date', 'number of services', time_values , user_area_services_period)
				graphic_area_services_per_time = FusionCharts("mscolumn3d", "ex6" , "525", "350", "chart-6", "json", data_source)
				services_stats_info ['graphic_area_services_per_time'] = graphic_area_services_per_time.render()

				services_stats_info['period_time']= period_of_time_selected
				return render (request, 'drylab/statsByServicesRequest.html', {'services_stats_info':services_stats_info})

			else:
				return render (request,'drylab/error_page.html', {'content':['There are no services created by ', 'For the time of period of between:',
																start_date , 'and', end_date]})
	else:
		form = ByServicesRequest()
	return render(request, 'drylab/statsByServicesRequest.html', {'form':form})





@login_required
def stats_by_samples_processed (request):
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
	if request.method == 'POST':
		form = ByDateUserStats(data=request.POST)
		if form.is_valid():
			# validate the input data in the form
			start_date = form['start_date'].data
			end_date = form['end_date'].data


	else:
		form = BySampleProcessed()
		return render(request, 'drylab/statsBySamplesProcessed.html', {'form':form})



@login_required
def stats_time_delivery (request):
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
	if request.method == 'POST':
		form = ByDateUserStats(data=request.POST)
		if form.is_valid():
			# validate the input data in the form
			start_date = form['start_date'].data
			end_date = form['end_date'].data


	else:
		form = TimeDelivery()
		return render(request, 'drylab/statsByDateUser.html', {'form':form})







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
		#import pdb ; pdb.set_trace()
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

