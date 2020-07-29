# -*- coding: utf-8 -*-
from django.shortcuts import get_object_or_404, render, redirect
from .models import *
from .forms import *
from .graphics import *
import os, re

from django.contrib.auth.decorators import login_required
from django.contrib.auth.models import Group
from django.conf import settings
# from django.template.loader import render_to_string
from django_utils.fusioncharts.fusioncharts import FusionCharts
# from django.core.mail import send_mail
from django.http import HttpResponse
from django.core.files.storage import FileSystemStorage

import datetime
import statistics
from iSkyLIMS_drylab import drylab_config
from iSkyLIMS_drylab.utils.testing_drylab_configuration import *
from iSkyLIMS_drylab.utils.drylab_common_functions import *
from iSkyLIMS_drylab.utils.handling_pipelines import *
from iSkyLIMS_drylab.utils.handling_forms import *
from iSkyLIMS_drylab.utils.configuration_functions import *

from smb.SMBConnection import SMBConnection
# from iSkyLIMS_drylab.models import RunProcess, Projects

####### Import libraries for static files
#from django.shortcuts import render_to_response
#from django.shortcuts import RequestContext
###### api to wetlab ############
try :
	from iSkyLIMS_wetlab.utils.api.wetlab_api import *
	wetlab_api_available = True
except:
	wetlab_api_available = False


#PNL
@login_required
def index(request):
    if Service.objects.all().exclude(serviceStatus = 'delivered').exclude(serviceStatus = 'rejected').exclude(serviceStatus = 'approved').exists():
        ongoing_services = Service.objects.all().exclude(serviceStatus = 'delivered').exclude(serviceStatus = 'rejected').exclude(serviceStatus = 'approved').order_by('serviceCreatedOnDate')
        service_list = []
        for service in ongoing_services:
            service_info = []
            service_info.append(service.serviceRequestNumber)
            service_info.append(service.serviceStatus)
            if service.serviceStatus == 'recorded':
                service_delivery_date = 'Not defined yet'
            else:
                if Resolution.objects.filter(resolutionServiceID = service).exists():

                    if Resolution.objects.filter(resolutionServiceID = service).last().resolutionEstimatedDate is not None:
                        service_delivery_date = Resolution.objects.filter(resolutionServiceID = service).last().resolutionEstimatedDate.strftime("%d %B, %Y")
                    else:
                        service_delivery_date = 'Not defined yet'
                else:
                    service_delivery_date = 'Not defined yet'
            service_info.append(service_delivery_date)
            service_list.append(service_info)

        return render(request, 'iSkyLIMS_drylab/index.html',{'service_list': service_list})
    else:
        return render(request, 'iSkyLIMS_drylab/index.html')

@login_required
def configuration_email(request):
    if request.user.username != 'admin':
        return redirect('')
    email_conf_data = get_email_data_from_file(__package__)
    if request.method == 'POST' and (request.POST['action']=='emailconfiguration'):
        email_user_field ={}
        for field in drylab_config.EMAIL_CONFIGURATION_FIELDS:
            email_user_field[field] = request.POST[field]

        if not create_email_conf_file (email_user_field, __package__) :
            error_message = drylab_config.ERROR_UNABLE_TO_SAVE_EMAIL_CONFIGURATION_SETTINGS
            return render(request, 'iSkyLIMS_drylab/configurationEmail.html',{'email_conf_data':email_user_field, 'error_message': error_message} )
        import importlib
        importlib.reload(drylab_config)
        return render(request, 'iSkyLIMS_drylab/configurationEmail.html',{'succesful_settings':True})
    return render(request, 'iSkyLIMS_drylab/configurationEmail.html',{'email_conf_data': email_conf_data})



@login_required
def service_request(request, serviceRequestType):
    request_type = {}
    if serviceRequestType == 'internal_sequencing':
        if request.method == "POST":
            if 'serviceProjects' in request.POST:
                from django.http import QueryDict
                #project_list = data.pop('serviceProjects', None)
                #data_dict = request.POST.dict()
                data = QueryDict('', mutable=True)
                data.update(request.POST)
                project_list = request.POST.getlist('serviceProjects')
                #data_dict['serviceAvailableService'] = request.POST.getlist('serviceAvailableService')
                dummy_value = data.pop('serviceProjects', None)
            else:
                data=request.POST

            form = ServiceRequestFormInternalSequencing(data = data, files = request.FILES)
            if form.is_valid():
                new_service = save_service_request_form(form, request.user, drylab_config.INTERNAL_SEQUENCING_UNIT)

                service_request_number = new_service.get_service_request_number()
                if 'serviceProjects' in request.POST:
                    stored_projects = store_projects_from_form(project_list, new_service)
                else:
                    stored_projects = ''

                ## Send mail to user and drylab admin group
                if drylab_config.EMAIL_USER_CONFIGURED :
                    email_data = {}
                    email_data['user_email'] = request.user.email
                    email_data['user_name'] = request.user.username
                    email_data['service_number'] = service_request_number
                    send_service_creation_confirmation_email(email_data)

                # PDF preparation file for confirmation of service request
                pdf_file = create_service_pdf_file(service_request_number, request.build_absolute_uri())

                # check if service allows to get data from external applications
                if len(services_allow_external_data(new_service, stored_projects)) > 0:
                    if drylab_config.EMAIL_USER_CONFIGURED :
                        send_required_preparation_pipeline_email(service_request_number)
                confirmation_result = {}
                confirmation_result['download_file'] = pdf_file
                confirmation_result['text'] = list(map(lambda st: str.replace(st, 'SERVICE_NUMBER', service_request_number), drylab_config.CONFIRMATION_TEXT_MESSAGE))
                return render(request,'iSkyLIMS_drylab/RequestForm.html',{'confirmation_result':confirmation_result})
            else:
                error_message = drylab_config.ERROR_UNABLE_TO_RECORD_YOUR_SERVICE
                request_type['type'] = 'Internal Sequencing'
                form = prepare_form_data_internal_sequencing(request.user)
                return render(request, 'iSkyLIMS_drylab/RequestForm.html',{'form': form, 'error_message': error_message , 'request_type': request_type})
        else:
            form = prepare_form_data_internal_sequencing(request.user)
            request_type['type'] = 'Internal Sequencing'
            return render(request, 'iSkyLIMS_drylab/RequestForm.html' , { 'form' : form , 'request_type': request_type})


    if serviceRequestType == 'external_sequencing':
        if request.method == "POST":
            form = ServiceRequestFormExternalSequencing(data=request.POST,files=request.FILES)
            if form.is_valid():
                new_service = save_service_request_form(form, request.user, request.POST['serviceSeqCenter'])
                service_request_number = new_service.get_service_request_number()

                ## Send email
                ## Send mail to user and drylab admin group
                if drylab_config.EMAIL_USER_CONFIGURED :
                    email_data = {}
                    email_data['user_email'] = request.user.email
                    email_data['user_name'] = request.user.username
                    email_data['service_number'] = service_request_number
                    send_service_creation_confirmation_email(email_data)

                # PDF preparation file for confirmation of service request
                #absolute_url = request.build_absolute_uri()
                pdf_file = create_service_pdf_file(service_request_number, request.build_absolute_uri())
                confirmation_result = {}
                confirmation_result['download_file'] = pdf_file
                confirmation_result['text'] = list(map(lambda st: str.replace(st, 'SERVICE_NUMBER', service_request_number), drylab_config.CONFIRMATION_TEXT_MESSAGE))
                return render(request,'iSkyLIMS_drylab/RequestForm.html',{'confirmation_result':confirmation_result})

        else:
            form = ServiceRequestFormExternalSequencing()
            form.fields['serviceAvailableService'].queryset = AvailableService.objects.filter(availServiceDescription__exact="Genomic data analysis").get_descendants(include_self=True)
            request_type['type'] = 'External Sequencing'
            return render(request, 'iSkyLIMS_drylab/RequestForm.html' , { 'form' : form ,  'request_external': 'request_external','request_type': request_type })


@login_required
def counseling_request(request):
    if request.method == "POST":
        form = ServiceRequestForm_extended(data=request.POST,files=request.FILES)

        if form.is_valid():
            new_service = form.save(commit=False)
            new_service.serviceStatus = "recorded"
            new_service.serviceUserId = User.objects.get(id=request.user.id)
            new_service.serviceRequestInt = increment_service_number(request.user)
            new_service.serviceRequestNumber = create_service_id(new_service.serviceRequestInt,request.user)
            new_service.save()
            form.save_m2m()
            ## Send email
            subject = 'Service ' + new_service.serviceRequestNumber + " has been recorded"
            body_message = 'Dear ' + request.user.username + "\n Your service " + new_service.serviceRequestNumber + " has been recorded. You will recieved the resolution of the request as soon as possible.\n Kind regards \n BU-ISCIII \n bioinformatica@isciii.es"
            from_user = 'bioinformatica@isciii.es'
            to_user = [request.user.email,'bioinformatica@isciii.es']
            send_mail (subject, body_message, from_user, to_user)
            # PDF preparation file for confirmation of service request
            information_to_include = get_data_for_service_confirmation(str(new_service.serviceRequestNumber))
            pdf_file_name = str(new_service.serviceRequestNumber) + '.pdf'
            pdf_file = create_pdf(request, information_to_include, drylab_config.REQUESTED_CONFIRMATION_SERVICE, pdf_file_name)

            pdf_url = pdf_file.replace(settings.BASE_DIR,'')
            download_file = '<a href="'+ pdf_url + '">Download the service request confirmation file</a>'
            return render(request,'django_utils/info_page.html',{'content':['Your service request has been successfully recorded.',
                                'The sequence number assigned for your request is: ', new_service.serviceRequestNumber,
                                'Keep this number safe for refering your request', download_file ,
                                'You will be contacted shortly.']})

        else:

            return render(request,'iSkyLIMS_drylab/error_page.html',{'content':['Your service request cannot be recorded.',
                                                'Check that all information is provided correctly.']})
    else:
        form = ServiceRequestForm_extended()

    form.fields['serviceAvailableService'].queryset = AvailableService.objects.filter(availServiceDescription__exact="Bioinformatics consulting and training").get_descendants(include_self=True)
    return render(request, 'iSkyLIMS_drylab/RequestForm.html' , { 'form' : form ,  'consulting_request': 'consulting_request'})


@login_required
def infrastructure_request(request):
    if request.method == "POST":
        form = ServiceRequestForm_extended(data=request.POST or None,files=request.FILES)

        if form.is_valid():
            new_service = form.save(commit=False)
            new_service.serviceStatus = "recorded"
            new_service.serviceUserId = User.objects.get(id=request.user.id)
            new_service.serviceRequestInt = increment_service_number(request.user)
            new_service.serviceRequestNumber = create_service_id(new_service.serviceRequestInt,request.user)

            new_service.save()
            form.save_m2m()
            ## Send email
            subject = 'Service ' + new_service.serviceRequestNumber + " has been recorded"
            body_message = 'Dear ' + request.user.username + "\n Your service " + new_service.serviceRequestNumber + " has been recorded. You will received the resolution of the request as soon as possible.\n Kind regards \n BU-ISCIII \n bioinformatica@isciii.es"
            from_user = 'bioinformatica@isciii.es'
            to_user = [request.user.email,'bioinformatica@isciii.es']
            send_mail (subject, body_message, from_user, to_user)

            information_to_include = get_data_for_service_confirmation(str(new_service.serviceRequestNumber))
            pdf_file_name = str(new_service.serviceRequestNumber) + '.pdf'
            pdf_file = create_pdf(request, information_to_include, drylab_config.REQUESTED_CONFIRMATION_SERVICE, pdf_file_name)

            pdf_url = pdf_file.replace(settings.BASE_DIR,'')
            download_file = '<a href="'+ pdf_url + '">Download the service request confirmation file</a>'
            return render(request,'django_utils/info_page.html',{'content':['Your service request has been successfully recorded.',
                                'The sequence number assigned for your request is: ', new_service.serviceRequestNumber,
                                'Keep this number safe for refering your request',download_file,'You will be contacted shortly.']})
        else:

            return render(request,'iSkyLIMS_drylab/error_page.html',{'content':['Your service request cannot be recorded.',
                                                'Check that all information is provided correctly.']})
    else:
        form = ServiceRequestForm_extended()

    form.fields['serviceAvailableService'].queryset = AvailableService.objects.filter(availServiceDescription__exact="User support").get_descendants(include_self=True)

    #form.helper[1].update_atrributes(hidden="true")
    return render(request, 'iSkyLIMS_drylab/RequestForm.html' , { 'form' : form , 'infrastructure_request': 'infrastructure_request'})


@login_required
def display_service (request, service_id):
    if request.user.is_authenticated:
        if not is_service_manager(request):
            return render (request,'iSkyLIMS_drylab/error_page.html', {'content':drylab_config.ERROR_USER_NOT_ALLOWED })
    else:
        #redirect to login webpage
        return redirect ('/accounts/login')
    if Service.objects.filter(pk=service_id).exists():
        # displays the service information with the latest changes done using the forms
        display_service_details = get_service_information(service_id)
        return render (request,'iSkyLIMS_drylab/display_service.html',{'display_service': display_service_details})
    else:
        return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['The service that you are trying to get does not exist ','Contact with your administrator .']})


@login_required
def search_service (request):
    if request.user.is_authenticated:
        if not is_service_manager(request):
            return render (request,'iSkyLIMS_drylab/error_page.html', {'content':drylab_config.ERROR_USER_NOT_ALLOWED })
    else:
        #redirect to login webpage
        return redirect ('/accounts/login')
    if request.method == 'POST' and request.POST['action'] == 'searchservice':

        service_number_request = request.POST['servicenumber']
        service_state = request.POST['servicestate']
        start_date=request.POST['startdate']
        end_date=request.POST['enddate']
        center = request.POST['center']
        user_name = request.POST['username']
        if service_number_request == '' and service_state == '' and start_date == '' and end_date == '' and center == '' and user_name =='':
            services_search_list = {}
            center_list_abbr = []
            center_availables = Center.objects.all().order_by ('centerAbbr')
            for center in center_availables:
                center_list_abbr.append (center.centerAbbr)
            services_search_list ['centers'] = center_list_abbr
            services_search_list ['status'] = STATUS_CHOICES

            return render( request,'iSkyLIMS_drylab/searchService.html',{'services_search_list': services_search_list })

        ### check the right format of start and end date
        if start_date != '':
            try:
                datetime.datetime.strptime(start_date, '%Y-%m-%d')
            except:
                return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['The format for the "Start Date Search" Field is incorrect ',
                                                                            'ADVICE:', 'Use the format  (DD-MM-YYYY)']})
        if end_date != '':
            try:
                datetime.datetime.strptime(end_date, '%Y-%m-%d')
            except:
                return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['The format for the "End Date Search" Field is incorrect ',
                                                                            'ADVICE:', 'Use the format  (DD-MM-YYYY)']})
        if service_number_request == '' and service_state == '':
            services_found = Service.objects.all()

        if service_number_request != '':
            # check if the requested service in the form matches exactly with the existing service in DB
            if Service.objects.filter(serviceRequestNumber__exact = service_number_request).exists():

                services_found = Service.objects.get(serviceRequestNumber__exact = service_number_request)
                redirect_page = '/drylab/display_service=' + str(services_found.id)
                return redirect (redirect_page)
            if Service.objects.filter(serviceRequestNumber__icontains = service_number_request).exists():
                services_found = Service.objects.filter(serviceRequestNumber__icontains = service_number_request)
            else:
                return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['No matches have been found for the service number ', service_number_request ]})

        if service_state != '':
            if service_number_request =='':
                services_found = Service.objects.all()
            if services_found.filter(serviceStatus__exact = service_state).exists():
                services_found = services_found.filter(serviceStatus__exact = service_state)
            else:
                return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['No matches have been found for the service number in state', service_state ]})
        if start_date !='' and end_date != '':
            if services_found.filter(serviceCreatedOnDate__range=(start_date, end_date)).exists():
                services_found = services_found.filter(serviceCreatedOnDate__range=(start_date, end_date))
            else:
                return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['There are no services containing ', service_number_request,
                                                        ' created between ', start_date, 'and the ', end_date]})
        if start_date !='' and end_date == '':
            if services_found.filter(serviceCreatedOnDate__gte = start_date).exists():
                services_found = services_found.filter(serviceCreatedOnDate__gte = start_date)
            else:
                return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['There are no services containing ', service_number_request,
                                                        ' created before ', start_date]})
        if start_date =='' and end_date != '':
            if services_found.filter(serviceCreatedOnDate__lte = end_date).exists():
                services_found = services_found.filter(serviceCreatedOnDate__lte = end_date)
            else:
                return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['There are no services containing ', service_number_request,
                                                        ' finish before ', end_date]})
        if center != '':
            if services_found.filter(serviceRequestNumber__icontains = center).exists():
                services_found = services_found.filter(serviceRequestNumber__icontains  = center)
            else:
                return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['There are no services related to the requested center', center]})

        if  user_name != '':
            if User.objects.filter (username__icontains = user_name).exists():
                user_id = User.objects.get (username__icontains = user_name).id
            else:
                return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['The user name  ', user_name, 'is not defined in iSkyLIMS']})
            if services_found.filter(serviceUserId = user_id).exists():
                services_found = services_found.filter(serviceUserId  = user_id)
            else:
                return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['There are no services requested by the user', center]})

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
            return render (request,'iSkyLIMS_drylab/searchService.html', {'display_multiple_services': display_multiple_services})
    services_search_list = {}

    center_list_abbr = []
    center_availables = Center.objects.all().order_by ('centerAbbr')
    for center in center_availables:
        center_list_abbr.append (center.centerAbbr)
    services_search_list ['centers'] = center_list_abbr
    services_search_list ['status'] = STATUS_CHOICES

    return render( request,'iSkyLIMS_drylab/searchService.html',{'services_search_list': services_search_list })


@login_required
def pending_services (request):
    if request.user.is_authenticated:
        if not is_service_manager(request):
            return render (request,'iSkyLIMS_drylab/error_page.html', {'content':drylab_config.ERROR_USER_NOT_ALLOWED })
    else:
        #redirect to login webpage
        return redirect ('/accounts/login')

    pending_services_details = {}
    recorded, queued, in_progress = {}, {}, {}
    if Service.objects.filter(serviceStatus__exact = 'recorded').exists():
        services_in_request = Service.objects.filter(serviceStatus__exact = 'recorded').order_by('-serviceCreatedOnDate')
        for services in services_in_request:
            recorded[services.id]= [services.get_service_information().split(';')]
        pending_services_details['recorded'] = recorded
    if Service.objects.filter(serviceStatus__exact = 'queued').exists():
        services_in_queued = Service.objects.filter(serviceStatus__exact = 'queued').order_by('-serviceCreatedOnDate')
        for services in services_in_queued:
            queued[services.id]= [services.get_service_information_with_service_name().split(';')]
        pending_services_details['queued'] = queued
    if Service.objects.filter(serviceStatus__exact = 'in_progress').exists():
        services_in_progress = Service.objects.filter(serviceStatus__exact = 'in_progress').order_by('-serviceCreatedOnDate')
        for services in services_in_progress:
            in_progress[services.id]= [services.get_service_information_with_service_name().split(';')]
        pending_services_details['in_progress'] = in_progress

    number_of_services = {}
    number_of_services ['RECORDED'] = len (recorded)
    number_of_services ['QUEUED'] = len (queued)
    number_of_services ['IN PROGRESS'] = len (in_progress)
    data_source = graphic_3D_pie('Number of Pending Services', '', '', '','fint',number_of_services)
    graphic_pending_services = FusionCharts("pie3d", "ex1" , "425", "350", "chart-1", "json", data_source)
    pending_services_details ['graphic_pending_services'] = graphic_pending_services.render()



    return render (request, 'iSkyLIMS_drylab/pendingServices.html', {'pending_services': pending_services_details})

@login_required
def add_resolution (request, service_id):
    if request.user.is_authenticated:
        if not is_service_manager(request):
            return render (request,'iSkyLIMS_drylab/error_page.html', {'content':drylab_config.ERROR_USER_NOT_ALLOWED })
    else:
        #redirect to login webpage
        return redirect ('/accounts/login')

    if request.method == "POST" and request.POST['action'] == 'addResolutionService' :

        #form = AddResolutionService(data=request.POST)
        if check_service_id_exists(request.POST['service_id']):

            resolution_data_form = get_add_resolution_data_form(request.POST)
            resolution_data_form['resolutionFullNumber'] = get_assign_resolution_full_number(resolution_data_form['service_id'], resolution_data_form['acronymName'])
            resolution_data_form['resolutionNumber'] = create_resolution_number(request.POST['service_id'])
            service_obj = get_service_obj_from_id(request.POST['service_id'])
            service_request_number = service_obj.get_service_request_number()
            if resolution_data_form['serviceAccepted'] == 'Accepted':
                service_obj.update_service_status("queued")
                service_obj.update_approved_date(datetime.date.today())
            else:
                service_obj.update_service_status("rejected")
                service_obj.update_rejected_date(datetime.date.today())

            new_resolution = Resolution.objects.create_resolution(resolution_data_form)

            if 'additional_parameters' in resolution_data_form:
                store_resolution_additional_parameter(resolution_data_form['additional_parameters'], new_resolution)

            # create a new resolution to be added to the service folder including the path where file is stored

            pdf_file = create_resolution_pdf_file(service_obj,new_resolution, request.build_absolute_uri())
            new_resolution.update_resolution_file(pdf_file)
            #pdf_name = resolution_data_form['resolutionNumber'] + ".pdf"
            #resolution_file = create_pdf(request,information, drylab_config.RESOLUTION_TEMPLATE, pdf_name)

            ## Send email
            if drylab_config.EMAIL_USER_CONFIGURED :
                email_data = {}
                email_data['user_email'] = request.user.email
                email_data['user_name'] = request.user.username
                email_data['service_number'] = service_request_number
                email_data['status'] = resolution_data_form['serviceAccepted']
                email_data['date'] = resolution_data_form['resolutionEstimatedDate']
                send_resolution_creation_email(email_data)
            created_resolution = {}
            created_resolution['resolution_number'] = resolution_data_form['resolutionNumber']
            return render(request,'iSkyLIMS_drylab/addResolution.html',{'created_resolution': created_resolution})
        else:
            return render (request, 'iSkyLIMS_drylab/error_page.html', {'content':drylab_config.ERROR_SERVICE_ID_NOT_FOUND})
    else:
        if check_service_id_exists(service_id):
            form_data = prepare_form_data_add_resolution(service_id)
            return render(request, 'iSkyLIMS_drylab/addResolution.html' , { 'form_data' : form_data})
        else:
            return render (request, 'iSkyLIMS_drylab/error_page.html', {'content':drylab_config.ERROR_SERVICE_ID_NOT_FOUND})


'''
def test (request):
    resolution_number = 'SRVIIER001.1'
    service_requested = 'SRVIIER001'
    from weasyprint import HTML, CSS
    from django.template.loader import get_template
    from django.core.files.storage import FileSystemStorage
    from django.http import HttpResponse
    from weasyprint.fonts import FontConfiguration



    information, user, resolution_data = {}, {}, {}
    # get service object
    service = Service.objects.get(serviceRequestNumber = service_requested)
    service_number ,run_specs, center, platform = service.get_service_information().split(';')
    # get resolution object
    resolution = Resolution.objects.get(resolutionNumber = resolution_number)
    resolution_info = resolution.get_resolution_information()
    # get profile object
    user_id = service.serviceUserId.id

    information['resolution_number'] = resolution_number
    information['requested_date'] = service.get_service_creation_time()
    information['resolution_date'] = resolution_info[4]
    information['nodes']= service.serviceAvailableService.all()
    user['name'] = service.serviceUserId.first_name
    user['surname'] = service.serviceUserId.last_name


    user_id = service.serviceUserId.id
    user['area'] = Profile.objects.get(profileUserID = user_id).profileArea
    user['center'] = Profile.objects.get(profileUserID = user_id).profileCenter
    user['position'] = Profile.objects.get(profileUserID = user_id).profilePosition
    user['phone'] = Profile.objects.get(profileUserID = user_id).profileExtension
    user['email'] = service.serviceUserId.email
    information['user'] = user
    resolution_data['folder'] = resolution_info[1]
    resolution_data['estimated_date'] = resolution_info[3]
    resolution_data['notes'] = resolution_info[6]
    resolution_data['decission'] = service.serviceStatus
    information['service_data'] = service.serviceNotes

    resolution_data['folder'] = resolution_info[1]
    information['resolution_data'] = resolution_data
    html_string = render_to_string('resolution_template.html', {'information': information})

    html = HTML(string=html_string, base_url=request.build_absolute_uri()).write_pdf('documents/drylab/res_pdf.pdf',stylesheets=[CSS(settings.STATIC_ROOT +
                                drylab_config.CSS_FOR_PDF)])

    fs = FileSystemStorage('documents/drylab')
    with fs.open('res_pdf.pdf') as pdf:
        response = HttpResponse(pdf, content_type='application/pdf')
        # save pdf file as attachment
        #response['Content-Disposition'] = 'attachment; filename="mypdf.pdf"'


        response['Content-Disposition'] = 'inline;filename=res_pdf.pdf'

    return response

'''

def add_new_resolution_file (conn, full_service_path,resolution_file,year):

    temp_file=resolution_file.split('/')
    resolution_name_file = temp_file[-1]
    resolution_remote_file = os.path.join(drylab_config.SAMBA_SERVICE_FOLDER,str(year),full_service_path,drylab_config.FOLDERS_FOR_SERVICES[1],resolution_name_file)

    try:
        with open(resolution_file ,'rb') as  res_samba_fp:
            conn.storeFile(drylab_config.SAMBA_SHARED_FOLDER_NAME, resolution_remote_file, res_samba_fp)
    except:
        return ( 'Unable to copy the resolution file ',resolution_remote_file,resolution_name_file)

    return True

'''
def create_service_structure (conn, service_request_file, service_file_uploaded, full_service_path, resolution_file):
    ## service_request_file and resolution_file contains the full path where these files
    ## are stored on iSkyLIMS. It means that OUTPUT_DIR_TEMPLATE value is added to thes variable
    ## to store the files on the remote system we need to have the full pathe where these files
    ## are located, but also the file name without including the path, in order to add only
    ## the file name to the remote path. To get only the file name we split the variable (containing
    ## path and file name ) to fetch only the file name


    # get the information for creating the subfolders
    time_now = datetime.datetime.now()
    year = str(time_now.year)
    # check if year directory already exists on remote server
    file_list = conn.listPath( drylab_config.SAMBA_SHARED_FOLDER_NAME, drylab_config.SAMBA_SERVICE_FOLDER)
    year_folder_exists = False
    for sh_file in file_list:
        if sh_file.filename == year:
            year_folder_exists = True
    year_folder = os.path.join(drylab_config.SAMBA_SERVICE_FOLDER, year)
    if not year_folder_exists :
        conn.createDirectory (drylab_config.SAMBA_SHARED_FOLDER_NAME, year_folder)

    service_path = os.path.join(year_folder, full_service_path)
    #create the directory for the new service
    conn.createDirectory (drylab_config.SAMBA_SHARED_FOLDER_NAME, service_path)
    for sub_folder in drylab_config.FOLDERS_FOR_SERVICES:
        sub_folder_path = os.path.join(service_path,sub_folder)
        conn.createDirectory(drylab_config.SAMBA_SHARED_FOLDER_NAME, sub_folder_path)

    #copy service confirmation file into request folder
    temp_file=resolution_file.split('/')
    resolution_name_file = temp_file[-1]
    resolution_remote_file = os.path.join(service_path,drylab_config.FOLDERS_FOR_SERVICES[1],resolution_name_file)

    try:
        with open(resolution_file ,'rb') as  res_samba_fp:
            conn.storeFile(drylab_config.SAMBA_SHARED_FOLDER_NAME, resolution_remote_file, res_samba_fp)
    except:
        return 'ERROR:: Unable to copy resolution file'
    temp_file=service_request_file.split('/')

    request_name_file = temp_file[-1]
    request_remote_file = os.path.join(service_path,drylab_config.FOLDERS_FOR_SERVICES[0],request_name_file)

    try:
        with open(service_request_file ,'rb') as  req_samba_fp:
            conn.storeFile(drylab_config.SAMBA_SHARED_FOLDER_NAME, request_remote_file, req_samba_fp)

    except:
            return 'ERROR:: Unable to copy service requested file '
    if service_file_uploaded != '':
        temp_file_name = service_file_uploaded.split('/')
        uploaded_name_file = temp_file_name [-1]
        uploaded_remote_file = os.path.join(service_path, drylab_config.FOLDERS_FOR_SERVICES[0],uploaded_name_file)

        try:
            with open(service_file_uploaded ,'rb') as  upload_samba_fp:
                conn.storeFile(drylab_config.SAMBA_SHARED_FOLDER_NAME, uploaded_remote_file, upload_samba_fp)
        except:
            return 'ERROR:: Unable to copy file uploaded by the investigator'

    # deleting the service_request_file and resolution_file from iSkyLIMS
    try:
        os.remove(service_request_file)
        os.remove(resolution_file)
    except:
        return 'ERROR:: Unable to delete the service_requested_file/ resolution_file'

    return True
'''

@login_required
def add_in_progress (request, resolution_id):
    if request.user.is_authenticated:
        try:
            groups = Group.objects.get(name = drylab_config.SERVICE_MANAGER)
            if groups not in request.user.groups.all():
                return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
        except:
            return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
    else:
        #redirect to login webpage
        return redirect ('/accounts/login')

    if Resolution.objects.filter(pk = resolution_id).exists():
        resolution = Resolution.objects.get(pk = resolution_id)
        service_to_update = resolution.resolutionServiceID
        # update the service status and in_porgress date

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
        return render (request,'django_utils/info_page.html',{'content':['Your resolution  request ', resolution.resolutionNumber,
                                'has been successfully upated to In Progress state']})
    else:

        return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['The resolution that you are trying to upadate does not exists ','Contact with your administrator .']})
    return

@login_required
def add_delivery (request , resolution_id):
    if request.user.is_authenticated:
        try:
            groups = Group.objects.get(name = drylab_config.SERVICE_MANAGER)
            if groups not in request.user.groups.all():
                return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
        except:
            return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
    else:
        #redirect to login webpage
        return redirect ('/accounts/login')
    if request.method == 'POST' :
        form = AddDeliveryService(data=request.POST)
        if form.is_valid():

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
            return render(request,'django_utils/info_page.html',{'content':['The service is now on Delivery status ']})
    else:
        if Resolution.objects.filter(pk = resolution_id).exists():

            form = AddDeliveryService()
            delivery_info = {}
            return render (request, 'iSkyLIMS_drylab/addDelivery.html', {'form':form, 'delivery_info': delivery_info})
        else:

            return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['The resolution that you are trying to upadate does not exists ','Contact with your administrator .']})
    return


@login_required
def stats_by_date_user (request):
    if request.user.is_authenticated:
        try:
            groups = Group.objects.get(name = drylab_config.SERVICE_MANAGER)
            if groups not in request.user.groups.all():
                return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
        except:
            return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
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
                    return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['Too many matches have been found for the user name field', user_name,
                                                                                'ADVICE:', 'Please write down one of the following user name and repeate again the search',
                                                                                name_string,]})
                else:
                    user_name_id = matched_names[0].id
                    user_name = matched_names[0].username
            else:
                return render (request,'iSkyLIMS_drylab/error_page.html', {'content':[user_name,'is not defined on database']})

            if start_date != '':
                try:
                    datetime.datetime.strptime(start_date, '%Y-%m-%d')
                except:
                    return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['The format for the "Start Date Search" Field is incorrect ',
                                                                                'ADVICE:', 'Use the format  (DD-MM-YYYY)']})
            if end_date != '':
                try:
                    datetime.datetime.strptime(end_date, '%Y-%m-%d')
                except:
                    return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['The format for the "End Date Search" Field is incorrect ',
                                                                                'ADVICE:', 'Use the format  (DD-MM-YYYY)']})

            services_user = Service.objects.filter(serviceUserId__exact = user_name_id).order_by('-serviceRequestNumber')
            if start_date != '' and end_date !='':
                if services_user.filter(serviceCreatedOnDate__range=(start_date,end_date)).exists():
                    services_user = services_user.filter(serviceCreatedOnDate__range=(start_date,end_date))
                else:
                    return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['There are no services created by ', user_name , 'For the time of period of between:',
                                                                start_date , 'and', end_date]})
            if start_date !='' and end_date == '':
                if services_user.filter(serviceCreatedOnDate__gte = start_date).exists():
                    services_user = services_user.filter(serviceCreatedOnDate__lte = start_date)
                else:
                    return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['There are no services created by ', user_name , 'Starting from ', start_date ]})
            if start_date =='' and end_date != '':
                if services_user.filter(serviceCreatedOnDate__lte = end_date).exists():
                    services_user = services_user.filter(serviceCreatedOnDate__lte = end_date)
                else:
                    return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['There are no services created by ', user_name , 'Finish before ', end_date ]})


            stats_info = {}
            service_by_user=[]
            for service_item in services_user:
                service_by_user.append(service_item.get_stats_information())


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

                service_time_tupla.append([key,service_time_dict[key]])
            data_source = column_graphic_tupla('Requested Services by:', user_name, '', '','fint',service_time_tupla)
            graphic_date_requested_services = FusionCharts("column3d", "ex3" , "600", "350", "chart-3", "json", data_source)
            stats_info ['graphic_date_requested_services'] = graphic_date_requested_services.render()


            return render (request, 'iSkyLIMS_drylab/statsByDateUser.html', {'stats_info':stats_info})
    else:
        form = ByDateUserStats()
        return render(request, 'iSkyLIMS_drylab/statsByDateUser.html', {'form':form})


@login_required
def stats_by_services_request (request):
    if request.user.is_authenticated:
        try:
            groups = Group.objects.get(name = drylab_config.SERVICE_MANAGER)
            if groups not in request.user.groups.all():
                return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
        except:
            return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
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
                    return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['The format for the "Start Date Search" Field is incorrect ',
                                                                                'ADVICE:', 'Use the format  (DD-MM-YYYY)']})
            if end_date != '':
                try:
                    end_date_format = datetime.datetime.strptime(end_date, '%Y-%m-%d')
                except:
                    return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['The format for the "End Date Search" Field is incorrect ',
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

                period_of_time_selected = str(' For the period between ' + start_date + ' and ' + end_date)
                #creating the graphic for requested services
                data_source = column_graphic_dict('Requested Services by:', period_of_time_selected , 'User names', 'Number of Services','fint',user_services)
                graphic_requested_services = FusionCharts("column3d", "ex1" , "525", "350", "chart-1", "json", data_source)
                services_stats_info ['graphic_requested_services_per_user'] = graphic_requested_services.render()
                #preparing stats for status of the services
                status_services ={}
                for service in services_found:
                    #user_id = service.serviceUserId.id

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
                    user_id = service.get_service_id()
                    if Profile.objects.filter(profileUserID = user_id).exists():
                        user_area = Profile.objects.get(profileUserID = user_id).profileArea
                    else:
                        user_area = 'No_user_area'

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
                    user_id = service.get_service_id()
                    if Profile.objects.filter(profileUserID = user_id).exists():
                        user_center = Profile.objects.get(profileUserID = user_id).profileCenter.centerAbbr
                    else:
                        user_center = 'Not defined'
                    if user_center in user_center_dict:
                        user_center_dict[user_center] +=1
                    else:
                        user_center_dict[user_center] = 1
                #creating the graphic for areas
                data_source = column_graphic_dict('Services requested per Center', period_of_time_selected, 'Center ', 'Number of Services','fint',user_center_dict)
                graphic_center_services = FusionCharts("column3d", "ex4" , "600", "350", "chart-4", "json", data_source)
                services_stats_info ['graphic_center_services'] = graphic_center_services.render()


                ################################################
                ## Preparing the statistics per period of time
                ################################################
                # calculating the period to be displayed the graphic (per month o per year)
                delta_dates = (end_date_format - start_date_format).days
                if delta_dates > 366 :
                    period_year_month = '%Y'
                else:
                    period_year_month = '%Y_%m'

                ## Preparing the statistics for Center on period of time
                user_services_period ={}
                center_period = {}
                time_values_dict = {}
                for service in services_found:
                    user_id = service.get_service_id()
                    date_service = service.serviceCreatedOnDate.strftime(period_year_month)
                    if Profile.objects.filter(profileUserID = user_id).exists():
                        user_center = Profile.objects.get(profileUserID = user_id).profileCenter.centerAbbr
                    else:
                        user_center = 'Not defined'
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


                ## Preparing the statistics for Area on period of time
                user_area_services_period ={}
                area_period = {}
                time_values_dict = {}
                for service in services_found:
                    user_id = service.serviceUserId.id
                    date_service = service.serviceCreatedOnDate.strftime(period_year_month)
                    if  Profile.objects.filter(profileUserID = user_id).exists():
                        user_area =  Profile.objects.get(profileUserID = user_id).profileArea
                    else:
                        user_center = 'Not defined'
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

                data_source = column_graphic_per_time ('Services requested by Area ',period_of_time_selected,  'date', 'number of services', time_values , user_area_services_period)
                graphic_area_services_per_time = FusionCharts("mscolumn3d", "ex6" , "525", "350", "chart-6", "json", data_source)
                services_stats_info ['graphic_area_services_per_time'] = graphic_area_services_per_time.render()

                services_stats_info['period_time']= period_of_time_selected


                # statistics on Requested Level 2 Services

                service_dict ={}
                for service in services_found :
                    service_request_list = service.serviceAvailableService.filter(level=2)
                    for service_requested in service_request_list:
                        service_name = service_requested.availServiceDescription
                        if service_name in service_dict:
                            service_dict [service_name] += 1
                        else:
                            service_dict [service_name] = 1

                #creating the graphic for requested services
                data_source = column_graphic_dict('Requested Services:', 'level 2 ', '', '','fint',service_dict)
                graphic_req_l2_services = FusionCharts("column3d", "ex7" , "800", "375", "chart-7", "json", data_source)
                services_stats_info ['graphic_req_l2_services'] = graphic_req_l2_services.render()

                # statistics on Requested Level 3 Services

                service_dict ={}
                for service in services_found :
                    service_request_list = service.serviceAvailableService.filter(level=3)
                    for service_requested in service_request_list:
                        service_name = service_requested.availServiceDescription
                        if service_name in service_dict:
                            service_dict [service_name] += 1
                        else:
                            service_dict [service_name] = 1

                #creating the graphic for requested services
                data_source = column_graphic_dict('Requested Services:', 'level 3 ', '', '','fint',service_dict)
                graphic_req_l3_services = FusionCharts("column3d", "ex8" , "800", "375", "chart-8", "json", data_source)
                services_stats_info ['graphic_req_l3_services'] = graphic_req_l3_services.render()





                return render (request, 'iSkyLIMS_drylab/statsByServicesRequest.html', {'services_stats_info':services_stats_info})

            else:
                return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['There are no services created by ', 'For the time of period of between:',
                                                                start_date , 'and', end_date]})
    else:
        form = ByServicesRequest()
    return render(request, 'iSkyLIMS_drylab/statsByServicesRequest.html', {'form':form})





@login_required
def stats_by_samples_processed (request):
    if request.user.is_authenticated:
        try:
            groups = Group.objects.get(name = drylab_config.SERVICE_MANAGER)
            if groups not in request.user.groups.all():
                return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
        except:
            return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
    else:
        #redirect to login webpage
        return redirect ('/accounts/login')
    if request.method == 'POST':
        form = ByDateUserStats(data=request.POST)
        if form.is_valid():
            # validate the input data in the form
            start_date = form['start_date'].data
            end_date = form['end_date'].data
        return render (request,'django_utils/info_page.html', {'content':['Under construction page']})

    else:
        form = BySampleProcessed()
        return render(request, 'iSkyLIMS_drylab/statsBySamplesProcessed.html', {'form':form})



@login_required
def stats_time_delivery (request):
    if request.user.is_authenticated:
        try:
            groups = Group.objects.get(name = drylab_config.SERVICE_MANAGER)
            if groups not in request.user.groups.all():
                return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
        except:
            return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
    else:
        #redirect to login webpage
        return redirect ('/accounts/login')
    if request.method == 'POST':
        form = ByDateUserStats(data=request.POST)
        if form.is_valid():
            # validate the input data in the form
            start_date = form['start_date'].data
            end_date = form['end_date'].data
        return render (request,'django_utils/info_page.html', {'content':['Under construction page']})

    else:
        form = TimeDelivery()
        return render(request, 'iSkyLIMS_drylab/statsByDateUser.html', {'form':form})


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
    if not request.user.is_authenticated:
        return redirect ('/accounts/login')
    if request.user.username != 'admin':
        return redirect('')

    user_connected = {}
    if get_current_users().exists():
        user_list_connected = get_current_users()
        user_data = []
        for user in user_list_connected:
            user_data.append([user.username, user.first_name, user.last_name, user.email])

        user_connected['user_data']= user_data

        user_connected['number_of_users'] = user_list_connected.count()
    return render (request, 'iSkyLIMS_drylab/openSessions.html', {'user_connected': user_connected })

@login_required
def user_login (request):
    if not request.user.is_authenticated:
        return redirect ('/accounts/login')
    if request.user.username != 'admin':
        return redirect('')

    user_data = []
    login_data = {}
    user_list = User.objects.all().order_by('-last_login')
    for user in user_list:
        user_data.append([user.username, user.first_name, user.last_name, user.email, user.last_login])
    login_data['user_data'] = user_data

    return render(request, 'iSkyLIMS_drylab/userLogin.html', {'login_data': login_data})

@login_required

def configuration_test (request):
    # check user privileges
    if request.user.is_authenticated:

        if not request.user.is_staff or not request.user.is_superuser:
            return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
    else:
        #redirect to login webpage
        return redirect ('/accounts/login')

    if request.method=='POST' and request.POST['action'] == 'basicTest':
        test_results = {}
        test_results['basic_checks_ok'] = 'OK'
        drylab_config_file = os.path.join(settings.BASE_DIR, 'iSkyLIMS_drylab', 'drylab_config.py')
        # check if access to databases are defined
        try:
            access_db = Service.objects.all()
            test_results['database_access'] = 'OK'
        except:
            test_results['database_access'] = 'NOK'

        # check if available services are defined

        list_available_services = AvailableService.objects.all()

        if len(list_available_services) == 0:
            test_results['services'] = ('Available services', 'NOK')
        else:
            test_results['services'] = ('Available services', 'OK')
        test_results['iSkyLIMS_settings'] = get_iSkyLIMS_settings()
        test_results['config_file'] = get_config_file(drylab_config_file)

        conn = open_samba_connection()
        if conn :
            test_results['samba_connection'] = 'OK'
        else:
            test_results['samba_connection'] = 'NOK'
        test_results['basic_checks_ok'] = 'OK'
        #if test_results['config_file']  and test_results['attr_files']  and test_results['database_access'] and test_results['samba_connection']:

        for result in test_results :
            if test_results[result] == 'NOK':
                test_results['basic_checks_ok'] = 'NOK'
                break

        return render (request,'iSkyLIMS_drylab/ConfigurationTest.html', {'test_results': test_results})
    elif request.method=='POST' and request.POST['action'] == 'resolutionTest':
        if 'Delete' in request.POST :
            delete_test_service ('SRVTEST-IIER001')
            return render(request,'iSkyLIMS_drylab/ConfigurationTest.html')

        resolution_results = {}
        service_requested = 'SRVTEST-IIER001'
        resolution_results['CreateService'], result = create_service_test(service_requested)

        if result == 'NOK' :
            resolution_results['create_service_ok'] = 'NOK'
            return render (request,'iSkyLIMS_drylab/ConfigurationTest.html', {'resolution_results': resolution_results})

        else:
            resolution_results['create_service_ok'] = 'OK'
            resolution_number = 'SRVTEST-IIER001.1'
            resolution_results ['resolution_test'] = create_resolution_test (resolution_number, service_requested)
            resolution_results['create_resolution_ok'] = 'OK'

            resolution_results['completed_ok']= 'OK'
            for result in resolution_results ['resolution_test']:
                if result[1] == 'NOK':
                    resolution_results['completed_ok']= 'NOK'
                break

            #service_request_file = os.path.join (settings.BASE_DIR, drylab_config.OUTPUT_DIR_TEMPLATE,str('test_resolution.pdf'))
            #service_file_uploaded = ''

            #create_service_structure (conn, service_request_file, service_file_uploaded, full_service_path, resolution_file)

            return render (request,'iSkyLIMS_drylab/ConfigurationTest.html', {'resolution_results': resolution_results})
    else:
        return render(request,'iSkyLIMS_drylab/ConfigurationTest.html')

@login_required
def define_pipeline_service(request):
    if request.user.is_authenticated:
        if not is_service_manager(request):
            return render (request,'iSkyLIMS_drylab/error_page.html', {'content':drylab_config.ERROR_USER_NOT_ALLOWED })
    else:
        return redirect ('/accounts/login')
    data_pipeline = get_data_form_pipeline()
    if request.method == 'POST' and request.POST['action'] == 'servicePipeline':
        pipeline_data_form = analyze_input_pipelines(request)
        if pipeline_version_exists(request.POST['pipelineName'], request.POST['pipelineVersion']):
            error_message = drylab_config.ERROR_PIPELINE_ALREADY_EXISTS
            data_pipeline.update(pipeline_data_form)
            return render(request,'iSkyLIMS_drylab/definePipelineService.html', {'data_pipeline': data_pipeline,'error_message': error_message})
        pipeline_data_form = analyze_input_pipelines(request)
        new_pipeline = Pipelines.objects.create_pipeline(pipeline_data_form)
        if 'additional_parameters' in  pipeline_data_form :
            store_pipeline_actions(new_pipeline, pipeline_data_form['additional_parameters'])
        defined_service_pipeline = get_pipeline_data_to_display(pipeline_data_form)

        set_default_service_pipeline(new_pipeline)
        return render(request,'iSkyLIMS_drylab/definePipelineService.html', {'defined_service_pipeline': defined_service_pipeline})

    return render(request,'iSkyLIMS_drylab/definePipelineService.html', {'data_pipeline': data_pipeline})

@login_required
def manage_pipelines(request):
    if request.user.is_authenticated:
        if not is_service_manager(request):
            return render (request,'iSkyLIMS_drylab/error_page.html', {'content':drylab_config.ERROR_USER_NOT_ALLOWED })
    else:
        #redirect to login webpage
        return redirect ('/accounts/login')
    pipelines_data = get_pipelines_for_manage()
    return render(request,'iSkyLIMS_drylab/managePipelines.html', {'pipelines_data': pipelines_data})

@login_required
def detail_pipeline(request,pipeline_id):
    if request.user.is_authenticated:
        if not is_service_manager(request):
            return render (request,'iSkyLIMS_drylab/error_page.html', {'content':drylab_config.ERROR_USER_NOT_ALLOWED })
    else:
        #redirect to login webpage
        return redirect ('/accounts/login')
    detail_pipelines_data = get_detail_pipeline_data(pipeline_id)
    return render(request,'iSkyLIMS_drylab/detailPipeline.html', {'detail_pipelines_data': detail_pipelines_data})
