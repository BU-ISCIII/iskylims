# -*- coding: utf-8 -*-
from django.shortcuts import get_object_or_404, render, redirect
from .models import *
#from .forms import *
from .utils.graphics import *
import os, re

from django.contrib.auth.decorators import login_required
from django.contrib.auth.models import Group
from django.conf import settings
# from django.template.loader import render_to_string
from django_utils.fusioncharts.fusioncharts import FusionCharts
# from django.core.mail import send_mail
from django.http import HttpResponse
from django.core.files.storage import FileSystemStorage

from datetime import date, datetime
import statistics
from iSkyLIMS_drylab import drylab_config
from iSkyLIMS_drylab.utils.testing_drylab_configuration import *
from iSkyLIMS_drylab.utils.drylab_common_functions import *
from iSkyLIMS_drylab.utils.handling_pipelines import *
from iSkyLIMS_drylab.utils.handling_request_services import *
from iSkyLIMS_drylab.utils.handling_resolutions import *
from iSkyLIMS_drylab.utils.handling_deliveries import *
#from iSkyLIMS_drylab.utils.handling_forms import *
from iSkyLIMS_drylab.utils.handling_multiple_files import *
from iSkyLIMS_core.utils.generic_functions import send_test_email, get_email_data

@login_required
def index(request):
    service_list = {}
    if Service.objects.filter(serviceStatus__exact = 'recorded').exists():
        r_service_objs =  Service.objects.filter(serviceStatus__exact = 'recorded').order_by('serviceCreatedOnDate')
        service_list['recorded'] = []
        for r_service_obj in r_service_objs:
            s_info = r_service_obj.get_service_information()
            s_info.append(r_service_obj.get_service_requested_user())
            service_list['recorded'].append(s_info)

    if Service.objects.all().exclude(serviceStatus__exact = 'delivered').exclude(serviceOnApprovedDate = None).exists():
        ongoing_services_objs = Service.objects.all().exclude(serviceStatus__exact = 'delivered').exclude(serviceOnApprovedDate = None).order_by('serviceOnApprovedDate')
        service_list['ongoing'] = []
        for ongoing_services_obj in ongoing_services_objs:
            s_info = ongoing_services_obj.get_service_information()
            s_info.append(ongoing_services_obj.get_delivery_date())
            service_list['ongoing'].append(s_info)

    return render(request, 'iSkyLIMS_drylab/index.html',{'service_list': service_list})



@login_required
def configuration_email(request):
    if request.user.username != 'admin':
        return redirect('/wetlab')
    email_conf_data = get_email_data()
    email_conf_data['EMAIL_ISKYLIMS'] = get_configuration_from_database('EMAIL_FOR_NOTIFICATIONS')
    if request.method == 'POST' and (request.POST['action']=='emailconfiguration'):
        result_email = send_test_email(request.POST)
        if result_email != 'OK':
            email_conf_data = get_email_data()
            email_conf_data['EMAIL_ISKYLIMS'] = request.POST['EMAIL_ISKYLIMS']
            email_conf_data['test_email'] = request.POST['test_email']
            return render(request, 'iSkyLIMS_drylab/configurationEmail.html',{'ERROR':result_email, 'email_conf_data': email_conf_data})
        save_database_configuration_value('EMAIL_FOR_NOTIFICATIONS', request.POST['EMAIL_ISKYLIMS'])
        return render(request, 'iSkyLIMS_drylab/configurationEmail.html',{'succesful_settings':True})
    return render(request, 'iSkyLIMS_drylab/configurationEmail.html',{'email_conf_data': email_conf_data})

@login_required
def request_sequencing_service(request):

	if request.POST and request.FILES :
		if 'file' in request.FILES:
			data = get_and_save_service_file(request)
			response = JSONResponse(data, mimetype='application/json')
			response['Content-Disposition'] = 'inline; filename=files.json'
			return response
		else:
			service_data_information = prepare_form_data_request_service_sequencing(request)
			return render(request,'iSkyLIMS_drylab/requestSequencingService.html',{'service_data_information':service_data_information})

	if request.method == 'POST' and request.POST['subAction'] == 'createservice':
		# check that at some services have been requested
		if len(request.POST.getlist('RequestedServices')) == 0 :
			service_data_information = prepare_form_data_request_service_sequencing(request)
			error_message = drylab_config.ERROR_NO_SERVICES_ARE_SELECTED
			return render(request,'iSkyLIMS_drylab/requestSequencingService.html',{'service_data_information':service_data_information,
									'error_message':error_message})

		new_service = create_new_save_sequencing_service_request(request)
		sample_stored = stored_samples_for_sequencing_request_service(request.POST, new_service)
		if 'files' in request.POST:
			add_files_to_service(request.POST.getlist('files'), new_service)
		## Send mail to user and drylab notification email
		email_data = {}
		if 'requestedForUserid' in request.POST:
			user_obj= User.objects.filter(pk__exact = request.POST['requestedForUserid']).last()
			email_data['user_name'] = user_obj.username
			email_data['user_email'] = user_obj.email
		else:
			email_data['user_email'] = request.user.email
			email_data['user_name'] = request.user.username
		email_data['service_number'] = new_service.get_service_request_number()
		email_result = send_service_creation_confirmation_email(email_data)
		# PDF preparation file for confirmation of service request
		confirmation_result = {}
		confirmation_result['download_file'] = create_service_pdf_file(new_service.get_service_request_number(), request.build_absolute_uri())

		service_request_number = new_service.get_service_request_number()
		confirmation_result['text'] = list(map(lambda st: str.replace(st, 'SERVICE_NUMBER', service_request_number), drylab_config.CONFIRMATION_TEXT_MESSAGE))
		if len(sample_stored) > 0 :
			confirmation_result['samples'] = sample_stored
		if email_result != 'OK':
			return render(request,'iSkyLIMS_drylab/requestSequencingService.html',{'confirmation_result':confirmation_result, 'error_message':email_result})
		return render(request,'iSkyLIMS_drylab/requestSequencingService.html',{'confirmation_result':confirmation_result})
	else:
		service_data_information = prepare_form_data_request_service_sequencing(request)
		return render(request,'iSkyLIMS_drylab/requestSequencingService.html',{'service_data_information':service_data_information})


@login_required
def counseling_request(request):

    if request.method == 'POST' and request.POST['action'] == 'createService':
        # check that at some services have been requested
        if len(request.POST.getlist('RequestedServices')) == 0 :
            service_data_information = prepare_form_data_request_counseling_service()
            error_message = drylab_config.ERROR_NO_SERVICES_ARE_SELECTED
            return render(request,'iSkyLIMS_drylab/requestCounselingService.html',{'service_data_information':service_data_information,
                           'error_message':error_message})

        new_service = create_new_save_counseling_service_request(request)

        email_data = {}
        email_data['user_email'] = request.user.email
        email_data['user_name'] = request.user.username
        email_data['service_number'] = new_service.get_service_request_number()
        send_service_creation_confirmation_email(email_data)
        confirmation_result = {}
        confirmation_result['download_file'] = create_service_pdf_file(new_service.get_service_request_number(), request.build_absolute_uri())

        service_request_number = new_service.get_service_request_number()
        confirmation_result['text'] = list(map(lambda st: str.replace(st, 'SERVICE_NUMBER', service_request_number), drylab_config.CONFIRMATION_TEXT_MESSAGE))

        return render(request,'iSkyLIMS_drylab/requestSequencingService.html',{'confirmation_result':confirmation_result})

    else:
        service_data_information = prepare_form_data_request_counseling_service()
        return render(request,'iSkyLIMS_drylab/requestCounselingService.html',{'service_data_information':service_data_information})


@login_required
def infrastructure_request(request):

    if request.method == 'POST' and request.POST['action'] == 'createService':
        # check that at some services have been requested
        if len(request.POST.getlist('RequestedServices')) == 0 :
            service_data_information = prepare_form_data_request_infrastructure_service()
            error_message = drylab_config.ERROR_NO_SERVICES_ARE_SELECTED
            return render(request,'iSkyLIMS_drylab/requestInfrastructureService.html',{'service_data_information':service_data_information,
                           'error_message':error_message})

        new_service = create_new_save_infrastructure_service_request(request)

        email_data = {}
        email_data['user_email'] = request.user.email
        email_data['user_name'] = request.user.username
        email_data['service_number'] = new_service.get_service_request_number()
        send_service_creation_confirmation_email(email_data)
        confirmation_result = {}
        confirmation_result['download_file'] = create_service_pdf_file(new_service.get_service_request_number(), request.build_absolute_uri())

        service_request_number = new_service.get_service_request_number()
        confirmation_result['text'] = list(map(lambda st: str.replace(st, 'SERVICE_NUMBER', service_request_number), drylab_config.CONFIRMATION_TEXT_MESSAGE))

        return render(request,'iSkyLIMS_drylab/requestInfrastructureService.html',{'confirmation_result':confirmation_result})

    else:
        service_data_information = prepare_form_data_request_counseling_service()
        return render(request,'iSkyLIMS_drylab/requestInfrastructureService.html',{'service_data_information':service_data_information})



@login_required
def display_service (request, service_id):
    if not request.user.is_authenticated:
        #redirect to login webpage
        return redirect ('/accounts/login')
    if Service.objects.filter(pk=service_id).exists():
        service_manager = is_service_manager(request)
        display_service_details = get_service_information(service_id, service_manager)
        return render (request,'iSkyLIMS_drylab/display_service.html',{'display_service': display_service_details})
    else:
	    return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['The service that you are trying to get does not exist ','Contact with your administrator .']})


@login_required
def search_service (request):
    if  not request.user.is_authenticated:
        #redirect to login webpage
        return redirect ('/accounts/login')
    services_search_list = {}

    center_list_abbr = []
    center_availables = Center.objects.all().order_by ('centerAbbr')
    for center in center_availables:
        center_list_abbr.append (center.centerAbbr)
    services_search_list ['centers'] = center_list_abbr
    services_search_list ['status'] = STATUS_CHOICES

    if 'iSkyLIMS_wetlab' in settings.INSTALLED_APPS:
        services_search_list['wetlab_app'] = True
    if not is_service_manager(request):
        services_search_list['username'] = request.user.username

    if request.method == 'POST' and request.POST['action'] == 'searchservice':

        service_number_request = request.POST['servicenumber']
        service_state = request.POST['servicestate']
        start_date=request.POST['startdate']
        end_date=request.POST['enddate']
        center = request.POST['center']

        user_name = request.POST['username']
        handeld_user = request.POST['userhandled']
        # Data related to samples,and run from wetlab application
        sample_name =  request.POST['sampleName'] if 'sampleName' in request.POST else ''
        project_name = request.POST['projectName'] if 'projectName' in request.POST else ''
        run_name = request.POST['runName'] if 'runName' in request.POST else ''

        if service_number_request == '' and service_state == '' and start_date == '' and end_date == '' and center == '' and sample_name == ''  and user_name =='' and project_name == '' and run_name == '' and handeld_user == '':
            return render( request,'iSkyLIMS_drylab/searchService.html',{'services_search_list': services_search_list })

        ### check the right format of start and end date
        if (request.POST['startdate'] != '' and not check_valid_date_format(request.POST['startdate'])) or  (request.POST['enddate'] != '' and not check_valid_date_format(request.POST['enddate'])):
            error_message = drylab_config.ERROR_INCORRECT_FORMAT_DATE
            return render( request,'iSkyLIMS_drylab/searchService.html',{'services_search_list': services_search_list , 'ERROR':error_message})

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
        else:
            services_found = Service.objects.all()

        if service_state != '':
            services_found = services_found.filter(serviceStatus__exact = service_state)
        if start_date !='' and end_date != '':
            services_found = services_found.filter(serviceCreatedOnDate__range=(start_date, end_date))
        if start_date !='' and end_date == '':
            services_found = services_found.filter(serviceCreatedOnDate__gte = start_date)
        if start_date =='' and end_date != '':
            services_found = services_found.filter(serviceCreatedOnDate__lte = end_date)
        if center != '':
            services_found = services_found.filter(serviceRequestNumber__icontains  = center)
        if  user_name != '':
            services_found = services_found.filter(serviceUserId__username__iexact  = user_name)
        if handeld_user != '':

            if not Resolution.objects.filter(resolutionAsignedUser__username__icontains = user_name).exists():
                error_message = drylab_config.ERROR_NO_MATCHES_FOUND_FOR_YOUR_SERVICE_SEARCH
                return render( request,'iSkyLIMS_drylab/searchService.html',{'services_search_list': services_search_list , 'ERROR':error_message})
            services_handled_by = []
            services_handled_by_objs = Resolution.objects.filter(resolutionAsignedUser__username__icontains = user_name)
            for services_handled_by_obj in services_handled_by_objs:
                services_handled_by.append(services_handled_by_obj.get_service_name())
            services_found = services_found.filter(serviceRequestNumber__in = services_handled_by)

        if project_name != '' or run_name != '' or sample_name != '':
            if project_name != '':
                samples_in_services = RequestedSamplesInServices.objects.filter(projectName__icontains = project_name, samplesInService__in = service_found)
            else:
                samples_in_services =  RequestedSamplesInServices.objects.filter(samplesInService__in = service_found)
            if run_name != '':
                samples_in_services = samples_in_services.filter(runName__icontains = run_name)
            if sample_name != '':
                samples_in_services = samples_in_services.filter(sampleNane__icontains = sample_name)
            service_list = []
            for project_in_service in project_in_services:
                service_list.append(project_in_service.samplesInService.pk)
            services_found = services_found.filter(pk__in = service_list)

        if len(services_found) == 0 :
            error_message = drylab_config.ERROR_NO_MATCHES_FOUND_FOR_YOUR_SERVICE_SEARCH
            return render( request,'iSkyLIMS_drylab/searchService.html',{'services_search_list': services_search_list , 'ERROR':error_message})
        #If only 1 service mathes the user conditions, then get the user information
        if len(services_found) == 1 :
            redirect_page = '/drylab/display_service=' + str(services_found[0].get_service_id())
            return redirect (redirect_page)
        else:
            display_multiple_services ={}
            s_list  = {}
            for service_item in services_found:
                data = []
                service_id = service_item.get_service_id()
                data.append(service_item.get_service_request_number())
                data.append(service_item.get_service_state())
                data.append(service_item.get_service_dates())
                data.append(service_item.get_service_request_center())
                data.append(get_projects_in_requested_samples(service_item))
                data.append(get_run_in_requested_samples(service_item))
                s_list [service_id]=[data]
            display_multiple_services['s_list'] = s_list
            return render (request,'iSkyLIMS_drylab/searchService.html', {'display_multiple_services': display_multiple_services})

    return render( request,'iSkyLIMS_drylab/searchService.html',{'services_search_list': services_search_list })


@login_required
def pending_services (request):
    if request.user.is_authenticated:
        if not is_service_manager(request):
            return render (request,'iSkyLIMS_drylab/error_page.html', {'content':drylab_config.ERROR_USER_NOT_ALLOWED })
    else:
        #redirect to login webpage
        return redirect ('/accounts/login')

    pending_services_details = get_pending_services_information()
    user_pending_services = get_user_pending_services_information(request.user.username)
    return render (request, 'iSkyLIMS_drylab/pendingServices.html', {'pending_services': pending_services_details, 'user_pending_services':user_pending_services})

@login_required
def add_resolution (request):
    if request.user.is_authenticated:
        if not is_service_manager(request):
            return render (request,'iSkyLIMS_drylab/error_page.html', {'content':drylab_config.ERROR_USER_NOT_ALLOWED })
    else:
        #redirect to login webpage
        return redirect ('/accounts/login')
    if request.method != 'POST' or not 'service_id' in request.POST:
        return render (request, 'iSkyLIMS_drylab/error_page.html', {'content':drylab_config.ERROR_SERVICE_ID_NOT_FOUND})

    if request.method == "POST" and request.POST['action'] == 'addResolutionService' :

        resolution_data_form = get_add_resolution_data_form(request.POST)
        new_resolution = create_new_resolution(resolution_data_form)
        if 'pipeline_ids' in resolution_data_form:
            add_pipelines_to_resolution(new_resolution, resolution_data_form['pipeline_ids'])
        # create a new resolution to be added to the service folder including the path where file is stored

        pdf_file = create_resolution_pdf_file(new_resolution, request.build_absolute_uri())
        new_resolution.update_resolution_file(pdf_file)
        #pdf_name = resolution_data_form['resolutionNumber'] + ".pdf"
        #resolution_file = create_pdf(request,information, drylab_config.RESOLUTION_TEMPLATE, pdf_name)

        ## Send email
        email_data = {}
        email_data['user_email'] = request.user.email
        email_data['user_name'] = request.user.username
        email_data['service_number'] = new_resolution.get_service_request_number()
        email_data['status'] = resolution_data_form['serviceAccepted']
        email_data['date'] = resolution_data_form['resolutionEstimatedDate']
        send_resolution_creation_email(email_data)
        created_resolution = {}
        created_resolution['resolution_number'] = resolution_data_form['resolutionNumber']
        ## Display pipeline parameters

        return render(request,'iSkyLIMS_drylab/addResolution.html',{'created_resolution': created_resolution})


    if request.method == "POST" and request.POST['action'] == 'formToaddResolutionService':
        resolution_form_data = prepare_form_data_add_resolution(request.POST)
        return render(request, 'iSkyLIMS_drylab/addResolution.html' , { 'resolution_form_data' : resolution_form_data})
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

@login_required
def add_in_progress (request):
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

    if request.method == "POST" and request.POST['action'] == 'inProgressResolutionService' :
        resolution_id = request.POST['resolution_id']
        if  not check_if_resolution_exists(resolution_id):
            error_message = drylab_config.ERROR_RESOLUTION_DOES_NOT_EXISTS
            return render (request,'iSkyLIMS_drylab/error_page.html', {'content':error_message})

        resolution_obj = get_resolution_obj_from_id(resolution_id)
        resolution_obj.update_resolution_in_progress_date()
        resolution_number = resolution_obj.get_resolution_number()
        service_obj = resolution_obj.get_service_obj()
        # check if services can change to "in progress"
        if allow_to_service_update_state (resolution_obj, 'in_progress'):
            # update the service status and in_porgress date
            service_obj.update_service_status('in_progress')
        email_data = {}
        email_data['user_email'] = service_obj.get_user_email()
        email_data['user_name'] = service_obj.get_username()
        email_data['resolution_number'] = resolution_number
        send_resolution_in_progress_email(email_data)
        in_progress_resolution = {}
        in_progress_resolution['resolution_number'] = resolution_number
        return render (request,'iSkyLIMS_drylab/addInProgress.html',{'in_progress_resolution':in_progress_resolution})

    error_message = drylab_config.ERROR_RESOLUTION_DOES_NOT_EXISTS
    return render (request,'iSkyLIMS_drylab/error_page.html', {'content':error_message})

@login_required
def add_delivery (request ):
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
    if request.method == 'POST' and request.POST['action'] == 'deliveryResolutionService':
        delivery_data = prepare_delivery_form (request.POST['resolution_id'])

        return render (request, 'iSkyLIMS_drylab/addDelivery.html', {'delivery_data':delivery_data})

    if request.method == 'POST' and request.POST['action'] == 'addDeliveryResolution':

        if (request.POST['startdate'] != '' and not check_valid_date_format(request.POST['startdate'])) or  (request.POST['enddate'] != '' and not check_valid_date_format(request.POST['enddate'])):
            delivery_data = prepare_delivery_form (request.POST['resolution_id'])
            error_message = drylab_config.ERROR_INCORRECT_FORMAT_DATE
            return render (request, 'iSkyLIMS_drylab/addDelivery.html', {'delivery_data':delivery_data, 'ERROR': error_message})

        delivery_recorded = store_resolution_delivery (request.POST)
        if delivery_recorded != None :
            email_data = {}
            email_data['user_email'] = request.user.email
            email_data['user_name'] = request.user.username
            email_data['resolution_number'] = delivery_recorded['resolution_number']
            send_delivery_service_email(email_data)
            resolution_obj = delivery_recorded['deliveryResolutionID']
            if allow_to_service_update_state (resolution_obj, 'delivered'):
                service_obj = resolution_obj.get_service_obj()
                service_obj.update_service_status('delivered')
            return render (request, 'iSkyLIMS_drylab/addDelivery.html', {'delivery_recorded': delivery_recorded})

    return render (request,'iSkyLIMS_drylab/error_page.html', {'content':['The resolution that you are trying to upadate does not exists ','Contact with your administrator .']})



@login_required
def stats_by_user (request):
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
    user_list = get_users_requested_services()
    if request.method == 'POST' and request.POST['action'] == 'userStatistics':
        # validate the input data in the form
        user_id = request.POST['userID']
        start_date = request.POST['start_date']
        end_date = request.POST['end_date']

        if start_date != '' and not check_valid_date_format(start_date):
            error_message = drylab_config.ERROR_INCORRECT_FORMAT_DATE
            return render(request, 'iSkyLIMS_drylab/statsByUser.html' , {'user_list': user_list,'ERROR':error_message })
        if end_date != '':
            if not check_valid_date_format(end_date):
                error_message = drylab_config.ERROR_INCORRECT_FORMAT_DATE
                return render(request, 'iSkyLIMS_drylab/statsByUser.html' , {'user_list': user_list,'ERROR':error_message })
        else:
            end_date  = date.today().strftime('%Y-%m-%d')

        if not User.objects.filter(pk__exact = user_id).exists():
            error_message = drylab_config.ERROR_USER_NOT_DEFINED
            return render(request, 'iSkyLIMS_drylab/statsByUser.html' , {'user_list': user_list,'ERROR':error_message })

        service_objs = Service.objects.filter(serviceUserId__exact = user_id).order_by('-serviceRequestNumber')
        if start_date != '':
            service_objs = service_objs.filter(serviceCreatedOnDate__gte = start_date)
        if end_date != '':
            service_objs = service_objs.filter(serviceCreatedOnDate__lte = end_date)
        if len(service_objs) == 0:
            error_message = drylab_config.ERROR_NO_MATCHES_FOUND_FOR_YOUR_SERVICE_SEARCH
            return render(request, 'iSkyLIMS_drylab/statsByUser.html' , {'user_list': user_list,'ERROR':error_message })

        stats_info = {}
        stats_info ['service_by_user'] = []
        stats_info ['user_name'] = User.objects.filter(pk__exact = user_id).last().username
        for service_item in service_objs:
            stats_info ['service_by_user'].append(service_item.get_stats_information())

        # perform calculation time media delivery for user
        if service_objs.filter(serviceStatus__exact = 'delivered').exists():
            delivery_services = service_objs.filter(serviceStatus__exact = 'delivered')
            delivery_time_in_days = []
            for service_item in delivery_services :
                delivery_time_in_days.append(int (service_item.get_time_to_delivery()))

            stats_info['time_mean_for_user']=  format(statistics.mean (delivery_time_in_days), '.2f')
        # preparing graphic for status of the services
        number_of_services = {}

        number_of_services ['RECORDED'] = service_objs.filter(serviceStatus__exact = 'recorded').count()
        number_of_services ['QUEUED'] = service_objs.filter(serviceStatus__exact = 'queued').count()
        number_of_services ['IN PROGRESS'] = service_objs.filter(serviceStatus__exact = 'in_progress').count()
        number_of_services ['DELIVERED'] = service_objs.filter(serviceStatus__exact = 'delivered').count()

        data_source = graphic_3D_pie('Status of Requested Services from:', stats_info ['user_name'], '', '','fint',number_of_services)
        graphic_by_user_date_services = FusionCharts("pie3d", "ex1" , "600", "350", "chart-1", "json", data_source)
        stats_info ['graphic_by_user_date_services'] = graphic_by_user_date_services.render()

        # getting statistics of the created services

        service_dict ={}
        for service_available in service_objs :
            service_list = service_available.serviceAvailableService.filter(level=3)
            for service in service_list:
                service_name = service.availServiceDescription
                if service_name in service_dict:
                    service_dict [service_name] += 1
                else:
                    service_dict [service_name] = 1
        #creating the graphic for requested services
        data_source = column_graphic_dict('Requested Services by:', stats_info ['user_name'], '', '','fint',service_dict)
        graphic_requested_services = FusionCharts("column3d", "ex2" , "550", "350", "chart-2", "json", data_source)
        stats_info ['graphic_requested_services'] = graphic_requested_services.render()

        # getting statistics for requested per time
        service_time_dict ={}
        for service_per_time in service_objs :
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
        data_source = column_graphic_tupla('Requested Services by:', stats_info ['user_name'], '', '','fint',service_time_tupla)
        graphic_date_requested_services = FusionCharts("column3d", "ex3" , "550", "350", "chart-3", "json", data_source)
        stats_info ['graphic_date_requested_services'] = graphic_date_requested_services.render()


        return render (request, 'iSkyLIMS_drylab/statsByUser.html', {'stats_info':stats_info})
    else:
        return render(request, 'iSkyLIMS_drylab/statsByUser.html' , {'user_list': user_list})


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
    if request.method == 'POST' and request.POST['action'] == 'serviceStatistics':
        start_date=request.POST['startdate']
        end_date=request.POST['enddate']
        if start_date != '' and not check_valid_date_format(start_date):
            return render(request, 'iSkyLIMS_drylab/statsByServicesRequest.html')
        if end_date != '':
            if not check_valid_date_format(end_date):
                return render(request, 'iSkyLIMS_drylab/statsByServicesRequest.html')
        else:
            end_date  = date.today().strftime('%Y-%m-%d')
        start_date_format = datetime.strptime(start_date, '%Y-%m-%d')
        end_date_format = datetime.strptime(end_date, '%Y-%m-%d')

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
        #form = ByServicesRequest()
        return render(request, 'iSkyLIMS_drylab/statsByServicesRequest.html')

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
    if request.method == 'POST' and request.POST['action'] == 'definePipeline':


        pipeline_data_form = analyze_input_pipelines(request)
        if pipeline_version_exists(request.POST['pipelineName'], request.POST['pipelineVersion']):
            error_message = drylab_config.ERROR_PIPELINE_ALREADY_EXISTS
            data_pipeline.update(pipeline_data_form)
            return render(request,'iSkyLIMS_drylab/definePipelineService.html', {'data_pipeline': data_pipeline,'error_message': error_message})
        #pipeline_data_form = analyze_input_pipelines(request)
        if request.FILES :
            fs = FileSystemStorage(location = os.path.join(settings.MEDIA_ROOT,drylab_config.PIPELINE_FILE_DIRECTORY))
            pipeline_data_form['filename'] = fs.save(str(pipeline_data_form['pipelineName'] + '_' + pipeline_data_form['pipelineVersion']),  request.FILES['pipelinefile'] )
        else:
            pipeline_data_form['filename'] =''
        new_pipeline = Pipelines.objects.create_pipeline(pipeline_data_form)
        if 'additional_parameters' in  pipeline_data_form :
            store_parameters_pipeline(new_pipeline, pipeline_data_form['additional_parameters'])
        defined_service_pipeline = get_defined_pipeline_data_to_display(new_pipeline)

        #set_default_service_pipeline(new_pipeline)
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



###################### Delete upload service file  #########################
def upload_service_file_delete(request,file_id):
    if request.method == 'DELETE':
        if not check_if_file_is_linked_to_service(file_id):
            delete_service_file(file_id)
            response = JSONResponse(True, mimetype='application/json')
            response['Content-Disposition'] = 'inline; filename=files.json'

            return response
        return HttpResponse(content='Not allowed', status=403, content_type='application/json')
    return HttpResponse(content='data not found', status=410, content_type='application/json')
