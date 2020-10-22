import json
from datetime import datetime
import os
from django.conf import settings
from django.contrib.auth.models import User

from iSkyLIMS_drylab import drylab_config
from iSkyLIMS_drylab.models import *
from iSkyLIMS_core.models import Samples

from iSkyLIMS_drylab.utils.drylab_common_functions import *


#### API from Wetlab ######
try :
	from iSkyLIMS_wetlab.utils.api.wetlab_api import *
	wetlab_api_available = True
except:
	wetlab_api_available = False


def check_service_id_exists(service_id):
    '''
    Description:
        The function check if service id exists
    Input:
        service_id      # id of the service
    Return:
        True if service id exists
    '''
    if Service.objects.filter(pk=service_id).exists():
        return True
    else:
        return False

def create_new_save_service_request(request):
    '''
    Description:
    	The function collect additional info and then save the form
    Input:
    	form      # form instance
    	user		# request user
    	unit 		# department who request the service
    Functions:
    	increment_service_number	# located at utils/drylab_common_functions
    	create_service_id			# located at utils/drylab_common_functions
    Return:
	new_service		# recorded instance of the form
    '''
    service_data = {}
    available_service_list = request.POST.getlist('RequestedServices')
    available_service_objs = []
    for av_service in available_service_list:
        available_service_objs.append(get_available_service_obj_from_id(av_service))

    if request.POST['center'] != '':
        service_data['serviceSeqCenter'] = request.POST['center']
    else:
        service_data['serviceSeqCenter'] = drylab_config.INTERNAL_SEQUENCING_UNIT
    service_data['serviceNotes'] = request.POST['description']
    #new_service = form.save(commit=False)
    service_data['serviceStatus'] = "recorded"
    service_data ['serviceUserId'] = request.user
    service_data['serviceRequestInt'] = increment_service_number(request.user.id)
    service_data['serviceRequestNumber'] = create_service_id(service_data['serviceRequestInt'],request.user.id)
	# Save the new service
    new_service = Service.objects.create_service(service_data)

	# Save files
    if 'uploadfile' in request.FILES :
        path = os.path.join(drylab_config.USER_REQUESTED_SERVICE_FILE_DIRECTORY)
        files = request.FILES.getlist('uploadfile')
        for file in files:
            file_name, full_path_file_name = store_file_from_form(file, path)
            new_file = RequestedServiceFile.objects.create_request_service_file(file_data)
			# service_data['serviceFile'] = full_path_file_name

    # Save the many-to-many data for the form
    for av_service_obj in available_service_objs:
        new_service.serviceAvailableService.add(av_service_obj)
    return new_service

def get_available_children_services_and_id(all_tree_services):
	'''
	Description:
		The function get the children available services from a query of service
	Input:
		all_tree_services  # queryset of available service
	Return:
		children_service
    '''
	children_services = []
	for t_services in all_tree_services:
		if t_services.get_children():
			continue
		children_services.append ([t_services.id, t_services.get_service_description()])
	return children_services

def get_available_service_obj_from_id(available_service_id):
    '''
	Description:
		The function get the available service obj  from the id
	Input:
		available_service_id  # id of the available service
	Return:
		avail_service_obj
    '''
    avail_service_obj = None
    if AvailableService.objects.filter(pk__exact = available_service_id).exists():
        avail_service_obj = AvailableService.objects.filter(pk__exact = available_service_id).last()
    return avail_service_obj

def get_service_obj_from_id(service_id):
    '''
	Description:
		The function get the  service obj  from the id
	Input:
		service_id  # id of the  service
	Return:
		service_obj
    '''
    service_obj = None
    if Service.objects.filter(pk__exact = service_id).exists():
        service_obj = Service.objects.filter(pk__exact = service_id).last()
    return service_obj

def get_service_information (service_id):
    service_obj = get_service_obj_from_id(service_id)
    display_service_details = {}

    #text_for_dates = ['Service Date Creation', 'Approval Service Date', 'Rejected Service Date']
    service_dates = []
    display_service_details['service_name'] = service_obj.get_service_request_number()
    display_service_details['service_id'] = service_id
    # get the list of projects
    #projects_in_service = {}
    if RequestedSamplesInServices.objects.filter(samplesInService = service_obj).exists():
        samples_in_service = RequestedSamplesInServices.objects.filter(samplesInService = service_obj)
        display_service_details['samples'] = []
        for sample in samples_in_service:
            display_service_details['samples'].append([sample.get_external_sample_id(), sample.get_external_sample_name()])

    #projects_class = service.serviceProjectNames.all()
    # for project in projects_class:
    #     project_id = project.id
    #     projects_in_service[project_id]=project.get_requested_project_name()
    #display_service_details['projects'] = projects_in_service
    display_service_details['user_name'] = service_obj.get_service_requested_user()
    user_input_file = service_obj.get_service_file()
    if user_input_file:
        display_service_details['file'] = os.path.join(settings.MEDIA_URL,user_input_file)
    display_service_details['state'] = service_obj.get_service_state()
    display_service_details['service_notes'] = service_obj.get_service_user_notes()
    #dates_for_services = service.get_service_dates()
    # for i in range(len(dates_for_services)):
    #     service_dates.append([text_for_dates[i],dates_for_services[i]])
    display_service_details['service_dates'] = zip (drylab_config.HEADING_SERVICE_DATES, service_obj.get_service_dates() )
    #display_service_details['service_dates'] = service_dates
    # if display_service_details['state'] != 'approved' and display_service_details['state'] != 'recorded':
        # get the proposal for the delivery date
    if Resolution.objects.filter(resolutionServiceID = service_obj).exists():
        last_resolution = Resolution.objects.filter(resolutionServiceID = service_obj).last()
        display_service_details['resolution_folder'] = last_resolution.get_resolution_number()
        #display_service_details['resolution_folder'] = resolution_folder
        resolution_estimated_date = last_resolution.get_resolution_estimated_date()


    # get all services
    display_service_details['nodes']= service_obj.serviceAvailableService.all()
    display_service_details['children_services'] = get_available_children_services_and_id(display_service_details['nodes'])
    # adding actions fields
    if service_obj.serviceStatus != 'rejected' or service_obj.serviceStatus != 'archived':
        display_service_details['add_resolution_action'] = service_id
        if len(display_service_details['children_services']) > 1:
            display_service_details['multiple_services'] = True
            if Resolution.objects.filter(resolutionServiceID = service_obj).exists():
                resolutions = Resolution.objects.filter(resolutionServiceID = service_obj)
                for resolution in resolutions:
                    pass
            else:
                display_service_details['first_resolution'] = True
    if service_obj.serviceStatus == 'queued':
        resolution_id = Resolution.objects.filter(resolutionServiceID = service_obj).last().id
        display_service_details['add_in_progress_action'] = resolution_id
    if service_obj.serviceStatus == 'in_progress':
        resolution_id = Resolution.objects.filter(resolutionServiceID = service_obj).last().id
        display_service_details['add_delivery_action'] = resolution_id


    if Resolution.objects.filter(resolutionServiceID = service_obj).exists():
        resolution_heading = drylab_config.HEADING_FOR_RESOLUTION_INFORMATION
        resolution_list = Resolution.objects.filter(resolutionServiceID = service_obj).order_by('resolutionState')
        resolution_info =[]
        for resolution_item in resolution_list :
            resolution_info.append([list(zip(resolution_heading,resolution_item.get_resolution_information()))])
        display_service_details['resolutions'] = resolution_info
        #import pdb; pdb.set_trace()
    if Resolution.objects.filter(resolutionServiceID = service_obj).exists():
        resolution_list = Resolution.objects.filter(resolutionServiceID = service_obj)
        delivery_info = []
        for resolution_id in resolution_list :
            if Delivery.objects.filter(deliveryResolutionID = resolution_id).exists():
                delivery = Delivery.objects.get(deliveryResolutionID = resolution_id)
                delivery_info.append([delivery.get_delivery_information()])
                display_service_details['delivery'] = delivery_info

    if service_obj.servicePipelines.all().exists():
        display_service_details['pipelines'] = {}
        display_service_details['pipelines']['heading'] = drylab_config.DISPLAY_NEW_DEFINED_PIPELINE
        display_service_details['pipelines']['services'] = []
        services_pipelines_objs = service_obj.servicePipelines.all()
        for service_pipeline in services_pipelines_objs:
            service_name = service_pipeline.get_pipleline_service()
            display_service_details['pipelines']['services'].append(service_pipeline.get_pipeline_basic())

    return display_service_details




def prepare_form_data_request_service_sequencing (request_user):
    '''
    Description:
    	The function get the information to display in the internal sequencing form
    Input:
    	request_user      # user instance who request the service
    Return:
    	display_service
    '''
    display_service = {}
    # getting projects from user sharing list
    sharing_list = []
    user_groups = request_user.groups.values_list('name',flat=True)
    for user in user_groups :
        if User.objects.filter(username__exact = user).exists():
            display_servicesharing_list.append(User.objects.get(username__exact = user).id)
    sharing_list.append(request_user.id)
    if wetlab_api_available :
        display_service['serviceProjects']= get_user_projects(sharing_list)
        display_service['serviceProjectsHeading']='User Projects'
    display_service['nodes'] = AvailableService.objects.filter(availServiceDescription__exact="Genomic data analysis").get_descendants(include_self=True)

    return display_service

def stored_samples_for_sequencing_request_service(sample_requested, new_service):
    '''
	Description:
		The function get the samples that were selected and store them on database
	Input:
		sample_requested      # table with all samples
        new_service     # service obj
	Functions:
		get_user_projects	# API from iSkyLIMS_wetlab located at file wetlab_api
	Return:
		display_service
	'''
    requested_sample_list = []
    requested_services_table = json.loads(sample_requested)
    heading = ['run_name', 'run_id', 'project_name', 'project_id','sample_name', 'sample_id']
    for row in requested_services_table:

        if row[-1] :
            data = {}
            for i in range(len(heading)):
                data[heading[i]] = row[i]
            data['samplesInService'] = new_service
            if Samples.objects.filter(sampleName__exact = data['sample_name']).exists():
                data['sample'] = Samples.objects.filter(sampleName__exact = data['sample']).last()
            else:
                data['sample'] = None

            req_samp_obj = RequestedSamplesInServices.objects.create_request_sample(data)
            requested_sample_list.append(data['sample_name'])

    return requested_sample_list
