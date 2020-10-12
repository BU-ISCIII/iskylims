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
    if 'uploadfile' in request.FILES :
        path = os.path.join(drylab_config.USER_REQUESTED_SERVICE_FILE_DIRECTORY)

        file_name, full_path_file_name = store_file_from_form(request.FILES['uploadfile'], path)
        service_data['serviceFile'] = full_path_file_name
    else:
        service_data['serviceFile'] = None
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
    # Save the many-to-many data for the form
    for av_service_obj in available_service_objs:
        new_service.serviceAvailableService.add(av_service_obj)
    return new_service

def get_available_service_obj_from_id(available_service_id):
    '''
	Description:
		The function get the information to display in the internal sequencing form
	Input:
		request_user      # user instance who request the service
	Return:
		avail_service_obj
    '''
    avail_service_obj = None
    if AvailableService.objects.filter(pk__exact = available_service_id).exists():
        avail_service_obj = AvailableService.objects.filter(pk__exact = available_service_id).last()
    return avail_service_obj

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
