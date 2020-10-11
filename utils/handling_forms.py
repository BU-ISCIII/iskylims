import json
from datetime import datetime
from django.conf import settings
from django.contrib.auth.models import Group, User

from iSkyLIMS_drylab import drylab_config
from iSkyLIMS_drylab.models import *
from iSkyLIMS_drylab.forms import *

#from iSkyLIMS_wetlab.models import Projects
from iSkyLIMS_drylab.utils.drylab_common_functions import *
#### API from Wetlab ######
try :
	from iSkyLIMS_wetlab.utils.api.wetlab_api import *
	wetlab_api_available = True
except:
	wetlab_api_available = False

def get_add_resolution_data_form(form_data):
	'''
	Description:
		The function extract the user form information and store it in a dictionary
	Input:
		form_data		# contains the user data form
	Constants:
		HEADING_ADDITIONAL_RESOLUTION_PARAMETERS
	Functions:
		rr
	Return:
		resolution_data_form
	'''
	resolution_data_form = {}

	resolution_data_form['service_id'] = form_data['service_id']
	resolution_data_form['resolutionEstimatedDate'] = datetime.datetime.strptime(form_data['resolutionEstimatedDate'],'%Y-%m-%d').date()
	resolution_data_form['acronymName'] = form_data['acronymName']
	resolution_data_form['resolutionAsignedUser'] = form_data['resolutionAsignedUser']
	resolution_data_form['serviceAccepted'] = form_data['serviceAccepted']
	resolution_data_form['resolutionNotes'] = form_data['resolutionNotes']

	# get additional parameters
	if 'parameters_data' in form_data:
		json_data = json.loads(form_data['parameters_data'])
		additional_parameters = []
		for row_index in range(len(json_data)) :
			if json_data[row_index][0] == '':
				continue
			parameter = {}
			for i in range(len(drylab_config.HEADING_ADDITIONAL_RESOLUTION_PARAMETERS)):
				parameter[drylab_config.HEADING_ADDITIONAL_RESOLUTION_PARAMETERS[i]] = json_data[row_index][i]
			additional_parameters.append(parameter)
		resolution_data_form['additional_parameters'] = additional_parameters
	return resolution_data_form


def prepare_form_data_service_internal_sequencing (request_user):
	'''
	Description:
		The function get the information to display in the internal sequencing form
	Input:
		request_user      # user instance who request the service
	Functions:
		get_user_projects	# API from iSkyLIMS_wetlab located at file wetlab_api
	Return:
		form
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
	#import pdb; pdb.set_trace()
	return display_service

def prepare_form_data_internal_sequencing (request_user):
	'''
	Description:
		The function get the information to display in the internal sequencing form
	Input:
		request_user      # user instance who request the service
	Functions:
		get_user_projects	# API from iSkyLIMS_wetlab located at file wetlab_api
	Return:
		form
	'''

	form = ServiceRequestFormInternalSequencing()
	# getting projects from user sharing list
	sharing_list = []
	user_groups = request_user.groups.values_list('name',flat=True)

	for user in user_groups :
		if User.objects.filter(username__exact = user).exists():
			sharing_list.append(User.objects.get(username__exact = user).id)
	sharing_list.append(request_user.id)
	if wetlab_api_available :
		form.fields['serviceProjects'].choices = get_user_projects(sharing_list)
		form.fields['serviceProjects'].label='User Projects'

	return form

def prepare_form_data_add_resolution(service_id):
	'''
	Description:
		The function collect additional info and then save the form
	Input:
		service_id		# id of the service
	Return:
	 	form_data
	'''
	form_data = {}
	service_obj= Service.objects.get(pk__exact = service_id)
	form_data['service_number'] = service_obj.get_service_request_number()

	if Resolution.objects.filter(resolutionServiceID = service_obj).exists():
		existing_resolution = Resolution.objects.filter(resolutionServiceID = service_obj).last()
		form_data['resolutionFullNumber'] = existing_resolution.get_resolution_number()
	users = User.objects.filter( groups__name = drylab_config.SERVICE_MANAGER)
	form_data['assigned_user'] =[]
	for user in users:
		form_data['assigned_user'].append([user.pk,user.username])
	form_data['heading'] = drylab_config.HEADING_ADDITIONAL_RESOLUTION_PARAMETERS
	form_data['service_id'] = service_id
	return form_data

def save_service_request_form(form, user, unit):
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
	# Create but dont save for following modif
	new_service = form.save(commit=False)
	new_service.serviceStatus = "recorded"
	new_service.serviceSeqCenter = unit

	new_service.serviceUserId = user
	new_service.serviceRequestInt = increment_service_number(user)
	new_service.serviceRequestNumber = create_service_id(new_service.serviceRequestInt,user)
	# Save the new instance
	new_service.save()
	# Save the many-to-many data for the form
	form.save_m2m()
	return new_service

def store_projects_from_form(project_list, service):
	'''
	Description:
		The function get the project name using the api to wetlab and store the
		service project on database
	Input:
		project_list      # List of project ids
		service		# service instance
	Functions:
		get_user_project_name	# API from iSkyLIMS_wetlab located at file wetlab_api
	Return:
		new_service		# recorded instance of the form
	'''

	project_names = get_user_project_name(project_list)
	created_projects = []
	for i in range(len(project_list)):
		data = {}
		data['projectService'] = service
		data['externalProjectKey'] = project_list[i]
		data['externalProjectName'] = project_names[i]
		created_projects.append(RequestedProjectInServices.objects.create_request_project_service(data))
	return created_projects
