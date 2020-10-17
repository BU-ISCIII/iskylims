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
