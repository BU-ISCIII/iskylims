from datetime import datetime
from iSkyLIMS_drylab.models import *
from iSkyLIMS_drylab.utils.handling_request_services import *


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


def prepare_form_data_add_resolution(form_data):
    '''
    Description:
    	The function collect additional info and then save the form
    Input:
    	form_data		# contains the user data form
    Return:
     	resolution_form_data
    '''
    resolution_form_data = {}
    service_obj= get_service_obj_from_id(form_data['service_id'])
    resolution_form_data['service_number'] = service_obj.get_service_request_number()
    all_tree_services = service_obj.serviceAvailableService.all()
    children_services = get_available_children_services_and_id(all_tree_services)
    import pdb; pdb.set_trace()
    if len(children_services) > 1 :
        if 'childrenServices' in form_data:
            list_of_ch_services = form_data.getlist('childrenServices')
            if len(list_of_ch_services) != len(children_services):
                selected_children_services = []
                for children in list_of_ch_services :
                    avail_serv_obj = get_available_service_obj_from_id(children)
                    selected_children_services.append([children , avail_serv_obj.get_service_description()])
                resolution_form_data['selected_children_services'] = selected_children_services
    if Resolution.objects.filter(resolutionServiceID = service_obj).exists():
        existing_resolution = Resolution.objects.filter(resolutionServiceID = service_obj).last()
        resolution_form_data['resolutionFullNumber'] = existing_resolution.get_resolution_number()
    users = User.objects.filter( groups__name = drylab_config.SERVICE_MANAGER)
    resolution_form_data['assigned_user'] =[]
    for user in users:
    	resolution_form_data['assigned_user'].append([user.pk,user.username])
    resolution_form_data['heading'] = drylab_config.HEADING_ADDITIONAL_RESOLUTION_PARAMETERS
    resolution_form_data['service_id'] = form_data['service_id']

    return resolution_form_data
