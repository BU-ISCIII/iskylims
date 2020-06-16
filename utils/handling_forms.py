from django.conf import settings
from django.contrib.auth.models import Group
from iSkyLIMS_drylab import drylab_config
from iSkyLIMS_drylab.models import *
from iSkyLIMS_drylab.forms import *

from iSkyLIMS_wetlab.models import Projects

def prepare_form_data_internal_sequencing (request_user):
	'''
	Description:
		The function get the information to display in the internal sequencing form
	Input:
		request_user      # user instance who request the service
	Return:
		form
	'''

	form = ServiceRequestFormInternalSequencing()
	# getting projects from user sharing list
	user_groups = request_user.groups.values_list('name',flat=True)
	if len (user_groups) > 0 :
		sharing_list = []
		for user in user_groups :
			if User.objects.filter(username__exact = user).exists():
				sharing_list.append(User.objects.get(username__exact = user).id)
		sharing_list.append(request_user.id)
		form.fields['serviceProjectNames'].queryset = Projects.objects.filter(user_id__in = sharing_list)
	else:
		form.fields['serviceProjectNames'].queryset = Projects.objects.filter(user_id__exact = request.user.id)
	form.fields['serviceAvailableService'].queryset = AvailableService.objects.filter(availServiceDescription__exact="Genomic data analysis").get_descendants(include_self=True)

	return form
