import datetime, time
import re, os
from django.conf import settings
from django.contrib.auth.models import Group, User
from iSkyLIMS_drylab import drylab_config
#from smb.SMBConnection import SMBConnection
from iSkyLIMS_drylab.models import *
from django.template.loader import render_to_string
from django.core.mail import send_mail
from django.core.files.storage import FileSystemStorage


#from iSkyLIMS_drylab.utils.handling_request_services import *
#from iSkyLIMS_drylab.utils.handling_resolutions import *

def check_valid_date_format (date):
    try:
        datetime.datetime.strptime(date, '%Y-%m-%d')
        return True
    except:
        return False

def create_pdf(absolute_url,information, template_file, pdf_file_name , out_dir):
    from weasyprint import HTML, CSS
    from django.template.loader import get_template
    from django.template.loader import render_to_string
    from weasyprint.fonts import FontConfiguration

    #font_config = FontConfiguration()
    html_string = render_to_string(template_file, {'information': information})
    pdf_dir =  os.path.join (settings.BASE_DIR, out_dir)
    if not os.path.exists(pdf_dir):
        os.makedirs(pdf_dir)
    pdf_file = os.path.join(pdf_dir, pdf_file_name)
    html = HTML(string=html_string, base_url=absolute_url).write_pdf(pdf_file,stylesheets=[CSS(settings.BASE_DIR + drylab_config.CSS_FOR_PDF)])

    return pdf_file


def create_service_id (service_number,user_name):
    '''
	Description:
		The function get the user center to build the service ID string
	Input:
		service_number      # number of the service
        user_name           # user name to get the center
	Constants:
		USER_CENTER_USED_WHEN_NOT_PROVIDED
        ABBREVIATION_USED_FOR_SERVICE_REQUEST
	Return:
		service_id
	'''
    try:
        user_center = Profile.objects.get(profileUserID = user_name).profileCenter.centerAbbr
    except:
        user_center = drylab_config.USER_CENTER_USED_WHEN_NOT_PROVIDED
    service_id = drylab_config.ABBREVIATION_USED_FOR_SERVICE_REQUEST + user_center + service_number
    return service_id



def is_service_manager (request):
    '''
    Description:
        The function will check if the logged user belongs to service
        manager group
    Input:
        request # contains the session information
    Variables:
        groups # drylab manager object group
    Return:
        Return True if the user belongs to service Manager, False if not
    '''
    try:
        groups = Group.objects.get(name = drylab_config.SERVICE_MANAGER)
        if groups not in request.user.groups.all():
            return False
    except:
        return False

    return True


def increment_service_number ( user_name):
    '''
    Description:
        The function will check if the logged user belongs to service
        manager group
    Input:
        user_name # contains the session information
    Return:
        service_number
    '''
    # check user center
    try:
        user_center = Profile.objects.get(profileUserID = user_name).profileCenter.centerAbbr
    except:
        user_center = drylab_config.USER_CENTER_USED_WHEN_NOT_PROVIDED
    # get latest service used for user's center
    if Service.objects.filter(serviceUserId__profile__profileCenter__centerAbbr=user_center).exists():
        number_request = Service.objects.filter(serviceUserId__profile__profileCenter__centerAbbr=user_center).last().serviceRequestInt
        if number_request == None:
            service_number = '001'
        else:
            service_number = str(int(number_request) + 1).zfill(3)
    else:
        service_number = '001'
    return service_number



def get_children_available_services():
    '''
    Description:
        The function collect the available services and the fields to present information
        in the form
    Return:
        children_services
    '''
    children_services = []
    if AvailableService.objects.filter(parent = None).exists():
        all_services = AvailableService.objects.filter(parent = None).get_descendants()
        for service in all_services :
            if not service.get_children().exists():
                children_services.append([service.pk, service.get_service_description()])
    return children_services

def get_user_sharing_lits(request_user):
    '''
    Description:
        The function get the primary key of the that are sharing their information
        If the request user is a service manager, the function return all user ids
    Input:
        request_user      # user obj
    Constant:

    Return:
        sharing_list
    '''
    # getting projects from user sharing list
    sharing_list = []
    user_groups = request_user.groups.values_list('name',flat=True)
    if drylab_config.SERVICE_MANAGER in user_groups:
        all_users = User.objects.all()
        for user in all_users:
            sharing_list.append(user.id)
    else:
        for user in user_groups :
            if User.objects.filter(username__exact = user).exists():
                sharing_list.append(User.objects.get(username__exact = user).id)
        sharing_list.append(request_user.id)
    return sharing_list


def store_file_from_form(file , path):
    '''
    Description:
        The function store the user input file and return the file_name.
    Input:
        file       # file to store
        path       # path to store the file
    Return:
        f_name      # contains the file name and the extension
        stored_file # contains the full path
    '''
    if '.' in file.name:
        split_filename=re.search('(.*)(\.\w+$)',file.name)
        file_name = split_filename.group(1)
        file_ext = split_filename.group(2)
    else:
        file_name = file.name
        file_ext = ''
    fs = FileSystemStorage()
    timestr = time.strftime("%Y%m%d%H%M%S")
    f_name = timestr + '_' + file_name  + file_ext
    full_path_file = os.path.join(path,f_name)

    fs.save(full_path_file, file )
    return f_name,  full_path_file
