import time
from datetime import datetime
import re, os
from django.conf import settings
from django.contrib.auth.models import Group, User
from iSkyLIMS_drylab import drylab_config
#from smb.SMBConnection import SMBConnection
from iSkyLIMS_drylab.models import *
from django.template.loader import render_to_string
from django.core.mail import send_mail
from django.core.files.storage import FileSystemStorage
from django.contrib.auth.models import User

#from iSkyLIMS_drylab.utils.handling_request_services import *
#from iSkyLIMS_drylab.utils.handling_resolutions import *

def check_valid_date_format (date):
    try:
        datetime.strptime(date, '%Y-%m-%d')
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

def get_email_data():
    '''
    Description:
        Fetch the email configuration file
    Return:
        email_data
    '''
    email_data = {}
    if EmailData.objects.all().exists():
        email_data_obj = EmailData.objects.last()
        email_data = email_data_obj.get_email_data()
    return email_data

def save_email_data(email_fields):
    '''
    Description:
        create the email configuration file . If exists the old information is deleted
    Input:
        email_fields    # Email fields settings
    Return:
        email_data_obj
    '''
    if EmailData.objects.all().exists():
        email_data_obj = EmailData.objects.last().update_data(email_fields)
    else:
        email_data_obj = EmailData.objects.create_email_data(email_fields)
    return email_data_obj

def increment_service_number ( request_user):
    '''
    Description:
        The function will check if the logged user belongs to service
        manager group
    Input:
        request_user # request user obj
    Return:
        service_number
    '''
    try:
        user_center = Profile.objects.get(profileUserID = request_user).profileCenter.centerAbbr
    except:
        user_center = drylab_config.USER_CENTER_USED_WHEN_NOT_PROVIDED
    # get latest service used for user's center
    import pdb; pdb.set_trace()
    if Service.objects.filter(serviceUserId__profile__profileCenter__centerAbbr__exact = user_center).exists():
        number_request = Service.objects.filter(serviceUserId__profile__profileCenter__centerAbbr__exact = user_center).last().get_service_request_integer()
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

def get_user_sharing_list(request_user):
    '''
    Description:
        The function get the primary key of the that are sharing their information
        If the request user is a service manager, the function return all user ids
    Input:
        request_user      # user obj
    Constant:
        SERVICE_MANAGER
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


def get_defined_username_and_ids ():
    '''
    Description:
        The function get the userid for all users defined in iSkyLIMS
    Return:
        userids_list
    '''
    userids_list = []
    if User.objects.all().exists():
        user_objs = User.objects.all().order_by('username')
        for user_obj in user_objs:
            userids_list.append([user_obj.username, user_obj.pk])
    return userids_list


def get_users_requested_services():
    '''
    Description:
        The function get the list of users that have requested any service.
    Return:
        user_list
    '''
    user_list = []
    if Service.objects.all().exists():
        user_ids = Service.objects.all().order_by('serviceUserId').values('serviceUserId').distinct()
        for user_id in user_ids:
            user_list.append([user_id['serviceUserId'], User.objects.filter(pk__exact = user_id['serviceUserId']).last().username])
    return user_list

def get_current_users():
    '''
    Description:
        The function returns the user that their session is active
    Return:
        User list
    '''
    from django.contrib.sessions.models import Session
    from django.utils import timezone
    active_sessions = Session.objects.filter(expire_date__gte=timezone.now())
    user_id_list = []
    for session in active_sessions:
        data = session.get_decoded()
        user_id_list.append(data.get('_auth_user_id', None))
    # Query all logged in users based on id list
    return User.objects.filter(id__in=user_id_list)


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
