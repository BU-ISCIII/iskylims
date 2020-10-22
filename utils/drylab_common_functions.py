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



def get_data_for_service_confirmation (service_requested):
    information = {}
    user = {}
    service_data ={}
    service = Service.objects.get(serviceRequestNumber = service_requested)
    service_number ,run_specs, center, platform = service.get_service_information().split(';')
    information['service_number'] = service_number
    information['requested_date'] = service.get_service_creation_time()
    information['nodes']= service.serviceAvailableService.all()
    user['name'] = service.serviceUserId.first_name
    user['surname'] = service.serviceUserId.last_name

    user_id = service.serviceUserId.id
    user['area'] = Profile.objects.get(profileUserID = user_id).profileArea
    user['center'] = Profile.objects.get(profileUserID = user_id).profileCenter
    user['phone'] = Profile.objects.get(profileUserID = user_id).profileExtension
    user['position'] = Profile.objects.get(profileUserID = user_id).profilePosition
    user['email'] = service.serviceUserId.email
    information['user'] = user
    projects_in_service = []
    projects_class = service.serviceProjectNames.all()
    for project in projects_class:
        projects_in_service.append(project.get_requested_project_name())
    service_data['projects'] = projects_in_service
    service_data['platform'] = platform
    service_data['run_specifications'] = run_specs
    service_data['center'] = center
    service_data['notes'] = service.serviceNotes
    if str(service.serviceFile) != '':
        full_path_file = str(service.serviceFile).split('/')
        stored_file = full_path_file[-1]
        temp_string_file = re.search('^\d+_(.*)', stored_file)
        service_data['file'] = temp_string_file.group(1)
    else:
        service_data['file'] = 'Not provided'
    information['service_data'] = service_data

    return information



def create_service_pdf_file (service_request_number, absolute_url):
    '''
    Description:
        The function collect the information to create the pdf file
    Input:
        request # contains the session information
    Functions:
        get_data_for_service_confirmation   # located at this file
        create_pdf                          # located at this file
    Constants:
        OUTPUT_DIR_SERVICE_REQUEST_PDF
    Return:
        pdf_file which contains the full path and name of the pdf file
    '''

    information_to_include = get_data_for_service_confirmation(service_request_number)
    pdf_file_name = service_request_number + '.pdf'
    full_path_pdf_file = create_pdf(absolute_url, information_to_include, drylab_config.REQUESTED_CONFIRMATION_SERVICE, pdf_file_name,  drylab_config.OUTPUT_DIR_SERVICE_REQUEST_PDF)
    pdf_file = full_path_pdf_file.replace(settings.BASE_DIR,'')
    return pdf_file



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

def send_service_creation_confirmation_email(email_data):
    '''
    Description:
        The function send the service email confirmation to user.
        Functions uses the send_email django core function to send the email
    Input:
        email_data      # Contains the information to include in the email
    Constant:
        SUBJECT_SERVICE_RECORDED
        BODY_SERVICE_RECORDED
        USER_EMAIL
    Return:
        None
    '''
    subject = drylab_config.SUBJECT_SERVICE_RECORDED.copy()
    subject.insert(1, email_data['service_number'])

    body_preparation = list(map(lambda st: str.replace(st, 'SERVICE_NUMBER', email_data['service_number']), drylab_config.BODY_SERVICE_RECORDED))
    body_preparation = list(map(lambda st: str.replace(st, 'USER_NAME', email_data['user_name']), body_preparation))
    body_message = '\n'.join(body_preparation)

    from_user = drylab_config.USER_EMAIL
    to_users = [email_data['user_email'], drylab_config.USER_EMAIL]
    send_mail (subject, body_message, from_user, to_users)
    return


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
