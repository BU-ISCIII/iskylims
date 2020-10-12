import datetime, time
import re, os
from django.conf import settings
from django.contrib.auth.models import Group, User
from iSkyLIMS_drylab import drylab_config
from smb.SMBConnection import SMBConnection
from iSkyLIMS_drylab.models import *
from django.template.loader import render_to_string
from django.core.mail import send_mail
from django.core.files.storage import FileSystemStorage



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

def create_resolution_pdf_file (service_obj,new_resolution, absolute_url):
    '''
    Description:
        The function collect the information to create the pdf file
    Input:
        request # contains the session information
    Functions:
        get_data_for_service_confirmation   # located at this file
        create_pdf                          # located at this file
    Constants:
        OUTPUT_DIR_RESOLUTION_PDF
    Return:
        pdf_file which contains the full path and name of the pdf file
    '''
    information_to_include = get_data_for_resolution(service_obj, new_resolution )

    pdf_file_name = new_resolution.get_resolution_number() + '.pdf'
    full_path_pdf_file = create_pdf(absolute_url, information_to_include, drylab_config.RESOLUTION_TEMPLATE, pdf_file_name, drylab_config.OUTPUT_DIR_RESOLUTION_PDF)
    pdf_file = full_path_pdf_file.replace(settings.BASE_DIR,'')
    return pdf_file

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


def get_assign_resolution_full_number(service_id, acronymName):
    '''
    Description:
        The function get the resolution full number if resolution already exists.
        Build the resolution  full number if it is the first resolution for the service
    Input:
        service_id # contains the service id
        acronymName # acronym name given to the service
    Functions:
        get_service_obj_from_id   # located at this file
    Return:
        resolution_full_number
    '''
    service_obj = get_service_obj_from_id(service_id)
    if Resolution.objects.filter(resolutionServiceID = service_id).exists():
        resolution_full_number = Resolution.objects.filter(resolutionServiceID = service_obj).last().get_resolution_number()
    else:
        resolution_full_number = ''
        resolution_full_number += service_obj.get_service_request_number() + '_'
        resolution_full_number += str(datetime.date.today()).replace('-','') + '_'
        resolution_full_number += acronymName + '_'
        resolution_full_number += service_obj.get_service_requested_user() + '_S'
    return resolution_full_number

def create_resolution_number(service_id):
    '''
    Description:
        The function create the resolution number and step it if more than 1 resolution
        have been created for the service.
    Input:
        service_id # contains the service id
    Functions:
        get_service_obj_from_id   # located at this file
    Return:
        resolution_number
    '''
    service_obj = get_service_obj_from_id(service_id)
    service_request_number = service_obj.get_service_request_number()
    if Resolution.objects.filter(resolutionServiceID = service_obj).exists():
        resolution_count =  Resolution.objects.filter(resolutionServiceID = service_obj).count()
        resolution_number = service_request_number + '.' + str(resolution_count +1 )
    else:
        resolution_number = service_request_number + '.1'
    return resolution_number


def get_data_for_resolution(service_obj, resolution_obj ):
    information, user, resolution_data = {}, {}, {}
    service_number ,run_specs, center, platform = service_obj.get_service_information().split(';')

    resolution_info = resolution_obj.get_resolution_information()
    # get profile object
    user_id = service_obj.serviceUserId.id

    information['resolution_number'] = resolution_obj.get_resolution_number()
    information['requested_date'] = service_obj.get_service_creation_time()
    information['resolution_date'] = resolution_info[4]
    information['nodes']= service_obj.serviceAvailableService.all()
    user['name'] = service_obj.serviceUserId.first_name
    user['surname'] = service_obj.serviceUserId.last_name

    user['area'] = Profile.objects.get(profileUserID = user_id).profileArea
    user['center'] = Profile.objects.get(profileUserID = user_id).profileCenter
    user['position'] = Profile.objects.get(profileUserID = user_id).profilePosition
    user['phone'] = Profile.objects.get(profileUserID = user_id).profileExtension
    user['email'] = service_obj.get_user_email()
    information['user'] = user
    resolution_info_split = resolution_info[1].split('_')
    resolution_data['acronym'] = resolution_info_split[2]
    resolution_data['estimated_date'] = resolution_info[3]
    resolution_data['notes'] = resolution_info[6]
    resolution_data['decission'] = service_obj.get_service_state()
    information['service_data'] = service_obj.get_service_user_notes()

    resolution_data['folder'] = resolution_info[1]
    information['resolution_data'] = resolution_data

    return information


def get_service_information (service_id):
    service_obj = Service.objects.get(pk=service_id)
    display_service_details = {}

    #text_for_dates = ['Service Date Creation', 'Approval Service Date', 'Rejected Service Date']
    service_dates = []
    display_service_details['service_name'] = service_obj.get_service_request_number()
    # get the list of projects
    #projects_in_service = {}
    if RequestedProjectInServices.objects.filter(projectService = service_obj).exists():
        projects_in_service = RequestedProjectInServices.objects.filter(projectService = service_obj)
        display_service_details['projects'] = []
        for project in projects_in_service:
            display_service_details['projects'].append([project.get_requested_external_project_id(), project.get_requested_project_name()])

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
    # adding actions fields
    if service_obj.serviceStatus != 'rejected' or service_obj.serviceStatus != 'archived':
        display_service_details['add_resolution_action'] = service_id
    if service_obj.serviceStatus == 'queued':
        resolution_id = Resolution.objects.filter(resolutionServiceID = service_obj).last().id
        display_service_details['add_in_progress_action'] = resolution_id
    if service_obj.serviceStatus == 'in_progress':
        resolution_id = Resolution.objects.filter(resolutionServiceID = service_obj).last().id
        display_service_details['add_delivery_action'] = resolution_id


    if Resolution.objects.filter(resolutionServiceID = service_obj).exists():
        resolution_list = Resolution.objects.filter(resolutionServiceID = service_obj)
        resolution_info =[]
        for resolution_item in resolution_list :
            resolution_info.append([resolution_item.get_resolution_information()])
        display_service_details['resolutions'] = resolution_info

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

def get_service_obj_from_id (service_request_id):
    '''
    Description:
        The function return the service instance from service id
    Input:
        service_request_id  # contains the id of the service request
    Return:
        service_request_obj
    '''
    if Service.objects.filter(pk__exact = service_request_id).exists():
        return Service.objects.get(pk__exact = service_request_id)
    return False

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


def get_latest_child_from_request_service(service):
    '''
    Description:
        The function collect the  child requested service
    Input:
        service
    Return:
        children_services
    '''
    return

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

def send_resolution_creation_email (email_data):
    '''
    Description:
        The function send the service email for resolution to user.
        Functions uses the send_email django core function to send the email
    Input:
        email_data      # Contains the information to include in the email
    Constant:
        SUBJECT_RESOLUTION_RECORDED
        BODY_RESOLUTION_RECORDED
        USER_EMAIL
    Return:
        None
    '''
    subject = drylab_config.SUBJECT_RESOLUTION_RECORDED.copy()
    subject.insert(1, email_data['service_number'])
    if email_data['status'] == 'Accepted':
        date = email_data['date'].strftime("%d %B, %Y")
        body_preparation = list(map(lambda st: str.replace(st, 'SERVICE_NUMBER', email_data['service_number']), drylab_config.BODY_RESOLUTION_ACCEPTED))
        body_preparation = list(map(lambda st: str.replace(st, 'USER_NAME', email_data['user_name']), body_preparation))
        body_preparation = list(map(lambda st: str.replace(st, 'STATUS', email_data['status']), body_preparation))
        body_preparation = list(map(lambda st: str.replace(st, 'DATE', date), body_preparation))
    else:
        body_preparation = list(map(lambda st: str.replace(st, 'SERVICE_NUMBER', email_data['service_number']), drylab_config.BODY_RESOLUTION_REJECTED))
        body_preparation = list(map(lambda st: str.replace(st, 'USER_NAME', email_data['user_name']), body_preparation))
        body_preparation = list(map(lambda st: str.replace(st, 'STATUS', email_data['status']), body_preparation))
    body_message = '\n'.join(body_preparation)

    from_user = drylab_config.USER_EMAIL
    to_users = [email_data['user_email'], drylab_config.USER_EMAIL]
    send_mail (subject, body_message, from_user, to_users)
    return

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

def store_resolution_additional_parameter(additional_parameters, resolution_obj):
    '''
    Description:
        The function store in database the additional resolution parameters.
    Input:
        additional_parameters       # Contains the list with the additional parameters
        resolution_obj              # resolution instance
    Return:
        None
    '''

    for additional_parameter in additional_parameters:

        parameter = {}
        parameter['resolution'] =resolution_obj
        for field in drylab_config.MAPPING_ADDITIONAL_RESOLUTION_PARAMETERS:
            parameter[field[0]] = additional_parameter[field[1]]
        new_parameter = ResolutionParameters.objects.create_resolution_parameters(parameter)
