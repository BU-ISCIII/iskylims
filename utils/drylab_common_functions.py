from django.conf import settings
from iSkyLIMS_drylab import drylab_config


def create_pdf(request,information, template_file, pdf_file_name):
    from weasyprint import HTML, CSS
    from django.template.loader import get_template

    from weasyprint.fonts import FontConfiguration


    #import pdb; pdb.set_trace()
    #font_config = FontConfiguration()
    html_string = render_to_string(template_file, {'information': information})
    pdf_file =  os.path.join (settings.BASE_DIR, drylab_config.OUTPUT_DIR_TEMPLATE , pdf_file_name)
    html = HTML(string=html_string, base_url=request.build_absolute_uri()).write_pdf(pdf_file,stylesheets=[CSS(settings.BASE_DIR + drylab_config.CSS_FOR_PDF)])

    return pdf_file


def create_service_id (service_number,user_name):
    user_center = Profile.objects.get(profileUserID = user_name).profileCenter.centerAbbr
    service_id = 'SRV' + user_center + service_number
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
        projects_in_service.append(project.get_project_name())
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

def get_service_information (service_id):
    service= Service.objects.get(pk=service_id)
    display_service_details = {}
    text_for_dates = ['Service Date Creation', 'Approval Service Date', 'Rejected Service Date']
    service_dates = []
    display_service_details['service_name'] = service.serviceRequestNumber
    # get the list of projects
    projects_in_service = {}
    projects_class = service.serviceProjectNames.all()
    for project in projects_class:
        project_id = project.id
        projects_in_service[project_id]=project.get_project_name()
    display_service_details['projects'] = projects_in_service
    display_service_details['user_name'] = service.serviceUserId.username
    display_service_details['file'] = os.path.join(settings.MEDIA_URL,str(service.serviceFile))
    display_service_details['state'] = service.serviceStatus
    display_service_details['service_notes'] = service.serviceNotes
    dates_for_services = service.get_service_dates()
    for i in range(len(dates_for_services)):
        service_dates.append([text_for_dates[i],dates_for_services[i]])
    display_service_details['service_dates'] = service_dates
    if service.serviceStatus != 'approved'and service.serviceStatus != 'recorded':
        # get the proposal for the delivery date
        #import pdb; pdb.set_trace()
        resolution_folder = Resolution.objects.filter(resolutionServiceID = service).last().resolutionFullNumber
        display_service_details['resolution_folder'] = resolution_folder
        resolution_estimated_date = Resolution.objects.filter(resolutionServiceID = service).last().resolutionEstimatedDate
        if resolution_estimated_date is None:
            resolution_estimated_date = "Not defined yet"
        display_service_details['estimated_delivery_date'] = resolution_estimated_date

    # get all services
    display_service_details['nodes']= service.serviceAvailableService.all()
    # adding actions fields
    if service.serviceStatus != 'rejected' or service.serviceStatus != 'archived':
        display_service_details['add_resolution_action'] = service_id
    if service.serviceStatus == 'queued':
        resolution_id = Resolution.objects.filter(resolutionServiceID = service).last().id
        display_service_details['add_in_progress_action'] = resolution_id
    if service.serviceStatus == 'in_progress':
        resolution_id = Resolution.objects.filter(resolutionServiceID = service).last().id
        display_service_details['add_delivery_action'] = resolution_id


    if Resolution.objects.filter(resolutionServiceID = service).exists():
        resolution_list = Resolution.objects.filter(resolutionServiceID = service)
        resolution_info =[]
        for resolution_item in resolution_list :
            resolution_info.append([resolution_item.get_resolution_information()])
        display_service_details['resolutions'] = resolution_info
    #import pdb; pdb.set_trace()
    if Resolution.objects.filter(resolutionServiceID = service).exists():
        resolution_list = Resolution.objects.filter(resolutionServiceID = service)
        delivery_info = []
        for resolution_id in resolution_list :
            if Delivery.objects.filter(deliveryResolutionID = resolution_id).exists():
                delivery = Delivery.objects.get(deliveryResolutionID = resolution_id)
                delivery_info.append([delivery.get_delivery_information()])
                display_service_details['delivery'] = delivery_info

    return display_service_details


def open_samba_connection():
    '''
    Description:
        The function open a samba connection with the parameter settings
        defined in frylab configuration file
    Return:
        conn object for the samba connection
    '''
    try:

        conn=SMBConnection(drylab_config.SAMBA_USER_ID, drylab_config.SAMBA_USER_PASSWORD, drylab_config.SAMBA_SHARED_FOLDER_NAME,
                            drylab_config.SAMBA_REMOTE_SERVER_NAME, use_ntlm_v2=drylab_config.SAMBA_NTLM_USED, domain = drylab_config.SAMBA_DOMAIN)
        conn.connect(drylab_config.SAMBA_IP_SERVER, int(drylab_config.SAMBA_PORT_SERVER))
    except:
        return False

    return conn



def increment_service_number ( user_name):
    # check user center
    #import pdb; pdb.set_trace()
    user_center = Profile.objects.get(profileUserID = user_name).profileCenter.centerAbbr
    # get latest service used for user's center
    if Service.objects.filter(serviceUserId__profile__profileCenter__centerAbbr=user_center).exists():
        number_request = Service.objects.filter(serviceUserId__profile__profileCenter__centerAbbr=user_center).last().serviceRequestInt
        service_number = str(int(number_request) + 1).zfill(3)
    else:
        service_number = '001'
    return service_number


