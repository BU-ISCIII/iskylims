import os
from django.conf import settings

from iSkyLIMS_drylab import drylab_config
import iSkyLIMS_drylab.models

from django.contrib.auth.models import User

def get_config_file (config_file):
    c_file = []
    try:
        with open (config_file ,'r') as fh:
            for line in fh:
                if 'PASSWORD' in line:
                    hide_passwd = line.split('=')
                    hide_passwd[1] = 'XXXXXXXXXXXXXXXXX'
                    line = ' = '.join(hide_passwd)
                line = line.replace('\n', '')
                c_file.append(line)
    except:
        return
    return c_file


def get_iSkyLIMS_settings():
    s_file = []
    settings_file = os.path.join(settings.BASE_DIR, 'iSkyLIMS' , 'settings.py')
    try:
        with open (settings_file ,'r') as fh:
            for line in fh:
                if 'PASSWORD' in line or 'SECRET_KEY' in line   :
                    if  not 'AUTH_PASSWORD_VALIDATORS' in line :
                        if '=' in line :
                            split_separator = '='
                        else :
                            split_separator = ':'
                        hide_passwd = line.split(split_separator)
                        hide_passwd[1] = 'XXXXXXXXXXXXXXXXX'
                        line = ' = '.join(hide_passwd)

                line = line.replace('\n', '')
                s_file.append(line)
    except:
        return

    return s_file

def create_service_test(service_requested):
    service_results = []
    # Check if test service exists
    if iSkyLIMS_drylab.models.Service.objects.filter(service_request_number__exact = service_requested).exists():
        delete_service =  iSkyLIMS_drylab.models.Service.objects.get(service_request_number__exact = service_requested)
        delete_service.delete()

    # Check user is defined in database
    if not User.objects.filter(username__exact = 'test_userDrylab').exists():
        user = User.objects.create_user(username='test_userDrylab',
                                 email='test_userDrylab@iSkyLIMS.com',
                                 password='test_userD')


    try:
        user_name = User.objects.get(username__exact = 'test_userDrylab')
        service_results.append(('User defined', 'OK'))
    except:
        service_results.append(('User defined', 'NOK'))
    try:
        service_platform = iSkyLIMS_drylab.models.Platform.objects.first()
        service_results.append(('Platform defined', 'OK'))
    except:
        service_results.append(('Platform defined', 'NOK'))
    try:
        service_file_ext = iSkyLIMS_drylab.models.FileExt.objects.first()
        service_results.append(('File extension defined', 'OK'))
    except:
        service_results.append(('File extension defined', 'NOK'))

    for i in range (len(service_results)):
        if 'NOK' in service_results[i]:
            return service_results, 'NOK'
    try:
        new_test_service = iSkyLIMS_drylab.models.Service(service_request_number = service_requested, serviceUserId = user_name,
                    servicePlatform = service_platform, serviceFileExt = service_file_ext,
                    serviceStatus = 'recorded')
        new_test_service.save()

        service_results.append(('Service Test creation', 'OK'))
        return service_results, 'OK'
    except:
        service_results.append(('Service Test creation', 'NOK'))
        return service_results, 'NOK'


def create_resolution_test (resolution_number, service_requested):
    from weasyprint import HTML, CSS

    from django.core.files.storage import FileSystemStorage
    from weasyprint.fonts import FontConfiguration

    resolution_test = []
    # get service object
    service = iSkyLIMS_drylab.models.Service.objects.get(service_request_number = service_requested)

    # Create resolution object
    try:
        test_resolution = iSkyLIMS_drylab.models.Resolution(resolution_serviceID = service, resolutionNumber = resolution_number, resolution_full_number = str('Test_'+ resolution_number))
        resolution_test.append(('Resolution creation', 'OK'))
    except:
        resolution_test.append(('Resolution creation', 'NOK'))
    #test_resolution.save()
    #resolution = Resolution.objects.get(resolutionNumber = resolution_number)
    resolution_info = test_resolution.get_resolution_information()

    # get profile object
    #user_id = service.serviceUserId.id

    resolution_test.append(('Folder structure creation', 'OK'))

    '''
    information['resolution_number'] = resolution_number
    information['requested_date'] = service.get_service_creation_time()
    information['resolution_date'] = resolution_info[4]
    information['nodes']= service.service_available_service.all()
    user['name'] = service.serviceUserId.first_name
    user['surname'] = service.serviceUserId.last_name


    user_id = service.serviceUserId.id
    user['area'] = Profile.objects.get(profileUserID = user_id).profileArea
    user['center'] = Profile.objects.get(profileUserID = user_id).profileCenter
    user['position'] = Profile.objects.get(profileUserID = user_id).profilePosition
    user['phone'] = Profile.objects.get(profileUserID = user_id).profileExtension
    user['email'] = service.serviceUserId.email
    information['user'] = user
    resolution_data['folder'] = resolution_info[1]
    resolution_data['estimated_date'] = resolution_info[3]
    resolution_data['notes'] = resolution_info[6]
    resolution_data['decission'] = service.serviceStatus
    information['service_data'] = service.service_notes

    resolution_data['folder'] = resolution_info[1]
    information['resolution_data'] = resolution_data
    html_string = render_to_string('resolution_template.html', {'information': information})

    html = HTML(string=html_string, base_url=request.build_absolute_uri()).write_pdf('documents/drylab/res_pdf.pdf',stylesheets=[CSS(settings.STATIC_ROOT +
                                drylab_config.CSS_FOR_PDF)])

    fs = FileSystemStorage('documents/drylab')
    with fs.open('res_pdf.pdf') as pdf:
        response = HttpResponse(pdf, content_type='application/pdf')
        # save pdf file as attachment
        #response['Content-Disposition'] = 'attachment; filename="mypdf.pdf"'


        response['Content-Disposition'] = 'inline;filename=res_pdf.pdf'
    '''
    return resolution_test

def delete_test_service(service_name):
    try:
        d_service = iSkyLIMS_drylab.models.Service.objects.get(service_request_number__exact = service_name)
        d_service.delete()
    except:
        return

    return
