from datetime import datetime
from django.core.mail import send_mail
from iSkyLIMS_drylab.models import *
from iSkyLIMS_drylab.utils.handling_request_services import get_available_service_obj_from_id
from iSkyLIMS_drylab.utils.handling_resolutions import get_resolution_obj_from_id
from iSkyLIMS_drylab.utils.handling_pipelines import get_pipeline_and_versions_for_available_service, get_pipeline_obj_from_id

def prepare_delivery_form(resolution_id):
    '''
    Description:
        The function get the services handled by the resolution to display them
        in the user form
    Input:
        resolution_id   # resolution id
    Functions:
        get_resolution_obj_from_id # located at iSkyLIMS_drylab.utils.handling_resolutions
        get_pipeline_and_versions_for_available_service     # located at iSkyLIMS_drylab.utils.handling_pipelines
        get_available_service_obj_from_id       # located at iSkyLIMS_drylab.utils.handling_request_services
    Return:
        delivery_data_form
    '''
    delivery_data_form ={}
    pipelines_data = []
    resolution_obj = get_resolution_obj_from_id(resolution_id)
    if resolution_obj != None :
        delivery_data_form['available_services'] = resolution_obj.get_available_services_ids()
        delivery_data_form['resolution_id'] = resolution_id
        delivery_data_form['resolution_number'] = resolution_obj.get_resolution_number()
        req_available_services_id = resolution_obj.get_available_services_ids()
        for avail_service in req_available_services_id:
            data = get_pipeline_and_versions_for_available_service(avail_service)
            if data :
                pipelines_data.append([ get_available_service_obj_from_id(avail_service).get_service_description() , data])
        delivery_data_form['pipelines_data'] = pipelines_data
    return delivery_data_form


def store_resolution_delivery(form_data):
    '''
    Description:
        The function get the information from the user from and create a new
        resolution delivery with the information
    Input:
        form_data   # user data form
    Functions:
        get_resolution_obj_from_id # located at iSkyLIMS_drylab.utils.handling_resolutions
        get_pipeline_and_versions_for_available_service     # located at iSkyLIMS_drylab.utils.handling_pipelines
        get_pipeline_obj_from_id       # located at iSkyLIMS_drylab.utils.handling_pipelines
    Return:
        delivery_data
    '''
    delivery_data = None
    resolution_obj = get_resolution_obj_from_id(form_data['resolution_id'])
    if resolution_obj != None :
        delivery_data = {}
        if form_data['startdate'] != '':
            delivery_data['executionStartDate'] = datetime.datetime.strptime(form_data['startdate'], '%Y-%m-%d')
        else:
            delivery_data['executionStartDate'] = None
        if form_data['startdate'] != '':
            delivery_data['executionEndDate'] = datetime.datetime.strptime(form_data['enddate'], '%Y-%m-%d')
        else:
            delivery_data['executionEndDate'] = None
        delivery_data['deliveryResolutionID'] = resolution_obj
        delivery_data['executionTime'] = form_data['time']
        delivery_data['permanentUsedSpace'] = form_data['pspace']
        delivery_data['temporaryUsedSpace'] = form_data['tspace']
        delivery_data['deliveryNotes'] = form_data['deliveryNotes']

        new_delivery = Delivery.objects.create_delivery(delivery_data)
        if 'pipelines' in form_data:
            pipelines_ids = form_data.getlist('pipelines')
            for pipeline_id in pipelines_ids:
                new_delivery.pipelinesInDelivery.add(get_pipeline_obj_from_id(pipeline_id))
        resolution_obj.update_resolution_in_delivered()
        service_obj = resolution_obj.get_service_obj()
        import pdb; pdb.set_trace()
        if service_obj.get_service_state() == 'In progress':
            service_obj.update_service_status('Delivered')
        '''
        if Resolution.objects.filter(resolutionServiceID = service_obj).exclude(resolutionState__resolutionStateName__exact = 'Delivery').exists():
            if Resolution.objects.filter(resolutionServiceID = service_obj, resolutionState__resolutionStateName__exact = 'In progress').exists():
                pending_resolution_obj = resolutionResolution.objects.filter(resolutionServiceID = service_obj, resolutionState__resolutionStateName__exact = 'In Progress').last()
                delivery_data['pending_resolution'] = pending_resolution_obj.get_resolution_number()
                service_obj.update_service_status('In progress')
            elif Resolution.objects.filter(resolutionServiceID = service_obj, resolutionState__resolutionStateName__exact = 'Recorded').exists():
                pending_resolution_obj = resolutionResolution.objects.filter(resolutionServiceID = service_obj, resolutionState__resolutionStateName__exact = 'Recorded').last()
                service_obj.update_service_status('Queued')

                delivery_data['pending_resolution'] = pending_resolution_obj.get_resolution_number()
            else:
                service_obj.update_service_status('Delivered')
        else:
            import pdb; pdb.set_trace()
            if service_obj.get_service_state() != 'Recorded' :
                service_obj.update_service_status('Delivered')
        '''
        
        delivery_data['resolution_number'] = resolution_obj.get_resolution_number()
    return delivery_data

def send_delivery_service_email (email_data):
    '''
    Description:
        The function send the email for delivery service to user.
        Functions uses the send_email django core function to send the email
    Input:
        email_data      # Contains the information to include in the email
    Constant:
        SUBJECT_RESOLUTION_DELIVERED
        BODY_RESOLUTION_DELIVERED
        USER_EMAIL
    Return:
        None
    '''
    subject_tmp = drylab_config.SUBJECT_RESOLUTION_DELIVERED.copy()
    subject_tmp.insert(1, email_data['resolution_number'])
    subject = ' '.join(subject_tmp)

    body_preparation = list(map(lambda st: str.replace(st, 'RESOLUTION_NUMBER', email_data['resolution_number']), drylab_config.BODY_RESOLUTION_DELIVERED))
    body_preparation = list(map(lambda st: str.replace(st, 'USER_NAME', email_data['user_name']), body_preparation))
    from_user = drylab_config.USER_EMAIL
    to_users = [email_data['user_email'], drylab_config.USER_EMAIL]
    body_message = '\n'.join(body_preparation)

    send_mail (subject, body_message, from_user, to_users)
    return
