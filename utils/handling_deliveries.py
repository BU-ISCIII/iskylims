from datetime import datetime
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
        new_delivery
    '''
    new_delivery = None
    resolution_obj = get_resolution_obj_from_id(form_data['resolution_id'])
    if resolution_obj != None :
        delivery_data = {}
        if form_data['startdate'] != '':
            delivery_data['executionStartDate'] = datetime.strptime(form_data['startdate'], '%Y-%m-%d')
        if form_data['startdate'] != '':
            delivery_data['executionEndDate'] = datetime.strptime(form_data['enddate'], '%Y-%m-%d')
        delivery_data['deliveryResolutionID'] = resolution_obj
        delivery_data['executionTime'] = data_form['time']
        delivery_data['permanentUsedSpace'] = data_form['pspace']
        delivery_data['temporaryUsedSpace'] = data_form['tspace']
        delivery_data['deliveryNotes'] = data_form['deliveryNotes']
        new_delivery = Delivery.objects.creat_delivery(delivery_data)
        if 'pipelines' in form_data:
            piplines_ids = form_data.getlist('pipelines')
            for pipeline_id in pipelines_ids:
                new_delivery.pipelinesInDelivery.add(get_pipeline_obj_from_id(pipeline_id))
    
    return new_delivery
