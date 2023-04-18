from datetime import datetime
from django.core.mail import send_mail
from iSkyLIMS_drylab.models import *
from iSkyLIMS_drylab.utils.handling_resolutions import get_resolution_obj_from_id
from iSkyLIMS_drylab.utils.handling_pipelines import get_pipelines_for_resolution


def prepare_delivery_form(resolution_id):
    """
    Description:
        The function get the services handled by the resolution to display them
        in the user form
    Input:
        resolution_id   # resolution id
    Functions:
        get_resolution_obj_from_id # located at iSkyLIMS_drylab.utils.handling_resolutions
        get_pipelines_for_service     # located at iSkyLIMS_drylab.utils.handling_pipelines
        get_available_service_obj_from_id       # located at iSkyLIMS_drylab.utils.handling_request_services
    Return:
        delivery_data_form
    """
    delivery_data_form = {}
    resolution_obj = get_resolution_obj_from_id(resolution_id)
    if resolution_obj != None:
        delivery_data_form[
            "available_services"
        ] = resolution_obj.get_available_services_ids()
        delivery_data_form["resolution_id"] = resolution_id
        delivery_data_form["resolution_number"] = resolution_obj.get_resolution_number()
        # req_available_services_id = resolution_obj.get_available_services_ids()

        delivery_data_form["pipelines_data"] = get_pipelines_for_resolution(
            resolution_obj
        )

    return delivery_data_form


def store_resolution_delivery(form_data):
    """
    Description:
        The function get the information from the user from and create a new
        resolution delivery with the information
    Input:
        form_data   # user data form
    Functions:
        get_resolution_obj_from_id # located at iSkyLIMS_drylab.utils.handling_resolutions
        get_pipelines_for_service     # located at iSkyLIMS_drylab.utils.handling_pipelines
        get_pipeline_obj_from_id       # located at iSkyLIMS_drylab.utils.handling_pipelines
    Return:
        delivery_data
    """
    delivery_data = None
    resolution_obj = get_resolution_obj_from_id(form_data["resolution_id"])
    if resolution_obj != None:
        delivery_data = {}
        if form_data["startdate"] != "":
            delivery_data["executionStartDate"] = datetime.strptime(
                form_data["startdate"], "%Y-%m-%d"
            )
        else:
            delivery_data["executionStartDate"] = None
        if form_data["startdate"] != "":
            delivery_data["executionEndDate"] = datetime.strptime(
                form_data["enddate"], "%Y-%m-%d"
            )
        else:
            delivery_data["executionEndDate"] = None
        delivery_data["deliveryResolutionID"] = resolution_obj
        delivery_data["executionTime"] = form_data["time"]
        delivery_data["permanentUsedSpace"] = form_data["pspace"]
        delivery_data["temporaryUsedSpace"] = form_data["tspace"]
        delivery_data["deliveryNotes"] = form_data["deliveryNotes"]

        Delivery.objects.create_delivery(delivery_data)

        resolution_obj.update_resolution_in_delivered()

        delivery_data["resolution_number"] = resolution_obj.get_resolution_number()
    return delivery_data


def send_delivery_service_email(email_data):
    """
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
    """
    subject_tmp = drylab_config.SUBJECT_RESOLUTION_DELIVERED.copy()
    subject_tmp.insert(1, email_data["resolution_number"])
    subject = " ".join(subject_tmp)

    body_preparation = list(
        map(
            lambda st: str.replace(
                st, "RESOLUTION_NUMBER", email_data["resolution_number"]
            ),
            drylab_config.BODY_RESOLUTION_DELIVERED,
        )
    )
    body_preparation = list(
        map(
            lambda st: str.replace(st, "USER_NAME", email_data["user_name"]),
            body_preparation,
        )
    )

    body_message = "\n".join(body_preparation)
    notification_user = (
        ConfigSetting.objects.filter(configurationName__exact="EMAIL_FOR_NOTIFICATIONS")
        .last()
        .get_configuration_value()
    )
    from_user = notification_user
    to_users = [email_data["user_email"], email_data["user_email"], notification_user]
    try:
        send_mail(subject, body_message, from_user, to_users)
    except:
        pass
    return
