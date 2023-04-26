# Generic imports
import datetime
import django.core.mail
import drylab.models
import drylab.utils
import drylab.config


def prepare_delivery_form(resolution_id):
    """
    Description:
        The function get the services handled by the resolution to display them
        in the user form
    Input:
        resolution_id  
    Functions:
        get_resolution_obj_from_id
        get_pipelines_for_service
        get_available_service_obj_from_id
    Return:
        delivery_data_form
    """
    delivery_data_form = {}
    resolution_obj = (
        drylab.utils.resolutions.get_resolution_obj_from_id(
            resolution_id
        )
    )
    if resolution_obj is not None:
        delivery_data_form[
            "available_services"
        ] = resolution_obj.get_available_services_ids()
        delivery_data_form["resolution_id"] = resolution_id
        delivery_data_form["resolution_number"] = resolution_obj.get_resolution_number()
        # req_available_services_id = resolution_obj.get_available_services_ids()

        delivery_data_form[
            "pipelines_data"
        ] = drylab.utils.pipelines.get_pipelines_for_resolution(
            resolution_obj
        )

    return delivery_data_form


def store_resolution_delivery(form_data):
    """
    Description:
        The function get the information from the user from and create a new
        resolution delivery with the information
    Input:
        form_data
    Functions:
        get_resolution_obj_from_id
        get_pipelines_for_service
        get_pipeline_obj_from_id
    Return:
        delivery_data
    """
    delivery_data = None
    resolution_obj = (
        drylab.utils.resolutions.get_resolution_obj_from_id(
            form_data["resolution_id"]
        )
    )
    if resolution_obj is not None:
        delivery_data = {}
        if form_data["startdate"] != "":
            delivery_data["executionStartDate"] = datetime.datetime.strptime(
                form_data["startdate"], "%Y-%m-%d"
            )
        else:
            delivery_data["executionStartDate"] = None
        if form_data["startdate"] != "":
            delivery_data["execution_end_date"] = datetime.datetime.strptime(
                form_data["enddate"], "%Y-%m-%d"
            )
        else:
            delivery_data["execution_end_date"] = None
        delivery_data["delivery_resolutionID"] = resolution_obj
        delivery_data["execution_time"] = form_data["time"]
        delivery_data["permanent_used_space"] = form_data["pspace"]
        delivery_data["temporary_used_space"] = form_data["tspace"]
        delivery_data["deliveryNotes"] = form_data["deliveryNotes"]

        drylab.models.Delivery.objects.create_delivery(delivery_data)

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
    subject_tmp = drylab.config.SUBJECT_RESOLUTION_DELIVERED.copy()
    subject_tmp.insert(1, email_data["resolution_number"])
    subject = " ".join(subject_tmp)

    body_preparation = list(
        map(
            lambda st: str.replace(
                st, "RESOLUTION_NUMBER", email_data["resolution_number"]
            ),
            drylab.config.BODY_RESOLUTION_DELIVERED,
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
        drylab.models.ConfigSetting.objects.filter(
            configuration_name__exact="EMAIL_FOR_NOTIFICATIONS"
        )
        .last()
        .get_configuration_value()
    )
    from_user = notification_user
    to_users = [email_data["user_email"], email_data["user_email"], notification_user]
    try:
        django.core.mail.send_mail(subject, body_message, from_user, to_users)
    except Exception:
        pass
    return
