# Generic imports
import datetime

import django.contrib.auth.models
import django.core.files.storage
from django.contrib.sessions.models import Session
from django.utils import timezone

# Local imports
import django_utils.models
import drylab.config
import drylab.models


def get_service_obj(service_id, input="pk"):
    """
    Description:
        The function get the  service obj  from the id
    Input:
        service_id  # id of the  service
        inpuf       # either pk or id. p.e SRVCNM123
    Return:
        service_obj
    """
    service_obj = None
    if input == "pk":
        if drylab.models.Service.objects.filter(pk__exact=service_id).exists():
            service_obj = drylab.models.Service.objects.filter(
                pk__exact=service_id
            ).last()
    elif input == "id":
        if drylab.models.Service.objects.filter(
            service_request_number__exact=service_id
        ).exists():
            service_obj = drylab.models.Service.objects.filter(
                service_request_number__exact=service_id
            ).last()

    return service_obj


def check_valid_date_format(date):
    try:
        datetime.datetime.strptime(date, "%Y-%m-%d")
        return True
    except Exception:
        return False


def create_service_id(service_number, user_name):
    """
    Description:
            The function get the user center to build the service ID string
    Input:
            service_number      # number of the service
    user_name           # user name to get the center
    Constants:
            USER_CENTER_USED_WHEN_NOT_PROVIDED
    Return:
            service_id
    """
    if get_configuration_from_database("USER_CENTER_USED_FOR_NAMING_SERVICE") == "True":
        try:
            user_center = django_utils.models.Profile.objects.get(
                profile_user_id=user_name
            ).profile_center.center_abbr
        except ValueError:
            user_center = drylab.config.USER_CENTER_USED_WHEN_NOT_PROVIDED
    else:
        user_center = ""
    abbr = get_configuration_from_database("ABBREVIATION_USED_FOR_SERVICE_REQUEST")

    service_id = abbr + user_center + service_number
    return service_id


def get_configuration_from_database(configuration_name):
    """
    Description:
        The function fetch from database the configuration setting value
    Input:
        configuration_name      # configuration settings name
    """
    configuration_value = ""
    if drylab.models.ConfigSetting.objects.filter(
        configuration_name__exact=configuration_name
    ).exists():
        configuration_settings_obj = drylab.models.ConfigSetting.objects.filter(
            configuration_name__exact=configuration_name
        ).last()
        configuration_value = configuration_settings_obj.get_configuration_value()
    return configuration_value


def is_service_manager(request):
    """The function will check if the logged user belongs to service
        manager group

    Parameters
    ----------
    request : user instance
        instance of the logged user

    Returns
    -------
    Boolean
        True if logged user is service manager. False if not
    """
    if request.user.groups.filter(name=drylab.config.SERVICE_MANAGER).exists():
        return True
    return False


def increment_service_number(request_user):
    """
    Description:
        The function will check if the logged user belongs to service
        manager group
    Input:
        request_user # request user obj
    Return:
        service_number
    """
    try:
        user_center = django_utils.models.Profile.objects.get(
            profile_user_id=request_user
        ).get_user_center_abbr()
    except Exception:
        user_center = drylab.config.USER_CENTER_USED_WHEN_NOT_PROVIDED
    # get latest service used for user's center

    if drylab.models.Service.objects.filter(service_center__exact=user_center).exists():
        number_request = (
            drylab.models.Service.objects.filter(service_center__exact=user_center)
            .last()
            .get_number()
        )
        if number_request is None:
            service_number = "001"
        else:
            service_number = str(int(number_request) + 1).zfill(3)
    else:
        service_number = "001"
    return service_number


def get_user_sharing_list(request_user):
    """
    Description:
        The function get the primary key of the that are sharing their information
        If the request user is a service manager, the function return all user ids
    Input:
        request_user      # user obj
    Constant:
        SERVICE_MANAGER
    Return:
        sharing_list
    """
    # getting projects from user sharing list
    sharing_list = []
    user_groups = request_user.groups.values_list("name", flat=True)
    if drylab.config.SERVICE_MANAGER in user_groups:
        all_users = django.contrib.auth.models.User.objects.all()
        for user in all_users:
            sharing_list.append(user.id)
    else:
        for user in user_groups:
            if django.contrib.auth.models.User.objects.filter(
                username__exact=user
            ).exists():
                sharing_list.append(
                    django.contrib.auth.models.User.objects.get(username__exact=user).id
                )
        sharing_list.append(request_user.id)
    return sharing_list


def get_defined_username_and_ids():
    """
    Description:
        The function get the userid for all users defined in iSkyLIMS
    Return:
        userids_list
    """
    userids_list = []
    if django.contrib.auth.models.User.objects.all().exists():
        user_objs = django.contrib.auth.models.User.objects.all().order_by("username")
        for user_obj in user_objs:
            userids_list.append([user_obj.username, user_obj.pk])
    return userids_list


def get_users_requested_services():
    """
    Description:
        The function get the list of users that have requested any service.
    Return:
        user_list
    """
    user_list = []
    if drylab.models.Service.objects.all().exists():
        user_ids = (
            drylab.models.Service.objects.all()
            .order_by("service_user_id")
            .values("service_user_id")
            .distinct()
        )
        for user_id in user_ids:
            user_list.append(
                [
                    user_id["service_user_id"],
                    django.contrib.auth.models.User.objects.filter(
                        pk__exact=user_id["service_user_id"]
                    )
                    .last()
                    .username,
                ]
            )
    return user_list


def get_current_users():
    """
    Description:
        The function returns the user that their session is active
    Return:
        User list
    """

    active_sessions = Session.objects.filter(expire_date__gte=timezone.now())
    user_id_list = []
    for session in active_sessions:
        data = session.get_decoded()
        user_id_list.append(data.get("_auth_user_id", None))
    # Query all logged in users based on id list
    return django.contrib.auth.models.User.objects.filter(id__in=user_id_list)


def save_database_configuration_value(configuration_name, configuration_value):
    """
    Description:
        The function saves configuration setting value. If not exists function create the configuration name
    Input:
        configuration_name       # configuration setting name
        configuration_value     # value for this configuration settings
    """
    if drylab.models.ConfigSetting.objects.filter(
        configuration_name__exact=configuration_name
    ).exists():
        config_settings_obj = drylab.models.ConfigSetting.objects.filter(
            configuration_name__exact=configuration_name
        ).last()
        config_settings_obj.set_configuration_value(configuration_value)
    else:
        config_settings_obj = drylab.models.ConfigSetting.objects.create_config_setting(
            configuration_name, configuration_value
        )
    return config_settings_obj
