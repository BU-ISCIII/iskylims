# Generic imports
import logging
import os
import re
import socket
import sys
import traceback
from datetime import datetime
from logging.config import fileConfig

from django.contrib.auth.models import User
from django.core.mail import send_mail
from smb.SMBConnection import SMBConnection

# Local imports
import core.models
import wetlab.config
import wetlab.models


def get_configuration_value(parameter_name):
    """
    Description:
        Function will get the parameter value defined in the configutration table
        if not exists return 'False'

    Input:
        parameter_name    #parameter name
    Return:
        parameter_value
    """
    parameter_value = "False"
    if wetlab.models.ConfigSetting.objects.filter(
        configuration_name__exact=parameter_name
    ).exists():
        parameter_obj = wetlab.models.ConfigSetting.objects.filter(
            configuration_name__exact=parameter_name
        ).last()
        parameter_value = parameter_obj.get_configuration_value()
    return parameter_value


def get_run_in_same_year_to_compare(run_object):
    """
    Description:
        The function return a list with all runs that have been created on the same year and they
        are based on the same chemistry.
    Input:
        run_object    # runProcess object
    Return:
        same_run_in_year run_object list
    """
    # get the chemistry type for the run, that will be used to compare runs with the same chemistry value
    chem_high_mid = wetlab.models.RunningParameters.objects.get(
        run_name_id__exact=run_object
    ).get_run_chemistry()
    run_different_chemistry = wetlab.models.RunningParameters.objects.all().exclude(
        chemistry__exact=chem_high_mid
    )
    run_year = run_object.get_run_year()

    start_date = str(run_year) + "-1-1"
    end_date = str(run_year) + "-12-31"
    same_run_in_year = wetlab.models.RunProcess.objects.filter(
        run_date__range=(start_date, end_date)
    ).exclude(run_name__in=run_different_chemistry)
    return same_run_in_year


def check_valid_date_format(date):
    try:
        datetime.strptime(date, "%Y-%m-%d")
        return True
    except ValueError:
        try:
            datetime.strptime(date, "%d-%m-%Y")
            return True
        except ValueError:
            return False


def get_samba_connection_data():
    """
    Description:
        Fetch the samba configuration from database
    Return:
        samba_data
    """
    samba_data = {}
    if wetlab.models.SambaConnectionData.objects.all().exists():
        samba_connection_obj = wetlab.models.SambaConnectionData.objects.all().last()
        samba_data = samba_connection_obj.get_samba_data()
    return samba_data


def save_samba_connection_data(data):
    """
    Description:
        store the information in database. If it is the first attempt to store data
        the instance is create if not the table is updated with the new data
    Input:
        data    # Dictionary having samba connection data
    Return:
        samba_connection_obj
    """
    if not wetlab.models.SambaConnectionData.objects.all().exists():
        samba_connection_obj = wetlab.models.SambaConnectionData.objects.create()
    else:
        samba_connection_obj = wetlab.models.SambaConnectionData.objects.all().last()
    samba_connection_obj.update_data(data)
    return samba_connection_obj


def open_samba_connection():
    """
    Description:
        The function open a samba connection with the parameter settings
        defined in wetlab configuration file
    Functions:
        get_samba_connection_data   # located at this file
    Return:
        conn object for the samba connection
    """
    logger = logging.getLogger(__name__)
    logger.debug("Starting function open_samba_connection")
    samba_data = get_samba_connection_data()

    try:
        if samba_data["samba_folder_name"] != "":
            samba_data["shared_folder_name"] = os.path.join(
                samba_data["shared_folder_name"],
                samba_data["samba_folder_name"],
            )
    except KeyError:
        string_message = "Samba connection data on database is empty"
        logging_errors(string_message, True, False)
        raise

    conn = SMBConnection(
        samba_data["user_id"],
        samba_data["user_password"],
        samba_data["shared_folder_name"],
        samba_data["remote_server_name"],
        use_ntlm_v2=samba_data["ntlm_used"],
        domain=samba_data["domain"],
        is_direct_tcp=samba_data["is_direct_tcp"],
    )

    if samba_data["host_name"]:
        conn.connect(
            socket.gethostbyname(samba_data["host_name"]),
            int(samba_data["port_server"]),
        )
    else:
        conn.connect(samba_data["ip_server"], int(samba_data["port_server"]))

    logger.debug("End function open_samba_connection")
    return conn


def find_xml_tag_text(input_file, search_tag):
    """
    Description:
        The function will look for the xml element tag in the file and it will return the text value
    Input:
        input_file  # file to find the tag
        search_tag  # xml tag to be found in the input_file
    Return:
        tag value
    """
    fh = open(input_file, "r")
    search_line = "<" + search_tag + ">(.*)</" + search_tag + ">"
    for line in fh:
        found_tag = re.search("^\s+ %s" % search_line, line)
        if found_tag:
            fh.close()
            return found_tag.group(1)
    fh.close()
    return "NOT FOUND"


def get_attributes_remote_file(conn, run_dir, remote_file):
    """
    Description:
        Function will fetch the file from remote server and copy on local
        directory
    Input:
        conn    # Samba connection object
        run_dir # run folder to fetch the file
        remote_file # file name to fetch on remote server
    Constants:
        SAMBA_SHARED_FOLDER_NAME
    Return:
        file_attributes if sucessful .
        Exception if attributes could not be fetched
    """
    logger = logging.getLogger(__name__)
    logger.debug("Starting function for getting remote attributes")
    try:
        file_attributes = conn.getAttributes(
            wetlab.config.SAMBA_SHARED_FOLDER_NAME, remote_file
        )
        logger.info("Got attributes from %s", remote_file)
    except Exception:
        string_message = "Unable to get attributes for " + remote_file
        logging_errors(string_message, True, False)
        raise Exception("Not get attributes")
    logger.debug("End function for  getting remote attributes")
    return file_attributes


def get_experiment_name_from_file(l_run_parameter):
    """
    Description:
        The function will get the experiment name  for the xml element tag in the
        file and it will return the experiment name value
    Input:
        l_run_parameter  # file to find the tag
    Functions:
        find_xml_tag_text # located at this file
    Variables:
        experiment_name # name of the experiment found in runParameter file
    Return:
        experiment_name
    """

    experiment_name = find_xml_tag_text(
        l_run_parameter, wetlab.config.EXPERIMENT_NAME_TAG
    )

    return experiment_name


def get_configuration_from_database(configuration_name):
    """
    Description:
        The function fetch from database the configuration setting value
    Input:
        configuration_name      # configuration settings name
    """
    configuration_value = ""
    if wetlab.models.ConfigSetting.objects.filter(
        configuration_name__exact=configuration_name
    ).exists():
        configuration_settings_obj = wetlab.models.ConfigSetting.objects.filter(
            configuration_name__exact=configuration_name
        ).last()
        configuration_value = configuration_settings_obj.get_configuration_value()
    return configuration_value


def get_allowed_user_for_sharing(request_user):
    """
    Description:
        The function get the primary key of the that are allowed view information
    Input:
        request_user      # user obj
    Return:
        sharing_list
    """
    # getting projects from user sharing list
    sharing_list = []
    user_groups = request_user.groups.values_list("name", flat=True)
    for user in user_groups:
        if User.objects.filter(username__exact=user).exists():
            sharing_list.append(User.objects.get(username__exact=user).id)
    sharing_list.append(request_user.id)
    return sharing_list


def get_userid_list():
    """
    Descripion:
        The function get the userid list defined in the system
    Return:
        user_id_list
    """
    user_id_list = []
    user_objs = User.objects.all()
    for user_obj in user_objs:
        user_id_list.append(user_obj.username)
    return user_id_list


def is_wetlab_manager(request):
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
    if request.user.groups.filter(name=wetlab.config.WETLAB_MANAGER).exists():
        return True
    return False


def normalized_data(set_data, all_data):
    """
    Description:
        The function is used to normalized data from diferent range of values
    Input:
        set_data    # contains a gruop of data
        all_data    # contains all data to be used for the normalization
    Variables:
        normalized_set_data # to keep the normalized set data
        normalized_all_data # to keep the normalized value of all data
    Return:
        normalized_set_data
        normalized_all_data.
    """
    normalized_set_data, normalized_all_data = [], []

    min_value = min(min(set_data), min(all_data))
    max_value = max(max(set_data), max(all_data))
    for value in set_data:
        normalized_set_data.append(format((value - min_value) / max_value, ".2f"))
    for value in all_data:
        normalized_all_data.append(format((value - min_value) / max_value, ".2f"))

    return normalized_set_data, normalized_all_data


def get_project_search_fields_form():
    project_form_data = {}

    project_form_data["run_states"] = []
    project_form_data["available_platforms"] = []
    project_form_data["available_sequencers"] = []
    run_states = wetlab.models.RunStates.objects.all()
    for r_state in run_states:
        project_form_data["run_states"].append(r_state.get_run_state_name())

    platforms = core.models.SequencingPlatform.objects.all()
    for platform in platforms:
        project_form_data["available_platforms"].append(platform.get_platform_name())
    sequencers = core.models.SequencerInLab.objects.all()
    for sequencer in sequencers:
        project_form_data["available_sequencers"].append(sequencer.get_sequencer_name())

    return project_form_data


def get_run_search_fields_form():
    run_form_data = {}

    run_form_data["run_states"] = []
    run_form_data["available_platforms"] = []
    run_form_data["available_sequencers"] = []
    run_states = wetlab.models.RunStates.objects.all()
    for r_state in run_states:
        run_form_data["run_states"].append(r_state.get_run_state_name())

    platforms = core.models.SequencingPlatform.objects.all()
    for platform in platforms:
        run_form_data["available_platforms"].append(platform.get_platform_name())
    machines = core.models.SequencerInLab.objects.all()
    for machine in machines:
        run_form_data["available_sequencers"].append(machine.get_sequencer_name())

    return run_form_data


def save_database_configuration_value(configuration_name, configuration_value):
    """
    Description:
        The function saves configuration setting value. If not exists function create the configuration name
    Input:
        configurationName       # configuration setting name
        configuration_value     # value for this configuration settings
    """
    if wetlab.models.ConfigSetting.objects.filter(
        configuration_name__exact=configuration_name
    ).exists():
        config_settings_obj = wetlab.models.ConfigSetting.objects.filter(
            configuration_name__exact=configuration_name
        ).last()
        config_settings_obj.set_configuration_value(configuration_value)
    else:
        config_settings_obj = wetlab.models.ConfigSetting.objects.create_config_setting(
            configuration_name, configuration_value
        )
    return config_settings_obj


def logging_errors(string_text, showing_traceback, print_on_screen):
    """
    Description:
        The function will log the error information to file.
        Optional can send an email to inform about the issue
    Input:
        logger # contains the logger object
        string_text # information text to include in the log
    Functions:

    Constant:
        SENT_EMAIL_ON_ERROR
        EMAIL_FOR_NOTIFICATIONS
    Variables:
        subject # text to include in the subject email
    """
    logger = logging.getLogger(__name__)
    logger.error("-----------------    ERROR   ------------------")
    logger.error(string_text)
    if wetlab.models.ConfigSetting.objects.filter(
        configuration_name__exact="SENT_EMAIL_ON_ERROR"
    ).exists():
        email_on_error_obj = wetlab.models.ConfigSetting.objects.filter(
            configuration_name__exact="SENT_EMAIL_ON_ERROR"
        ).last()
        if email_on_error_obj.get_configuration_value() == "TRUE":
            if wetlab.models.ConfigSetting.objects.filter(
                configuration_name__exact="EMAIL_FOR_NOTIFICATIONS"
            ).exists():
                email_on_notification_obj = wetlab.models.ConfigSetting.objects.filter(
                    configuration_name__exact="EMAIL_FOR_NOTIFICATIONS"
                ).last()
                email_notification = email_on_notification_obj.get_configuration_value()
                if "@" in email_notification:
                    subject = "Error found on wetlab when running crontab"
                    try:
                        send_mail(
                            subject,
                            string_text,
                            email_notification,
                            [email_notification],
                        )
                    except Exception:
                        logger.error(
                            "*************UNABLE TO SEND ERROR EMAIL TO USER *****************"
                        )
                else:
                    logger.error(
                        "****** INVALID EMAIL FORMAT.  EMAIL IS NOT SENT ***************"
                    )
    if showing_traceback:
        logger.error("################################")
        logger.error(traceback.format_exc())
        logger.error("################################")
    logger.error("-----------------    END ERROR   --------------")
    if print_on_screen:
        from datetime import datetime

        print("********* ERROR **********")
        print(string_text)
        print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        print("Check log for detail information")
        print("******* END ERROR ********")
    return ""


def logging_warnings(string_text, print_on_screen):
    """
    Description:
        The function will log the error information to file.
        Optional can send an email to inform about the issue
    Input:
        logger # contains the logger object
        string_text # information text to include in the log
    """
    logger = logging.getLogger(__name__)
    logger.warning("-----------------    WARNING   ------------------")
    logger.warning(string_text)
    logger.warning("-----------------    END WARNING   --------------")
    if print_on_screen:
        from datetime import datetime

        print("******* WARNING ********")
        print(string_text)
        print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        print("Check log for detail information")
        print("**** END WARNING *******")
    return ""


def open_log(config_file):
    """
    Description:
        The function will create the log object to write all logging information
    Input:
        logger_name    # contains the logger name that will be included
                        in the log file
    Constant:
        LOGGING_CONFIG_FILE
    Return:
        logger object
    """
    fileConfig(config_file)
    logger = logging.getLogger(__name__)
    return logger
