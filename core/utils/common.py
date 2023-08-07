# Generic imports
import smtplib
import datetime
from django.conf import settings
from django.contrib.auth.models import User
from django.core.mail import send_mail

# local imports
import core.core_config
import core.models


def get_installed_apps():
    return settings.APPS_NAMES


def get_friend_list(user_name):
    """Function get the user names that are in their friend list

    Parameters
    ----------
    user_name : user object
        instance of the user

    Returns
    -------
    List
        The list contains user objects that are in the friend list
    """
    friend_list = []
    user_groups = user_name.groups.values_list("name", flat=True)
    if len(user_groups) > 0:
        for user in user_groups:
            if User.objects.filter(username__exact=user).exists():
                friend_list.append(User.objects.get(username__exact=user))

    friend_list.append(user_name)
    return friend_list


def get_inital_sample_settings_values(apps_name):
    """
    Description:
        The function get the defined values for species, origin of the sample,
        type of molecules, and type of protocol.
    Input:
        apps_name    # application which the data belongs to
    Return:
        initial_data
    """
    initial_data = {}
    if core.models.Species.objects.filter(apps_name__exact=apps_name).exists():
        species_objs = core.models.Species.objects.filter(apps_name__exact=apps_name)
        initial_data["species_data"] = []
        for species_obj in species_objs:
            initial_data["species_data"].append(
                [species_obj.get_name(), species_obj.get_id()]
            )
    if core.models.LabRequest.objects.filter(apps_name__exact=apps_name).exists():
        lab_request_objs = core.models.LabRequest.objects.filter(
            apps_name__exact=apps_name
        )
        initial_data["lab_request_data"] = []
        for lab_request_obj in lab_request_objs:
            data = lab_request_obj.get_all_data()
            data.append(lab_request_obj.get_id())
            initial_data["lab_request_data"].append(data)
    if core.models.MoleculeType.objects.filter(apps_name__exact=apps_name).exists():
        initial_data["molecule_type_data"] = []
        molecule_type_objs = core.models.MoleculeType.objects.filter(
            apps_name__exact=apps_name
        )
        for molecule_type_obj in molecule_type_objs:
            initial_data["molecule_type_data"].append(
                [molecule_type_obj.get_name(), molecule_type_obj.get_id()]
            )
    if core.models.ProtocolType.objects.filter(apps_name__exact=apps_name).exists():
        protocol_type_objs = core.models.ProtocolType.objects.filter(
            apps_name__exact=apps_name
        ).order_by("molecule")
        initial_data["protocol_type_data"] = []
        for protocol_type_obj in protocol_type_objs:
            initial_data["protocol_type_data"].append(
                [
                    protocol_type_obj.get_name(),
                    protocol_type_obj.get_id(),
                    protocol_type_obj.get_molecule_type(),
                ]
            )
    if core.models.StateInCountry.objects.filter(apps_name__exact=apps_name).exists():
        state_objs = core.models.StateInCountry.objects.filter(
            apps_name__exact=apps_name
        ).order_by("state_name")
        initial_data["states_data"] = []
        for state_obj in state_objs:
            initial_data["states_data"].append(
                [state_obj.get_state_name(), state_obj.get_state_id()]
            )
    if core.models.City.objects.filter(apps_name__exact=apps_name).exists():
        city_objs = core.models.City.objects.filter(
            apps_name__exact=apps_name
        ).order_by("city_name")
        initial_data["cities_data"] = []
        for city_obj in city_objs:
            initial_data["cities_data"].append(
                [city_obj.get_city_name(), city_obj.get_city_id()]
            )

    return initial_data


def get_email_data():
    """
    Description:
        Fetch the email configuration file
    Constant:
        EMAIL_HOST
        EMAIL_PORT
        EMAIL_HOST_USER
        EMAIL_HOST_PASSWORD
        EMAIL_USE_TLS
    Return:
        email_data
    """
    email_data = {}
    email_data["EMAIL_HOST"] = settings.EMAIL_HOST
    email_data["EMAIL_PORT"] = settings.EMAIL_PORT
    email_data["USER_EMAIL"] = settings.EMAIL_HOST_USER
    email_data["USER_PASSWORD"] = settings.EMAIL_HOST_PASSWORD
    email_data["USE_TLS"] = settings.EMAIL_USE_TLS
    return email_data


def sheet_header_to_field_name(header, field_info):
    """function that changes header/verbose names to 
    field names as defined in models.py

    Parameters
    ----------
    header
        list with the header names in the web excel table
    field_info
        dictionary with the map between the verbose names and the db field names

    Returns
    -------
        field_names: same header list but changing the verbose names to the field names.
        IMPORTANT: IN THE SAME ORDER!!
    """
    field_names = [
        name
        for item in header
        for name, verbose_name in field_info.items()
        if item == verbose_name
    ]
    return field_names


def jspreadsheet_to_dict(heading, data):
    """Convert the list of item list into a list where each item in the list
    is a dictionary and keys are the heading values

    Args:
        heading (list): the list which contains the key values
        data (list): List of data to be mapped to dictianary

    Returns:
        (list): List which each item has a dictionary with the heading as key values
    """
    c_data = [
        {heading[idx]: item[idx] for idx in range(len(heading))}
        for item in data
        if any(item)
    ]
    return c_data


def send_test_email(form_data):
    """
    Description:
        Get the configuration data from the user form and send a test email
    Constant:
        EMAIL_HOST
        EMAIL_PORT
        EMAIL_HOST_USER
        EMAIL_HOST_PASSWORD
        EMAIL_USE_TLS
    Return:
        email_data
    """
    settings.EMAIL_HOST = form_data["EMAIL_HOST"]
    settings.EMAIL_PORT = form_data["EMAIL_PORT"]
    settings.EMAIL_HOST_USER = form_data["USER_EMAIL"]
    settings.EMAIL_HOST_PASSWORD = form_data["USER_PASSWORD"]
    settings.EMAIL_USE_TLS = True if "USE_TLS" in form_data == "True" else False

    from_user = form_data["EMAIL_ISKYLIMS"]
    to_users = [form_data["test_email"]]
    subject = "testing email from iSlyLIMS"
    body_message = "This is a email test to verify iSkyLIMS"
    try:
        send_mail(subject, body_message, from_user, to_users)
        return "OK"
    except smtplib.SMTPException as e:
        return str(e)
    except Exception as e:
        return str(e)


def save_inital_sample_setting_value(apps_name, data):
    """
    Description:
        The function get the defined values for species, origin of the sample,
        type of molecules, and type of protocol.
    Input:
        apps_name    # application which the data belongs to
        data        # information to store on database
    Constants:
        ERROR_SPECIES_ALREADY_DEFINED
        ERROR_LABORATORY_REQUEST_ALREADY_DEFINED
        ERROR_MOLECULE_TYPE_ALREADY_DEFINED
        ERROR_PROTOCOL_TYPE_ALREADY_DEFINED
    Return:
        setting_defined
    """
    setting_defined = {}
    if "species" in data:
        species_data = {}
        species_data["apps_name"] = apps_name
        species_data["name"] = data["species"]
        if core.models.Species.objects.filter(
            species_name__iexact=species_data["name"],
            apps_name__exact=species_data["apps_name"],
        ).exists():
            setting_defined["ERROR"] = [
                core.core_config.ERROR_SPECIES_ALREADY_DEFINED,
                species_data["name"],
            ]
            return setting_defined

        core.models.Species.objects.create_new_specie(species_data)
        setting_defined["settings"] = "Species"
        setting_defined["value"] = species_data["name"]

    if "lab_request" in data:
        lab_request_data = {}
        lab_request_data["apps_name"] = apps_name
        lab_request_data["labName"] = data["lab_request"]["labRequestName"]
        lab_request_data["labNameCoding"] = data["lab_request"]["labRequesCoding"]
        lab_request_data["labUnit"] = data["lab_request"]["department"]
        lab_request_data["labContactName"] = data["lab_request"]["contact"]
        lab_request_data["labPhone"] = data["lab_request"]["phone"]
        lab_request_data["labEmail"] = data["lab_request"]["email"]
        lab_request_data["address"] = data["lab_request"]["address"]
        lab_request_data["city"] = data["lab_request"]["city"]
        if core.models.LabRequest.objects.filter(
            lab_name_coding__iexact=lab_request_data["labNameCoding"],
            apps_name__exact=lab_request_data["apps_name"],
        ).exists():
            setting_defined["ERROR"] = [
                core.core_config.ERROR_LABORATORY_REQUEST_ALREADY_DEFINED,
                lab_request_data["labNameCoding"],
            ]
            return setting_defined
        core.models.LabRequest.objects.create_lab_request(lab_request_data)
        setting_defined["settings"] = "Lab Request"
        setting_defined["value"] = lab_request_data["labName"]
    if "molecule_type" in data:
        molecule_type_data = {}
        molecule_type_data["apps_name"] = apps_name
        molecule_type_data["moleculeType"] = data["molecule_type"]
        if core.models.MoleculeType.objects.filter(
            molecule_type__iexact=molecule_type_data["moleculeType"],
            apps_name__exact=molecule_type_data["apps_name"],
        ).exists():
            setting_defined["ERROR"] = [
                core.core_config.ERROR_MOLECULE_TYPE_ALREADY_DEFINED,
                molecule_type_data["moleculeType"],
            ]
            return setting_defined
        core.models.MoleculeType.objects.create_molecule_type(molecule_type_data)
        setting_defined["settings"] = "Molecule Type"
        setting_defined["value"] = molecule_type_data["moleculeType"]

    if "protocol_type" in data:
        protocol_type_data = {}
        protocol_type_data["apps_name"] = apps_name
        protocol_type_data["protocol_type"] = data["protocol_type"][0]
        if data["protocol_type"][1] == "None":
            protocol_type_data["molecule"] = None
        else:
            protocol_type_data["molecule"] = data["protocol_type"][1]
        if core.models.ProtocolType.objects.filter(
            protocol_type__iexact=protocol_type_data["protocol_type"],
            molecule__molecule_type__iexact=protocol_type_data["molecule"],
            apps_name__exact=protocol_type_data["apps_name"],
        ).exists():
            setting_defined["ERROR"] = [
                core.core_config.ERROR_PROTOCOL_TYPE_ALREADY_DEFINED,
                protocol_type_data["protocol_type"],
            ]
            return setting_defined
        core.models.ProtocolType.objects.create_protocol_type(protocol_type_data)
        setting_defined["settings"] = "Protocol Type"
        setting_defined["value"] = protocol_type_data["protocol_type"]

    if "state" in data:
        state_data = {}
        if core.models.StateInCountry.objects.filter(
            state_name__iexact=data["state"], apps_name__exact=apps_name
        ).exists():
            setting_defined["ERROR"] = [
                core.core_config.ERROR_STATE_ALREADY_DEFINED,
                data["state"],
            ]
            return setting_defined
        state_data["apps_name"] = apps_name
        state_data["state"] = data["state"]
        core.models.StateInCountry.objects.create_new_state(state_data)
        setting_defined["settings"] = "State"
        setting_defined["value"] = data["state"]

    if "city" in data:
        city_data = {}
        if core.models.City.objects.filter(
            city_name__iexact=data["city"]["cityName"], apps_name__exact=apps_name
        ).exists():
            setting_defined["ERROR"] = [
                core.core_config.ERROR_CITY_ALREADY_DEFINED,
                data["city"]["cityName"],
            ]
            return setting_defined
        city_data["apps_name"] = apps_name
        city_data["state"] = data["city"]["state"]
        city_data["cityName"] = data["city"]["cityName"]
        city_data["latitude"] = data["city"]["latitude"]
        city_data["longitude"] = data["city"]["longitude"]
        core.models.City.objects.create_new_city(city_data)
        setting_defined["settings"] = "core.models.City"
        setting_defined["value"] = data["city"]["cityName"]
    return setting_defined


def week_month_number_to_date(input_data, type_number, value_param=None, format=None):
    """Convert year + number of week or year + number of month into a complete date format.
    If value_param is set the converted date is used as key and value_param as value.

    Parameters
    ----------
    input_data : QuerySet list
        List of dictionnaries, with keys "year" and "week" used to map, and
        optional the value_param to create a dictionnary.
    type_number : string
        Type of number recived, number of week or number of month
    value_param: str, optional
        key in the dictionnary where the value is located.
    format : str, optional
        output of the date, in string or as datetime object if format is None,
        by default None
    Return:
        List of Dictionnaries if value_param is set
        List if value_param is Nome
    """
    output_date = []
    for record in input_data:
        if type_number == "week":
            week = "{year}-W{week}-1".format(year=record["year"], week=record["week"])
            timestamp = datetime.datetime.strptime(week, "%Y-W%W-%w")
        else:
            month = "{year}-M{month}-1".format(
                year=record["year"], month=record["month"]
            )
            timestamp = datetime.datetime.strptime(month, "%Y-M%m-%M")
        if format is not None:
            timestamp = datetime.datetime.strftime(timestamp, format)
        if value_param is not None:
            date_value = {}
            date_value[timestamp] = record[value_param]
            output_date.append(date_value)
        else:
            output_date.append(timestamp)
    return output_date
