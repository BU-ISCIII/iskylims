import json, re
from iSkyLIMS_core.core_config import *
from iSkyLIMS_core.models import *
from django.contrib.auth.models import User


def create_new_protocol(new_protocol, protocol_type, description, apps_name):
    """
    Description:
        The function create a new item in database.
    Input:
        new_protocol # new protocol name
        protocol_type # protocol type for the new protocol
    Return:
        ID of the new created object.
    """
    protocol_type_obj = ProtocolType.objects.get(
        protocol_type__exact=protocol_type, apps_name__exact=apps_name
    )
    new_protocol_object = Protocols(
        type=protocol_type_obj, name=new_protocol, description=description
    )
    new_protocol_object.save()
    return new_protocol_object.pk


def check_if_protocol_exists(protocol, app_name):
    """
    Description:
        The function return True if protocol exists. False if not
    Input:
        protocol # "protocol name" or "protocol id" to be checked
    Return:
        True/False.
    """
    if type(protocol) is int:
        if Protocols.objects.filter(
            pk__exact=protocol, type__apps_name__exact=app_name
        ).exists():
            return True
        else:
            return False
    else:
        if Protocols.objects.filter(
            name__exact=protocol, type__apps_name__exact=app_name
        ).exists():
            return True
        else:
            return False


def check_if_protocol_parameters_exists(protocol):
    if ProtocolParameters.objects.filter(protocol_id__name__exact=protocol).exists():
        return True
    else:
        return False


def define_table_for_prot_parameters(protocol_id):
    """
    Description:
        The function return a dictionary with the information to create the table
        for defining the parameters used in the protocol
    Input:
        protocol_id # protocol id  to get protocol information
    Return:
        prot_parameters
    """
    prot_parameters = {}
    protocol_obj = Protocols.objects.get(pk__exact=protocol_id)

    prot_parameters["protocol_name"] = protocol_obj.get_name()
    prot_parameters["protocol_id"] = protocol_id
    prot_parameters["heading"] = HEADING_FOR_DEFINING_PROTOCOL_PARAMETERS
    return prot_parameters


def display_available_protocols(app_name):
    """
    Description:
        The function return a list with all defined protocols that contains
        molecule definition. This means to exclude any other protocol that their
        parameters are not stored in iSkyLIMS_core.
    Return:
        molecule_protocol_list and other_protocol_list .
    """

    molecule_protocol_list = []
    if (
        ProtocolType.objects.filter(apps_name__exact=app_name)
        .exclude(molecule=None)
        .exists()
    ):
        protocol_types = (
            ProtocolType.objects.filter(apps_name__exact=app_name)
            .exclude(molecule=None)
            .order_by("molecule")
        )
        for protocol_type in protocol_types:
            molecule_type = protocol_type.get_molecule_type()
            prot_type_str = protocol_type.get_name()
            if Protocols.objects.filter(type=protocol_type).exists():
                protocols = Protocols.objects.filter(type=protocol_type).order_by(
                    "type"
                )
                for protocol in protocols:
                    data_prot = []
                    data_prot.append(molecule_type)
                    data_prot.append(prot_type_str)
                    data_prot.append(protocol.get_name())
                    data_prot.append(protocol.pk)
                    if ProtocolParameters.objects.filter(protocol_id=protocol).exists():
                        data_prot.append(True)
                    else:
                        data_prot.append(False)
                    molecule_protocol_list.append(data_prot)
    other_protocol_list = []
    if ProtocolType.objects.filter(molecule=None, apps_name__exact=app_name).exists():
        protocol_types = ProtocolType.objects.filter(
            molecule=None, apps_name__exact=app_name
        )
        for protocol_type in protocol_types:
            prot_type_str = protocol_type.get_name()
            if Protocols.objects.filter(type=protocol_type).exists():
                protocols = Protocols.objects.filter(type=protocol_type).order_by(
                    "type"
                )
                for protocol in protocols:
                    data_prot = []
                    data_prot.append(prot_type_str)
                    data_prot.append(protocol.get_name())
                    data_prot.append(protocol.pk)
                    if ProtocolParameters.objects.filter(protocol_id=protocol).exists():
                        data_prot.append(True)
                    else:
                        data_prot.append(False)
                    other_protocol_list.append(data_prot)

    return molecule_protocol_list, other_protocol_list


def display_protocol_list():
    """
    Description:
        The function return the protocol list defined.
    Return:
        protocol_list
    """
    protocol_list = []
    if Protocols.objects.all().exists():
        protocols_objs = Protocols.objects.all().order_by("type").order_by("name")
        for protocol_obj in protocols_objs:
            protocol_list.append(
                [protocol_obj.get_protocol_id(), protocol_obj.get_name()]
            )
    return protocol_list


def get_defined_protocols(app_name, exclude_non_molecule):
    """
    Description:
        The function return a list with all protocol defined.
        If exclude_non_molecule = True, then only protocols which involves molecules are returned
        If false return all protocols defined
    Input:
        app_name : Name of the application to return only the protocols defined under this application
        exclude_non_molecule : True/False if exclude  non molecules are requested
    Variables:
        defined_protocols: list with all protocols
    Return:
        defined_protocols.
    """
    defined_protocols = []
    if exclude_non_molecule:
        if (
            Protocols.objects.filter(type__apps_name__exact=app_name)
            .exclude(type__molecule=None)
            .exists()
        ):
            protocols_obj = (
                Protocols.objects.filter(type__apps_name__exact=app_name)
                .exclude(type__molecule=None)
                .order_by("type")
            )
            for protocol_obj in protocols_obj:
                defined_protocols.append(protocol_obj.get_name())
    else:
        if Protocols.objects.filter(type__apps_name__exact=app_name).exists():
            protocols_obj = Protocols.objects.filter(
                type__apps_name__exact=app_name
            ).order_by("type")
            for protocol_obj in protocols_obj:
                defined_protocols.append(protocol_obj.get_name())
    return defined_protocols


def get_protocol_from_prot_types(prot_types):
    protocols = {}

    for prot_type in prot_types:
        prot_names = []
        prots = Protocols.objects.filter(type__protocol_type__exact=prot_type)
        for prot in prots:
            prot_names.append(prot.get_name())
        protocols[prot_type] = prot_names
    return protocols


def display_protocol_types(app_name):
    """
    Description:
        The function return a list with all protocol types defined.
    Return:
        protocol_types_list.
    """
    protocol_types_list = []
    if ProtocolType.objects.filter(apps_name__exact=app_name).exists():
        protocol_types = ProtocolType.objects.filter(
            apps_name__exact=app_name
        ).order_by("molecule")
        for protocol_type in protocol_types:
            protocol_types_list.append(protocol_type.get_name())
    return protocol_types_list


def get_all_protocol_info(protocol_id):
    """
    Description:
        The function return a dictionary with all definition parameters fro the given protocol.
    Return:
        protocol_data.
    """
    protocol_data = {}
    protocol_data["parameters"] = []
    protocol_obj = Protocols.objects.get(pk__exact=protocol_id)

    if ProtocolParameters.objects.filter(protocol_id=protocol_obj).exists():
        protocol_data["parameter_heading"] = HEADING_FOR_DEFINING_PROTOCOL_PARAMETERS
        protocol_data["protocol_name"] = protocol_obj.get_name()
        protocol_parameters = ProtocolParameters.objects.filter(
            protocol_id=protocol_obj
        ).order_by("parameterOrder")
        for parameter in protocol_parameters:
            protocol_data["parameters"].append(parameter.get_all_parameter_info())
        protocol_data["protocol_id"] = protocol_id
    return protocol_data


def get_protocol_obj_from_name(protocol_name):
    if Protocols.objects.filter(name__exact=protocol_name).exists():
        return Protocols.objects.get(name__exact=protocol_name)
    else:
        return None


def get_protocol_obj_from_id(protocol_id):
    if Protocols.objects.filter(pk__exact=protocol_id).exists():
        return Protocols.objects.get(pk__exact=protocol_id)
    else:
        return None


def get_protocol_fields(protocol_id):
    """
    Description:
        The function return the parameters definded for the protocol id.
    Input:
        protocol_id       # id of the protocol
    Return:
        info_protocol
    """
    parameters_protocol = {}
    protocol_obj = Protocols.objects.get(pk__exact=protocol_id)
    if ProtocolParameters.objects.filter(protocol_id__exact=protocol_obj).exists():
        protocol_parameter_objs = ProtocolParameters.objects.filter(
            protocol_id__exact=protocol_obj
        ).order_by("parameterOrder")
        parameter_list = []
        for protocol_parameter_obj in protocol_parameter_objs:
            parameter_data = protocol_parameter_obj.get_protocol_fields_for_javascript()
            parameter_data.insert(1, "")
            parameter_data.append(protocol_parameter_obj.get_parameter_protocol_id())
            parameter_list.append(parameter_data)

        parameters_protocol["heading"] = HEADING_FOR_MODIFY_PROTOCOL_FIELDS
        parameters_protocol["protocol_id"] = protocol_id
        parameters_protocol["protocol_name"] = protocol_obj.get_name()
        parameters_protocol["fields"] = parameter_list

    return parameters_protocol


def get_protocol_parameters(protocol_obj):
    """
    Description:
        The function return a list of the used parameters .
    Return:
        protocol_parameter_list.
    """
    protocol_parameter_list = []
    if ProtocolParameters.objects.filter(protocol_id=protocol_obj).exists():
        protocol_parameters = ProtocolParameters.objects.filter(
            protocol_id=protocol_obj, parameterUsed=True
        ).order_by("parameterOrder")
        for protocol_parameter in protocol_parameters:
            protocol_parameter_list.append(protocol_parameter.get_parameter_name())
    return protocol_parameter_list


def get_protocol_parameters_and_type(protocol_obj):
    """
    Description:
        The function return a list of the used parameters and types.
    Return:
        protocol_parameter_type_list.
    """
    protocol_parameter_type_list = []
    if ProtocolParameters.objects.filter(protocol_id=protocol_obj).exists():
        protocol_parameters = ProtocolParameters.objects.filter(
            protocol_id=protocol_obj, parameterUsed=True
        ).order_by("parameterOrder")
        for protocol_parameter in protocol_parameters:
            heading_item = []
            heading_item.append(protocol_parameter.get_parameter_name())
            heading_item.append(protocol_parameter.get_parameter_type())
            heading_item.append(
                protocol_parameter.get_parameter_option_values().split(",")
            )
            protocol_parameter_type_list.append(heading_item)
    return protocol_parameter_type_list


def get_protocol_parameter_obj_from_id(protocol_parameter_id):
    """
    Description:
        The function return the protocol parameter obj
    Input:
        protocol_parameter_id       # id of the protocol parameter
    Return:
        protocol_parameter_obj
    """
    if ProtocolParameters.objects.filter(pk__exact=protocol_parameter_id).exists():
        return ProtocolParameters.objects.get(pk__exact=protocol_parameter_id)
    return None


def get_project_name_by_id(protocol_id):
    """
    Description:
        The function return a protocol_name for a ginven protocol id.
    Return:
        protocol_name.
    """
    if Protocols.objects.filter(pk__exact=protocol_id).exists():
        protocol_name = Protocols.objects.get(pk__exact=protocol_id).get_name()
        return protocol_name
    else:
        "None"


def modify_fields_in_protocol(form_data):
    """
    Description:    The function get the protocol field value and check if there is
        some changes. If change then replace the old values by thenew ones
    Input:
        form_data     # form data from user
    Return:
        saved_fields #
    """
    saved_fields = {}
    protocol_id = form_data["protocol_id"]
    protocol_obj = get_protocol_obj_from_id(protocol_id)
    saved_fields["fields"] = []
    saved_fields["heading"] = HEADING_FOR_DEFINING_PROTOCOL_PARAMETERS
    saved_fields["protocol_name"] = protocol_obj.get_name()

    json_data = json.loads(form_data["table_data1"])
    fields = HEADING_FOR_MODIFY_PROTOCOL_FIELDS
    for row_data in json_data:
        if row_data[0] == "" and row_data[1] == "":
            continue
        p_fields = {}

        for i in range(len(fields)):
            p_fields[fields[i]] = row_data[i]

        if row_data[fields.index("Field type")] == "Option List":
            option_list_values = row_data[fields.index("Option Values")].split(",")
            clean_value_list = []
            for opt_value in option_list_values:
                value = opt_value.strip()
                if value != "":
                    clean_value_list.append(value)

            p_fields["Option Values"] = ",".join(clean_value_list)
        else:
            p_fields["Option Values"] = ""
        if row_data[0] == "" and row_data[1] != "":
            # new field
            p_fields["Parameter name"] = row_data[1]
            p_fields["protocol_id"] = protocol_obj
            saved_fields["fields"].append(
                ProtocolParameters.objects.create_protocol_parameter(
                    p_fields
                ).get_all_parameter_info()
            )
            continue
        if row_data[0] != "" and row_data[1] != "":
            # rename field name
            p_fields["Field name"] = row_data[1]
        else:
            p_fields["Field name"] = row_data[0]
        # Update  Field
        protocol_parameter_obj = get_protocol_parameter_obj_from_id(row_data[-1])
        if not protocol_parameter_obj:
            # Unable to find the object class. Skipping this change
            continue
        protocol_parameter_obj.update_protocol_fields(p_fields)
        saved_fields["fields"].append(protocol_parameter_obj.get_all_parameter_info())
    return saved_fields


def set_protocol_parameters(request):
    protocol_id = request.POST["protocol_id"]
    json_data = json.loads(request.POST["table_data1"])
    parameters = HEADING_FOR_DEFINING_PROTOCOL_PARAMETERS
    protocol_id_obj = Protocols.objects.get(pk__exact=protocol_id)

    saved_parameters = []
    stored_parameters = {}
    for row_data in json_data:

        if row_data[0] == "":
            continue
        prot_parameters = {}

        prot_parameters["protocol_id"] = protocol_id_obj
        for i in range(len(parameters)):
            prot_parameters[parameters[i]] = row_data[i]

        if row_data[parameters.index("Parameter Type")] == "Option List":
            option_list_values = row_data[parameters.index("Option Values")].split(",")
            clean_value_list = []
            for opt_value in option_list_values:
                value = opt_value.strip()
                if value != "":
                    clean_value_list.append(value)

            prot_parameters["Option Values"] = ",".join(clean_value_list)
        else:
            prot_parameters["Option Values"] = ""

        saved_parameters.append(
            ProtocolParameters.objects.create_protocol_parameter(
                prot_parameters
            ).get_parameter_name()
        )
    stored_parameters["parameters"] = saved_parameters
    stored_parameters["protocol_name"] = protocol_id_obj.get_name()

    return stored_parameters
