import json

from iSkyLIMS_clinic.clinic_config import *
from iSkyLIMS_clinic.models import *
from iSkyLIMS_clinic.utils.managed_samples import *
from iSkyLIMS_core.utils.handling_protocols import *


def add_result_protocol_parameters(request):
    c_sample_updated_list = []
    result_parameter_value = {}
    json_data = json.loads(request.POST["parameters_data"])
    c_samples = request.POST["c_samples"].split(",")
    parameter_heading = request.POST["heading_in_excel"].split(",")
    parameters_length = len(json_data[0])
    protocol_used_obj = get_sample_clinic_obj_from_id(c_samples[0]).get_protocol_obj()
    fixed_heading_length = parameters_length - len(parameter_heading)
    for row_index in range(len(json_data)):
        c_sample_obj = get_sample_clinic_obj_from_id(c_samples[row_index])

        empty_values = True
        for p_index in range(fixed_heading_length, parameters_length):
            result_parameter_value["parameter_id"] = ProtocolParameters.objects.get(
                protocol_id=protocol_used_obj,
                parameterName__exact=parameter_heading[p_index - fixed_heading_length],
            )
            result_parameter_value["clinicSample_id"] = c_sample_obj
            value = json_data[row_index][p_index]
            result_parameter_value["parameterValue"] = value
            ResultParameterValue.objects.create_result_parameter_value(
                    result_parameter_value
                )
            if value != "":
                empty_values = False
        if not empty_values:
            c_sample_obj.set_state("Completed")
            c_sample_updated_list.append(c_sample_obj.get_id())
    return c_sample_updated_list


def define_table_for_result_parameters(c_sample_ids):
    """
    Description:

    Input:

    Return:
        prot_parameters
    """
    result_parameters = {}
    prot_used_in_table = ""
    showed_c_samples = []
    pending_c_samples = []

    for c_sample in c_sample_ids:
        c_sample_obj = get_sample_clinic_obj_from_id(c_sample)
        prot_in_c_sample = c_sample_obj.get_protocol()
        if prot_used_in_table == "":
            prot_used_in_table = prot_in_c_sample
            protocol_used_obj = c_sample_obj.get_protocol_obj()
            if ProtocolParameters.objects.filter(
                protocol_id__exact=protocol_used_obj
            ).exists():
                protocol_parameters = ProtocolParameters.objects.filter(
                    protocol_id__exact=protocol_used_obj, paramete_used=True
                ).order_by("parameter_order")
                parameter_list = []
                for parameter in protocol_parameters:
                    parameter_list.append(parameter.get_parameter_name())
                length_param_heading = len(parameter_list)
                result_parameters["fix_heading"] = HEADING_FOR_RESULT_TO_PROTOCOL
                result_parameters["param_heading"] = parameter_list
                result_parameters["data"] = []
        if prot_used_in_table == prot_in_c_sample:
            showed_c_samples.append(c_sample)
            data_param = [""] * length_param_heading
            data_c_sample = c_sample_obj.get_sample_info_for_protocol()
            data = data_c_sample + data_param
            result_parameters["data"].append(data)

        else:
            pending_c_samples.append(c_sample_obj.get_id())

    if len(pending_c_samples) > 0:
        result_parameters["pending_id"] = ",".join(pending_c_samples)
    result_parameters["c_sample_id"] = ",".join(showed_c_samples)
    result_parameters["heading_in_excel"] = ",".join(parameter_list)
    return result_parameters


def get_table_result_to_protocol(c_samples_pending_protocol):
    result_protocol = {}
    clinic_samples = []
    result_protocol["data"] = []
    result_protocol["heading"] = HEADING_FOR_RESULT_TO_PROTOCOL
    app_name = __package__.split(".")[0]
    # defined_protocols, other_protocol_list = display_available_protocols(app_name)
    defined_protocol_types = display_protocol_types(app_name)

    result_protocol["protocols_dict"] = get_protocol_from_prot_types(
        defined_protocol_types
    )

    result_protocol["protocol_list"] = []
    result_protocol["protocol_type"] = defined_protocol_types
    result_protocol["protocol_filter_selection"] = []
    for key, value in result_protocol["protocols_dict"].items():
        for protocol_value in value:
            if check_if_protocol_parameters_exists(protocol_value):
                result_protocol["protocol_list"].append(protocol_value)
        result_protocol["protocol_filter_selection"].append(
            [key, result_protocol["protocol_list"]]
        )

    result_protocol["table_length"] = len(HEADING_FOR_RESULT_TO_PROTOCOL)
    result_protocol["number_of_samples"] = len(c_samples_pending_protocol)
    for c_sample in c_samples_pending_protocol:
        data = c_sample.get_sample_info_for_protocol()
        for x in range(len(data), len(HEADING_FOR_RESULT_TO_PROTOCOL)):
            data.append("")
        result_protocol["data"].append(data)
        clinic_samples.append(c_sample.get_id())
    result_protocol["clinic_samples"] = ",".join(clinic_samples)
    return result_protocol


def record_result_protocol(request):
    json_data = json.loads(request.POST["table_data"])
    if "samples_in_list" in request.POST:
        c_samples = request.POST["c_samples"]
    else:
        c_samples = request.POST["c_samples"].split(",")
    # invalid_c_samples = []
    invalid_c_samples_ids = []
    invalid = {}
    # valid_c_samples = []
    valid_c_samples_ids = []
    for i in range(len(c_samples)):
        c_sample_obj = get_sample_clinic_obj_from_id(c_samples[i])
        if json_data[i][-1] == "" or not c_sample_obj:
            # invalid_c_samples.append(json_data[i])
            invalid_c_samples_ids.append(c_samples[i])
            continue
        c_sample_obj.set_protocol(json_data[i][-1])
        c_sample_obj.set_state("Pending results")
        valid_c_samples_ids.append(c_samples[i])
    if len(invalid_c_samples_ids) > 0:
        # invalid['invalid_c_samples'] = invalid_c_samples
        invalid["invalid_c_samples_ids"] = ",".join(invalid_c_samples_ids)

    return invalid, valid_c_samples_ids
