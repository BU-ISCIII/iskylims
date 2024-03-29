# Generic imports
import json

# Local imports
import clinic.clinic_config
import clinic.models
import core.utils.common
import core.utils.samples


def analyze_and_store_patient_data(user_post, user):
    stored_samples = []
    analyze_data = {}
    not_match = []
    incomplete_clinic_samples = []
    incomplete_clinic_samples_ids = []
    json_data = json.loads(user_post["patient_data"])
    clinic_samples = user_post["clinic_samples"].split(",")
    heading = clinic.clinic_config.ADDITIONAL_HEADING_FOR_RECORDING_SAMPLES
    for c_samples_id in range(len(clinic_samples)):
        # check if patient history number matches

        if core.utils.samples.check_empty_fields(json_data[c_samples_id]):
            incomplete_clinic_samples.append(json_data[c_samples_id])
            incomplete_clinic_samples_ids.append(clinic_samples[c_samples_id])
            continue

        history_number = json_data[c_samples_id][heading.index("History Number")]

        patient_obj = get_patient_obj(history_number)

        if not patient_obj:
            not_match.append(json_data[c_samples_id][0:6])
            incomplete_clinic_samples.append(json_data[c_samples_id])
            incomplete_clinic_samples_ids.append(clinic_samples[c_samples_id])
            continue
        else:
            patient_data = {}
            patient_data["patient_id"] = patient_obj
            for map_column in clinic.clinic_config.MAP_ADDITIONAL_HEADING_TO_DATABASE:
                patient_data[map_column[1]] = json_data[c_samples_id][
                    heading.index(map_column[0])
                ]
                # match_ids.append(clinic_samples[c_samples_id])
            patient_data["doctor_id"] = get_doctor_obj(
                json_data[c_samples_id][heading.index("Doctor")]
            )
            patient_data["serviceUnit_id"] = get_service_unit_obj(
                (json_data[c_samples_id][heading.index("Requested Service by")])
            )

            clinic_obj = get_clinic_sample_obj_from_id(clinic_samples[c_samples_id])
            clinic_obj.update(patient_data)
            stored_samples.append(
                [
                    json_data[c_samples_id][0],
                    history_number,
                    json_data[c_samples_id][heading.index("Requested Service by")],
                ]
            )
            # create suspicion data for this clinic sample
            suspicion_data = {}
            suspicion_data["patient_id"] = patient_obj
            suspicion_data["clinicSample_id"] = clinic_obj
            suspicion_data["description"] = json_data[c_samples_id][
                heading.index("Suspicion")
            ]
            clinic.models.SuspicionHistory.objects.create_suspicion_history(
                suspicion_data
            )
    if stored_samples:
        analyze_data[
            "stored_samples_heading"
        ] = clinic.clinic_config.HEADING_FOR_STORED_PATIENT_DATA
    if incomplete_clinic_samples_ids:
        analyze_data[
            "heading"
        ] = clinic.clinic_config.ADDITIONAL_HEADING_FOR_RECORDING_SAMPLES
        analyze_data["heading_length"] = len(
            clinic.clinic_config.ADDITIONAL_HEADING_FOR_RECORDING_SAMPLES
        )
        analyze_data["incomplete_clinic_samples"] = incomplete_clinic_samples
        analyze_data["data_length"] = len(incomplete_clinic_samples)
        analyze_data["incomplete_clinic_samples_ids"] = ",".join(
            incomplete_clinic_samples_ids
        )
        # analyze_data['data'] = not_match
        analyze_data["s_request_by"] = get_service_units()
        analyze_data["doctor"] = get_available_doctor()

    if not_match:
        analyze_data["not_match"] = not_match
        analyze_data["heading_not_match"] = clinic.clinic_config.HEADING_FOR_NOT_MATCH
    # analyze_data['not_match_samples_ids'] =not_match_samples_ids
    if stored_samples:
        analyze_data["stored_samples"] = stored_samples

    return analyze_data


def check_if_need_update(c_sample_state):
    c_samples_in_patient = get_clinic_sample_in_state(c_sample_state)

    if c_samples_in_patient:
        for c_sample in c_samples_in_patient:
            sample_core_state = c_sample.get_sample_core_state()
            mapped_sample_state = clinic.clinic_config.MAPPING_SAMPLES_CORE_VS_CLINIC[
                sample_core_state
            ]
            if mapped_sample_state == "Patient update":
                continue
            # Update the clinic sample with new State
            c_sample.set_state(mapped_sample_state)

    return


def check_if_sample_c_exists(sample_c_id):
    if clinic.models.ClinicSampleRequest.objects.filter(pk__exact=sample_c_id).exists():
        return True
    else:
        return False


def collect_sample_data_for_search():
    """
    Description:
        The function will get the doctor names and unit services definded to inclued in
        the option values to display in the select menus.
    Functions:
        get_available_doctor
        get_service_units
    Return:
        search_sample_data.
    """
    search_sample_data = {}
    search_sample_data["doctors"] = get_available_doctor()
    search_sample_data["requested_service_by"] = get_service_units()
    return search_sample_data


def get_clinic_sample_in_state(state):
    if clinic.models.ClinicSampleRequest.objects.filter(
        clinic_sample_state__clinic_state__exact=state
    ).exists():
        return clinic.models.ClinicSampleRequest.objects.filter(
            clinic_sample_state__clinic_state__exact=state
        )
    return


def define_clinic_samples(sample_list, user, state):
    clinic_sample_list = []
    for sample_id in sample_list:
        c_sample_data = {}
        sample_obj = core.utils.samplesget_sample_obj_from_id(sample_id)
        c_sample_data["sampleCore"] = sample_obj
        c_sample_data["patientCore"] = sample_obj.get_sample_patient_obj()
        c_sample_data["user"] = user
        c_sample_data["state"] = "Defined"
        new_clinic_sample = (
            clinic.models.ClinicSampleRequest.objects.create_clinic_sample(
                c_sample_data
            )
        )
        clinic_sample_list.append(new_clinic_sample.get_id())
    return clinic_sample_list


def display_one_sample_info(id):
    sample_info = {}
    clinic_sample_obj = get_clinic_sample_obj_from_id(id)
    core_sample_obj = clinic_sample_obj.get_core_sample_obj()
    p_core_obj = core_sample_obj.get_sample_patient_obj()

    sample_info["s_name"] = core_sample_obj.get_sample_name()
    sample_info["sample_core_info"] = clinic_sample_obj.get_sample_core_info()
    sample_info[
        "sample_core_heading"
    ] = clinic.clinic_config.HEADING_FOR_DISPLAY_SAMPLE_MAIN_INFORMATION

    p_main_info = clinic_sample_obj.get_patient_information()

    # collect patient name and code
    sample_info["patient_basic_info"] = list(
        zip(
            clinic.clinic_config.HEADING_FOR_DISPLAY_PATIENT_BASIC_INFORMATION,
            p_main_info,
        )
    )
    # get patient data
    if clinic.models.PatientData.objects.filter(patien_core__exact=p_core_obj).exists():
        patient_data_obj = clinic.models.PatientData.objects.get(
            patienCore__exact=p_core_obj
        )
        p_data_info = patient_data_obj.get_patient_full_data()
        sample_info["patient_data"] = list(
            zip(
                clinic.clinic_config.HEADING_FOR_DISPLAY_PATIENT_ADDITIONAL_INFORMATION,
                p_data_info,
            )
        )
    sample_info["not_massive"] = core.utils.samples.get_all_sample_information(
        core_sample_obj.get_sample_id(), False
    )
    sample_info["massive"] = core.utils.samples.get_all_sample_information(
        core_sample_obj.get_sample_id(), True
    )
    return sample_info


def display_sample_list(sample_c_list):
    display_sample_list_info = {}
    display_sample_list_info[
        "heading"
    ] = clinic.clinic_config.HEADING_SEARCH_LIST_SAMPLES_CLINIC
    sample_c_data = []
    for sample_c in sample_c_list:
        sample_c_obj = clinic.models.ClinicSampleRequest.objects.get(pk__exact=sample_c)
        sample_c_data.append(sample_c_obj.get_sample_info_for_list())

    display_sample_list_info["sample_c_data"] = sample_c_data
    return display_sample_list_info


def get_clinic_sample_obj_from_id(id):
    clinic_sample_obj = clinic.models.ClinicSampleRequest.objects.get(pk__exact=id)
    return clinic_sample_obj


def get_clinic_samples_by_state(state):
    if clinic.models.ClinicSampleRequest.objects.filter(
        clinic_sample_state__clinic_state__exact=state
    ).exists():
        c_samples = clinic.models.ClinicSampleRequest.objects.filter(
            clinic_sample_state__clinic_state__exact=state
        )
        c_samples_ids = []
        for c_sample in c_samples:
            c_samples_ids.append(c_sample.get_id())
        return c_samples_ids
    else:
        return None


def get_available_doctor():
    doctor_list = []
    if clinic.models.Doctor.objects.all().exists():
        doctors = clinic.models.Doctor.objects.all().order_by("doctorName")
        for doctor in doctors:
            doctor_list.append(doctor.get_name())
    return doctor_list


def get_doctor_obj(doctor_name):
    if clinic.models.Doctor.objects.filter(doctor_name__exact=doctor_name).exists():
        return clinic.models.Doctor.objects.get(doctor_name__exact=doctor_name)
    else:
        return None


def get_patient_obj(history_number):
    if core.models.PatientCore.objects.filter(
        history_number__exact=history_number
    ).exists():
        return core.models.PatientCore.objects.get(history_number__exact=history_number)
    else:
        return None


def get_clinic_samples_in_state(user, state):
    """
    Description:
        The function will return an object list with samples which are in teh variable state.
    Input:
        user  # pending samples limited to the logged user and the friend list
        state #
    Functions:
        get_friend_list # located at core.utils.common
    Return:
        samples or False
    """
    user_friend_list = core.utils.common.get_friend_list(user)
    if clinic.models.ClinicSampleRequest.objects.filter(
        clinic_sample_state__clinic_state__exact=state,
        sample_request_user__in=user_friend_list,
    ).exists():
        samples = clinic.models.ClinicSampleRequest.objects.filter(
            clinic_sample_state__clinic_state__exact=state,
            sample_request_user__in=user_friend_list,
        ).order_by("priority")
        return samples
    else:
        return False


def get_clinic_samples_defined_state(user):
    """
    Description:
        The function will return a list with samples which are in defined state.
    Input:
        user  # pending samples limited to the logged user and the friend list
    Variables:
        sample_type_names # list containing all sample types names
    Functions:
        get_friend_list # located at core.utils.common
    Return:
        samples_in_state.
    """
    sample_information = []
    samples_in_state = {}
    samples_obj = get_clinic_samples_in_state(user, "Defined")
    if samples_obj:
        for sample_obj in samples_obj:
            sample_information.append(sample_obj.get_info_for_defined_state())
        samples_in_state["sample_information"] = sample_information
        samples_in_state[
            "sample_heading"
        ] = clinic.clinic_config.HEADING_FOR_DISPLAY_SAMPLE_DEFINED_STATE
        samples_in_state["length"] = len(sample_information)
        return samples_in_state
    else:
        samples_in_state["length"] = 0
        return samples_in_state


def get_clinic_samples_patient_sequencing_state(user, state):
    """
    Description:
        The function will return a list with samples which are in Patient upadate state.
    Input:
        user  # pending samples limited to the logged user and the friend list
    Variables:
        sample_type_names # list containing all sample types names
    Functions:
        get_clinic_samples_in_state # located at this file
    Return:
        samples_in_state.
    """
    sample_information = []
    samples_in_state = {}
    samples_obj = get_clinic_samples_in_state(user, state)
    if samples_obj:
        for sample_obj in samples_obj:
            sample_information.append(
                sample_obj.get_info_for_patient_sequencing_state()
            )
        samples_in_state["sample_information"] = sample_information
        samples_in_state[
            "sample_heading"
        ] = clinic.clinic_config.HEADING_FOR_DISPLAY_SAMPLE_PATIENT_SEQUENCING_STATE
        samples_in_state["length"] = len(sample_information)
        return samples_in_state
    else:
        samples_in_state["length"] = 0
        return samples_in_state


def get_clinic_sample_obj_from_sample_id(sample_core_id):
    if clinic.models.ClinicSampleRequest.objects.filter(
        sample_core__exact=sample_core_id
    ).exists():
        c_sample_obj = clinic.models.ClinicSampleRequest.objects.get(
            sample_core__exact=sample_core_id
        )
        return c_sample_obj
    else:
        return


def get_clinic_samples_pending_results(user, state):
    """
    Description:
        The function will return a list with samples which are in pending protocol state.
    Input:
        user  # pending samples limited to the logged user and the friend list
    Functions:
        get_friend_list # located at core.utils.common
    Return:
        samples_in_state.
    """
    sample_information = []
    samples_in_state = {}
    samples_obj = get_clinic_samples_in_state(user, state)
    if samples_obj:
        for sample_obj in samples_obj:
            c_sample = sample_obj.get_info_for_patient_sequencing_state()
            c_sample.append(sample_obj.get_id())
            sample_information.append(c_sample)
        samples_in_state["sample_information"] = sample_information
        samples_in_state[
            "sample_heading"
        ] = clinic.clinic_config.HEADING_FOR_DISPLAY_SAMPLE_PENDING_RESULT_STATE
        samples_in_state["length"] = len(sample_information)
        return samples_in_state
    else:
        samples_in_state["length"] = 0
        return samples_in_state


def get_service_unit_obj(unit_name):
    if clinic.models.ServiceUnits.objects.filter(
        service_unit_name__exact=unit_name
    ).exists():
        return clinic.models.ServiceUnits.objects.get(
            service_unit_name__exact=unit_name
        )
    else:
        return None


def get_service_units():
    service_unit_list = []
    if clinic.models.ServiceUnits.objects.all().exists():
        services = clinic.models.ServiceUnits.objects.all().order_by("serviceUnitName")
        for service in services:
            service_unit_list.append(service.get_name())
    return service_unit_list


def get_samples_clinic_in_search(data_request):
    clinic_s_list = []
    if clinic.models.ClinicSampleRequest.objects.all().exists():
        clinic_s_found = clinic.models.ClinicSampleRequest.objects.all()
    else:
        return clinic_s_list
    if data_request["patientCode"] != "":
        clinic_s_found = clinic_s_found.filter(
            patient_core__patient_code__icontains=data_request["patientCode"]
        )

    if data_request["patientName"] != "":
        clinic_s_found = clinic_s_found.filter(
            patientCore__patientName__icontains=data_request["patientName"]
        )
    if data_request["patientSurname"] != "":
        clinic_s_found = clinic_s_found.filter(
            patient_core__patient_surname__icontains=data_request["patientSurname"]
        )

    if data_request["sample_name"] != "":
        clinic_s_found = clinic_s_found.filter(
            sample_core__sample_name__icontains=data_request["sample_name"]
        )
        if len(clinic_s_found) == 1:
            clinic_s_list.append(clinic_s_found[0].pk)
            return clinic_s_list
    if data_request["doctor"] != "":
        clinic_s_found = clinic_s_found.filter(
            doctor_id__doctor_name__exact=data_request["doctor"]
        )
    if data_request["requestedby"] != "":
        clinic_s_found = clinic_s_found.filter(
            service_unit_id__service_unit_name__exact=data_request["requestedby"]
        )

    if data_request["start_date"] != "" and data_request["end_date"] != "":
        clinic_s_found = clinic_s_found.filter(
            generated_at___range=(data_request["start_date"], data_request["end_date"])
        )

    if data_request["start_date"] != "" and data_request["end_date"] == "":
        clinic_s_found = clinic_s_found.filter(
            generated_at__gte=data_request["start_date"]
        )

    if data_request["start_date"] == "" and data_request["end_date"] != "":
        clinic_s_found = clinic_s_found.filter(run_date__lte=data_request["end_date"])

    for clinic_s in clinic_s_found:
        clinic_s_list.append(clinic_s.pk)

    return clinic_s_list


def prepare_patient_form(clinic_samples_ids):
    patient_info = {}
    patient_info["data"] = []
    patient_info[
        "heading"
    ] = clinic.clinic_config.ADDITIONAL_HEADING_FOR_RECORDING_SAMPLES
    heading_length = len(clinic.clinic_config.ADDITIONAL_HEADING_FOR_RECORDING_SAMPLES)
    for clinic_s_id in clinic_samples_ids:
        clinic_sample_obj = get_clinic_sample_obj_from_id(clinic_s_id)
        data_sample = [""] * heading_length
        data_sample[0] = clinic_sample_obj.get_sample_name()
        patient_info["data"].append(data_sample)
    patient_info["clinic_samples_ids"] = ",".join(clinic_samples_ids)
    patient_info["heading_length"] = heading_length
    patient_info["data_length"] = len(patient_info["data"])
    patient_info["doctor"] = get_available_doctor()
    patient_info["s_request_by"] = get_service_units()
    patient_info["clinic_samples"] = ",".join(clinic_samples_ids)
    return patient_info


def update_clinic_sample_state_from_core_sample_id(sample_list, state):
    for core_sample_id in sample_list:
        core_sample_obj = core.utils.samples.get_sample_obj_from_id(core_sample_id)
        clinic_sample_obj = clinic.models.ClinicSampleRequest.objects.get(
            sampleCore=core_sample_obj
        )
        clinic_sample_obj.set_state(state)
    return
