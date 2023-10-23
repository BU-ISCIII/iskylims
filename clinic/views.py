# Generic imports
import datetime
import json
from django.conf import settings
from django.contrib.auth.decorators import login_required
from django.shortcuts import redirect, render

# Local imports
import clinic.clinic_config
import clinic.models
import clinic.utils.patient
import clinic.utils.projects
import clinic.utils.samples

import core.utils.commercial_kits
import core.utils.patient_projects
import core.utils.protocols
import core.utils.samples
import clinic.utils.common
import core.models


def index(request):
    org_name = clinic.utils.common.get_configuration_from_database("ORGANIZATION_NAME")
    return render(request, "clinic/index.html", {"organization_name": org_name})


@login_required
def add_commercial_kit(request):
    app_name = __package__.split(".")[0]
    # exclude the protocols that have no molecule involved
    defined_protocols = core.utils.protocols.get_defined_protocols(app_name, True)

    commercial_kits_data = core.utils.commercial_kits.get_data_for_commercial_kits(
        "Non_NGS"
    )

    if request.method == "POST" and request.POST["action"] == "addCommercialKit":
        if core.utils.commercial_kits.get_commercial_kit_id(request.POST["kitName"]):
            return render(
                request,
                "clinic/addCommercialKit.html",
                {
                    "defined_protocols": defined_protocols,
                    "invalid_name": request.POST["kitName"],
                },
            )
        new_kit = core.utils.commercial_kits.store_commercial_kit(request.POST)
        new_kit_data = core.utils.commercial_kits.get_commercial_kit_basic_data(new_kit)
        return render(
            request,
            "clinic/addCommercialKit.html",
            {"new_kit_data": new_kit_data},
        )
    else:
        return render(
            request,
            "clinic/addCommercialKit.html",
            {
                "defined_protocols": defined_protocols,
                "commercial_kits_data": commercial_kits_data,
            },
        )


@login_required
def add_user_lot_commercial_kit(request):
    if request.method == "POST" and request.POST["action"] == "addUserLotKit":
        if core.utils.commercial_kits.get_lot_user_commercial_kit_obj(
            request.POST["nickName"]
        ):
            defined_kits = core.utils.commercial_kits.get_defined_commercial_kits()
            return render(
                request,
                "clinic/addUserLotCommercialKit.html",
                {
                    "defined_kits": defined_kits,
                    "invalid_name": request.POST["nickName"],
                },
            )
        new_lot_kit = core.utils.commercial_kits.store_lot_user_commercial_kit(
            request.POST, request.user
        )
        new_lot_kit_data = (
            core.utils.commercial_kits.get_lot_user_commercial_kit_basic_data(
                new_lot_kit
            )
        )
        return render(
            request,
            "clinic/addUserLotCommercialKit.html",
            {"new_lot_kit_data": new_lot_kit_data},
        )
    else:
        defined_kits = core.utils.commercial_kits.get_defined_commercial_kits()
        return render(
            request,
            "clinic/addUserLotCommercialKit.html",
            {"defined_kits": defined_kits},
        )


@login_required
def assign_project(request):
    if request.method == "POST" and request.POST["action"] == "addPatientProject":
        defined_project = {}
        (
            defined_project["fields"],
            defined_project["project_id"],
        ) = clinic.utils.projects.assign_project_patient(request.POST, __package__)
        defined_project["patient_id"] = request.POST["patient_id"]

        return render(
            request,
            "clinic/addPatientProject.html",
            {"defined_project": defined_project},
        )
    elif request.method == "POST" and request.POST["action"] == "defineProjectFields":
        project_fields_added = clinic.utils.projects.add_project_fields(request.POST)
        return render(
            request,
            "clinic/addPatientProject.html",
            {"project_fields_added": project_fields_added},
        )

    return render(request, "clinic/addPatientProject.html")


@login_required
def create_new_patient_project(request):
    defined_projects = core.utils.patient_projects.get_defined_patient_projects(
        __package__
    )
    if request.method == "POST" and request.POST["action"] == "addNewProject":
        new_project = core.utils.patient_projects.create_patient_project(
            request.POST, __package__
        )
        if "ERROR" in new_project:
            return render(
                request,
                "clinic/createNewProject.html",
                {
                    "defined_projects": defined_projects,
                    "error": clinic.clinic_config.ERROR_MESSAGE_FOR_PROJECT_NAME_EXISTS,
                },
            )
        return render(
            request,
            "clinic/createNewProject.html",
            {"new_project": new_project},
        )

    return render(
        request,
        "clinic/createNewProject.html",
        {"defined_projects": defined_projects},
    )


@login_required
def create_protocol(request):
    # get the list of defined protocols
    (
        defined_protocols,
        other_protocol_list,
    ) = core.utils.protocols.display_available_protocols(__package__)
    defined_protocol_types = core.utils.protocols.display_protocol_types(__package__)

    if request.method == "POST" and request.POST["action"] == "addNewProtocol":
        new_protocol = request.POST["newProtocolName"]
        protocol_type = request.POST["protocolType"]
        description = request.POST["description"]

        if core.utils.protocols.check_if_protocol_exists(new_protocol, __package__):
            return render(
                request,
                "clinic/createProtocol.html",
                {"content": ["Protocol Name ", new_protocol, "Already exists."]},
            )
        new_protocol_id = core.utils.protocols.create_new_protocol(
            new_protocol, protocol_type, description, __package__
        )

        return render(
            request,
            "clinic/createProtocol.html",
            {
                "defined_protocols": defined_protocols,
                "defined_protocol_types": defined_protocol_types,
                "new_defined_protocol": new_protocol,
                "new_protocol_id": new_protocol_id,
                "other_protocol_list": other_protocol_list,
            },
        )

    return render(
        request,
        "clinic/createProtocol.html",
        {
            "defined_protocols": defined_protocols,
            "defined_protocol_types": defined_protocol_types,
            "other_protocol_list": other_protocol_list,
        },
    )


@login_required
def create_sample_projects(request):
    # get the information of defined sample Projects
    defined_samples_projects = core.utils.samples.get_info_for_defined_sample_projects(
        __package__
    )

    if request.method == "POST" and request.POST["action"] == "addNewSampleProject":
        sample_project_name = request.POST["sampleProyectName"]
        # description = request.POST['description']

        if core.models.SampleProjects.objects.filter(
            sample_project_name__iexact=sample_project_name,
            apps_name__exact=__package__,
        ).exists():
            error_message = clinic.clinic_config.ERROR_SAMPLE_PROJECT_ALREADY_EXISTS
            return render(
                request,
                "clinic/createSampleProjects.html",
                {
                    "defined_samples_projects": defined_samples_projects,
                    "error_message": error_message,
                },
            )
        new_sample_project_id = core.utils.samples.create_new_sample_project(
            request.POST, __package__
        )
        new_defined_sample_project = sample_project_name
        return render(
            request,
            "clinic/createSampleProjects.html",
            {
                "defined_samples_projects": defined_samples_projects,
                "new_sample_project_id": new_sample_project_id,
                "new_defined_sample_project": new_defined_sample_project,
            },
        )

    return render(
        request,
        "clinic/createSampleProjects.html",
        {"defined_samples_projects": defined_samples_projects},
    )


@login_required
def define_extraction_molecules(request):
    extraction_molecules = {}

    extraction_molecules[
        "extract_molecule"
    ] = core.utils.samples.get_sample_objs_in_state("defined", request.user)
    if request.method == "POST" and request.POST["action"] == "continueWithMolecule":
        # processing the samples selected by user
        pass
        return render(
            request,
            "clinic/defineExtractionMolecules.html",
            {"extraction_molecules": extraction_molecules},
        )
        # return the samples that match user request
    elif request.method == "POST" and request.POST["action"] == "searchSamples":
        pass
        return render(
            request,
            "clinic/defineExtractionMolecules.html",
            {"extraction_molecules": extraction_molecules},
        )
        # group the samples selected by user
    elif request.method == "POST" and request.POST["action"] == "addingMoreSamples":
        pass
        return render(
            request,
            "clinic/defineExtractionMolecules.html",
            {"extraction_molecules": extraction_molecules},
        )
    else:
        pass
    return render(
        request,
        "clinic/defineExtractionMolecules.html",
        {"extraction_molecules": extraction_molecules},
    )


@login_required
def define_new_patient(request):
    if request.method == "POST" and request.POST["action"] == "defineNewPatient":
        defined_patient = clinic.utils.patient.create_new_patient(
            request.POST, __package__
        )
        if "ERROR" in defined_patient:
            patient_definition_data = clinic.utils.patient.fields_for_new_patient(
                __package__
            )
            return render(
                request,
                "clinic/defineNewPatient.html",
                {
                    "patient_definition_data": patient_definition_data,
                    "error": clinic.clinic_config.ERROR_MESSAGE_FOR_PATIENT_CODE_EXISTS,
                },
            )
        return render(
            request,
            "clinic/defineNewPatient.html",
            {"defined_patient": defined_patient},
        )
    elif request.method == "POST" and request.POST["action"] == "defineBatchPatient":
        if "patientExcel" in request.FILES:
            patient_batch_data = clinic.utils.patient.read_batch_patient_file(
                request.FILES["patientExcel"]
            )
            if "ERROR" in patient_batch_data:
                patient_definition_data = clinic.utils.patient.fields_for_new_patient(
                    __package__
                )
                return render(
                    request,
                    "clinic/defineNewPatient.html",
                    {
                        "patient_definition_data": patient_definition_data,
                        "error": patient_batch_data["ERROR"],
                    },
                )
            defined_batch_patient = clinic.utils.patient.store_batch_patient(
                patient_batch_data
            )
        return render(
            request,
            "clinic/defineNewPatient.html",
            {"defined_batch_patient": defined_batch_patient},
        )
    elif (
        request.method == "POST"
        and request.POST["action"] == "addAdditionalInformation"
    ):
        clinic.utils.patient.add_additional_information(request.POST)
        return redirect(
            "display_patient_information", patient_id=request.POST["patient_id"]
        )
        # return render(request, 'clinic/defineNewPatient.html' ,{'defined_patient': defined_patient})
    elif request.method == "POST" and request.POST["action"] == "defineProjectFields":
        project_fields_added = clinic.utils.projects.add_project_fields(request.POST)
        return render(
            request,
            "clinic/defineNewPatient.html",
            {"project_fields_added": project_fields_added},
        )
    else:
        patient_definition_data = clinic.utils.patient.fields_for_new_patient(
            __package__
        )
        return render(
            request,
            "clinic/defineNewPatient.html",
            {"patient_definition_data": patient_definition_data},
        )


@login_required
def define_new_patient_history(request):
    return


@login_required
def record_samples(request):
    """
    Functions :
        analyze_input_samples
        analyze_input_sample_project_fields
        prepare_sample_input_table
        prepare_sample_project_input_table
        analyze_reprocess_data
    """
    # Record new samples
    if request.method == "POST" and request.POST["action"] == "recordsample":
        req_user = request.user.username
        excel_json_data = json.loads(request.POST["table_data"])

        sample_recorded = core.utils.samples.analyze_input_samples(
            req_user, excel_json_data
        )

        # if no samples are in any of the options, displays the inital page
        if (
            "defined_samples" not in sample_recorded
            and "pre_defined_samples" not in sample_recorded
            and "invalid_samples" not in sample_recorded
            and "incomplete_samples" not in sample_recorded
        ):
            sample_information = core.utils.samples.prepare_sample_input_table(
                __package__
            )
            return render(
                request,
                "clinic/defineNewSamples.html",
                {"sample_information": sample_information},
            )

        if "valid_samples" in sample_recorded:
            # create the clinic sample  in Define state
            clinic_sample_list = clinic.utils.samples.define_clinic_samples(
                sample_recorded["valid_samples"], request.user, "Defined"
            )
            """
            for sample_id in sample_recorded['valid_samples_ids']:
                c_sample_data = {}
                sample_obj = get_sample_obj_from_id(sample_id)
                c_sample_data['sampleCore'] = sample_obj
                c_sample_data['patientCore'] = sample_obj.get_sample_patient_obj()
                c_sample_data['user'] = request.user
                c_sample_data['state'] = 'Defined'
                new_clinic_sample = ClinicSampleRequest.objects.create_clinic_sample(c_sample_data)

                clinic_sample_list.append(new_clinic_sample.get_id())
            """
            clinic_samples_ids = ",".join(clinic_sample_list)
            sample_recorded["clinic_samples_ids"] = clinic_samples_ids
        if "incomplete_samples" in sample_recorded:
            sample_recorded.update(
                core.utils.samples.prepare_sample_input_table(__package__)
            )
            sample_recorded["number_of_samples"] = len(
                sample_recorded["incomplete_samples"]
            )

        if "pre_defined_samples_id" in sample_recorded:
            # create the clinic sample  in Define state
            clinic_sample_list = clinic.utils.samples.define_clinic_samples(
                sample_recorded["pre_defined_samples_id"], request.user, "Pre-Defined"
            )

            sample_recorded.update(
                core.utils.samples.prepare_sample_project_input_table(
                    sample_recorded["pre_defined_samples_id"]
                )
            )

        return render(
            request,
            "clinic/defineNewSamples.html",
            {"sample_recorded": sample_recorded},
        )

    # display the form to show the samples in pre-defined state that user requested to complete
    elif (
        request.method == "POST"
        and request.POST["action"] == "select_samples_pre_defined"
    ):
        if "samples_in_list" in request.POST:
            pre_defined_samples_id = request.POST.getlist("samples")
        sample_recorded = core.utils.samples.prepare_sample_project_input_table(
            pre_defined_samples_id
        )

        return render(
            request,
            "wetlab/recordSample.html",
            {"sample_recorded": sample_recorded},
        )

    # Add the additional information related to the project
    elif request.method == "POST" and request.POST["action"] == "sampleprojectdata":
        sample_recorded = core.utils.samples.analyze_input_sample_project_fields(
            request.POST
        )
        clinic_sample_for_update = request.POST["pre_defined_id"]
        clinic.utils.samples.update_clinic_sample_state_from_core_sample_id(
            clinic_sample_for_update, "Defined"
        )
        if request.POST["pending_pre_defined"] != "":
            sample_recorded.update(
                core.utils.samples.prepare_sample_project_input_table(
                    request.POST["pending_pre_defined"]
                )
            )
            return render(
                request,
                "wetlab/recordSample.html",
                {"sample_recorded": sample_recorded},
            )
        else:
            return render(
                request,
                "wetlab/recordSample.html",
                {"sample_recorded": sample_recorded},
            )

    else:
        sample_information = core.utils.samples.prepare_sample_input_table(__package__)
        return render(
            request,
            "clinic/defineNewSamples.html",
            {"sample_information": sample_information},
        )


@login_required
def define_patient_information(request):
    if request.method == "POST" and request.POST["action"] == "continueWithPatient":
        if "samples_in_list" in request.POST:
            patient_information = clinic.utils.samples.prepare_patient_form(
                request.POST.getlist("c_samples")
            )
        else:
            patient_information = clinic.utils.samples.prepare_patient_form(
                request.POST["clinic_samples"].split(",")
            )
        return render(
            request,
            "clinic/definePatientInformation.html",
            {"patient_information": patient_information},
        )
    elif request.method == "POST" and request.POST["action"] == "storePatientInfo":
        updated_information = clinic.utils.samples.analyze_and_store_patient_data(
            request.POST, request.user
        )

        return render(
            request,
            "clinic/definePatientInformation.html",
            {"updated_information": updated_information},
        )
    else:
        clinic_samples = clinic.utils.samples.get_clinic_samples_by_state("Defined")
        if not clinic_samples:
            return render(
                request,
                "clinic/definePatientInformation.html",
                {"no_samples": True},
            )
        else:
            patient_information = clinic.utils.samples.prepare_patient_form(
                clinic_samples
            )

        return render(
            request,
            "clinic/definePatientInformation.html",
            {"patient_information": patient_information},
        )


@login_required
def define_project_fields(request, project_id):
    if request.method == "POST" and request.POST["action"] == "defineProjectFields":
        recorded_project_fields = core.utils.patient_projects.set_project_fields(
            request.POST
        )
        return render(
            request,
            "clinic/defineProjectFields.html",
            {"recorded_project_fields": recorded_project_fields},
        )
    else:
        project_fields = core.utils.patient_projects.define_table_for_project_fields(
            project_id
        )
    return render(
        request,
        "clinic/defineProjectFields.html",
        {"project_fields": project_fields},
    )


@login_required
def define_protocol_parameters(request, protocol_id):
    if (
        request.method == "POST"
        and request.POST["action"] == "define_protocol_parameters"
    ):
        recorded_prot_parameters = core.utils.protocols.set_protocol_parameters(request)
        return render(
            request,
            "clinic/defineProtocolParameters.html",
            {"recorded_prot_parameters": recorded_prot_parameters},
        )
    else:
        if not core.utils.protocols.check_if_protocol_exists(protocol_id, __package__):
            return render(
                request,
                "clinic/error_page.html",
                {
                    "content": [
                        "The requested Protocol does not exist",
                        "Create the protocol name before assigning custom protocol parameters.",
                    ]
                },
            )
        prot_parameters = core.utils.protocols.define_table_for_prot_parameters(
            protocol_id
        )
        return render(
            request,
            "clinic/defineProtocolParameters.html",
            {"prot_parameters": prot_parameters},
        )


@login_required
def define_sample_projects_fields(request, sample_project_id):
    # get the list of defined sample Projects

    if (
        request.method == "POST"
        and request.POST["action"] == "defineSampleProjectFields"
    ):
        sample_project_field_data = core.utils.samples.set_sample_project_fields(
            request.POST
        )

        return render(
            request,
            "clinic/defineSampleProjectFields.html",
            {"sample_project_field_data": sample_project_field_data},
        )

    else:
        if not core.utils.samples.check_if_sample_project_id_exists(sample_project_id):
            return render(
                request,
                "clinic/error_page.html",
                {
                    "content": [
                        "The requested Protocol does not exist",
                        "Create the protocol name before assigning custom protocol parameters.",
                    ]
                },
            )

        sample_project_data = core.utils.samples.define_table_for_sample_project_fields(
            sample_project_id
        )
        return render(
            request,
            "clinic/defineSampleProjectFields.html",
            {"sample_project_data": sample_project_data},
        )


@login_required
def display_patient_information(request, patient_id):
    display_patient_info = clinic.utils.patient.display_one_patient_info(
        patient_id, __package__
    )
    if "ERROR" in display_patient_info:
        return render(
            request,
            "clinic/displayPatientInformation.html",
            {"ERROR": "ERROR"},
        )
    else:
        return render(
            request,
            "clinic/displayPatientInformation.html",
            {"display_patient_info": display_patient_info},
        )
    return


@login_required
def display_patient_project(request, project_id):
    if not core.utils.patient_projects.check_if_project_exists(project_id, __package__):
        return render(
            request,
            "clinic/error_page.html",
            {
                "content": [
                    "The project that you are trying to get ",
                    "DOES NOT exists .",
                ]
            },
        )
    project_data = core.utils.patient_projects.get_all_project_info(project_id)

    return render(request, "clinic/displayProject.html", {"project_data": project_data})


@login_required
def display_protocol(request, protocol_id):
    if not core.utils.protocols.check_if_protocol_exists(protocol_id, __package__):
        return render(
            request,
            "clinic/error_page.html",
            {
                "content": [
                    "The protocol that you are trying to get ",
                    "DOES NOT exists .",
                ]
            },
        )
    protocol_data = core.utils.protocols.get_all_protocol_info(protocol_id)

    return render(
        request,
        "clinic/displayProtocol.html",
        {"protocol_data": protocol_data},
    )


@login_required
def display_result_protocol(request, result_protocol_id):
    if not core.utils.protocols.check_if_protocol_exists(
        result_protocol_id, __package__
    ):
        return render(
            request,
            "clinic/error_page.html",
            {
                "content": [
                    "The result protocol that you are trying to get ",
                    "DOES NOT exists .",
                ]
            },
        )
    result_protocol_data = core.utils.protocols.get_all_protocol_info(
        result_protocol_id
    )

    return render(
        request,
        "clinic/displayResultProtocol.html",
        {"result_protocol_data": result_protocol_data},
    )


@login_required
def display_sample_clinic_info(request, sample_c_id):
    """
    Description:
        The function will get the option values to display in the select menus.
    Functions:
        display_one_sample_info         # located at clinic/utils/managed_samples.py
        collect_sample_data_for_search  # located at clinic/utils/managed_samples.py
    Return:
        patient_definition_data.
    """
    if clinic.utils.samples.check_if_sample_c_exists(sample_c_id):
        display_sample_info = clinic.utils.samples.display_one_sample_info(sample_c_id)

        return render(
            request,
            "clinic/displaySampleClinicInfo.html",
            {"display_sample_info": display_sample_info},
        )
    else:
        search_sample_data = clinic.utils.samples.collect_sample_data_for_search()
        error_message = [
            "The clinic sample that you are trying to get",
            "DOES NOT exists .",
        ]
        return render(
            request,
            "clinic/searchSample.html",
            {"search_sample_data": search_sample_data, "error_message": error_message},
        )


@login_required
def display_sample_project(request, sample_project_id):
    samples_project_data = core.utils.samples.get_info_to_display_sample_project(
        sample_project_id
    )
    if "ERROR" in samples_project_data:
        error_message = samples_project_data["ERROR"]
        return render(request, "clinic/error_page.html", {"content": error_message})

    return render(
        request,
        "clinic/displaySampleProject.html",
        {"samples_project_data": samples_project_data},
    )


def pending_to_update(request):
    if "wetlab" in settings.INSTALLED_APPS:
        for c_sample_state in ["Sequencing", "Patient update"]:
            clinic.utils.samples.check_if_need_update(c_sample_state)
    pending = {}

    # get the samples in defined state
    pending["defined"] = clinic.utils.samples.get_clinic_samples_defined_state(
        request.user
    )

    pending[
        "patient_update"
    ] = clinic.utils.samples.get_clinic_samples_patient_sequencing_state(
        request.user, "Patient update"
    )

    pending[
        "sequencing"
    ] = clinic.utils.samples.get_clinic_samples_patient_sequencing_state(
        request.user, "Sequencing"
    )
    pending[
        "pending_protocol"
    ] = clinic.utils.samples.get_clinic_samples_pending_results(
        request.user, "Pending protocol"
    )
    pending[
        "pending_results"
    ] = clinic.utils.samples.get_clinic_samples_pending_results(
        request.user, "Pending results"
    )

    pending[
        "graphic_pending_samples"
    ] = clinic.utils.samples.pending_clinic_samples_for_grafic(pending).render()

    return render(request, "clinic/pendingToUpdate.html", {"pending": pending})


@login_required
def search_sample(request):
    search_sample_data = clinic.utils.samples.collect_sample_data_for_search()
    if request.method == "POST" and (request.POST["action"] == "searchSample"):
        data_request = {}
        data_request["sampleName"] = request.POST["sampleName"]
        data_request["patientName"] = request.POST["patientName"]
        data_request["patientSurname"] = request.POST["patientSurname"]
        data_request["patientCode"] = request.POST["patientCode"]
        data_request["doctor"] = request.POST["doctor"]
        data_request["requestedby"] = request.POST["requestedby"]
        data_request["start_date"] = request.POST["startdate"]
        data_request["end_date"] = request.POST["enddate"]

        # check that some values are in the request if not return the form
        if len([v for v in data_request.values() if v != ""]) == 0:
            return render(
                request,
                "clinic/searchSample.html",
                {"search_sample_data": search_sample_data},
            )
        # Check for valid date format
        if data_request["start_date"] != "":
            try:
                datetime.datetime.strptime(data_request["start_date"], "%Y-%m-%d")
            except Exception:
                return render(
                    request,
                    "clinic/searchSample.html",
                    {
                        "search_sample_data": search_sample_data,
                        "Error": clinic.clinic_config.ERROR_MESSAGE_FOR_INCORRECT_START_SEARCH_DATE,
                    },
                )
        if data_request["end_date"] != "":
            try:
                datetime.datetime.strptime(data_request["end_date"], "%Y-%m-%d")
            except Exception:
                return render(
                    request,
                    "clinic/searchSample.html",
                    {
                        "search_sample_data": search_sample_data,
                        "Error": clinic.clinic_config.ERROR_MESSAGE_FOR_INCORRECT_END_SEARCH_DATE,
                    },
                )
        # Patient name length must be longer than 5 characters
        if data_request["patientName"] != "" and len(data_request["patientName"]) < 4:
            return render(
                request,
                "clinic/searchSample.html",
                {
                    "search_sample_data": search_sample_data,
                    "Error": clinic.clinic_config.ERROR_MESSAGE_FOR_SORT_PATIENT_NAME,
                },
            )
        sample_c_list = clinic.utils.samples.get_samples_clinic_in_search(data_request)

        if len(sample_c_list) == 0:
            return render(
                request,
                "clinic/searchSample.html",
                {
                    "search_sample_data": search_sample_data,
                    "Error": clinic.clinic_config.ERROR_MESSAGE_FOR_NO_MATCH_IN_SEARCH,
                },
            )
        if len(sample_c_list) == 1:
            display_sample_info = clinic.utils.samples.display_one_sample_info(
                sample_c_list[0]
            )

            return render(
                request,
                "clinic/displaySampleClinicInfo.html",
                {"display_sample_info": display_sample_info},
            )
        else:
            display_sample_list_info = clinic.utils.samples.display_sample_list(
                sample_c_list
            )
            return render(
                request,
                "clinic/displaySampleClinicInfo.html",
                {"display_sample_list_info": display_sample_list_info},
            )

    else:
        return render(
            request,
            "clinic/searchSample.html",
            {"search_sample_data": search_sample_data},
        )


@login_required
def search_patient(request):
    s_patient_data = clinic.utils.patient.from_data_for_search_patient()
    if request.method == "POST" and (request.POST["action"] == "searchPatient"):
        data_request = {}
        data_request["p_name"] = request.POST["patientname"]
        data_request["p_surname"] = request.POST["patientsurname"]
        data_request["p_code"] = request.POST["patientcode"]

        len_search_values = len(set(data_request.values()))

        if len_search_values == 1:
            return render(request, "clinic/searchPatient.html")

        if data_request["p_name"] != "" and len(data_request["p_name"]) < 4:
            return render(
                request,
                "clinic/searchPatient.html",
                {"Error": clinic.clinic_config.ERROR_MESSAGE_FOR_SORT_PATIENT_NAME},
            )
        patient_list = clinic.utils.patient.get_patients_in_search(data_request)

        if len(patient_list) == 0:
            return render(
                request,
                "clinic/searchPatient.html",
                {
                    "Error": clinic.clinic_config.ERROR_MESSAGE_FOR_NO_MATCH_IN_SEARCH,
                    "s_patient_data": s_patient_data,
                },
            )
        if len(patient_list) == 1:
            return redirect("display_patient_information", patient_id=patient_list[0])
        else:
            display_patient_list_info = clinic.utils.patient.display_patient_list(
                patient_list
            )
            return render(
                request,
                "clinic/searchPatient.html",
                {"display_patient_list_info": display_patient_list_info},
            )

    else:
        return render(
            request,
            "clinic/searchPatient.html",
            {"s_patient_data": s_patient_data},
        )


@login_required
def set_molecule_values(request):
    if request.method == "POST" and request.POST["action"] == "continueWithMolecule":
        if request.POST["c_samples"] == "":
            return render(
                request,
                "clinic/error_page.html",
                {"content": ["There was no sample selected "]},
            )
        if "samples_in_list" in request.POST:
            c_samples = request.POST.getlist("c_samples")
        else:
            c_samples = request.POST["c_samples"].split(",")

        molecule_protocol = core.utils.samples.get_table_record_molecule(
            c_samples, __package__
        )
        if "ERROR" in molecule_protocol:
            return render(
                request,
                "clinic/error_page.html",
                {"content": ["There was no valid sample selected "]},
            )

        molecule_protocol["samples"] = ",".join(c_samples)

        return render(
            request,
            "clinic/setMoleculeValues.html",
            {"molecule_protocol": molecule_protocol},
        )

    elif (
        request.method == "POST" and request.POST["action"] == "updateMoleculeProtocol"
    ):
        molecule_recorded = core.utils.samples.record_molecules(request)

        if "heading" not in molecule_recorded:
            samples = request.POST["samples"].split(",")
            molecule_protocol = core.utils.samples.get_table_record_molecule(
                samples, __package__
            )
            molecule_protocol["data"] = molecule_recorded["incomplete_molecules"]
            molecule_protocol["samples"] = ",".join(samples)
            return render(
                request,
                "clinic/setMoleculeValues.html",
                {"molecule_protocol": molecule_protocol},
            )
        else:
            if "incomplete_molecules" in molecule_recorded:
                samples = molecule_recorded["incomplete_molecules_ids"].split(",")
                molecule_recorded.update(
                    core.utils.samples.get_table_record_molecule(samples, __package__)
                )

            return render(
                request,
                "clinic/setMoleculeValues.html",
                {"molecule_recorded": molecule_recorded},
            )

    elif (
        request.method == "POST"
        and request.POST["action"] == "displayMoleculeParameters"
    ):
        if "samples_in_list" in request.POST:
            molecules = request.POST.getlist("molecules")
        else:
            molecules = request.POST["molecules"].split(",")
        show_molecule_parameters = (
            core.utils.samples.display_molecule_protocol_parameters(
                molecules, request.user
            )
        )
        return render(
            request,
            "clinic/setMoleculeValues.html",
            {"show_molecule_parameters": show_molecule_parameters},
        )

    elif request.method == "POST" and request.POST["action"] == "addMoleculeParameters":
        (
            added_molecule_protocol_parameters,
            sample_updated_list,
        ) = core.utils.samples.add_molecule_protocol_parameters(request)
        # Update the clinic sample request state

        for sample_updated in sample_updated_list:
            clinic.utils.samples.get_clinic_sample_obj_from_sample_id(
                sample_updated
            ).set_state("Pending results")

        if "pending" in request.POST:
            molecules = request.POST["pending"].split(",")
            show_molecule_parameters = (
                core.utils.samples.display_molecule_protocol_parameters(
                    molecules, request.user
                )
            )
            return render(
                request,
                "clinic/setMoleculeValues.html",
                {
                    "added_molecule_protocol_parameters": added_molecule_protocol_parameters,
                    "show_molecule_parameters": show_molecule_parameters,
                },
            )
        else:
            return render(
                request,
                "clinic/setMoleculeValues.html",
                {
                    "added_molecule_protocol_parameters": added_molecule_protocol_parameters
                },
            )

    else:
        register_user = request.user.username
        display_list = core.utils.samples.get_defined_samples(register_user)

        return render(
            request,
            "clinic/setMoleculeValues.html",
            {"display_list": display_list},
        )
