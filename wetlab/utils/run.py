# Generic imports
import json
import os
import string
import datetime

# from Bio.Seq import Seq
from django.conf import settings

# Local imports
# import core.models import SequencingConfiguration
import core.utils.commercial_kits
import core.utils.protocols
import core.utils.samples
import django_utils.models
import wetlab.config
import wetlab.models
import wetlab.utils.library
import wetlab.utils.pool
import wetlab.utils.samplesheet


def check_run_already_defined_by_crontab(exp_name, pool_ids):
    """Check if user tries to defined a run that it is already
    Args:
        exp_name (string): name of the experiment name
        pool_id (string): number of the pool id
    Returns:
        dictionay: detailing error message or in "defined" key boolean value.
        True if already defined False if not.
    """
    if not wetlab.models.RunProcess.objects.filter(run_name__iexact=exp_name).exists():
        return {"defined": False}
    # experiment name already exists, check if samples are the same
    # to confirm the match
    run_obj = wetlab.models.RunProcess.objects.filter(run_name__iexact=exp_name).last()
    sample_sheet = run_obj.get_sample_file()
    f_name = os.path.join(settings.MEDIA_ROOT, sample_sheet)
    try:
        with open(f_name, "r") as fh:
            file_data = fh.readlines()
    except FileNotFoundError:
        error_message = str(
            wetlab.config.ERROR_RUN_NAME_BY_CRONTAB_ALREADY_CREATED
            + " but  "
            + wetlab.config.ERROR_SAMPLE_SHEET_NOT_FOUND_WHEN_CREATED_BY_CRONTAB
        )
        return {"ERROR": error_message}
    sample_in_s_sheet = wetlab.utils.samplesheet.get_samples_in_sample_sheet(file_data)
    sample_in_pools = wetlab.utils.pool.get_sample_name_in_pools(pool_ids)
    if len(sample_in_pools) != len(sample_in_s_sheet["samples"]):
        return {"ERROR": wetlab.config.ERROR_EXISTING_RUN_WITH_DIF_SAMPLES_AS_IN_CRON}
    for sample in sample_in_pools:
        if sample not in sample_in_s_sheet["samples"]:
            return {
                "ERROR": wetlab.config.ERROR_EXISTING_RUN_WITH_DIF_SAMPLES_AS_IN_CRON
            }
    return {"defined": True}


def check_valid_data_for_creation_run(form_data, user_obj):
    """
    Description:
        Function checks if polls are compatible and if run name is not defined
    Input:
        form_data    # user form
        user_obj    # user who is requesting the creation run
    Functions:
        check_pools_compatible   # located at this file
        get_library_prep_in_pools    # located at this file
    Return:
        True if all ckecks are ok, ERROR and error message if not
    """
    error = {}
    pool_ids = form_data.getlist("poolID")
    result_compatibility = check_pools_compatible(form_data)

    if "ERROR" in result_compatibility:
        return result_compatibility

    lib_prep_ids = get_library_prep_in_pools(pool_ids)
    if len(lib_prep_ids) == 0:
        error["ERROR"] = wetlab.config.ERROR_POOLS_WITH_NO_LIBRARY

    return "OK"


def create_run_in_pre_recorded_and_get_data_for_confirmation(form_data, user_obj):
    """
    Description:
        Function get the user data and create a new run instance. Function also gets
        the pool information for data confirmation
    Input:
        form_data   # user form
        user_obj    # user who is requesting the creation run
    Functions:
        get_library_preparation_data_in_run     # located at this file
        get_stored_user_sample_sheet            # located at this file
        fetch_reagent_kits_used_in_run          # located at this file
    Return:
        display_sample_information
    """
    display_sample_information = {}
    pool_ids = form_data.getlist("poolID")
    lib_prep_ids = get_library_prep_in_pools(pool_ids)
    try:
        center_requested_id = (
            django_utils.models.Profile.objects.filter(profile_user_id=user_obj)
            .last()
            .profile_center.id
        )
        center_requested_by = django_utils.models.Center.objects.get(
            pk__exact=center_requested_id
        )
    except Exception:
        center_requested_by = None
    reagent_kit_objs = fetch_reagent_kits_used_in_run(form_data)

    display_sample_information = get_library_preparation_data_in_run(
        lib_prep_ids, pool_ids
    )
    display_sample_information.update(get_stored_user_sample_sheet(lib_prep_ids))
    # update Reagents kits
    new_run_obj = wetlab.models.RunProcess(
        run_name=form_data["experimentName"],
        sample_sheet="",
        state=wetlab.models.RunStates.objects.get(run_state_name__exact="pre_recorded"),
        center_requested_by=center_requested_by,
    )
    new_run_obj.save()
    for reagent_kit_obj in reagent_kit_objs:
        new_run_obj.reagent_kit.add(reagent_kit_obj)

    for pool in pool_ids:
        # changed from version 3.1.0 the relation in pools
        new_run_obj.set_pool(pool)
    display_sample_information["experiment_name"] = form_data["experimentName"]
    display_sample_information["run_process_id"] = new_run_obj.get_run_id()
    return display_sample_information


def collect_data_and_update_library_preparation_samples_for_run(data_form, user):
    """
    Description:
        Function collect the information in the form and update the library preparation
        data with the value confirmed bu user.
        Return a dictionary with the user form
    Input:
        data_form    # user data form
    Constants:
        HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_MISEQ_SINGLE_READ
        HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_MISEQ_PAIRED_END
        HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_NEXTSEQ_SINGLE_READ
        HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_NEXTSEQ_PAIRED_END
        MAP_USER_SAMPLE_SHEET_TO_DATABASE_MISEQ_SINGLE_READ
        MAP_USER_SAMPLE_SHEET_TO_DATABASE_MISEQ_PAIRED_END
        MAP_USER_SAMPLE_SHEET_TO_DATABASE_NEXTSEQ_SINGLE_READ
        MAP_USER_SAMPLE_SHEET_TO_DATABASE_NEXTSEQ_PAIRED_END
    Return:
        record_data with the information collected from the user form
    """
    record_data = {}
    lib_prep_ids = data_form["lib_prep_ids"].split(",")
    lib_prep_unique_ids = data_form["lib_prep_unique_ids"].split(",")
    json_data = json.loads(data_form["s_sheet_data"])
    record_data["single_read"] = data_form["single_read"]
    projects = []
    record_data["version"] = get_iem_version_from_user_sample_sheet(lib_prep_ids)

    if data_form["platform_type"] == "MiSeq":
        record_data["platform"] = "MiSeq"
        if record_data["single_read"] == "TRUE":
            if record_data["version"] == "4":
                heading = (
                    wetlab.config.HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_MISEQ_SINGLE_READ_VERSION_4.copy()
                )
                mapping_fields = (
                    wetlab.config.MAP_USER_SAMPLE_SHEET_TO_DATABASE_MISEQ_SINGLE_READ_VERSION_4
                )
            else:
                heading = (
                    wetlab.config.HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_MISEQ_SINGLE_READ_VERSION_5.copy()
                )
                mapping_fields = (
                    wetlab.config.MAP_USER_SAMPLE_SHEET_TO_DATABASE_MISEQ_SINGLE_READ_VERSION_5
                )
        else:
            if record_data["version"] == "4":
                heading = (
                    wetlab.config.HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_MISEQ_PAIRED_END_VERSION_4.copy()
                )
                mapping_fields = (
                    wetlab.config.MAP_USER_SAMPLE_SHEET_TO_DATABASE_MISEQ_PAiRED_END_VERSION_4
                )
            else:
                heading = (
                    wetlab.config.HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_MISEQ_PAIRED_END_VERSION_5.copy()
                )
                mapping_fields = (
                    wetlab.config.MAP_USER_SAMPLE_SHEET_TO_DATABASE_MISEQ_PAiRED_END_VERSION_5
                )
    else:
        record_data["platform"] = "NextSeq"
        if record_data["single_read"] == "TRUE":
            heading = (
                wetlab.config.HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_NEXTSEQ_SINGLE_READ
            )
            mapping_fields = (
                wetlab.config.MAP_USER_SAMPLE_SHEET_TO_DATABASE_NEXTSEQ_SINGLE_READ
            )
        else:
            heading = (
                wetlab.config.HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_NEXTSEQ_PAIRED_END.copy()
            )
            mapping_fields = (
                wetlab.config.MAP_USER_SAMPLE_SHEET_TO_DATABASE_NEXTSEQ_PAIRED_END
            )
    sample_sheet_data_field = []
    # Store the confirmation information form user form
    for row_index in range(len(json_data)):
        confirmation_data = {}
        right_id = lib_prep_ids[lib_prep_unique_ids.index(json_data[row_index][0])]
        lib_prep_obj = wetlab.utils.library.get_lib_prep_obj_from_id(right_id)

        for mapping_index in range(len(mapping_fields)):
            confirmation_data[mapping_fields[mapping_index][1]] = json_data[row_index][
                heading.index(mapping_fields[mapping_index][0])
            ]
        projects.append(json_data[row_index][heading.index("Sample_Project")])
        # keep the old value of the user sample sheet and update the library prepation
        confirmation_data["user_sample_sheet"] = (
            lib_prep_obj.get_user_sample_sheet_obj()
        )
        lib_prep_obj.update_library_preparation_with_indexes(confirmation_data)
        sample_sheet_data_field.append(json_data[row_index])
    if record_data["platform"] == "NextSeq":
        record_data["index_well"] = False
        requires_remove_index_plate_well = False
        if "Index_Plate_Well" in heading:
            index = heading.index("Index_Plate_Well")
            record_data["index_well"] = True
            requires_remove_index_plate_well = True
            for data_line in sample_sheet_data_field:
                if data_line[index] != "":
                    requires_remove_index_plate_well = False
                    record_data["index_well"] = True
                    break
        if requires_remove_index_plate_well:
            heading.remove("Index_Plate_Well")
            for data_line in sample_sheet_data_field:
                del data_line[index]
    record_data["sample_sheet_data_field"] = sample_sheet_data_field
    record_data["sample_sheet_data_heading"] = heading
    record_data["run_obj"] = get_run_obj_from_id(data_form["run_process_id"])
    # get data for sample Sheet

    record_data["investigator"] = user.username
    record_data["application"] = data_form["application"]
    if "instrument" in data_form:
        record_data["instrument"] = data_form["instrument"]
    else:
        record_data["instrument"] = ""
    record_data["assay"] = data_form["assay"]
    if "collection_index" in data_form:
        record_data["collection_index"] = data_form["collection_index"]
    else:
        record_data["collection_index"] = ""
    record_data["reads"] = data_form["reads"]
    record_data["adapter"] = data_form["adapter"]
    if "adapter2" in data_form:
        record_data["adapter2"] = data_form["adapter2"]

    record_data["exp_name"] = record_data["run_obj"].get_run_name()
    record_data["lib_prep_ids"] = lib_prep_ids
    record_data["projects"] = list(set(projects))

    return record_data


def create_new_projects_added_to_run(project_list, run_obj, user_obj):
    """
    Description:
        The function create the project instance and join them to the run obj
    Input:
        project_list    # List of project to create
        run_obj     # run object
        user_obj    # user object
    Ouput:
        project_objs
    """
    projects_objs = []
    duplicated_projects = []
    if wetlab.models.ConfigSetting.objects.filter(
        configuration_name__exact="PROJECTS_ALLOWED_IN_MULTIPLE_RUNS"
    ).exists():
        allow_projects_in_multi_run = (
            wetlab.models.ConfigSetting.objects.filter(
                configuration_name__exact="PROJECTS_ALLOWED_IN_MULTIPLE_RUNS"
            )
            .last()
            .get_configuration_value()
        )
    else:
        allow_projects_in_multi_run = "TRUE"
    for project in project_list:
        if wetlab.models.Projects.objects.filter(project_name__iexact=project).exists():
            if allow_projects_in_multi_run == "FALSE":
                duplicated_projects.append(project)
                continue
            project_obj = wetlab.models.Projects.objects.filter(
                project_name__iexact=project
            ).last()
            projects_objs.append(project_obj)
        else:
            project_data = {}
            project_data["user_id"] = user_obj
            project_data["projectName"] = project
            project_obj = wetlab.models.Projects.objects.create_new_empty_project(
                project_data
            )
            projects_objs.append(project_obj)
    if len(duplicated_projects) > 0:
        # delete defined projects
        for project_obj in projects_objs:
            project_obj.delete()
        error_message = wetlab.config.ERROR_NOT_ALLOWED_REPEATED_PROJECTS
        error_message += ". Check: " + ", ".join(duplicated_projects)
        return {"ERROR": error_message}

    # join project with run
    for project_obj in projects_objs:
        project_obj.run_process.add(run_obj)

    return projects_objs


def get_pool_objs_from_ids(pool_ids):
    """
    Description:
        The function get the pool instance from the list of pool ids, and return a lisf of pool objs
    Return:
        pool_obj
    """
    pool_objs = []
    for pool in pool_ids:
        pool_objs.append(wetlab.models.LibraryPool.objects.get(pk__exact=pool))
    return pool_objs


def get_pool_adapters(pool_objs):
    """
    Description:
        The function get the adapters used for each pool.
        return a dictionary with adapter as a key and the pool name list as value
    Return:
        adapters
    """
    adapters = {}
    for pool in pool_objs:
        adapter = pool.get_adapter()
        if adapter not in adapters:
            adapters[adapter] = []
        adapters[adapter].append(pool.get_pool_name())

    return adapters


def get_pool_duplicated_index(pool_objs):
    """
    Description:
        The function get the single read and paired ened  used for each pool.
        return a dictionary with single_paired as a key and the pool name list as value
    Functions:
        check_if_duplicated_index   # located at wetlab/utils/pool_preparation.py
    Return:
        False if no duplicated found or the result_index cotaining the duplicated samples index
    """
    library_preparation_ids = []
    for pool in pool_objs:
        lib_prep_objs = wetlab.models.LibPrepare.objects.filter(pools=pool)
        for lib_prep_obj in lib_prep_objs:
            library_preparation_ids.append(lib_prep_obj.get_id())

        result_index = wetlab.utils.pool.check_if_duplicated_index(
            library_preparation_ids
        )
    if "True" in result_index:
        return "False"
    else:
        return result_index


def check_pools_compatible(data_form):
    """
    Description:
        The function in case that more than one pools are used for
        a new run it will check if, single_paired, adapters are the same.
        It checks that there are not duplicate_index.
    Functions:
        get_pool_objs_from_ids       # located at this file
        get_pool_adapters            # located at this file

        get_pool_duplicated_index   # located at this file
    Return:
        True if all cheks are ok, or error message
    """
    error = {}
    pool_ids = data_form.getlist("poolID")
    if len(pool_ids) == 1:
        return "True"
    pool_objs = get_pool_objs_from_ids(pool_ids)
    # get adapters used in the pools
    adapters = get_pool_adapters(pool_objs)
    if len(adapters) > 1:
        error_message = wetlab.config.ERROR_DIFFERENT_ADAPTERS_USED_IN_POOL.copy()
        error["ERROR"] = error_message.insert(1, ",".join(adapters))
        return error

    duplicated_index = get_pool_duplicated_index(pool_objs)
    if "False" not in duplicated_index:
        error_message = (
            wetlab.config.ERROR_DUPLICATED_INDEXES_FOUND_IN_DIFFERENT_POOLS.copy()
        )
        for duplicated in duplicated_index["incompatible_index"]:
            error_message.insert(1, ",".join(duplicated))
        error["ERROR"] = error_message
        return error
    return "True"


def store_confirmation_sample_sheet(fields):
    """
    Description:
        The function store the sample sheet for the run, using
        the template and returning the file name including the relative path
    Input:
        fields          # dictionary having all information for creating the sample sheet
    Constants:
        SAMPLE_SHEET
        RUN_TEMP_DIRECTORY
        RUN_SAMPLE_SHEET_DIRECTORY
        MEDIA_ROOT
        TEMPLATE_FILES_DIRECTORY
        SAMPLE_SHEET_TWO_ADAPTERS_TEMPLATE_NAME
        SAMPLE_SHEET_ONE_ADAPTER_TEMPLATE_NAME
    Return:
        sample sheet name including the relative path
    """
    exp_name_in_file = fields["exp_name"].replace(" ", "̣_")
    today_date = datetime.datetime.today().strftime("%Y%m%d_%H%M%S")
    file_name = str(
        exp_name_in_file + "_" + today_date + "_" + wetlab.config.SAMPLE_SHEET
    )
    tmp_file_relative_path = os.path.join(wetlab.config.RUN_TEMP_DIRECTORY, file_name)
    ss_file_relative_path = os.path.join(
        wetlab.config.RUN_SAMPLE_SHEET_DIRECTORY, file_name
    )
    ss_file_full_path = os.path.join(settings.MEDIA_ROOT, tmp_file_relative_path)

    today_date = datetime.datetime.today().strftime("%d/%m/%Y")

    if "adapter2" in fields:
        fields["adapter"] = str(
            fields["adapter"] + "\nAdapterRead2," + fields["adapter2"]
        )
    d = {
        "investigator": fields["investigator"],
        "exp_name": fields["exp_name"],
        "date": today_date,
        "application": fields["application"],
        "instrument": fields["instrument"],
        "assay": fields["assay"],
        "collection_index": fields["collection_index"],
        "reads": fields["reads"],
        "adapter": fields["adapter"],
    }

    if fields["platform"] == "MiSeq":
        if fields["version"] == "4":
            template_file = os.path.join(
                settings.MEDIA_ROOT,
                wetlab.config.TEMPLATE_FILES_DIRECTORY,
                wetlab.config.SAMPLE_SHEET_MISEQ_VERSION_4_TEMPLATE_NAME,
            )
        else:
            template_file = os.path.join(
                settings.MEDIA_ROOT,
                wetlab.config.TEMPLATE_FILES_DIRECTORY,
                wetlab.config.SAMPLE_SHEET_MISEQ_VERSION_5_TEMPLATE_NAME,
            )
    else:
        template_file = os.path.join(
            settings.MEDIA_ROOT,
            wetlab.config.TEMPLATE_FILES_DIRECTORY,
            wetlab.config.SAMPLE_SHEET_NEXTSEQ_VERSION_5_TEMPLATE_NAME,
        )
    with open(template_file, "r") as filein:
        # filein = open(template_file, 'r')
        ss_template = string.Template(filein.read())

    updated_info = ss_template.substitute(d)
    fh = open(ss_file_full_path, "w")
    fh.write(updated_info)
    fh.write(",".join(fields["sample_sheet_data_heading"]) + "\n")
    for sample in fields["sample_sheet_data_field"]:
        fh.write(",".join(sample) + "\n")
    fh.close()
    # store sample sheet in database
    run_obj = wetlab.models.RunProcess.objects.get(run_name=fields["run_obj"])
    run_obj.update_sample_sheet(ss_file_full_path, file_name)

    os.remove(ss_file_full_path)
    return ss_file_relative_path


def get_experiment_name(run_id):
    """
    Description:
        The function returns the experiment name
    Input:
        run_id # run process id
    Return:
        experiment_name
    """
    run_obj = wetlab.models.RunProcess.objects.get(pk__exact=run_id)
    experiment_name = run_obj.get_run_name()
    return experiment_name


def get_library_preparation_data_in_run(lib_prep_ids, pool_ids):
    """
    Description:
        The function returns the information for the library preparation
    Input:
        lib_prep_ids    # list having the library preparation id
        pool_ids        # pool id list
    Constant:
        HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_PAIREDEND
        HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_SINGLEREAD
    Function:
        get_type_read_sequencing             # located at this file
        collect_lib_prep_data_for_new_run    # located at this file
    Return:
        display_sample_information
    """
    display_sample_information = {}

    pool_objs = get_pool_objs_from_ids(pool_ids)
    sequencers_platform = []
    for pool_obj in pool_objs:
        sequencers_platform.append(pool_obj.get_platform_name())
    if "MiSeq" in sequencers_platform:
        platform_in_pool = "MiSeq"
    else:
        platform_in_pool = "NextSeq"
    # single_paired = get_type_read_sequencing(pool_ids)
    lib_prep_data = collect_lib_prep_data_for_new_run(lib_prep_ids, platform_in_pool)
    display_sample_information["data"] = lib_prep_data["data"]
    display_sample_information["heading"] = lib_prep_data["heading"]
    display_sample_information["lib_prep_ids"] = ",".join(lib_prep_ids)
    display_sample_information["lib_prep_unique_ids"] = ",".join(
        lib_prep_data["uniqueID_list"]
    )

    display_sample_information["date"] = datetime.datetime.today().strftime("%Y%m%d")
    display_sample_information["single_read"] = lib_prep_data["single_read"]
    display_sample_information["platform_type"] = platform_in_pool
    return display_sample_information


def get_iem_version_from_user_sample_sheet(lib_prep_id):
    """
    Description:
        The function returns the latest IEM version used for the requested library
        preparation
    Input:
        lib_prep_id   # id of the library preparation
    Return:
        iem_version
    """
    iem_versions = wetlab.utils.library.get_iem_version_for_library_prep_ids(
        lib_prep_id
    )
    if len(iem_versions) == 1:
        iem_version = iem_versions[0]
    else:
        # multiple version used in IEM, select the latest version
        latest_version = 0
        for version in iem_versions:
            version = int(version)
            if version > latest_version:
                latest_version = version
        iem_version = str(latest_version)
    return iem_version


def get_stored_user_sample_sheet(lib_prep_ids):
    """
    Description:
        The function returns the stored values of the user sample sheet
    Input:
        lib_prep_ids   # id of the library preparation
    Return:
        sample_sheet_data
    """
    sample_sheet_data = {}
    iem_version = get_iem_version_from_user_sample_sheet(lib_prep_ids)
    if iem_version == "0":
        lib_prep_obj = wetlab.utils.library.get_lib_prep_obj_from_id(lib_prep_ids[0])
    else:
        for lib_prep_id in lib_prep_ids:
            lib_prep_obj = wetlab.utils.library.get_lib_prep_obj_from_id(lib_prep_id)
            if lib_prep_obj.get_iem_version() == iem_version:
                break
    u_sample_sheet_obj = lib_prep_obj.get_user_sample_sheet_obj()
    fields = [
        "collection_index",
        "application",
        "instrument",
        "adapter1",
        "adapter2",
        "assay",
        "reads",
    ]
    not_index_fields_in_version_4 = [0, 2]
    data = u_sample_sheet_obj.get_all_data()

    for i in range(len(fields)):
        if iem_version == "4" and i in not_index_fields_in_version_4:
            continue
        else:
            sample_sheet_data[fields[i]] = data[i]
    return sample_sheet_data


def collect_lib_prep_data_for_new_run(lib_prep_ids, platform_in_pool):
    """
    Description:
        The function returns the library preparation data for each one in the lib_prep_ids.
        Sample_uniqueID is modified by adding '-' and the number of reused for the library preparation
    Input:
        lib_prep_ids        # list of library preparations
        platform_in_pool    # platfrom used to add new fields if needed
    Return:
        data
    """
    lib_data = {}
    data = []
    single_read = True
    uniqueID_list = []
    iem_version = get_iem_version_from_user_sample_sheet(lib_prep_ids)
    for lib_prep_id in lib_prep_ids:
        lib_prep_obj = wetlab.models.LibPrepare.objects.get(pk__exact=lib_prep_id)
        # get index 7 and index 5. then index 5 is deleted if all samples are single end
        row_data = lib_prep_obj.get_info_for_run_paired_end()
        if single_read:
            if lib_prep_obj.get_i5_index() != "":
                single_read = False
        if platform_in_pool == "MiSeq":
            if iem_version == "5":
                row_data.insert(8, lib_prep_obj.get_manifest())
                row_data.insert(9, lib_prep_obj.get_genome_folder())
            else:
                row_data.insert(8, lib_prep_obj.get_genome_folder())
        else:
            row_data.insert(4, lib_prep_obj.get_index_plate_well())
        row_data[0] = row_data[0] + "-" + lib_prep_obj.get_reused_value()
        uniqueID_list.append(row_data[0])
        data.append(row_data)
    lib_data["data"] = data
    lib_data["uniqueID_list"] = uniqueID_list
    lib_data["single_read"] = "TRUE" if single_read else "FALSE"

    # if there is not index 5 in any of the samples. Then delete the I5 index colunm for all samples
    if single_read:
        for item in data:
            if platform_in_pool == "MiSeq":
                del item[7]
                del item[6]
            else:
                # sample well was inserted in position 3, so delete index for i5 are stepped in one
                del item[8]
                del item[7]
    if platform_in_pool == "MiSeq":
        if single_read:
            if iem_version == "4":
                lib_data["heading"] = (
                    wetlab.config.HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_MISEQ_SINGLE_READ_VERSION_4
                )
            else:
                lib_data["heading"] = (
                    wetlab.config.HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_MISEQ_SINGLE_READ_VERSION_5
                )
        else:
            if iem_version == "4":
                lib_data["heading"] = (
                    wetlab.config.HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_MISEQ_PAIRED_END_VERSION_4
                )
            else:
                lib_data["heading"] = (
                    wetlab.config.HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_MISEQ_PAIRED_END_VERSION_5
                )
    else:
        if single_read:
            lib_data["heading"] = (
                wetlab.config.HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_NEXTSEQ_SINGLE_READ
            )
        else:
            lib_data["heading"] = (
                wetlab.config.HEADING_FOR_COLLECT_INFO_FOR_SAMPLE_SHEET_NEXTSEQ_PAIRED_END
            )
    return lib_data


def increase_reuse_if_samples_exists(sample_list):
    """
    Description:
        The function check if sample was exists (only this function is requested when
        run is reused, by uploading the sample sheet
    Input:
        sample_list        # list of samples fetched from sample sheet
    Return:
        samples_reused
    """
    samples_reused = []
    for sample in sample_list:
        # PROBLEM GET SAMPLE OBJ FROM SAMPLE NAME IS INCORRECT AS SAMPLE NAMES CAN BE REPEATED
        sample_obj = core.utils.samples.get_sample_obj_from_sample_name(sample)
        if sample_obj:
            core.utils.samples.update_sample_reused(sample_obj.get_sample_id())
            molecules = core.utils.samples.get_molecule_objs_from_sample(sample_obj)
            last_molecule_obj = molecules.reverse()[0]
            core.utils.samples.update_molecule_reused(
                sample_obj.get_sample_id(), last_molecule_obj.get_molecule_code_id()
            )
            wetlab.utils.library.update_library_preparation_for_reuse(samples_reused)
            samples_reused.append(sample)

    return samples_reused


def get_library_prep_in_pools(pool_ids):
    """
    Description:
        The function get the pool id list and returns the the ids for the LibraryPreparation
    Return:
        lib_prep_ids
    """
    lib_prep_ids = []

    for pool_id in pool_ids:
        if wetlab.models.LibPrepare.objects.filter(pools__exact=pool_id).exists():
            lib_prep_objs = wetlab.models.LibPrepare.objects.filter(
                pools__exact=pool_id
            ).order_by("register_user")
            for lib_prep_obj in lib_prep_objs:
                lib_prep_ids.append(lib_prep_obj.get_id())
    return lib_prep_ids


def get_available_pools_for_run():
    """
    Description:
        The function get the pool which are not used in a run. It splits in 2 keys. Pool which are not assigned
        yet to a run and the pools that are associated but there is information missing to be filled.
    Return:
        pools_to_update
    """
    pools_to_update = {}

    # get the pools that were selected
    if wetlab.models.LibraryPool.objects.filter(
        pool_state__pool_state__exact="Selected"
    ).exists():
        pool_objs = wetlab.models.LibraryPool.objects.filter(
            pool_state__pool_state__exact="Selected"
        ).order_by("platform")
        pools_to_update["pools_available"] = {}
        for pool_obj in pool_objs:
            platform = pool_obj.get_platform_name()
            if platform not in pools_to_update["pools_available"]:
                pools_to_update["pools_available"][platform] = []
            pools_to_update["pools_available"][platform].append(pool_obj)
    # get the pools that are associated to a run but not yet completed
    if (
        wetlab.models.LibraryPool.objects.filter(
            pool_state__pool_state__exact="Selected"
        )
        .exclude(runprocess=None)
        .exists()
    ):
        pools_to_update["defined_runs"] = (
            wetlab.models.LibraryPool.objects.filter(
                pool_state__pool_state__exact="Selected"
            )
            .exclude(runprocess=None)
            .order_by("run_process_id")
        )

    return pools_to_update


def get_pool_info(pools_to_update):
    """
    Description:
        The function get the information for pool which are not used in a run.
        And information for pools that are only defined the run name
    Input:
        pools_to_update     # contains the pool objects split in pools_available and defined_runs
                        pools_available is a dictionary which contains as key the
                        platform and value a list of the pool objects
    Constant:
        HEADING_FOR_SELECTING_POOLS
        HEADING_FOR_INCOMPLETED_SELECTION_POOLS
    Functions:
        get_lot_reagent_commercial_kits # located at core/utils/handling_commercial_kits
    Return:
        pool_info
    """
    pool_info = {}
    if "pools_available" in pools_to_update:
        pool_data = {}
        pool_data["heading"] = wetlab.config.HEADING_FOR_SELECTING_POOLS
        pool_data["platform"] = {}
        reagents_kits = {}
        commercial_kits = {}
        for platform, protocol_objs in pools_to_update["pools_available"].items():
            pool_data["platform"][platform] = []
            for pool in protocol_objs:
                data = pool.get_info()
                # compare the number of samples to check that no samples are deleted
                if int(pool.get_number_of_samples()) == len(
                    wetlab.models.LibPrepare.objects.filter(pools=pool)
                ):
                    data.append(pool.get_id())
                    pool_data["platform"][platform].append(data)

                    # get the reagents kits used for the platform
                    (
                        reagents_kits[platform],
                        commercial_kits[platform],
                    ) = core.utils.commercial_kits.get_lot_reagent_commercial_kits(
                        platform
                    )
                    configutation_name = (
                        wetlab.models.LibPrepare.objects.filter(pools=pool)
                        .last()
                        .get_user_sample_sheet_obj()
                        .get_sequencing_configuration_name()
                    )
                    user_lot_configuration_kit = (
                        core.utils.commercial_kits.get_lot_reagent_from_comercial_kit(
                            configutation_name
                        )
                    )
                    if len(user_lot_configuration_kit) > 0:
                        reagents_kits[platform].append(user_lot_configuration_kit)
                        commercial_kits[platform] = str(
                            commercial_kits[platform] + "," + configutation_name
                        )
                else:
                    if "invalid_run_data" not in pool_info:
                        pool_info["invalid_run_data"] = {}
                        pool_info["invalid_run_data"]["data"] = []
                    pool_info["invalid_run_data"]["data"].append(data)

        pool_info["pool_data"] = pool_data
        pool_info["reagents_kits"] = reagents_kits
        pool_info["commercial_kits"] = commercial_kits
    if "defined_runs" in pools_to_update:
        run_data = {}
        tmp_data = {}
        # run_data['r_name'] = {}
        for pool in pools_to_update["defined_runs"]:
            run_name = pool.get_run_name()
            if int(pool.get_number_of_samples()) == len(
                wetlab.models.LibPrepare.objects.filter(pools=pool)
            ):
                if run_name not in tmp_data:
                    tmp_data[run_name] = {}
                    tmp_data[run_name]["data"] = []
                tmp_data[run_name]["run_id"] = pool.get_run_id()
                pool_name = pool.get_pool_name()
                pool_code = pool.get_pool_code_id()
                pool_numbers = pool.get_number_of_samples()
                tmp_data[run_name]["data"].append([pool_name, pool_code, pool_numbers])
            else:
                data = pool.get_info()
                data.append(pool.get_id())

                if "invalid_run_data" not in pool_info:
                    pool_info["invalid_run_data"] = {}
                    pool_info["invalid_run_data"]["data"] = []
                pool_info["invalid_run_data"]["data"].append(data)
        run_info_data = []
        for r_name, values in tmp_data.items():
            run_info_data.append([r_name, values["data"], values["run_id"]])
        if len(run_info_data) > 0:
            run_data["heading"] = wetlab.config.HEADING_FOR_INCOMPLETED_SELECTION_POOLS
            run_data["data"] = run_info_data
            pool_info["run_data"] = run_data
        if "invalid_run_data" in pool_info:
            pool_info["invalid_run_data"][
                "heading"
            ] = wetlab.config.HEADING_FOR_INCOMPLETED_SELECTION_POOLS
    return pool_info


def get_run_obj_from_id(run_id):
    """
    Description:
        The function return the runprocess object for the run_id
    Input:
        run_id     # run_process id
    Return:
        run_obj
    """
    run_obj = wetlab.models.RunProcess.objects.get(pk__exact=run_id)
    return run_obj


def get_run_user_lot_kit_used_in_sample(sample_id):
    """
    Description:
        The function return the runprocess object for the run_id
    Input:
        sample_id   # sample id
    Constamt:
        HEADING_FOR_DISPLAY_ADDITIONAL_KIT_LIBRARY_PREPARATION
    Return:
        kit_data
    """
    kit_data = {}
    kit_data["run_kits_from_sample"] = {}
    if wetlab.models.LibPrepare.objects.filter(sample_id__pk__exact=sample_id).exists():
        kit_data["heading_run_kits"] = (
            wetlab.config.HEADING_FOR_DISPLAY_KIT_IN_RUN_PREPARATION
        )
        library_preparation_items = wetlab.models.LibPrepare.objects.filter(
            sample_id__pk__exact=sample_id
        ).order_by("protocol_id")
        for lib_prep in library_preparation_items:
            pool_objs = lib_prep.pools.all()

            for pool_obj in pool_objs:
                run_obj = pool_obj.get_run_obj()
                if run_obj is not None:
                    run_name = run_obj.get_run_name()
                    kit_data["run_kits_from_sample"][run_name] = []
                    run_kit_objs = run_obj.reagent_kit.all()
                    for run_kit_obj in run_kit_objs:
                        data = [lib_prep.get_lib_prep_code()]
                        data.append(run_kit_obj.get_lot_number())
                        data.append(run_kit_obj.get_commercial_kit())
                        data.append(run_kit_obj.get_expiration_date())
                        data.append(run_obj.get_run_generated_date())

                        kit_data["run_kits_from_sample"][run_name].append(data)

    return kit_data


def display_available_pools():
    """
    Description:
        The function call 2 functions:
        - get_available_pools_for_run to get pool instances that are available
        - get_pool_info for adding the information to display
    Functions:
        get_available_pools_for_run     # located at this file
        get_pool_info                    # located at this file
    Return:
        display_pools_for_run
    """
    display_pools_for_run = {}
    pools_to_update = get_available_pools_for_run()
    if pools_to_update:
        display_pools_for_run = get_pool_info(pools_to_update)
    return display_pools_for_run


def get_pool_instance_from_id(pool_id):
    pool_obj = wetlab.models.LibraryPool.objects.get(pk__exact=pool_id)
    return pool_obj


def fetch_reagent_kits_used_in_run(form_data):
    """
    Description:
        The function fetch the reagent kits in the form
        Return an object list with the reagents user kits.
    Input:
        form_data    # data from the user form
    Fucntion:
        update_usage_user_lot_kit       # located at core.utils.handoling_commercial_kits
    Return:
        user_reagents_kit_objs
    """
    user_reagents_kit_objs = []
    if form_data["commercialKits"] != "":
        commercial_kit_names = form_data["commercialKits"].split(",")
        for kit_name in commercial_kit_names:
            if form_data[kit_name] == "":
                continue
            user_reagents_kit_objs.append(
                core.utils.commercial_kits.update_usage_user_lot_kit(
                    form_data[kit_name]
                )
            )

    return user_reagents_kit_objs


def link_pool_with_existing_run(exp_name, pool_ids):
    run_obj = wetlab.models.RunProcess.objects.filter(run_name__iexact=exp_name).last()
    for pool in pool_ids:
        pool_obj = get_pool_instance_from_id(pool)
        run_obj.set_pool(pool_obj)
    return None
