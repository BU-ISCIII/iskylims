import logging
import os
import re
import shutil
import xml.etree.ElementTree as ET
from datetime import datetime

from django.conf import settings
from django.contrib.auth.models import User
from interop import py_interop_run, py_interop_run_metrics, py_interop_summary

from wetlab.models import *
from wetlab.config import *

from .common import *
from .samplesheet import *

############################
# CRONTAB COMMON FUNCTIONS #
############################


def get_run_disk_utilization(conn, run_folder, experiment_name):
    """
    Description:
        Function to get the size of the run directory on the  remote server.
        It will use the function get_size_dir to get the size utilization
        for each subfolder
    Input:
        conn                # Connection samba object
        run_folder          # root folder to start the checking file size
        experiment_name     # Experiment name
    Functions:
        get_samba_application_shared_folder
        get_size_dir
    Return:
        disk_utilization    # in the last iteraction will return the total size of the folder
    """
    logger = logging.getLogger(__name__)
    logger.debug("%s : Starting function get_run_disk_utilization", experiment_name)
    shared_folder = get_samba_shared_folder()
    application_folder = get_samba_application_shared_folder()
    full_path_run = os.path.join(application_folder, run_folder)

    try:
        get_full_list = conn.listPath(shared_folder, run_folder)
    except Exception:
        string_message = experiment_name + " : Unable to get the folder " + run_folder
        logging_errors(string_message, True, False)
        logger.debug(
            "%s : End function get_run_disk_utilization with Exception", experiment_name
        )
        raise

    rest_of_dir_size = 0
    data_dir_size = 0
    images_dir_size = 0
    MEGA_BYTES = 1024 * 1024
    disk_utilization = {}
    logger.info("Start folder iteration ")
    for item_list in get_full_list:
        if item_list.filename == "." or item_list.filename == "..":
            continue
        if item_list.filename == "Data":
            logger.info(
                "%s : Starting getting disk space utilization for Data Folder",
                experiment_name,
            )
            dir_data = os.path.join(run_folder, "Data")
            data_dir_size = get_size_dir(dir_data, conn, shared_folder)

        elif item_list.filename == "Images":
            logger.info(
                "%s : Starting getting disk space utilization for Images Folder",
                experiment_name,
            )
            dir_images = os.path.join(full_path_run, "Images")
            images_dir_size = get_size_dir(dir_images, conn, shared_folder)

        if item_list.isDirectory:
            item_dir = os.path.join(full_path_run, item_list.filename)
            rest_of_dir_size += get_size_dir(item_dir, conn, shared_folder)
        else:
            rest_of_dir_size += item_list.file_size

    disk_utilization["useSpaceFastaMb"] = "{0:,}".format(
        round(data_dir_size / MEGA_BYTES)
    )
    disk_utilization["useSpaceImgMb"] = "{0:,}".format(
        round(images_dir_size / MEGA_BYTES)
    )
    disk_utilization["useSpaceOtherMb"] = "{0:,}".format(
        round(rest_of_dir_size / MEGA_BYTES)
    )

    logger.info(
        "%s : End  disk space utilization for runID  %s", experiment_name, run_folder
    )
    logger.debug("%s : End function get_run_disk_utilization", experiment_name)
    return disk_utilization


def get_size_dir(directory, conn, shared_folder):
    """
    Description:
        Recursive function to get the size of the run directory on the remote server.
    Input:
        conn            # Connection samba object
        directory       # root folder to start the checking file size
        shared_folder   # shared remote folder
    Return:
        count_file_size # in the last iteraction will return the total folder size
    """
    count_file_size = 0
    file_list = conn.listPath(shared_folder, directory)
    for sh_file in file_list:
        if sh_file.isDirectory:
            if sh_file.filename == "." or sh_file.filename == "..":
                continue
            sub_directory = os.path.join(directory, sh_file.filename)
            count_file_size += get_size_dir(sub_directory, conn, shared_folder)
        else:
            count_file_size += sh_file.file_size

    return count_file_size


def assign_used_library_in_run(
    run_process_obj,
    l_sample_sheet_path,
    experiment_name,
):
    """
    Description:
        The function first check if the run has already assinged to projects.
        If not , it reads the sample sheet, fetch the projects and asign them to the run.
    Input:
        run_process_obj     # runProcess object
        l_sample_sheet_path   # path the sample sheet
        experiment_name     # experiment name
    Functions:
        get_projects_in_run # located at utils/sample_sheet_utils.py
    Return:
        None
    """
    logger = logging.getLogger(__name__)
    logger.debug("%s : Starting function assign_used_library_in_run", experiment_name)
    if run_process_obj.get_index_library() == "None":
        index_library_name = get_index_library_name(l_sample_sheet_path)
        run_process_obj.update_index_library(index_library_name)
        logger.info("%s : Defined index library name", experiment_name)
    else:
        logger.info("%s : Already defined index library name", experiment_name)
    logger.debug("%s : End function assign_used_library_in_run", experiment_name)
    return


def assign_projects_to_run(run_process_obj, sample_sheet_file, experiment_name):
    """
    Description:
        The function first check if the run has already assinged to projects.
        If not , it reads the sample sheet, fetch the projects and asign them to the run.
    Input:
        run_process_obj     # runProcess object
        sample_sheet_file   # path the sample sheet
        experiment_name     # experiment name
    Functions:
        get_projects_in_run # located at utils/sample_sheet_utils.py
    Return:
        None
    """
    logger = logging.getLogger(__name__)
    logger.debug("%s : Starting function assign_projects_to_run", experiment_name)
    if run_process_obj.projects_set.all().exists():
        logger.info("%s  : Projects already defined in the run", experiment_name)
        logger.debug("%s : End function assign_projects_to_run", experiment_name)
        return
    # fetch the project from sample sheet
    projects_objs = []
    project_with_users = get_projects_in_run(sample_sheet_file)
    for project, user in project_with_users.items():
        if not Projects.objects.filter(project_name__exact=project).exists():
            project_data = {}
            project_data["projectName"] = project
            if User.objects.filter(username__iexact=user).exists():
                project_data["user_id"] = User.objects.filter(
                    username__iexact=user
                ).last()
            else:
                project_data["user_id"] = None
            project_obj = Projects.objects.create_new_empty_project(project_data)
            projects_objs.append(project_obj)
            # assign project to run
            project_obj.run_process.add(run_process_obj)
            logger.info("%s : Project name  %s added ", experiment_name, project)
        else:
            project_obj = Projects.objects.filter(project_name__exact=project).last()
            if project_obj not in projects_objs:
                projects_objs.append(project_obj)
                project_obj.run_process.add(run_process_obj)
                logger.info("%s : Project name  %s added ", experiment_name, project)
    logger.debug("%s : End function assign_projects_to_run", experiment_name)
    return


def check_sequencer_status_from_log_file(
    log_file_content, log_cycles, number_of_cycles, experiment_name
):
    """
    Description:
        The function checks in the logs if run was canceled. If not canceled, then check
        the number of created log files are the same as they expected to have when
        run processing is completed
    Input:
        log_file_content    # content of the log file
        log_cycles          # number of cycles executed on sequencer
        number_of_cycles    # number of cycle to check in logs
        experiment_name     # experiment name
    Return:
        completed and run_completion_date in case the sequencer ends sucessfully.
        If not return status of the run (cancelled/still_running)
        ERROR returns if not able to fecth log files
    """
    logger = logging.getLogger(__name__)
    logger.debug(
        "%s : Starting function check_sequencer_status_from_log_file", experiment_name
    )
    run_completion_date = ""

    if "Cancel" in log_file_content:
        run_process_obj = RunProcess.objects.filter(
            run_name__exact=experiment_name
        ).last()
        if run_process_obj.get_forced_continue_on_error():
            logger.warning(
                "%s : Forced to continue on execution that was canceled",
                experiment_name,
            )
        else:
            logger.warning("%s : Sequencer execution was canceled", experiment_name)
            status = "cancelled"
            logger.debug(
                "%s : End function check_sequencer_status_from_log_file",
                experiment_name,
            )
            return status, run_completion_date
    elif log_cycles != number_of_cycles:
        logger.info("%s : run sequencer is still running", experiment_name)
        status = "still_running"
        logger.debug(
            "%s : End function check_sequencer_status_from_log_file", experiment_name
        )
        return status, run_completion_date
    # Fetch the run completcion time
    status = "completed"
    logger.info("%s : run in sequencer is completed", experiment_name)
    last_line_in_file = log_file_content.split("\n")[-2]
    last_log_time = last_line_in_file.split(" ")[0:2]
    last_log_time[0] = str("20" + last_log_time[0])
    string_completion_date = " ".join(last_log_time)
    run_completion_date = datetime.strptime(
        string_completion_date, "%Y-%m-%d %H:%M:%S.%f"
    )
    logger.debug(
        "%s : End function check_sequencer_status_from_log_file", experiment_name
    )
    return status, run_completion_date


def check_sequencer_status_from_completion_file(l_run_completion, experiment_name):
    """
    Description:
        The function will check if the run in sequencer was successful
    Input:
        l_run_completion  # local path for the run completion file
        experiment_name   # Contains the experiment name
    Functions:
        find_xml_tag_text # located at utils.common
    Constant:
        COMPLETION_TAG
    Return
        True if successfuly completed
    """
    logger = logging.getLogger(__name__)
    logger.debug(
        "%s : Starting function for check_sequencer_status_from_completion_file",
        experiment_name,
    )
    # check if NextSEq run have been successful completed
    status_run = find_xml_tag_text(l_run_completion, COMPLETION_TAG)
    if status_run not in COMPLETION_SUCCESS:
        logger.info(
            "%s : Run in sequencer was not completed but %s",
            experiment_name,
            status_run,
        )
        string_message = (
            experiment_name
            + " : Sequencer Run was not completed. Reason was "
            + status_run
        )
        logging_warnings(string_message, False)
        logger.debug(
            "%s : End function for check_sequencer_status_from_completion_file",
            experiment_name,
        )
        return False
    else:
        logger.info("%s : Sequencer Run successfuly completed ", experiment_name)
        logger.debug(
            "%s : End function for check_sequencer_status_from_completion_file",
            experiment_name,
        )
        return True


def check_sequencer_run_is_completed(
    conn, run_folder, platform, number_of_cycles, experiment_name
):
    """
    Description:

    Input:
        conn                # Connection Samba object
        run_folder          # run folder on the remote server
        platform            # platform name
        number_of_cycles    # number of cycle to check in logs
        experiment_name     # experiment name
    Functions:
        get_samba_application_shared_folder     # located at this file
        get_latest_run_procesing_log            # located at this file
        check_sequencer_status_from_log_file    # located at this file
    Constants:
        PLATFORM_WAY_TO_CHECK_RUN_COMPLETION
        RUN_LOG_FOLDER
        RUN_COMPLETION_XML_FILE
        RUN_COMPLETION_TXT_FILE
    Return:
        status and run_completion_date
        ERROR returns if not able to fecth log files or not defined the way to
        check the termination of the run in sequencer
    """
    logger = logging.getLogger(__name__)
    logger.debug(
        "%s : Starting function check_sequencer_run_is_completed", experiment_name
    )
    way_to_check = ""
    for method in PLATFORM_WAY_TO_CHECK_RUN_COMPLETION:
        if method[0] not in platform:
            continue
        way_to_check = method[1]

    if way_to_check == "logs":
        log_folder = os.path.join(
            get_samba_application_shared_folder(), run_folder, RUN_LOG_FOLDER
        )
        try:
            log_cycles, log_file_content = get_latest_run_procesing_log(
                conn, log_folder, experiment_name
            )
        except Exception:
            string_message = (
                experiment_name
                + " : Unable to fetch the log files on run folder "
                + run_folder
                + "/"
                + log_folder
            )
            logging_errors(string_message, True, False)
            logger.debug(
                "%s : End function check_sequencer_run_is_completed with exeception",
                experiment_name,
            )
            return {"ERROR": 18}, ""
        status, run_completion_date = check_sequencer_status_from_log_file(
            log_file_content, log_cycles, number_of_cycles, experiment_name
        )
        logger.debug(
            "%s : End function check_sequencer_run_is_completed", experiment_name
        )
        return status, run_completion_date

    elif way_to_check == "xml_file":
        l_run_completion = os.path.join(RUN_TEMP_DIRECTORY, RUN_COMPLETION_XML_FILE)
        s_run_completion = os.path.join(
            get_samba_application_shared_folder(), run_folder, RUN_COMPLETION_XML_FILE
        )

        try:
            l_run_completion = fetch_remote_file(
                conn, run_folder, s_run_completion, l_run_completion
            )
            logger.info(
                "%s : Sucessfully fetch of Completion status file", experiment_name
            )
        except Exception:
            logger.warning(
                "%s : Completion status file is not present on the run folder %s",
                experiment_name,
                run_folder,
            )
            logger.debug(
                "%s : End function check_sequencer_run_is_completed", experiment_name
            )
            return "still_running", ""

        if check_sequencer_status_from_completion_file(
            l_run_completion, experiment_name
        ):
            os.remove(l_run_completion)
            logger.info("%s : Deleted Run Completion file", experiment_name)
            shared_folder = get_samba_shared_folder()
            conversion_attributes = conn.getAttributes(shared_folder, s_run_completion)
            run_completion_date = datetime.fromtimestamp(
                int(conversion_attributes.create_time)
            ).strftime("%Y-%m-%d %H:%M:%S")
            logger.debug(
                "%s : End function for handling NextSeq run with exception",
                experiment_name,
            )
            return "completed", run_completion_date
        os.remove(l_run_completion)
        logger.info("%s : Deleted Run Completion file", experiment_name)
        logger.debug(
            "%s : End function check_sequencer_run_is_completed with exception",
            experiment_name,
        )
        return "cancelled", ""
    elif way_to_check == "txt_file":
        # l_run_completion = os.path.join(RUN_TEMP_DIRECTORY, RUN_COMPLETION_TXT_FILE)
        s_run_completion = os.path.join(
            get_samba_application_shared_folder(), run_folder, RUN_COMPLETION_TXT_FILE
        )

        try:
            shared_folder = get_samba_shared_folder()
            conversion_attributes = conn.getAttributes(shared_folder, s_run_completion)
            logger.info(
                "%s : Sucessfully fetch of Completion status file", experiment_name
            )
        except Exception:
            logger.warning(
                "%s : Completion status file is not present on the run folder %s",
                experiment_name,
                run_folder,
            )
            logger.debug(
                "%s : End function check_sequencer_run_is_completed", experiment_name
            )
            return "still_running", ""
        try:
            run_completion_date = datetime.fromtimestamp(
                int(conversion_attributes.create_time)
            ).strftime("%Y-%m-%d %H:%M:%S")
        except Exception:
            run_completion_date = ""
            logger.warning(
                "%s : Unable to collect date for completion file ", experiment_name
            )
            logger.debug(
                "%s : End function check_sequencer_run_is_completed", experiment_name
            )
        return "completed", run_completion_date
    else:
        string_message = (
            experiment_name
            + " : way to check the completion run is not defined for  "
            + platform
        )
        logging_errors(string_message, True, False)
        logger.debug(
            "%s : End function check_sequencer_run_is_completed with exeception",
            experiment_name,
        )
        return {"ERROR": 25}, ""


def copy_sample_sheet_to_remote_folder(
    conn, sample_sheet_path, run_folder, experiment_name
):
    """
    Description:
        The functions copy the sample sheet to the remote run folder
    Input:
        conn                # Connection Samba object
        sample_sheet_path   # path of the sample sheet
        run_folder          # run folder on the remote server
        experiment_name     # experiment name
    Functions:
        get_samba_application_shared_folder     # located at this file
        copy_to_remote_file                     # located at this file
    Constants:
        EMPTY_FIELDS_IN_SEQUENCER
    Return:
        new_sequencer_obj
    """
    logger = logging.getLogger(__name__)
    logger.debug(
        "%s : Starting function copy_sample_sheet_to_remote_folder", experiment_name
    )

    logger.info(
        "%s : Copy sample sheet to remote folder %s", experiment_name, run_folder
    )
    s_sample = os.path.join(
        get_samba_application_shared_folder(), run_folder, SAMPLE_SHEET
    )

    try:
        copy_to_remote_file(
            conn, run_folder, s_sample, sample_sheet_path
        )
        logger.info(
            "%s : Sucessfully copy Sample sheet to remote folder", experiment_name
        )
    except Exception:
        string_message = (
            experiment_name
            + ": Unable to copy the Sample Sheet to remote folder "
            + run_folder
        )
        logging_errors(string_message, True, False)
        handling_errors_in_run(experiment_name, "23")
        logger.debug(
            "%s : End function for copy_sample_sheet_to_remote_folder with exception",
            experiment_name,
        )
        raise Exception
    logger.debug(
        "%s : End function copy_sample_sheet_to_remote_folder", experiment_name
    )
    return


def copy_to_remote_file(conn, run_dir, remote_file, local_file):
    """
    Description:
        Function will fetch the file from remote server and copy on local  directory
    Input:
        conn    # Samba connection object
        run_dir # run folder to fetch the file
        remote_file # file name to fetch on remote server
        local_file # local copy of the file fetched
    Constants:
        SAMBA_SHARED_FOLDER_NAME
    Return:
        True if file was successfuly copy.
        Exception if file could not be fetched
    """
    logger = logging.getLogger(__name__)
    logger.debug("Starting function for copy file to remote")
    with open(local_file, "rb") as r_par_fp:
        try:
            samba_folder = get_samba_shared_folder()
            conn.storeFile(samba_folder, remote_file, r_par_fp)
            logger.info("Saving the file %s to remote server", local_file)
        except Exception:
            string_message = (
                "Unable to copy the " + local_file + "file on folder " + run_dir
            )
            logging_errors(string_message, True, False)
            raise Exception("File not copied")
    logger.debug("End function for copy file to remote")
    return True


def create_new_sequencer_lab_not_defined(sequencer_name, num_of_lanes, experiment_name):
    """
    Description:

        creates a new entry in database wit only the sequencer name and the lane numbers
    Input:
        sequencer_name    # sequencer name
        num_of_lanes        # number of lanes
        experiment_name     # experiment name
    Functions:
        get_sequencer_lanes_number_from_file # located at this file
    Constants:
        EMPTY_FIELDS_IN_SEQUENCER
    Return:
        new_sequencer_obj
    """
    logger = logging.getLogger(__name__)
    logger.debug(
        "%s : Starting function create_new_sequencer_lab_not_defined", experiment_name
    )
    seq_data = {}
    for item in EMPTY_FIELDS_IN_SEQUENCER:
        seq_data[item] = None
    seq_data["sequencerNumberLanes"] = num_of_lanes
    seq_data["sequencerName"] = sequencer_name
    new_sequencer_obj = SequencerInLab.objects.create_sequencer_in_lab(seq_data)
    logger.info("%s : Created the new sequencer in database", experiment_name)
    logger.debug(
        "%s : End function create_new_sequencer_lab_not_defined", experiment_name
    )
    return new_sequencer_obj


def fetch_remote_file(conn, run_dir, remote_file, local_file):
    """
    Description:
        Function will fetch the file from remote server and copy on local
        directory
    Input:
        conn    # Samba connection object
        run_dir # run folder to fetch the file
        remote_file # file name to fetch on remote server
        local_file # local copy of the file fetched
    Return:
        local_file if the file was successfuly copy.
        Exception if file could not be fetched
    """
    logger = logging.getLogger(__name__)
    logger.debug("%s : Starting function for fetching remote file", run_dir)
    with open(local_file, "wb") as r_par_fp:
        try:
            samba_folder = get_samba_shared_folder()
            conn.retrieveFile(samba_folder, remote_file, r_par_fp)
            logger.info(
                "%s : Retrieving the remote %s file for %s ",
                run_dir,
                remote_file,
                local_file,
            )
        except Exception:
            string_message = (
                "Unable to fetch the " + local_file + " file on folder : " + run_dir
            )
            logging_errors(string_message, True, False)
            os.remove(local_file)
            logger.debug("%s : End function for fetching remote file", run_dir)
            raise Exception("File not found")
    logger.debug("%s : End function for fetching remote file", run_dir)
    return local_file


def get_latest_run_procesing_log(conn, log_folder, experiment_name):
    """
    Description:
        The function will find the latest log file for the input folder
    Input:
        conn        # samba connection object
        log_folder  # folder name
        experiment_name     # experiment name
    Constant:
        RUN_TEMP_DIRECTORY
    Return:
        number_of_cycles
        file_content
    """
    logger = logging.getLogger(__name__)
    logger.debug("%s : Starting function get_latest_run_procesing_log", experiment_name)
    shared_folder = get_samba_shared_folder()
    folder_logs = os.path.join("/", get_samba_application_shared_folder(), log_folder)

    remote_file_list = conn.listPath(shared_folder, folder_logs)
    max_cycle = -1
    logger.info("%s : Succesful connection to fetch logs files", experiment_name)
    for sfh in remote_file_list:
        if sfh.isDirectory:
            continue
        file_remote = sfh.filename
        if file_remote.endswith(".log"):
            log_file = re.search(".*_Cycle(\d+)_.*", file_remote)

            cycle_number = int(log_file.group(1))
            if cycle_number > max_cycle:
                max_cycle = cycle_number
                latest_log = file_remote
    logger.info("%s : Fetching the latest log file  %s  ", experiment_name, latest_log)
    temporary_log = os.path.join(RUN_TEMP_DIRECTORY, "miseq_cycle.log")
    s_latest_log = os.path.join(log_folder, latest_log)

    temporary_log = fetch_remote_file(conn, log_folder, s_latest_log, temporary_log)
    logger.info(
        "%s : copied to tmp folder the log is : %s", experiment_name, s_latest_log
    )
    with open(temporary_log, "r", encoding="utf8") as fh:
        log_file_content = fh.read()

    os.remove(temporary_log)
    logger.debug("%s : End function get_latest_run_procesing_log", experiment_name)
    return max_cycle, log_file_content


def get_new_runs_from_remote_server(processed_runs, conn, shared_folder):
    """
    Description:
        The function fetch the folder names from the remote server and
        returns a list containing the folder names that have not been
        processed yet.
    Input:
        processed_runs  # full path and name of the file
        conn # samba connection object
        shared_folder   # shared folder in the remote server
    Functions:
        get_samba_application_shared_folder     # located at this file
    Return:
        new runs
    """
    logger = logging.getLogger(__name__)
    logger.debug("Starting function get_new_runs_on_remote_server")
    new_runs = []
    run_data_root_folder = os.path.join("/", get_samba_application_shared_folder())
    logger.debug("Shared folder  on remote server is : %s", run_data_root_folder)
    run_folder_list = conn.listPath(shared_folder, run_data_root_folder)
    for sfh in run_folder_list:
        if sfh.isDirectory:
            folder_run = sfh.filename
            if folder_run == "." or folder_run == "..":
                continue
            # if the run folder has been already process continue searching
            if folder_run in processed_runs:
                continue
            else:
                logger.info(" %s  : Found new folder run ", folder_run)
                new_runs.append(folder_run)
    logger.debug("End function get_new_runs_on_remote_server")
    return new_runs


def get_remote_sample_sheet(conn, new_run, experiment_name):
    """
    Description:
        The function will get the sample sheet from remote server, and it return the
        path of the sample sheet
    Input:
        conn                # samba connection instance
        new_run             # run folder on the remote server
        experiment_name     # experiment name
    Constants:
        RUN_TEMP_DIRECTORY
        SAMPLE_SHEET
        SAMBA_APPLICATION_FOLDER_NAME
        SAMPLE_SHEET
    Functions:
        get_samba_application_shared_folder     # located at this file
    Return:
        l_sample_sheet_path
    """
    logger = logging.getLogger(__name__)
    logger.debug("%s  : Starting function get_remote_sample_sheet", experiment_name)

    l_sample_sheet_path = os.path.join(RUN_TEMP_DIRECTORY, SAMPLE_SHEET)
    s_sample_sheet_path = os.path.join(
        get_samba_application_shared_folder(), new_run, SAMPLE_SHEET
    )
    try:
        fetch_remote_file(
            conn, new_run, s_sample_sheet_path, l_sample_sheet_path
        )
        logger.info("%s : Sucessfully fetch of Sample Sheet file", experiment_name)
    except Exception:
        error_message = "Unable to fetch Sample Sheet file for folder :" + new_run
        logging_errors(error_message, True, False)
        logger.debug("%s  : End function get_remote_sample_sheet", experiment_name)
        return None

    logger.debug("%s  : End function get_remote_sample_sheet", experiment_name)
    return l_sample_sheet_path


def get_run_platform_from_file(l_run_parameter):
    """
    Description:
        The function will get the run platform for the xml element tag in the
        file and it will return the platform used
    Input:
        l_run_parameter  # file to find the tag
    Return:
        platform
    """
    platform = find_xml_tag_text(l_run_parameter, APPLICATION_NAME_TAG)

    return platform


def get_run_process_obj_or_create_if_not_exists(experiment_name):
    """
    Description:
        The function get the run_proces obj or it is created if does not exists
    Input:
        experiment_name     # experiment name
    Return:
        run_process_obj
    """
    logger = logging.getLogger(__name__)
    logger.debug("Starting function get_run_process_obj_or_create_if_not_exists")
    if RunProcess.objects.filter(run_name__exact=experiment_name).exists():
        run_process_obj = RunProcess.objects.filter(
            run_name__exact=experiment_name
        ).last()
    else:
        run_data = {}
        run_data["experiment_name"] = experiment_name
        run_process_obj = RunProcess.objects.create_new_run_from_crontab(run_data)
        logger.info("%s  : New RunProcess instance created", experiment_name)
    logger.debug("End function get_run_process_obj_or_create_if_not_exists")
    return run_process_obj


def get_sequencer_obj_or_create_if_no_exists(running_parameters, experiment_name):
    """
    Description:
        The function get the sequencer obj or it is created if does not exists
    Input:
        running_parameters      # information in case a new sequencer must be defined
        experiment_name         # experiment name
    Functions:
        create_new_sequencer_lab_not_defined    # located at this file
    Return:
        sequencer_obj
    """
    logger = logging.getLogger(__name__)
    logger.debug("Starting function get_sequencer_obj_or_create_if_no_exists")
    if SequencerInLab.objects.filter(
        sequencer_name__exact=running_parameters["instrument"]
    ).exists():
        sequencer_obj = SequencerInLab.objects.filter(
            sequencer_name__exact=running_parameters["instrument"]
        ).last()

    else:
        string_message = (
            experiment_name
            + " : "
            + running_parameters["instrument"]
            + " no sequencer defined "
        )
        logging_errors(string_message, True, False)
        sequencer_obj = create_new_sequencer_lab_not_defined(
            running_parameters["instrument"],
            running_parameters["running_data"]["NumLanes"],
            experiment_name,
        )
        logger.info(
            "%s : Continue the proccess after creating the new sequencer",
            experiment_name,
        )

    logger.debug("End function get_sequencer_obj_or_create_if_no_exists")
    return sequencer_obj


def get_samba_application_shared_folder():
    """
    Description:
        The function get in database the shared application folder used for samba connections
        Application folder is a sub_directory of the shared folder. In many cases this value is empty
    Return:
        samba_folder_name
    """
    return SambaConnectionData.objects.last().get_samba_application_folder_name()


def get_samba_shared_folder():
    """
    Description:
        The function get in database the shared folder used for samba connections
    Return:
        samba_folder_name
    """
    return SambaConnectionData.objects.last().get_samba_shared_folder_name()


def handling_errors_in_run(experiment_name, error_code):
    """
    Description:
        Function will manage the error situation where the run must be
        set to run state ERROR
    Input:
        experiment_name # name of the run to be updated
        error_code      # Error code
    Return:
        True
    """
    logger = logging.getLogger(__name__)
    logger.debug("%s : Starting function handling_errors_in_run", experiment_name)
    logger.info("%s : Set run to ERROR state", experiment_name)
    if RunProcess.objects.filter(run_name__exact=experiment_name).exists():
        run_process_obj = RunProcess.objects.filter(
            run_name__exact=experiment_name
        ).last()
        run_process_obj.set_run_error_code(error_code)
        logger.info("%s : is now on ERROR state", experiment_name)
    else:
        logger.info(
            "%s : experiment name is not defined yet in database", experiment_name
        )
    logger.debug("%s : End function handling_errors_in_run", experiment_name)
    return True


def parsing_run_info_and_parameter_information(
    l_run_info, l_run_parameter, experiment_name
):
    """
    Description:
        The function is called for parsing the RunInfo and RunParameter  files.
        After parsing the RunningParameters database table will be
        updated with a new row having the parsed data
        Empty values will be set for MiSeq runs that exist on NextSeq
        but not in MiSeq runs
    Input:
        run_info    # contains the path for RunInfo.xml file
        run_parameter # contains the path for RunParameter.xml file
        experiment_name      # contains the experiment name
    CONSTANTS:
        FIELDS_TO_COLLECT_FROM_RUN_INFO_FILE
        SETUP_TAG
        NUMBER_OF_LANES_TAG
        APPLICATION_NAME_TAG
     Return:
        parsing_data
    """
    logger = logging.getLogger(__name__)
    logger.debug("%s : Starting function parsing_run_information", experiment_name)
    running_data = {}
    parsing_data = {}
    image_channel = []

    ############################
    # parsing RunInfo.xml file #
    ############################
    run_data = ET.parse(l_run_info)
    run_root = run_data.getroot()
    logger.info("%s : parsing the runInfo.xml file ", experiment_name)
    p_run = run_root[0]
    # getting the common values NextSeq and MiSeq
    logger.info("%s  : Fetching Flowcell and FlowcellLayout data ", experiment_name)
    running_data["Flowcell"] = p_run.find("Flowcell").text
    try:
        running_data["FlowcellLayout"] = p_run.find("FlowcellLayout").attrib
    except Exception:
        running_data["FlowcellLayout"] = ""
        string_message = (
            experiment_name
            + " : Parameter  FlowcellLayout  not found in RunParameter.xml"
        )
        logging_warnings(string_message, False)

    for i in run_root.iter("Name"):
        image_channel.append(i.text)

    running_data["ImageChannel"] = image_channel
    try:
        running_data["ImageDimensions"] = p_run.find("ImageDimensions").attrib
    except Exception:
        running_data["ImageDimensions"] = ""
        logger.debug(
            "%s : There is no image dimesions on runInfo file", experiment_name
        )
    # get the instrument for NextSeq run
    parsing_data["instrument"] = p_run.find("Instrument").text

    # parsing RunParameter.xml file
    
    logger.info("%s : Parsing the runParameter.xml file  ", experiment_name)
    parameter_data = ET.parse(l_run_parameter)
    parameter_data_root = parameter_data.getroot()
    # getting the common values NextSeq and MiSeq
    for field in FIELDS_TO_COLLECT_FROM_RUN_INFO_FILE:
        try:
            running_data[field] = parameter_data_root.find(field).text
        except Exception:
            # get the tags item for searching when tagas are in different Caps and lower combination
            # because of new sintax in NovaSeq
            try:
                tag_found_in_case_insensitive = False
                for element in parameter_data_root.iter():
                    if field.lower() == element.tag.lower():
                        running_data[field] = parameter_data_root.find(element.tag).text
                        tag_found_in_case_insensitive = True
                        break
                if not tag_found_in_case_insensitive:
                    running_data[field] = ""
                    string_message = (
                        experiment_name
                        + " : Parameter "
                        + field
                        + " not found looking for case insensitive in RunParameter.xml"
                    )
                    logging_warnings(string_message, False)
            except Exception:
                running_data[field] = ""
                string_message = (
                    experiment_name
                    + " : Parameter "
                    + field
                    + " unable to fetch in RunParameter.xml"
                )
                logging_warnings(string_message, False)

    # get the nuber of lanes in case sequencer lab is not defined
    if parameter_data_root.find(SETUP_TAG):
        param_in_setup = ["ApplicationVersion", "NumTilesPerSwath"]
        for i in range(len(param_in_setup)):
            try:
                running_data[param_in_setup[i]] = (
                    parameter_data_root.find("Setup").find(param_in_setup[i]).text
                )
            except Exception:
                string_message = (
                    experiment_name
                    + " : Parameter "
                    + param_in_setup[i]
                    + " not found in RunParameter.xml"
                )
                logging_warnings(string_message, False)
                continue
        # collect information for MiSeq and NextSeq
        for setup_field in FIELDS_TO_FETCH_FROM_SETUP_TAG:
            try:
                running_data[setup_field] = (
                    parameter_data_root.find(SETUP_TAG).find(setup_field).text
                )
            except Exception:
                running_data[setup_field] = ""
                string_message = (
                    experiment_name
                    + " : Parameter in Setup -- "
                    + setup_field
                    + " unable to fetch in RunParameter.xml"
                )
                logging_warnings(string_message, False)

        if "MiSeq" in running_data[APPLICATION_NAME_TAG]:
            # initialize paramters in case there are not exists on runParameter file
            for i in range(len(READ_NUMBER_OF_CYCLES)):
                running_data[READ_NUMBER_OF_CYCLES[i]] = ""
            # get the length index number for reads and indexes for MiSeq Runs
            for run_info_read in parameter_data_root.iter(RUN_INFO_READ_TAG):
                try:
                    index_number = int(run_info_read.attrib[NUMBER_TAG]) - 1
                    running_data[
                        READ_NUMBER_OF_CYCLES[index_number]
                    ] = run_info_read.attrib[NUMBER_CYCLES_TAG]
                except Exception:
                    string_message = (
                        experiment_name
                        + " : Parameter RunInfoRead: Read Number not found in RunParameter.xml"
                    )
                    logging_warnings(string_message, False)
                    continue
    else:
        # Collect information for NovaSeq
        for novaseq_field in FIELDS_NOVASEQ_TO_FETCH_TAG:
            try:
                running_data[novaseq_field] = parameter_data_root.find(
                    novaseq_field
                ).text
            except Exception:
                running_data[novaseq_field] = ""
                string_message = (
                    experiment_name
                    + " : Parameter in Setup -- "
                    + novaseq_field
                    + " unable to fetch in RunParameter.xml"
                )
                logging_warnings(string_message, False)
    # get date for miSeq and NextSeq with the format yymmdd
    date = p_run.find("Date").text
    try:
        run_date = datetime.strptime(date, "%y%m%d")
    except Exception:
        # get date for novaseq sequencer
        date = p_run.find("Date").text.split(" ")[0]
        try:
            run_date = datetime.strptime(date, "%m/%d/%Y")
        except Exception:
            run_date = ""

    # updating the date fetched from the Date tag for run and project
    logger.debug("%s : Found date that was recorded the Run %s", experiment_name, date)
    parsing_data["running_data"] = running_data
    parsing_data["run_date"] = run_date
    logger.debug("%s : End function nextseq_parsing_run_information", experiment_name)
    return parsing_data


def save_run_parameters_data_to_database(
    run_parameters, run_process_obj, experiment_name
):
    """
    Description:
        The function save the run parameters if they are not store yet
    Input:
        run_parameters      # dictionnary with the running parameters
        run_process_obj      # run process object
        experiment_name     # name of the run
    Return:
        run_parameter_obj

    """
    logger = logging.getLogger(__name__)
    logger.debug(
        "%s : Starting function save_run_parameters_data_to_database", experiment_name
    )

    if RunningParameters.objects.filter(run_name_id=run_process_obj).exists():
        run_parameter_objs = RunningParameters.objects.filter(
            run_name_id=run_process_obj
        )
        for run_parameter_obj in run_parameter_objs:
            logger.info(
                "%s  : Deleting RunParameters object from database", experiment_name
            )
            run_parameter_obj.delete()
    run_parameter_obj = RunningParameters.objects.create_running_parameters(
        run_parameters, run_process_obj
    )
    logger.info("%s  : Created RunParameters object on database", experiment_name)
    logger.debug(
        "%s : End function save_run_parameters_data_to_database", experiment_name
    )
    return run_parameter_obj


def store_sample_sheet_if_not_defined_in_run(
    run_process_obj, l_sample_sheet_path, experiment_name
):
    """
    Description:
        The function will move the sample sheet from the local temporary
        folder to the folder destination defined in RUN_SAMPLE_SHEET_DIRECTORY
        It will update the field sampleSheet in database
    Input:
        l_sample_sheet_path  # local copy of sample sheet file
        experiment_name  # name of the run
    Import:
        datetime
        os
    Variable:
        run_updated     # RunProcess object
        new_sample_sheet_name  # renamed sample sheet including time to have
                                a unique file name
        now     # Present time of the server
        run_updated # RunProcess object for miseq run
        sample_sheet_on_database # sample sheet path used in database
        timestr     # Present time including 3 digits for miliseconds
    Return:
        sample_sheet_on_database
    """
    logger = logging.getLogger(__name__)
    logger.debug(
        "%s : Starting the function store_sample_sheet_in_run", experiment_name
    )
    # Get the present time in miliseconds to add to have an unique file name
    now = datetime.now()
    timestr = now.strftime("%Y%m%d-%H%M%S.%f")[:-3]
    new_sample_sheet_name = "SampleSheet" + timestr + ".csv"

    new_sample_sheet_file = os.path.join(
        settings.MEDIA_ROOT, RUN_SAMPLE_SHEET_DIRECTORY, new_sample_sheet_name
    )
    logger.info("%s : new sample sheet name %s", experiment_name, new_sample_sheet_file)
    # Path to be included in database
    sample_sheet_on_database = os.path.join(
        RUN_SAMPLE_SHEET_DIRECTORY, new_sample_sheet_name
    )
    # Move sample sheet to final folder
    os.rename(l_sample_sheet_path, new_sample_sheet_file)
    # Update the run with the sample sheet information  (full_path, relative_path, file_name)
    run_process_obj.update_sample_sheet(new_sample_sheet_file, new_sample_sheet_name)

    logger.info("%s : Updated runProccess table with the sample sheet", experiment_name)
    logger.debug("%s : End function store_sample_sheet_in_run", experiment_name)
    return sample_sheet_on_database


def waiting_time_expired(run_process_obj, time_to_check, maximun_time, experiment_name):
    """
    Description:
        The function get the time run was recorded to compare  with the present time.
        If the value is less that the allowed time to wait  will return False.
        True is returned if the time is bigger
    Input:
        run_process_obj     # run process object
        time_to_check       # reference time to be checked
        maximun_time        # maximm number of days to wait
        experiment_name     # experiment name to be checked
    Return:
        True if the number of days is bigger that the maximum number of days to wait
    """
    logger = logging.getLogger(__name__)
    logger.debug("%s : Starting function waiting_time_expired", experiment_name)
    today = datetime.now().date()
    number_of_days = abs((today - time_to_check).days)
    if number_of_days > int(maximun_time):
        logger.info("%s  : Waiting time already exceeded", experiment_name)
        logger.debug("%s  : End function waiting_time_expired", experiment_name)
        return True
    else:
        logger.info("%s  : It is allowed to waiting more time", experiment_name)
        logger.debug("%s  : End function waiting_time_expired", experiment_name)
        return False


#########################
# RUN METRICS FUNCTIONS #
#########################


def create_run_metric_graphics(
    run_metric_folder, run_process_obj, run_folder, experiment_name
):
    """
    Description:
        The function create an entry on database with the run graphics
        by using the run metrics files
    Input:
        run_metric_folder   # local folder with the run metric files
        run_process_obj     # RunProcess object for this run
        run_folder          # run folder to store the figures
        experiment_name     # Experiment name
    Constants:
        RUN_METRIC_GRAPHIC_COMMANDS
        INTEROP_PATH
        RUN_IMAGES_DIRECTORY
        MEDIA_ROOT
    Return:
        graphic_stats_obj
    """
    logger = logging.getLogger(__name__)
    logger.debug("%s : Starting create_run_metric_graphics", experiment_name)
    present_working_dir = os.getcwd()
    run_graphic_dir = os.path.join(
        settings.MEDIA_ROOT, RUN_IMAGES_DIRECTORY, run_folder
    )
    try:
        if os.path.exists(run_graphic_dir):
            shutil.rmtree(run_graphic_dir)
        os.mkdir(run_graphic_dir)
        logger.info("%s : created new directory %s", experiment_name, run_graphic_dir)
    except Exception:
        string_message = (
            experiment_name
            + " : Unable to create folder to store graphics on "
            + run_graphic_dir
        )
        logging_errors(string_message, True, False)
        logger.debug(
            "%s : End create_run_metric_graphics with exception", experiment_name
        )
        return {"ERROR": 28}
    os.chdir(run_graphic_dir)
    logger.info(
        "%s : Changed working direcory to copy run metric graphics", experiment_name
    )
    # create the graphics
    logger.info("%s : Creating plot graphics for run id ", experiment_name)

    full_path_run_processing_tmp = os.path.join(
        present_working_dir, RUN_TEMP_DIRECTORY_PROCESSING
    )
    for graphic in RUN_METRIC_GRAPHIC_COMMANDS:
        graphic_command = os.path.join(INTEROP_PATH, graphic)
        plot_command = graphic_command + full_path_run_processing_tmp + "  | gnuplot"
        logger.debug(
            "%s : command used to create graphic is : %s", experiment_name, plot_command
        )
        os.system(plot_command)

    os.chdir(present_working_dir)
    logger.info("%s : Returning back the working directory", experiment_name)

    # removing the processing_ character in the file names
    graphic_files = os.listdir(run_graphic_dir)
    logger.info("%s : Renaming the graphic files", experiment_name)

    for graphic_file in graphic_files:
        old_file_name = os.path.join(run_graphic_dir, graphic_file)
        split_file_name = graphic_file.split("_")
        if not split_file_name[1].endswith(PLOT_EXTENSION):
            split_file_name[1] = split_file_name[1] + PLOT_EXTENSION
        new_file_name = os.path.join(run_graphic_dir, split_file_name[1])
        os.rename(old_file_name, new_file_name)
        logger.debug(
            "%s : Renamed file from %s to %s",
            experiment_name,
            old_file_name,
            new_file_name,
        )

    # saving the graphic location in database
    graphic_stats_obj = GraphicsStats.objects.create_graphic_run_metrics(
        run_process_obj, run_folder
    )

    logger.info("%s : Store Graphic plots in database", experiment_name)
    logger.debug("%s : End function create_graphics", experiment_name)
    return {"graph_stats_obj": graphic_stats_obj}


def delete_existing_run_metrics_table_processed(run_process_obj, experiment_name):
    """
    Description:
        The function will check if exists data stored on StatsRunSummary  and/or
        StatsRunRead for the run
    Input:
        run_process_obj     # runProcess object
        experiment_name     # experiment name
    Return:
        None
    """
    logger = logging.getLogger(__name__)
    logger.debug("%s : Starting function check_run_metrics_processed", experiment_name)
    if StatsRunSummary.objects.filter(runprocess_id=run_process_obj).exists():
        run_summary_objs = StatsRunSummary.objects.filter(runprocess_id=run_process_obj)
        for run_summary_obj in run_summary_objs:
            run_summary_obj.delete()
        logger.info("%s : Deleted rows on StatsRunSummary table", experiment_name)

    if StatsRunRead.objects.filter(runprocess_id=run_process_obj).exists():
        run_read_objs = StatsRunRead.objects.filter(runprocess_id=run_process_obj)
        for run_read_obj in run_read_objs:
            run_read_obj.delete()
        logger.info("%s : Deleted rows on StatsRunSummary table", experiment_name)

    if GraphicsStats.objects.filter(runprocess_id=run_process_obj).exists():
        graph_stats_objs = GraphicsStats.objects.filter(runprocess_id=run_process_obj)
        for graph_stats_obj in graph_stats_objs:
            graph_stats_obj.delete()
        logger.info("%s : Deleted rows on GraphicsStats table", experiment_name)
    logger.debug("%s : End function check_run_metrics_processed", experiment_name)
    return None


def delete_run_metric_files(experiment_name):
    """
    Description:
        The function delete the files used for collecting the run metrics
    Input:
        experiment_name     # experiment name
    Return:
        None
    """
    logger = logging.getLogger(__name__)
    logger.debug("%s : Starting function delete_run_metric_files", experiment_name)
    local_metric_folder = os.path.join(RUN_TEMP_DIRECTORY_PROCESSING, RUN_METRIC_FOLDER)
    l_run_parameter = os.path.join(RUN_TEMP_DIRECTORY_PROCESSING, RUN_PARAMETER_FILE)
    l_run_info = os.path.join(RUN_TEMP_DIRECTORY_PROCESSING, RUN_INFO)
    local_files = [l_run_parameter, l_run_info]
    for local_file in local_files:
        if os.path.exists(local_file):
            try:
                os.remove(local_file)
            except Exception:
                string_message = experiment_name + " : Unable to delete " + local_file
                logging_errors(string_message, True, False)
                continue
    logger.info("%s : Deleted temporary files", experiment_name)
    if os.path.exists(local_metric_folder):
        try:
            shutil.rmtree(local_metric_folder)
            logger.info("%s : Deleted Folder %s", experiment_name, local_metric_folder)
        except Exception:
            string_message = (
                experiment_name + " : Unable to delete  folder " + local_metric_folder
            )
            logging_errors(string_message, True, False)

    logger.debug("%s : End function delete_run_metric_files", experiment_name)
    return


def get_run_metric_files(conn, run_folder, experiment_name):
    """
    Description:
        The function will collect the run metric files created by sequencer as part of the run process.
    Input:
        conn # Connection samba object
        run_folder   # folder run to fetch the remote files
        experiment_name # experiment name
    Constant:
        RUN_INFO
        RUN_METRIC_FOLDER
        RUN_TEMP_DIRECTORY
        RUN_PARAMETER_FILE
        STATISTICS_FOLDER
    Functions:
        get_samba_application_shared_folder
        fetch_remote_file
    Return:
        copied_files
    """
    logger = logging.getLogger(__name__)
    logger.debug("%s : Starting function get_run_metric_files", experiment_name)

    # runInfo needed for run metrics stats
    l_run_info = os.path.join(RUN_TEMP_DIRECTORY_PROCESSING, RUN_INFO)
    s_run_info = os.path.join(
        get_samba_application_shared_folder(), run_folder, RUN_INFO
    )
    # runParameters needed for run metrics stats
    l_run_parameter = os.path.join(RUN_TEMP_DIRECTORY_PROCESSING, RUN_PARAMETER_FILE)
    s_run_parameter = os.path.join(
        get_samba_application_shared_folder(), run_folder, RUN_PARAMETER_FILE
    )
    l_metric_folder = os.path.join(RUN_TEMP_DIRECTORY_PROCESSING, RUN_METRIC_FOLDER)
    s_metric_folder = os.path.join(
        get_samba_application_shared_folder(), run_folder, RUN_METRIC_FOLDER
    )
    copied_files = {}

    if not os.path.exists(l_metric_folder):
        try:
            os.makedirs(l_metric_folder)
            logger.info("%s : Created folder %s", experiment_name, l_metric_folder)
        except Exception:
            string_message = (
                experiment_name + " : cannot create folder on " + l_metric_folder
            )
            logging_errors(string_message, True, False)
            logger.debug(
                "%s : End function get_run_metric_files with error", experiment_name
            )
            return {"ERROR": 26}

    try:
        l_run_info = fetch_remote_file(conn, run_folder, s_run_info, l_run_info)
        logger.info("%s : Sucessfully fetch of RunInfo file", experiment_name)
    except Exception:
        string_message = experiment_name + " : Unable to fetch " + s_run_info
        logging_errors(string_message, True, False)
        logger.debug(
            "%s : End function get_run_metric_files with error", experiment_name
        )
        return {"ERROR": 20}
    copied_files[RUN_INFO] = l_run_info

    try:
        l_run_parameter = fetch_remote_file(
            conn, run_folder, s_run_parameter, l_run_parameter
        )
        logger.info("%s : Sucessfully fetch of RunParameter file", experiment_name)
    except Exception:
        string_message = experiment_name + " : Unable to fetch " + s_run_parameter
        logging_errors(string_message, True, False)
        logger.debug(
            "%s : End function get_run_metric_files with error", experiment_name
        )
        return {"ERROR": 21}
    copied_files[RUN_PARAMETER_FILE] = l_run_parameter

    try:
        file_list = conn.listPath(get_samba_shared_folder(), s_metric_folder)
        logger.info(
            "%s : InterOp folder found at  %s", experiment_name, s_metric_folder
        )
    except Exception:
        string_message = experiment_name + " : Unable to fetch " + s_run_parameter
        logging_errors(string_message, True, False)
        shutil.rmtree(l_metric_folder)
        logger.debug(
            "%s : End function get_run_metric_files with error", experiment_name
        )
        return {"ERROR": 27}
    # copy all binary files in interop folder to local  documents/wetlab/tmp/processing/interop
    copied_files[RUN_METRIC_FOLDER] = []
    try:
        for sh in file_list:
            if sh.isDirectory:
                continue
            else:
                run_metrics_file_name = sh.filename
                s_run_metric_file = os.path.join(s_metric_folder, run_metrics_file_name)
                l_run_metric_file = os.path.join(l_metric_folder, run_metrics_file_name)
                l_run_metric_file = fetch_remote_file(
                    conn, run_folder, s_run_metric_file, l_run_metric_file
                )
                # copied_files[RUN_METRIC_FOLDER].append(l_run_metric_file)
    except Exception:
        string_message = experiment_name + " : Unable to fetch " + s_run_metric_file
        logging_errors(string_message, True, False)
        shutil.rmtree(l_metric_folder)
        logger.debug(
            "%s : End function get_run_metric_files with error", experiment_name
        )
        return {"ERROR": 27}

    logger.debug("%s : End function get_run_metric_files", experiment_name)
    return copied_files


def parsing_run_metrics_files(
    local_run_metric_folder, run_process_obj, experiment_name
):
    """
    Description:
        The function parse the information from the run metric files
    Input:
        local_run_metric_folder   # local folder with the run metric files
        run_process_obj           # RunProcess object for this run
        experiment_name           # experiment name
    Import:
        py_interop_run
        py_interop_run_metrics
    Variables:
        bin_run_stats_summary_list # list of dictionnary with the summary
                                    information
        run_stats_read_list  # list of dictionnary with the read information
    Return:
        bin_run_stats_summary_list, run_stats_read_list
    """
    logger = logging.getLogger(__name__)
    logger.debug("%s : Starting function parsing_run_metrics", experiment_name)
    run_param_obj = RunningParameters.objects.get(run_name_id=run_process_obj)
    # get the number of lanes for the run
    number_of_lanes = int(run_param_obj.get_number_of_lanes())
    # get number of reads for the run
    num_of_reads = run_param_obj.get_number_of_reads()
    logger.info(
        "%s : Fetched run information  needed for running metrics", experiment_name
    )

    run_metrics = py_interop_run_metrics.run_metrics()
    valid_to_load = py_interop_run.uchar_vector(py_interop_run.MetricCount, 0)
    py_interop_run_metrics.list_summary_metrics_to_load(valid_to_load)
    run_metrics.read(local_run_metric_folder)
    summary = py_interop_summary.run_summary()
    py_interop_summary.summarize_run_metrics(run_metrics, summary)

    bin_run_stats_summary_list = []
    logger.info("%s : Starts collecting data for run metric ", experiment_name)
    # get the Run Summary for each Read
    for read_level in range(num_of_reads):
        run_summary_stats_level = {}
        # summary yield total
        run_summary_stats_level["yieldTotal"] = format(
            summary.at(read_level).summary().yield_g(), ".3f"
        )
        # summary projected total yield
        run_summary_stats_level["projectedTotalYield"] = format(
            summary.at(read_level).summary().projected_yield_g(), ".3f"
        )

        # percent yield
        run_summary_stats_level["aligned"] = format(
            summary.at(read_level).summary().percent_aligned(), ".3f"
        )
        # Error rate
        run_summary_stats_level["errorRate"] = format(
            summary.at(read_level).summary().error_rate(), ".3f"
        )
        # intensity cycle 1
        run_summary_stats_level["intensityCycle"] = str(
            round(summary.at(read_level).summary().first_cycle_intensity())
        )
        # Q30
        run_summary_stats_level["biggerQ30"] = format(
            summary.at(read_level).summary().percent_gt_q30(), ".3f"
        )

        run_summary_stats_level["level"] = str(read_level + 1)

        bin_run_stats_summary_list.append(run_summary_stats_level)
    logger.info("%s : Parsed run Metrics on summary level ", experiment_name)

    # get the run summary for Total
    run_summary_stats_level = {}
    # total summary
    run_summary_stats_level["yieldTotal"] = format(
        summary.total_summary().yield_g(), ".3f"
    )
    # total projected_yield_g
    run_summary_stats_level["projectedTotalYield"] = format(
        summary.total_summary().projected_yield_g(), ".3f"
    )
    # total percent aligned
    run_summary_stats_level["aligned"] = format(
        summary.total_summary().percent_aligned(), ".3f"
    )
    # total error rate
    run_summary_stats_level["errorRate"] = format(
        summary.total_summary().error_rate(), ".3f"
    )
    # total intensity cycle
    run_summary_stats_level["intensityCycle"] = str(
        round(summary.total_summary().first_cycle_intensity())
    )
    # total Q 30
    run_summary_stats_level["biggerQ30"] = format(
        summary.total_summary().percent_gt_q30(), ".3f"
    )

    run_summary_stats_level["level"] = "Total"

    logger.info("%s : Parsed run Metrics on Total lane", experiment_name)

    bin_run_stats_summary_list.append(run_summary_stats_level)

    # get the run summary for non index
    run_summary_stats_level = {}
    # non index yield
    run_summary_stats_level["yieldTotal"] = format(
        summary.nonindex_summary().yield_g(), ".3f"
    )
    #  non index projected yield
    run_summary_stats_level["projectedTotalYield"] = format(
        summary.nonindex_summary().projected_yield_g(), ".3f"
    )

    # non index percent aligned
    run_summary_stats_level["aligned"] = format(
        summary.nonindex_summary().percent_aligned(), ".3f"
    )
    # non index percent error rate
    run_summary_stats_level["errorRate"] = format(
        summary.nonindex_summary().error_rate(), ".3f"
    )
    # non index intensity cycle
    run_summary_stats_level["intensityCycle"] = str(
        round(summary.nonindex_summary().first_cycle_intensity())
    )
    # non index Q 30
    run_summary_stats_level["biggerQ30"] = format(
        summary.nonindex_summary().percent_gt_q30(), ".3f"
    )

    run_summary_stats_level["level"] = "Non Index"
    logger.info("%s : Parsed run metric for Non Index lane", experiment_name)

    bin_run_stats_summary_list.append(run_summary_stats_level)

    # information per reads
    run_stats_read_list = []
    # lan_summary= py_interop_summary.lane_summary()
    # Tiles
    for read_number in range(num_of_reads):
        for lane_number in range(number_of_lanes):
            logger.info(
                "%s : Processing run metrics stats on Read %s and on Lane %s",
                experiment_name,
                read_number,
                lane_number,
            )
            run_read_stats_level = {}
            run_read_stats_level["tiles"] = str(
                int(summary.at(read_number).at(lane_number).tile_count()) * 2
            )
            # Density (k/mm2) divide the value by 1000 to have it K/mm2
            # get the +/- with the steddev
            try:
                read_lane_density_mean = str(
                    round(
                        float(summary.at(read_number).at(lane_number).density().mean())
                        / 1000
                    )
                )
                read_lane_density_stddev = str(
                    round(
                        float(
                            summary.at(read_number).at(lane_number).density().stddev()
                        )
                        / 1000
                    )
                )
            except Exception:
                read_lane_density_mean = "NaN"
                read_lane_density_stddev = "NaN"
                string_message = experiment_name + " : Unable to convert to float "
                logging_warnings(string_message, False)
            run_read_stats_level["density"] = (
                read_lane_density_mean
                + "  "
                + chr(177)
                + "  "
                + read_lane_density_stddev
            )
            # cluster _pf  in %
            try:
                read_lane_percent_pf_mean = format(
                    summary.at(read_number).at(lane_number).percent_pf().mean(), ".3f"
                )
                read_lane_percent_pf_stddev = format(
                    summary.at(read_number).at(lane_number).percent_pf().stddev(), ".3f"
                )
            except Exception:
                read_lane_percent_pf_mean = "NaN"
                read_lane_percent_pf_stddev = "NaN"
                string_message = (
                    experiment_name
                    + " : Unable to format to float read_lane_percent_pf"
                )
                logging_warnings(string_message, False)
            run_read_stats_level["cluster_PF"] = (
                read_lane_percent_pf_mean
                + "  "
                + chr(177)
                + "  "
                + read_lane_percent_pf_stddev
            )
            # phas/ prepas in %
            try:
                read_lane_phasing_mean = format(
                    summary.at(read_number).at(lane_number).phasing().mean(), ".3f"
                )
                read_lane_phasing_dev = format(
                    summary.at(read_number).at(lane_number).phasing().stddev(), ".1f"
                )
                read_lane_prephasing_mean = format(
                    summary.at(read_number).at(lane_number).prephasing().mean(), ".3f"
                )
                read_lane_prephasing_stddev = format(
                    summary.at(read_number).at(lane_number).prephasing().stddev(), ".3f"
                )
            except Exception:
                (
                    read_lane_phasing_mean,
                    read_lane_phasing_dev,
                    read_lane_prephasing_mean,
                    read_lane_prephasing_stddev,
                ) = ("NaN", "NaN", "NaN", "NaN")
                string_message = (
                    experiment_name + " : Unable to format to float read_lane_phasing"
                )
                logging_warnings(string_message, False)
            run_read_stats_level["phas_prephas"] = (
                read_lane_phasing_mean
                + "  "
                + chr(177)
                + "  "
                + read_lane_phasing_dev
                + "  /  "
                + read_lane_prephasing_mean
                + "  "
                + chr(177)
                + "  "
                + read_lane_prephasing_stddev
            )
            # reads (M)
            try:
                run_read_stats_level["reads"] = format(
                    float(summary.at(read_number).at(lane_number).reads()) / 1000000,
                    ".3f",
                )
            except Exception:
                run_read_stats_level["reads"] = "NaN"
                string_message = (
                    experiment_name
                    + " : Unable to format to float run_read_stats_level[reads]"
                )
                logging_warnings(string_message, False)
            # reads PF (M)
            try:
                run_read_stats_level["reads_PF"] = format(
                    float(summary.at(read_number).at(lane_number).reads_pf()) / 1000000,
                    ".3f",
                )
            except Exception:
                run_read_stats_level["reads_PF"] = "NaN"
                string_message = (
                    experiment_name
                    + " : Unable to format to float run_read_stats_level[reads_PF]"
                )
                logging_warnings(string_message, False)
            # percent q30
            try:
                run_read_stats_level["q30"] = format(
                    summary.at(read_number).at(lane_number).percent_gt_q30(), ".3f"
                )
            except Exception:
                run_read_stats_level["q30"] = "NaN"
                string_message = (
                    experiment_name
                    + " : Unable to format to float run_read_stats_level[q30]"
                )
                logging_warnings(string_message, False)
            # yield _g
            try:
                run_read_stats_level["yields"] = format(
                    summary.at(read_number).at(lane_number).yield_g(), ".3f"
                )
            except Exception:
                run_read_stats_level["yields"] = "NaN"
                string_message = (
                    experiment_name
                    + " : Unable to format to float run_read_stats_level[yields]"
                )
                logging_warnings(string_message, False)
            # cycles err Rate
            try:
                run_read_stats_level["cyclesErrRated"] = str(
                    summary.at(read_number)
                    .at(lane_number)
                    .cycle_state()
                    .error_cycle_range()
                    .first_cycle()
                )
            except Exception:
                run_read_stats_level["cyclesErrRated"] = "NaN"
                string_message = (
                    experiment_name
                    + " : Unable to format to float run_read_stats_level[cyclesErrRated]"
                )
                logging_warnings(string_message, False)
            # percent_aligned
            try:
                read_lane_percent_aligned_mean = format(
                    summary.at(read_number).at(lane_number).percent_aligned().mean(),
                    ".3f",
                )
                read_lane_percent_aligned_stddev = format(
                    summary.at(read_number).at(lane_number).percent_aligned().stddev(),
                    "3f",
                )
            except Exception:
                read_lane_percent_aligned_mean, read_lane_percent_aligned_stddev = (
                    "NaN",
                    "NaN",
                )
                string_message = (
                    experiment_name
                    + " : Unable to format to float read_lane_percent_aligned_mean"
                )
                logging_warnings(string_message, False)
            run_read_stats_level["aligned"] = (
                read_lane_percent_aligned_mean
                + "  "
                + chr(177)
                + "  "
                + read_lane_percent_aligned_stddev
            )
            # error rate
            try:
                read_lane_error_rate_mean = format(
                    summary.at(read_number).at(lane_number).error_rate().mean(), ".3f"
                )
                read_lane_error_rate_stddev = format(
                    summary.at(read_number).at(lane_number).error_rate().stddev(), ".3f"
                )
            except Exception:
                read_lane_error_rate_mean, read_lane_error_rate_stddev = "NaN", "NaN"
                string_message = (
                    experiment_name
                    + " : Unable to format to float read_lane_error_rate_mean"
                )
                logging_warnings(string_message, False)
            run_read_stats_level["errorRate"] = (
                read_lane_error_rate_mean
                + "  "
                + chr(177)
                + "  "
                + read_lane_error_rate_stddev
            )
            # error rate_35
            try:
                read_lane_error_rate_35_mean = format(
                    summary.at(read_number).at(lane_number).error_rate_35().mean(),
                    ".3f",
                )
                read_lane_error_rate_35_stddev = format(
                    summary.at(read_number).at(lane_number).error_rate_35().stddev(),
                    ".3f",
                )
            except Exception:
                read_lane_error_rate_35_mean, read_lane_error_rate_35_stddev = (
                    "NaN",
                    "NaN",
                )
                string_message = (
                    experiment_name
                    + " : Unable to format to float read_lane_error_rate_35_mean"
                )
                logging_warnings(string_message, False)
            run_read_stats_level["errorRate35"] = (
                read_lane_error_rate_35_mean
                + "  "
                + chr(177)
                + "  "
                + read_lane_error_rate_35_stddev
            )
            # error rate 50
            try:
                read_lane_error_rate_50_mean = format(
                    summary.at(read_number).at(lane_number).error_rate_50().mean(),
                    ".3f",
                )
                read_lane_error_rate_50_stddev = format(
                    summary.at(read_number).at(lane_number).error_rate_50().stddev(),
                    ".3f",
                )
            except Exception:
                read_lane_error_rate_50_mean, read_lane_error_rate_50_stddev = (
                    "NaN",
                    "NaN",
                )
                string_message = (
                    experiment_name
                    + " : Unable to format to float read_lane_error_rate_50_mean"
                )
                logging_warnings(string_message, False)
            run_read_stats_level["errorRate50"] = (
                read_lane_error_rate_50_mean
                + "  "
                + chr(177)
                + "  "
                + read_lane_error_rate_50_stddev
            )
            # error rate 75
            try:
                read_lane_error_rate_75_mean = format(
                    summary.at(read_number).at(lane_number).error_rate_75().mean(),
                    ".3f",
                )
                read_lane_error_rate_75_stddev = format(
                    summary.at(read_number).at(lane_number).error_rate_75().stddev(),
                    ".3f",
                )
            except Exception:
                read_lane_error_rate_75_mean, read_lane_error_rate_75_stddev = (
                    "NaN",
                    "NaN",
                )
                string_message = (
                    experiment_name
                    + " : Unable to format to float read_lane_error_rate_75_mean"
                )
                logging_warnings(string_message, False)
            run_read_stats_level["errorRate75"] = (
                read_lane_error_rate_75_mean
                + "  "
                + chr(177)
                + "  "
                + read_lane_error_rate_75_stddev
            )
            # error rate 100
            try:
                read_lane_error_rate_100_mean = format(
                    summary.at(read_number).at(lane_number).error_rate_100().mean(),
                    ".3f",
                )
                read_lane_error_rate_100_stddev = format(
                    summary.at(read_number).at(lane_number).error_rate_100().stddev(),
                    ".3f",
                )
            except Exception:
                read_lane_error_rate_100_mean, read_lane_error_rate_100_stddev = (
                    "NaN",
                    "NaN",
                )
                string_message = (
                    experiment_name
                    + " : Unable to format to float read_lane_error_rate_100_mean"
                )
                logging_warnings(string_message, False)
            run_read_stats_level["errorRate100"] = (
                read_lane_error_rate_100_mean
                + "  "
                + chr(177)
                + "  "
                + read_lane_error_rate_100_stddev
            )
            # intensity cycle 1
            try:
                read_lane_intensity_cycle_mean = format(
                    summary.at(read_number)
                    .at(lane_number)
                    .first_cycle_intensity()
                    .mean(),
                    ".3f",
                )  # get tiles for read 1 and lane 1
                read_lane_intensity_cycle_stddev = format(
                    summary.at(read_number)
                    .at(lane_number)
                    .first_cycle_intensity()
                    .stddev(),
                    ".3f",
                )
            except Exception:
                read_lane_intensity_cycle_mean, read_lane_intensity_cycle_stddev = (
                    "NaN",
                    "NaN",
                )
                string_message = (
                    experiment_name
                    + " : Unable to format to float read_lane_intensity_cycle_mean"
                )
                logging_warnings(string_message, False)
            run_read_stats_level["intensityCycle"] = (
                read_lane_intensity_cycle_mean
                + "  "
                + chr(177)
                + "  "
                + read_lane_intensity_cycle_stddev
            )

            run_read_stats_level["read"] = str(read_number + 1)
            run_read_stats_level["lane"] = str(lane_number + 1)
            # append run_read_stats_level information to run_stats_read_list
            run_stats_read_list.append(run_read_stats_level)

    logger.debug("%s : End function parsing_run_metrics", experiment_name)
    return bin_run_stats_summary_list, run_stats_read_list


#######################
# BCL2FASTQ FUNCTIONS #
#######################


def check_demultiplexing_folder_exists(conn, run_folder, experiment_name):
    """
    Description:
        The function check if the folder of demultiplexing files exists
    Input:
        conn                # samba connection instance
        run_folder          # run folder on the remote server
        experiment_name     # Experiment name
    Constants:
        STATS_FILE_PATH
        CONVERSION_STATS_FILE
    Functions:
        get_samba_application_shared_folder
        get_samba_shared_folder
    Return:
        bcl2fastq_finish_date
    """
    logger = logging.getLogger(__name__)
    logger.debug(
        "%s : Starting function check_demultiplexing_folder_exists", experiment_name
    )
    bcl2fastq_finish_date = ""
    statistics_folder = os.path.join(
        get_samba_application_shared_folder(),
        run_folder,
        STATS_FILE_PATH,
    )

    try:
        conn.listPath(get_samba_shared_folder(), statistics_folder)
    except Exception:
        string_message = (
            experiment_name
            + " : Unable to fetch folder demultiplexing at  "
            + statistics_folder
        )
        logging_warnings(string_message, True)
        logger.debug(
            "%s : End function check_demultiplexing_folder_exists with warning",
            experiment_name,
        )
        return {"ERROR": 29}

    logger.info(
        "%s : bcl2fastq has been completed . Collecting date when finish the bcl2fastq",
        experiment_name,
    )

    s_conversion_stats = os.path.join(statistics_folder, CONVERSION_STATS_FILE)
    try:
        conversion_attributes = conn.getAttributes(
            get_samba_shared_folder(), s_conversion_stats
        )
    except Exception:
        string_message = experiment_name + " : Unable to fetch " + CONVERSION_STATS_FILE
        logging_errors(string_message, True, False)
        logger.debug(
            "%s : End function check_demultiplexing_folder_exists with error",
            experiment_name,
        )
        return {"ERROR": 31}
    bcl2fastq_finish_date = datetime.fromtimestamp(
        int(conversion_attributes.create_time)
    ).strftime("%Y-%m-%d %H:%M:%S")
    logger.info("%s : Collected Bcl2Fastq time ", experiment_name)

    logger.debug("%s : End function waiting_time_expired", experiment_name)
    return bcl2fastq_finish_date


def delete_existing_bcl2fastq_table_processed(run_process_obj, experiment_name):
    """
    Description:
        The function check if the exists tables created for bcl2fatsq state are already defined.
        If so they are deleted to avoid wrong duplication information
    Input:
        run_process_obj     # run process instance
        experiment_name     # Experiment name
    Return:
        True
    """
    logger = logging.getLogger(__name__)
    logger.debug(
        "%s : Starting function delete_existing_bcl2fastq_table_processed",
        experiment_name,
    )
    if RawDemuxStats.objects.filter(runprocess_id=run_process_obj).exists():
        raw_demux_objs = RawDemuxStats.objects.filter(runprocess_id=run_process_obj)
        for raw_demux_obj in raw_demux_objs:
            raw_demux_obj.delete()
        logger.info("%s : Deleted RawDemuxStats tables ", experiment_name)

    if StatsFlSummary.objects.filter(runprocess_id=run_process_obj).exists():
        stats_fl_summary_objs = StatsFlSummary.objects.filter(
            runprocess_id=run_process_obj
        )
        for stats_fl_summary_obj in stats_fl_summary_objs:
            stats_fl_summary_obj.delete()
        logger.info("%s : Deleted StatsFlSummary tables ", experiment_name)

    if StatsLaneSummary.objects.filter(runprocess_id=run_process_obj).exists():
        stats_lane_summary_objs = StatsLaneSummary.objects.filter(
            runprocess_id=run_process_obj
        )
        for stats_lane_summary_obj in stats_lane_summary_objs:
            stats_lane_summary_obj.delete()
        logger.info("%s : Deleted StatsLaneSummary tables ", experiment_name)
    if SamplesInProject.objects.filter(run_process_id=run_process_obj).exists():
        sample_in_project_objs = SamplesInProject.objects.filter(
            run_process_id=run_process_obj
        )
        for sample_in_project_obj in sample_in_project_objs:
            sample_in_project_obj.delete()
        logger.info("%s : Deleted SamplesInProject tables ", experiment_name)
    logger.debug(
        "%s : End function delete_existing_bcl2fastq_table_processed", experiment_name
    )
    return True


def get_demultiplexing_files(conn, run_folder, experiment_name):
    """
    Description:
        The function check if the folder of demultiplexing files exists
    Input:
        conn                # samba connection instance
        run_folder          # run folder on the remote server
        experiment_name     # Experiment name
    Constants:
        STATS_FILE_PATH
        CONVERSION_STATS_FILE
    Functions:
        get_samba_application_shared_folder
        get_samba_shared_folder
        fetch_remote_file
    Return:
        demux_files
    """
    logger = logging.getLogger(__name__)
    logger.debug("%s : Starting function get_demultiplexing_files", experiment_name)
    statistics_folder = os.path.join(
        get_samba_application_shared_folder(), run_folder, STATS_FILE_PATH
    )
    try:
        conn.listPath(get_samba_shared_folder(), statistics_folder)
    except Exception:
        string_message = (
            experiment_name
            + " : Unable to fetch folder demultiplexing at  "
            + statistics_folder
        )
        logging_errors(string_message, True, False)
        logger.debug(
            "%s : End function get_demultiplexing_files with exception", experiment_name
        )
        return {"ERROR": 29}
    # conversion stats file
    l_conversion_stats = os.path.join(RUN_TEMP_DIRECTORY, CONVERSION_STATS_FILE)
    s_conversion_stats = os.path.join(statistics_folder, CONVERSION_STATS_FILE)
    # demultiplexion stats file
    l_demux_stats = os.path.join(RUN_TEMP_DIRECTORY, DEMULTIPLEXION_STATS_FILE)
    s_demux_stats = os.path.join(statistics_folder, DEMULTIPLEXION_STATS_FILE)
    demux_files = {}
    try:
        demux_files["conversion_stats"] = fetch_remote_file(
            conn, run_folder, s_conversion_stats, l_conversion_stats
        )
        logger.info("%s : Sucessfully fetch of ConversionStats file", experiment_name)
    except Exception:
        string_message = (
            experiment_name
            + " : cannot copy or fetch the "
            + CONVERSION_STATS_FILE
            + " file"
        )
        logging_errors(string_message, True, False)
        logger.debug(
            "%s : End function get_demultiplexing_files with execption", experiment_name
        )
        return {"ERROR": 31}
    try:
        demux_files["demux_stats"] = fetch_remote_file(
            conn, run_folder, s_demux_stats, l_demux_stats
        )
        logger.info("%s : Sucessfully fetch of demultiplexing file", experiment_name)
    except Exception:
        string_message = (
            experiment_name
            + " : cannot copy or fetch the "
            + DEMULTIPLEXION_STATS_FILE
            + " file"
        )
        logging_errors(string_message, True, False)
        os.remove(l_conversion_stats)
        logger.info("%s : deleted %s file", experiment_name, l_conversion_stats)
        logger.debug(
            "%s : End function get_demultiplexing_files with execption", experiment_name
        )
        return {"ERROR": 31}
    logger.debug("%s : End function get_demultiplexing_files ", experiment_name)
    return demux_files


def parsing_demux_and_conversion_files(demux_files, number_of_lanes, experiment_name):
    """
    Description:
        The function will parse demultiplexing and the conversion files that are created during the bcl2fastq process
    Input:
        demux_files             # local copy of the demultiplexting  and conversion files
        conversion_file         # local copy of the conversion file
        experiment_name         # Experiment name
    Return:
        parsed_result
    """
    logger = logging.getLogger(__name__)
    logger.debug(
        "%s : Starting function parsing_demux_and_conversion_files", experiment_name
    )

    parsed_result = {}
    projects = []
    demux_parse = ET.parse(demux_files["demux_stats"])
    root = demux_parse.getroot()

    logger.info("%s : Start parsing demux file", experiment_name)
    for child in root.iter("Project"):
        projects.append(child.attrib["name"])
    total_samples = 0

    for i in range(len(projects)):
        p_temp = root[0][i]
        barcodeCount = []
        p_b_count, one_mismatch_count = [], []
        project_parsed_information = {}

        samples = p_temp.findall("Sample")
        sample_all_index = len(samples) - 1

        for c in p_temp[sample_all_index].iter("BarcodeCount"):
            barcodeCount.append(c.text)

        for c in p_temp[sample_all_index].iter("PerfectBarcodeCount"):
            p_b_count.append(c.text)

        # Fill with zeroes if not PerfectBarcodeCount is written in file
        if len(barcodeCount) != len(p_b_count):
            for i_bar in range(len(barcodeCount)):
                if barcodeCount[i_bar] == "0":
                    p_b_count.insert(i_bar, "0")

        # look for One mismatch barcode
        if p_temp[sample_all_index].find("OneMismatchBarcodeCount") is None:
            for fill in range(number_of_lanes):
                one_mismatch_count.append("NaN")
        else:
            for c in p_temp[sample_all_index].iter("OneMismatchBarcodeCount"):
                one_mismatch_count.append(c.text)

        # one_mismatch_count.append(one_m_count)
        project_parsed_information["BarcodeCount"] = barcodeCount
        project_parsed_information["PerfectBarcodeCount"] = p_b_count
        project_parsed_information["sampleNumber"] = len(samples) - 1
        project_parsed_information["OneMismatchBarcodeCount"] = one_mismatch_count

        parsed_result[projects[i]] = project_parsed_information

        if projects[i] != "default" and projects[i] != "all":
            total_samples += len(samples) - 1
        logger.info(
            "%s : Completed parsing from demux file for project %s",
            experiment_name,
            projects[i],
        )
    # overwrite the value for total samples
    logger.info(
        "%s : Completed parsed information for all projects in stats files",
        experiment_name,
    )
    parsed_result["all"]["sampleNumber"] = total_samples

    conversion_stat = ET.parse(demux_files["conversion_stats"])
    root_conv = conversion_stat.getroot()
    projects = []
    logger.info("%s : Starting parsing for conversion file", experiment_name)
    for child in root_conv.iter("Project"):
        projects.append(child.attrib["name"])
    for i in range(len(projects)):
        p_temp = root_conv[0][i]
        # samples=p_temp.findall('Sample')
        # sample_all_index=len(samples)-1

        list_raw_yield, list_raw_yield_q30 = [], []
        list_raw_qualityscore, list_pf_yield = [], []
        list_pf_yield_q30, list_pf_qualityscore = [], []
        for sample_index in range(len(p_temp.findall("Sample"))):
            if p_temp[sample_index].attrib["name"] != "all":
                continue
            for lane_index in range(len(p_temp[sample_index][0].findall("Lane"))):
                raw_yield_value = 0
                raw_yield_q30_value = 0
                raw_quality_value = 0
                pf_yield_value = 0
                pf_yield_q30_value = 0
                pf_quality_value = 0

                for tile_index in range(
                    len(p_temp[sample_index][0][lane_index].findall("Tile"))
                ):
                    # get the yield value for RAW and for read 1 and 2
                    try:
                        for c in p_temp[sample_index][0][lane_index][tile_index][
                            0
                        ].iter("Yield"):
                            raw_yield_value += int(c.text)
                    except Exception:
                        pass

                    # get the yield Q30 value for RAW  and for read 1 and 2
                    try:
                        for c in p_temp[sample_index][0][lane_index][tile_index][
                            0
                        ].iter("YieldQ30"):
                            raw_yield_q30_value += int(c.text)
                    except Exception:
                        pass
                    try:
                        for c in p_temp[sample_index][0][lane_index][tile_index][
                            0
                        ].iter("QualityScoreSum"):
                            raw_quality_value += int(c.text)
                    except Exception:
                        pass
                    # get the yield value for PF and for read 1 and 2
                    try:
                        for c in p_temp[sample_index][0][lane_index][tile_index][
                            1
                        ].iter("Yield"):
                            pf_yield_value += int(c.text)
                    except Exception:
                        pass
                    # get the yield Q30 value for PF and for read 1 and 2
                    try:
                        for c in p_temp[sample_index][0][lane_index][tile_index][
                            1
                        ].iter("YieldQ30"):
                            pf_yield_q30_value += int(c.text)
                    except Exception:
                        pass
                    try:
                        for c in p_temp[sample_index][0][lane_index][tile_index][
                            1
                        ].iter("QualityScoreSum"):
                            pf_quality_value += int(c.text)
                    except Exception:
                        pass

                list_raw_yield.append(str(raw_yield_value))
                list_raw_yield_q30.append(str(raw_yield_q30_value))
                list_raw_qualityscore.append(str(raw_quality_value))
                list_pf_yield.append(str(pf_yield_value))
                list_pf_yield_q30.append(str(pf_yield_q30_value))
                list_pf_qualityscore.append(str(pf_quality_value))

        parsed_result[projects[i]]["RAW_Yield"] = list_raw_yield
        parsed_result[projects[i]]["RAW_YieldQ30"] = list_raw_yield_q30
        parsed_result[projects[i]]["RAW_QualityScore"] = list_raw_qualityscore
        parsed_result[projects[i]]["PF_Yield"] = list_pf_yield
        parsed_result[projects[i]]["PF_YieldQ30"] = list_pf_yield_q30
        parsed_result[projects[i]]["PF_QualityScore"] = list_pf_qualityscore
        logger.info(
            "%s : Completed parsing for conversion stats for project %s",
            experiment_name,
            projects[i],
        )

    unknow_lanes = []
    unknow_barcode_start_index = len(projects)
    counter = 0
    logger.info("%s : Collecting the Top Unknow Barcodes", experiment_name)
    for un_child in root_conv.iter("TopUnknownBarcodes"):
        un_index = unknow_barcode_start_index + counter
        p_temp = root_conv[0][un_index][0]
        unknow_barcode_lines = p_temp.findall("Barcode")
        unknow_bc_count = []
        for lanes in unknow_barcode_lines:
            unknow_bc_count.append(lanes.attrib)

        unknow_lanes.append(unknow_bc_count)
        counter += 1

    if len(unknow_lanes) != number_of_lanes:
        for index_bar in range(len(barcodeCount)):
            if barcodeCount[index_bar] == "0":
                empty_data = {"count": "0", "sequence": "Not Applicable"}
                fill_data = [empty_data] * 10
                unknow_lanes.insert(index_bar, fill_data)
    parsed_result["TopUnknownBarcodes"] = unknow_lanes
    logger.debug(
        "%s : End function parsing_demux_and_conversion_files", experiment_name
    )

    return parsed_result


def parsing_demux_sample_project(demux_files, number_of_lanes, experiment_name):
    """
    Description:
        The function will parse demultiplexing and the conversion files
        to get the samples inside the project
    Input:
        demux_files         # local copy of the demultiplexting  and conversion files
        number_of_lanes     # numer of lane to be fetched, based on the sequencer type
        experiment_name     # Experiment name
    Return:
        parsed_result
    """
    logger = logging.getLogger(__name__)
    logger.debug("%s : Starting function parsing_demux_sample_project", experiment_name)

    parsed_result = {}
    demux_stat = ET.parse(demux_files["demux_stats"])
    root = demux_stat.getroot()
    projects = []

    logger.info(
        "%s : Starting parsing DemultiplexingStats.XML for getting Sample information",
        experiment_name,
    )
    for child in root.iter("Project"):
        projects.append(child.attrib["name"])

    for i in range(len(projects)):
        if projects[i] == "default" or projects[i] == "all":
            continue
        p_temp = root[0][i]
        samples = p_temp.findall("Sample")
        sample_dict = {}
        for index in range(len(samples)):
            sample_name = samples[index].attrib["name"]
            if sample_name == "all":
                continue
            barcodeCount, perfectBarcodeCount = 0, 0

            sample_stats = {}
            sample_stats["barcodeName"] = samples[index].find("Barcode").attrib["name"]

            for bar_count in p_temp[index][0].iter("BarcodeCount"):
                barcodeCount += int(bar_count.text)

            for p_bar_count in p_temp[index][0].iter("PerfectBarcodeCount"):
                perfectBarcodeCount += int(p_bar_count.text)

            sample_stats["BarcodeCount"] = barcodeCount
            sample_stats["PerfectBarcodeCount"] = perfectBarcodeCount
            sample_dict[sample_name] = sample_stats

            parsed_result[projects[i]] = sample_dict
    logger.info(
        "%s : Complete parsing from demux file for sample and for project %s",
        experiment_name,
        projects[i],
    )

    conversion_stat = ET.parse(demux_files["conversion_stats"])
    root_conv = conversion_stat.getroot()
    projects = []
    logger.info("%s : Starting conversion for conversion file", experiment_name)
    for child in root_conv.iter("Project"):
        projects.append(child.attrib["name"])
    for i in range(len(projects)):
        if projects[i] == "default" or projects[i] == "all":
            continue
        p_temp = root_conv[0][i]
        samples = p_temp.findall("Sample")

        for s_index in range(len(samples)):
            sample_name = samples[s_index].attrib["name"]
            if sample_name == "all":
                continue
            raw_yield_value = 0
            raw_yield_q30_value = 0
            raw_quality_value = 0
            pf_yield_value = 0
            pf_yield_q30_value = 0
            pf_quality_value = 0

            for l_index in range(number_of_lanes):
                tiles_index = len(p_temp[s_index][0][l_index].findall("Tile"))
                for t_index in range(tiles_index):
                    # get the yield value for RAW and for read 1 and 2
                    for c in p_temp[s_index][0][l_index][t_index][0].iter("Yield"):
                        raw_yield_value += int(c.text)
                        # get the yield Q30 value for RAW  and for read 1 and 2
                    for c in p_temp[s_index][0][l_index][t_index][0].iter("YieldQ30"):
                        raw_yield_q30_value += int(c.text)
                    for c in p_temp[s_index][0][l_index][t_index][0].iter(
                        "QualityScoreSum"
                    ):
                        raw_quality_value += int(c.text)
                    # get the yield value for PF and for read 1 and 2
                    for c in p_temp[s_index][0][l_index][t_index][1].iter("Yield"):
                        pf_yield_value += int(c.text)
                    # get the yield Q30 value for PF and for read 1 and 2
                    for c in p_temp[s_index][0][l_index][t_index][1].iter("YieldQ30"):
                        pf_yield_q30_value += int(c.text)
                    for c in p_temp[s_index][0][l_index][t_index][1].iter(
                        "QualityScoreSum"
                    ):
                        pf_quality_value += int(c.text)

            parsed_result[projects[i]][sample_name]["RAW_Yield"] = raw_yield_value
            parsed_result[projects[i]][sample_name][
                "RAW_YieldQ30"
            ] = raw_yield_q30_value
            parsed_result[projects[i]][sample_name][
                "RAW_QualityScore"
            ] = raw_quality_value
            parsed_result[projects[i]][sample_name]["PF_Yield"] = pf_yield_value
            parsed_result[projects[i]][sample_name]["PF_YieldQ30"] = pf_yield_q30_value
            parsed_result[projects[i]][sample_name][
                "PF_QualityScore"
            ] = pf_quality_value
        logger.info(
            "%s : Completed parsing for xml stats for project %s",
            experiment_name,
            projects[i],
        )

    logger.debug("%s : End function parsing_demux_sample_project", experiment_name)

    return parsed_result


def process_and_store_fl_summary_data(
    parsed_data, run_process_obj, number_of_lanes, experiment_name
):
    """
    Description:
        The function manages the parsed data to save them in StatsFlSummary
    Input:
        parsed_data         # dictionnary with parsed data
        run_process_obj     # runProcess object
        number_of_lanes     # number of lanes handled by Sequencer
        experiment_name     # Experiment name
    Return:
        None
    """
    logger = logging.getLogger(__name__)
    logger.debug(
        "%s : Starting function process_and_store_fl_summary_data", experiment_name
    )
    M_BASE = 1.004361 / 1000000
    run_param_obj = RunningParameters.objects.get(run_name_id=run_process_obj)
    number_of_lanes = int(run_param_obj.get_number_of_lanes())

    for project in parsed_data.keys():
        project_flowcell = {}
        if project == "TopUnknownBarcodes":
            continue

        logger.info(
            "%s : Start processing flow Summary for project %s",
            experiment_name,
            project,
        )
        flow_raw_cluster, flow_pf_cluster, flow_yield_mb = 0, 0, 0

        for fl_item in range(number_of_lanes):
            # make the calculation for Flowcell

            flow_raw_cluster += int(parsed_data[project]["BarcodeCount"][fl_item])
            flow_pf_cluster += int(parsed_data[project]["PerfectBarcodeCount"][fl_item])
            flow_yield_mb += float(parsed_data[project]["PF_Yield"][fl_item]) * M_BASE
        logger.debug(
            "%s : flow_pf_cluster value is %s", experiment_name, flow_pf_cluster
        )
        project_flowcell["flowRawCluster"] = "{0:,}".format(flow_raw_cluster)
        project_flowcell["flowPfCluster"] = "{0:,}".format(flow_pf_cluster)
        project_flowcell["flowYieldMb"] = "{0:,}".format(round(flow_yield_mb))
        project_flowcell["sampleNumber"] = parsed_data[project]["sampleNumber"]

        if project == "all" or project == "default":
            project_flowcell["project_id"] = None
            project_flowcell["defaultAll"] = project
        else:
            project_flowcell["project_id"] = Projects.objects.filter(
                project_name__iexact=project
            ).last()
            project_flowcell["defaultAll"] = None
        project_flowcell["runprocess_id"] = run_process_obj
        # store in database

        logger.info(
            "%s : End processing flow Summary for project %s", experiment_name, project
        )
        StatsFlSummary.objects.create_fl_summary(project_flowcell)

    logger.debug("%s : End function process_and_store_fl_summary_data", experiment_name)
    return


def process_and_store_lane_summary_data(
    parsed_data, run_process_obj, number_of_lanes, experiment_name
):
    """
    Description:
        The function manages the parsed data to save them in StatsLaneSummary
    Input:
        parsed_data         # dictionnary with parsed data
        run_process_obj     # runProcess object
        number_of_lanes     # number of lanes handled by Sequencer
        experiment_name     # Experiment name
    Return:
        None
    """
    logger = logging.getLogger(__name__)
    logger.debug(
        "%s : Starting function process_and_store_lane_summary_data", experiment_name
    )

    M_BASE = 1.004361 / 1000000
    total_cluster_lane = parsed_data["all"]["PerfectBarcodeCount"]
    for project in parsed_data.keys():
        if project == "TopUnknownBarcodes":
            continue
        logger.info(
            "%s : Start processing Lane Summary for project %s",
            experiment_name,
            project,
        )
        for i in range(number_of_lanes):
            # get the lane information
            project_lane = {}
            project_lane["lane"] = str(i + 1)
            pf_cluster_int = int(parsed_data[project]["PerfectBarcodeCount"][i])

            project_lane["pfCluster"] = "{0:,}".format(pf_cluster_int)
            try:
                project_lane["perfectBarcode"] = format(
                    int(parsed_data[project]["PerfectBarcodeCount"][i])
                    * 100
                    / int(parsed_data[project]["BarcodeCount"][i]),
                    ".3f",
                )
            except Exception:
                project_lane["perfectBarcode"] = "0"
            try:
                project_lane["percentLane"] = format(
                    float(int(pf_cluster_int) / int(total_cluster_lane[i])) * 100, ".3f"
                )
            except Exception:
                project_lane["percentLane"] = "0"
            project_lane["oneMismatch"] = parsed_data[project][
                "OneMismatchBarcodeCount"
            ][i]
            project_lane["yieldMb"] = "{0:,}".format(
                round(float(parsed_data[project]["PF_Yield"][i]) * M_BASE)
            )
            try:
                project_lane["biggerQ30"] = format(
                    float(parsed_data[project]["PF_YieldQ30"][i])
                    * 100
                    / float(parsed_data[project]["PF_Yield"][i]),
                    ".3f",
                )
            except Exception:
                project_lane["biggerQ30"] = "0"
            try:
                project_lane["meanQuality"] = format(
                    float(parsed_data[project]["PF_QualityScore"][i])
                    / float(parsed_data[project]["PF_Yield"][i]),
                    ".3f",
                )
            except Exception:
                project_lane["meanQuality"] = "0"
            if project == "all" or project == "default":
                project_lane["project_id"] = None
                project_lane["defaultAll"] = project
            else:
                project_lane["project_id"] = Projects.objects.filter(
                    project_name__exact=project
                ).last()
                project_lane["defaultAll"] = None

            project_lane["runprocess_id"] = run_process_obj

        logger.info(
            "%s : Processed information for project %s", experiment_name, project
        )
        StatsLaneSummary.objects.create_lane_summary(project_lane)
        logger.info(
            "%s : Saved information to StatsLaneSummary for project %s ",
            experiment_name,
            project,
        )
    logger.debug(
        "%s : End function process_and_store_lane_summary_data", experiment_name
    )
    return


def process_and_store_raw_demux_project_data(
    parsed_data, run_process_obj, experiment_name
):
    """
    Description:
        The function manages the parsed data to save them in RawDemuxStats
    Input:
        parsed_data         # dictionnary with parsed data
        run_process_obj     # runProcess object
        experiment_name     # Experiment name
    Functions:
        create_new_empty_project
    Return:
        None
    """
    logger = logging.getLogger(__name__)
    logger.debug(
        "%s : Starting function process_and_store_raw_demux_project_data",
        experiment_name,
    )
    project_list = []
    # check first if all project found in the parsing are already defined.
    if run_process_obj.projects_set.all().exists():
        project_objs = run_process_obj.projects_set.all()
        for project_obj in project_objs:
            project_list.append(project_obj.get_project_name())

    for project in parsed_data.keys():
        if project == "TopUnknownBarcodes" or project == "all" or project == "default":
            continue
        if project not in project_list:
            # create project and link to the run
            # check if project exists yet
            if not Projects.objects.filter(project_name__exact=project).exists():
                new_project_obj = Projects.objects.create_new_empty_project(
                    {"user_id": None, "projectName": project}
                )
            else:
                new_project_obj = Projects.objects.filter(
                    project_name__exact=project
                ).last()
            new_project_obj.runProcess.add(run_process_obj)
            string_message = (
                experiment_name
                + " : Created  project name "
                + project
                + "Because it was not store"
            )
            logging_warnings(string_message, True)

    logger.info("%s : Processing demultiplexing raw project data", experiment_name)
    for project in parsed_data.keys():
        project_raw_data = {}
        # ignore on this step the unknowbarcode parsed information
        if project == "TopUnknownBarcodes":
            continue
        logger.info(
            "%s : Start processing project %s for raw parsed data",
            experiment_name,
            project,
        )
        if project == "all" or project == "default":
            project_raw_data["project_id"] = None
            project_raw_data["defaultAll"] = project
        else:
            project_raw_data["project_id"] = Projects.objects.filter(
                project_name__exact=project
            ).last()
            project_raw_data["defaultAll"] = None

        project_raw_data["rawYield"] = parsed_data[project]["RAW_Yield"]
        project_raw_data["rawYieldQ30"] = parsed_data[project]["RAW_YieldQ30"]
        project_raw_data["rawQuality"] = parsed_data[project]["RAW_QualityScore"]
        project_raw_data["PF_Yield"] = parsed_data[project]["PF_Yield"]
        project_raw_data["PF_YieldQ30"] = parsed_data[project]["PF_YieldQ30"]
        project_raw_data["PF_QualityScore"] = parsed_data[project]["PF_QualityScore"]

        RawDemuxStats.objects.create_stats_run_read(
            project_raw_data, run_process_obj
        )
        logger.info(
            "%s : End processing project %s for raw parsed data",
            experiment_name,
            project,
        )
    logger.info(
        "%s : Completed processing demultiplexing raw project data", experiment_name
    )
    logger.debug(
        "%s : End function process_and_store_raw_demux_project_data", experiment_name
    )
    return


def process_and_store_samples_projects_data(
    parsed_data, run_process_obj, experiment_name
):
    """
    Description:
        The function will parse demultiplexing and the conversion files
        that are created during the bcl2fastq process
    Input:
        parsed_data         # dictionnary with parsed data
        run_process_obj     # runProcess object
        experiment_name     # Experiment name
    Functions:
        get_sample_with_user_owner  # Located at wetlab/utils/sample_sheet_utils.py
    Return:
        None
    """
    logger = logging.getLogger(__name__)
    logger.debug(
        "%s : Starting function process_and_store_samples_projects_data",
        experiment_name,
    )
    # Read sample sheet to get the user id for each smaple
    sample_sheet = run_process_obj.get_sample_file()
    samples_with_user_ids = get_sample_with_user_owner(sample_sheet)
    # get the total number of read per lane
    M_BASE = 1.004361 / 1000000

    for project in parsed_data:
        # find the total number of PerfectBarcodeCount in the procjec to make percent calculations
        logger.info(
            "%s : Starting fetching sample project %s", experiment_name, project
        )
        total_perfect_barcode_count = 0
        for sample in parsed_data[project]:
            total_perfect_barcode_count += parsed_data[project][sample][
                "PerfectBarcodeCount"
            ]

        for sample in parsed_data[project]:
            project_sample_data = {}
            project_sample_data["sampleName"] = sample
            project_sample_data["barcodeName"] = parsed_data[project][sample][
                "barcodeName"
            ]
            perfect_barcode = int(parsed_data[project][sample]["PerfectBarcodeCount"])
            project_sample_data["pfClusters"] = "{0:,}".format(perfect_barcode)
            try:
                project_sample_data["percentInProject"] = format(
                    float(perfect_barcode) * 100 / total_perfect_barcode_count, ".2f"
                )
            except Exception:
                project_sample_data["percentInProject"] = "0"
            project_sample_data["yieldMb"] = "{0:,}".format(
                round(float(parsed_data[project][sample]["PF_Yield"]) * M_BASE)
            )
            if parsed_data[project][sample]["PF_Yield"] > 0:
                bigger_q30 = format(
                    float(parsed_data[project][sample]["PF_YieldQ30"])
                    * 100
                    / float(parsed_data[project][sample]["PF_Yield"]),
                    ".3f",
                )
                mean_quality = format(
                    float(parsed_data[project][sample]["PF_QualityScore"])
                    / float(parsed_data[project][sample]["PF_Yield"]),
                    ".3f",
                )
            else:
                bigger_q30 = 0
                mean_quality = 0

            project_sample_data["qualityQ30"] = bigger_q30
            project_sample_data["meanQuality"] = mean_quality
            project_sample_data["project_id"] = Projects.objects.filter(
                project_name__exact=project
            ).last()
            project_sample_data["runProcess_id"] = run_process_obj
            try:
                project_sample_data["user_id"] = samples_with_user_ids[sample]
            except KeyError:
                raise KeyError(33)

            SamplesInProject.objects.create_sample_project(
                project_sample_data
            )

        logger.info(
            "%s : Collected  sample data for the project %s", experiment_name, project
        )
    logger.debug(
        "%s : End function process_and_store_samples_projects_data", experiment_name
    )
    return


def process_and_store_unknown_barcode_data(
    parsed_data, run_process_obj, number_of_lanes, experiment_name
):
    """
    Description:
        The function manages the parsed data to save them in RawTopUnknowBarcodes
    Input:
        parsed_data         # dictionnary with parsed data
        run_process_obj     # runProcess object
        number_of_lanes     # number of lanes handled by Sequencer
        experiment_name     # Experiment name
    Return:
        None
    """
    logger = logging.getLogger(__name__)
    logger.debug(
        "%s : Starting function process_and_store_unknown_barcode_data", experiment_name
    )
    processed_barcode_data = []
    for project in parsed_data:
        if project == "TopUnknownBarcodes":
            logger.info("%s : Collecting TopUnknownBarcodes data ", experiment_name)
            for un_lane in range(number_of_lanes):
                logger.info("Processing lane %s for TopUnknownBarcodes", un_lane)
                top_number = 1

                for barcode_line in parsed_data[project][un_lane]:
                    unknow_barcode = {}
                    unknow_barcode["runprocess_id"] = run_process_obj
                    unknow_barcode["lane_number"] = str(un_lane + 1)
                    unknow_barcode["top_number"] = str(top_number)
                    unknow_barcode["count"] = "{0:,}".format(int(barcode_line["count"]))
                    unknow_barcode["sequence"] = barcode_line["sequence"]
                    top_number += 1

                RawTopUnknowBarcodes.objects.create_unknow_barcode(
                    unknow_barcode
                )

    logger.debug(
        "%s : End function process_and_store_unknown_barcode_data", experiment_name
    )
    return processed_barcode_data
