# Generic imports
import datetime
import logging
import os

from django.conf import settings

# Local imports
import wetlab.config
import wetlab.models
import wetlab.utils.common
import wetlab.utils.crontab_process
import wetlab.utils.samplesheet


def get_list_processed_runs():
    """
    Description:
        The function get the run folder id from the Running Parameter
        table. This list will be used to compare agains the folder on
        remote server.
    Return:
        processed_runs
    """
    logger = logging.getLogger(__name__)
    logger.debug("Starting function get_list_processed_runs")
    processed_runs = []
    r_parameters_objs = wetlab.models.RunningParameters.objects.all()

    for r_parameter in r_parameters_objs:
        processed_runs.append(r_parameter.get_run_folder())

    run_objs = wetlab.models.RunProcess.objects.all()
    for run in run_objs:
        if r_parameters_objs.filter(run_name_id=run).exists():
            run_folder = r_parameters_objs.get(run_name_id=run).get_run_folder()
        else:
            run_folder = ""

        if run_folder not in processed_runs:
            processed_runs.append(run.run_name)

    logger.info("run processed list is filled")
    logger.debug("End function get_list_processed_runs")
    return processed_runs


def search_update_new_runs(request_reason):
    """
    Description:
        The function will check if there are new run folders in the remote
        server.
        Get the runParameter file to identify if run is NextSeq or miSeq
        to execute its dedicate handler process.
    Input:
        request_reason      # define if crontab process request the search or was for testing
    Functions:
        assign_projects_to_run
        assign_used_library_in_run
        copy_sample_sheet_to_remote_folder
        fetch_remote_file
        open_samba_connection
        get_list_processed_runs
        get_new_runs_on_remote_server
        get_new_runs_from_remote_server
        get_experiment_name_from_file
        get_remote_sample_sheet
        get_samba_application_shared_folder
        get_samba_shared_folder
        parsing_run_info_and_parameter_information
        store_sample_sheet_if_not_defined_in_run
    Constants:
        PROCESSED_RUN_FILE
        RUN_TEMP_DIRECTORY
        SAMBA_SHARED_FOLDER_NAME
        SAMPLE_SHEET
    Return:
        new_processed_runs # List with all run successfully processed
    """
    logger = logging.getLogger(__name__)
    logger.debug("Starting function for searching new runs")
    processed_runs = get_list_processed_runs()

    try:
        conn = wetlab.utils.common.open_samba_connection()
        logger.info("Sucessfully  SAMBA connection for search_update_new_runs")
    except Exception:
        string_message = (
            "Unable to open SAMBA connection for the process search update runs"
        )
        # raising the exception to stop crontab
        wetlab.utils.common.logging_errors(string_message, True, True)
        raise

    new_runs = wetlab.utils.crontab_process.get_new_runs_from_remote_server(
        processed_runs, conn, wetlab.utils.crontab_process.get_samba_shared_folder()
    )
    if len(new_runs) > 0:
        for new_run in new_runs:
            l_run_parameter_path = os.path.join(
                wetlab.config.RUN_TEMP_DIRECTORY, wetlab.config.RUN_PARAMETER_FILE
            )
            s_run_parameter_path = os.path.join(
                wetlab.utils.crontab_process.get_samba_application_shared_folder(),
                new_run,
                wetlab.config.RUN_PARAMETER_FILE,
            )
            try:
                l_run_parameter = wetlab.utils.crontab_process.fetch_remote_file(
                    conn, new_run, s_run_parameter_path, l_run_parameter_path
                )
                logger.info("%s : Sucessfully fetch of RunParameter file", new_run)
                experiment_name = wetlab.utils.common.get_experiment_name_from_file(
                    l_run_parameter
                )
            except Exception:
                error_message = (
                    "Unable to fetch RunParameter file for folder :" + new_run
                )
                wetlab.utils.common.logging_errors(error_message, True, False)
                experiment_name = "Experiment name NOT FOUND"
                wetlab.utils.common.logging_errors(error_message, True, False)
                # check the run folder creation date to allow more time before
                # setting the run to error
                f_created_date = wetlab.utils.common.get_samba_atribute_data(
                    conn,
                    wetlab.utils.crontab_process.get_samba_shared_folder(),
                    new_run,
                    "create_time",
                )
                time_to_check = datetime.datetime.utcfromtimestamp(
                    f_created_date
                ).date()
                max_time_for_run_parameters = (
                    wetlab.models.ConfigSetting.objects.filter(
                        configuration_name__exact="MAXIMUM_TIME_WAIT_RUN_PARAMETERS"
                    )
                    .last()
                    .get_configuration_value()
                )
                if wetlab.utils.crontab_process.waiting_time_expired(
                    time_to_check, max_time_for_run_parameters, experiment_name
                ):
                    # Maximum time for waiting run parameter in folder has expired
                    # we don't have experiment name when there is no run_parameters file.
                    # We used the run folder name instead.
                    run_process_obj = wetlab.utils.crontab_process.get_run_process_obj_or_create_if_not_exists(
                        new_run
                    )
                    wetlab.utils.crontab_process.handling_errors_in_run(new_run, "21")
                else:
                    logger.debug(
                        "%s : RunParameter not in folder run. Allowing more time",
                        experiment_name,
                    )
                logger.debug("%s : Deleting RunParameter file", experiment_name)
                logger.debug(
                    "%s : Aborting the process for this run. Continue with the next.",
                    experiment_name,
                )
                continue

            if request_reason == "crontab_request":
                if (
                    experiment_name == ""
                    or experiment_name == "NOT FOUND"
                    or "test" in experiment_name.lower()
                ):
                    if experiment_name == "":
                        string_message = new_run + " : Experiment name is empty"
                        wetlab.utils.common.logging_errors(string_message, False, False)
                    elif experiment_name == "NOT FOUND":
                        string_message = (
                            new_run + " : Experiment name field was not found in file"
                        )
                        wetlab.utils.common.logging_errors(string_message, False, False)
                    else:
                        string_message = (
                            new_run + " : Ignoring test folder " + experiment_name
                        )
                        logger.info(string_message)
                    os.remove(l_run_parameter)
                    logger.info(" %s  : Deleted temporary run parameter file", new_run)
                    continue
            else:
                if experiment_name != request_reason:
                    logger.info("ignoring test folder %s", experiment_name)
                    os.remove(l_run_parameter)
                    logger.info(" %s  : Deleted temporary run parameter file", new_run)
                    continue

            # Fetch run info
            l_run_info_path = os.path.join(
                wetlab.config.RUN_TEMP_DIRECTORY, wetlab.config.RUN_INFO
            )
            s_run_info_path = os.path.join(
                wetlab.utils.crontab_process.get_samba_application_shared_folder(),
                new_run,
                wetlab.config.RUN_INFO,
            )
            try:
                l_run_info = wetlab.utils.crontab_process.fetch_remote_file(
                    conn, new_run, s_run_info_path, l_run_info_path
                )
                logger.info("%s : Sucessfully fetch of RunInfo file", experiment_name)
            except Exception:
                string_message = (
                    experiment_name
                    + " : Unable to fetch the RunInfo file on folder "
                    + new_run
                )
                wetlab.utils.common.logging_errors(string_message, True, False)
                run_process_obj = wetlab.utils.crontab_process.get_run_process_obj_or_create_if_not_exists(
                    experiment_name
                )
                wetlab.utils.crontab_process.handling_errors_in_run(
                    experiment_name, "20"
                )
                # cleaning up the RunParameter in local temporaty file
                logger.debug("%s : Deleting RunParameter file", experiment_name)
                os.remove(l_run_parameter)
                logger.debug(
                    "%s : Aborting the process for this run. Continue with the next.",
                    experiment_name,
                )
                continue

            logger.debug(
                "%s : Found the experiment name called : %s", new_run, experiment_name
            )
            exclude_states = ["Error", "Recorded"]
            if (
                wetlab.models.RunProcess.objects.filter(run_name__exact=experiment_name)
                .exclude(state__run_state_name__in=exclude_states)
                .exists()
            ):
                # This situation should not occurr.
                # The run_processed file should  have this value. To avoid new iterations with this run
                # we update the run process file with this run and continue  with the next item
                run_state = wetlab.models.RunProcess.objects.get(
                    run_name__exact=experiment_name
                ).get_state()
                string_message = (
                    new_run
                    + " :  experiment name  state "
                    + experiment_name
                    + " in incorrect state. Run state is "
                    + run_state
                )
                wetlab.utils.common.logging_errors(string_message, False, False)
                logger.info(
                    "%s : Deleting temporary runParameter file", experiment_name
                )
                os.remove(l_run_parameter)
                logger.info(
                    "%s : RunParameter file. Local copy deleted", experiment_name
                )
                continue

            running_parameters = (
                wetlab.utils.crontab_process.parsing_run_info_and_parameter_information(
                    l_run_info, l_run_parameter, experiment_name
                )
            )
            logger.info("%s  : Deleting runParameter file", experiment_name)
            os.remove(l_run_parameter)
            logger.info("%s  : Deleting runInfo file", experiment_name)
            os.remove(l_run_info)

            run_process_obj = wetlab.utils.crontab_process.get_run_process_obj_or_create_if_not_exists(
                experiment_name
            )
            sequencer_obj = (
                wetlab.utils.crontab_process.get_sequencer_obj_or_create_if_no_exists(
                    running_parameters, experiment_name
                )
            )
            run_process_obj = run_process_obj.set_used_sequencer(sequencer_obj)
            if isinstance(running_parameters["run_date"], datetime.datetime):
                run_process_obj.set_run_date(running_parameters["run_date"])
            logger.info("%s : Sequencer  stored on database", experiment_name)

            if run_process_obj.get_sample_file() == "":
                # Fetch sample Sheet from remote server
                l_sample_sheet_path = (
                    wetlab.utils.crontab_process.get_remote_sample_sheet(
                        conn, new_run, experiment_name
                    )
                )
                if not l_sample_sheet_path:
                    logger.debug(
                        "%s : End the process. Waiting more time to get Sample Sheet file",
                        experiment_name,
                    )
                    # Delete the previous collected data to create again when file is availabñe
                    run_process_obj.delete()
                    continue
                # check if sampleSheet contains userID in description
                if (
                    wetlab.utils.common.get_configuration_value(
                        "DESCRIPTION_IN_SAMPLE_SHEET_MUST_HAVE_USERNAME"
                    )
                    == "TRUE"
                ):
                    user_id_list = wetlab.utils.common.get_userid_list()
                    file_read = wetlab.utils.samplesheet.read_user_iem_file(
                        l_sample_sheet_path
                    )
                    users = wetlab.utils.samplesheet.validate_userid_in_user_iem_file(
                        file_read, user_id_list
                    )

                    if "ERROR" in users:
                        string_message = (
                            experiment_name
                            + " : Description field does not contains userid or userid is not defined in iskylims."
                        )
                        logger.info(string_message)
                        wetlab.utils.crontab_process.handling_errors_in_run(
                            experiment_name, "1"
                        )
                        continue

                wetlab.utils.crontab_process.assign_projects_to_run(
                    run_process_obj, l_sample_sheet_path, experiment_name
                )
                wetlab.utils.crontab_process.assign_used_library_in_run(
                    run_process_obj, l_sample_sheet_path, experiment_name
                )
                wetlab.utils.crontab_process.store_sample_sheet_if_not_defined_in_run(
                    run_process_obj, l_sample_sheet_path, experiment_name
                )
            else:
                if (
                    wetlab.utils.common.get_configuration_value(
                        "COPY_SAMPLE_SHEET_TO_REMOTE"
                    )
                    == "TRUE"
                    and "NextSeq"
                    in running_parameters["running_data"][
                        wetlab.config.APPLICATION_NAME_TAG
                    ]
                ):
                    sample_sheet = run_process_obj.get_sample_file()
                    sample_sheet_path = os.path.join(settings.MEDIA_ROOT, sample_sheet)
                    try:
                        wetlab.utils.crontab_process.copy_sample_sheet_to_remote_folder(
                            conn, sample_sheet_path, new_run, experiment_name
                        )
                    except Exception:
                        string_message = (
                            experiment_name
                            + " : Unable to copy Sample Sheet to Remote folder"
                            + new_run
                        )
                        wetlab.utils.common.logging_errors(string_message, True, False)
                        wetlab.utils.crontab_process.handling_errors_in_run(
                            experiment_name, "23"
                        )
                        logger.debug(
                            "%s : Aborting the process. Exiting with exception",
                            experiment_name,
                        )
                        continue  # returning to handle next run folder
            wetlab.utils.crontab_process.save_run_parameters_data_to_database(
                running_parameters["running_data"], run_process_obj, experiment_name
            )
            logger.info(
                "%s : RunParameters information  stored on database", experiment_name
            )
            run_process_obj.set_run_state("Sample Sent")

    logger.info("Clossing SAMBA connection")
    conn.close()
    logger.debug("End function searching new runs. Returning handle runs ")
    return


def handle_not_completed_run():
    """
    Description:
        The function will search for run that are not completed yet.
        It will use a dictionary to group the runs in their state.
        Then it with the  check if there are new run folder in the remote
        server.
        Get the runParameter file to identify if run is NextSeq or miSeq
        to execute its dedicate handler process.

    Input:
        logger # log object for logging
    Functions:
        open_samba_connection   # located in utils.common.py

        save_new_miseq_run      # located at this file
    Constants:
        PROCESSED_RUN_FILE
        RUN_TEMP_DIRECTORY
        SAMBA_SHARED_FOLDER_NAME
        SAMPLE_SHEET
    """
    logger = logging.getLogger(__name__)
    logger.debug("Starting function for search_not_completed_run")
    try:
        conn = wetlab.utils.common.open_samba_connection()
        logger.info(
            "Sucessfully  SAMBA connection for the process_run_in_recorded_state"
        )
    except Exception:
        string_message = (
            "Unable to open SAMBA connection for the process search update runs"
        )
        # raising the exception to stop crontab
        wetlab.utils.common.logging_errors(string_message, True, False)
        logger.debug("End function for search_not_completed_run")
        raise Exception

    runs_to_handle = {}

    # state_list_be_processed = ['Sample Sent','Processing Run','Processed Run', 'Processing Bcl2fastq',
    #                                'Processed Bcl2fastq', 'Recorded']
    state_list_be_processed = [
        "Recorded",
        "Sample Sent",
        "Processing Run",
        "Processed Run",
        "Processing Bcl2fastq",
        "Processed Bcl2fastq",
    ]
    # get the list for all runs that are not completed
    for state in state_list_be_processed:
        # run_state_obj = wetlab.models.RunStates.objects.filter(runStateName__exact = state).last()

        if wetlab.models.RunProcess.objects.filter(
            state__run_state_name__exact=state
        ).exists():
            runs_to_handle[state] = []
            runs_in_state_objs = wetlab.models.RunProcess.objects.filter(
                state__run_state_name__exact=state
            )
            for run_in_state_obj in runs_in_state_objs:
                runs_to_handle[state].append(run_in_state_obj)

    for state in runs_to_handle.keys():
        logger.info("Start processing the run found for state %s", state)
        if state == "Recorded":
            manage_run_in_recorded_state(conn, runs_to_handle[state])

        elif state == "Sample Sent":
            manage_run_in_sample_sent_processing_state(conn, runs_to_handle[state])
        elif state == "Processing Run":
            manage_run_in_sample_sent_processing_state(conn, runs_to_handle[state])
        elif state == "Processed Run":
            manage_run_in_processed_run_state(conn, runs_to_handle[state])
        elif state == "Processing Bcl2fastq":
            manage_run_in_processing_bcl2fastq_state(conn, runs_to_handle[state])
        elif state == "Processed Bcl2fastq":
            manage_run_in_processed_bcl2fastq_state(conn, runs_to_handle[state])
        else:
            for run_obj in runs_to_handle[state]:
                experiment_name = run_obj.get_run_name()
                string_message = (
                    experiment_name + " : Is in state not supported by crontab process"
                )
                wetlab.utils.common.logging_errors(string_message, False, False)

    logger.debug("End function for search_not_completed_run")
    return


def manage_run_in_recorded_state(conn, run_process_objs):
    """
    Description:
        The funtion get the runs in recorded state. If the RunningParameter table does not exists
        the run is not handled. Handle by the search_not_completed_run function.
        The function get the sample sheet and extract the project, and library information
        to update the run information.
        If time pass and not Sample Sheet is not available before the Expiration time
        the run state is set to error.
    Input:
        run_process_objs    # list of runProcess objects that are in recorded
    Constants:
        MAXIMUM_TIME_WAIT_SAMPLE_SHEET
        COPY_SAMPLE_SHEET_TO_REMOTE
    Functions:
        get_remote_sample_sheet
        waiting_time_expired
        logging_errors
        handling_errors_in_run
        assign_projects_to_run
        assign_used_library_in_run
        store_sample_sheet_if_not_defined_in_run
        copy_sample_sheet_to_remote_folder
    Return:
        None
    """
    logger = logging.getLogger(__name__)
    logger.debug(" Starting function manage_run_in_recorded_state")
    for run_process_obj in run_process_objs:
        experiment_name = run_process_obj.get_run_name()
        if not wetlab.models.RunningParameters.objects.filter(
            run_name_id=run_process_obj
        ).exists():
            logger.info(
                "%s : Ignore this run. Waiting for existing folder on remote server",
                experiment_name,
            )
            continue
        run_folder = (
            wetlab.models.RunningParameters.objects.filter(run_name_id=run_process_obj)
            .last()
            .get_run_folder()
        )
        if run_process_obj.get_sample_file() == "":
            # sample sheet does not included yet in run process. Try to get it from remote server
            l_sample_sheet_path = wetlab.utils.crontab_process.get_remote_sample_sheet(
                conn, run_folder, experiment_name
            )
            if not l_sample_sheet_path:
                maximun_time = (
                    wetlab.models.ConfigSetting.objects.filter(
                        configuration_name__exact="MAXIMUM_TIME_WAIT_SAMPLE_SHEET"
                    )
                    .last()
                    .get_configuration_value()
                )
                time_to_check = (
                    run_process_obj.get_run_generated_date_no_format().date()
                )
                if not wetlab.utils.crontab_process.waiting_time_expired(
                    time_to_check, maximun_time, experiment_name
                ):
                    logger.debug(
                        "%s : End the process. Waiting more time to get Sample Sheet file",
                        experiment_name,
                    )
                    continue
                else:
                    string_message = (
                        experiment_name
                        + " : Expired time for waiting for Sample Sheet file on folder "
                        + run_folder
                    )
                    wetlab.utils.common.logging_errors(string_message, False, False)
                    wetlab.utils.crontab_process.handling_errors_in_run(
                        experiment_name, 19
                    )
                    logger.debug(
                        "%s : Aborting the process. Exceeded the waiting time for fetching ths sample sheet",
                        experiment_name,
                    )
                    continue
            wetlab.utils.crontab_process.assign_projects_to_run(
                run_process_obj, l_sample_sheet_path, experiment_name
            )
            wetlab.utils.crontab_process.assign_used_library_in_run(
                run_process_obj, l_sample_sheet_path, experiment_name
            )
            wetlab.utils.crontab_process.store_sample_sheet_if_not_defined_in_run(
                run_process_obj, l_sample_sheet_path, experiment_name
            )

        if (
            wetlab.config.COPY_SAMPLE_SHEET_TO_REMOTE
            and "NextSeq" in run_process_obj.get_run_platform()
        ):
            sample_sheet_path = run_process_obj.get_sample_file()

            try:
                wetlab.utils.crontab_process.copy_sample_sheet_to_remote_folder(
                    conn, sample_sheet_path, run_folder, experiment_name
                )
            except Exception:
                logger.info(
                    "%s : Aborting process, Unable to copy sammple sheet to remote server",
                    experiment_name,
                )

                continue

        run_process_obj.set_run_state("Sample Sent")
        logger.info("%s  : is now on Sample Sent state", experiment_name)
    logger.debug(" End function manage_run_in_recorded_state")
    return


def manage_run_in_sample_sent_processing_state(conn, run_process_objs):
    """
    Description:
        The funtion get the runs in sample sent state.
    Input:
        conn                # samba connection instance
        run_process_objs    # list of runProcess objects that are in recorded
    Functions:
        check_log_for_run_completions
        waiting_time_expired
        handling_errors_in_run
    Return:
        None
    """
    logger = logging.getLogger(__name__)
    logger.debug(" Starting function manage_run_in_sample_sent_processing_state")
    for run_process_obj in run_process_objs:
        experiment_name = run_process_obj.get_run_name()
        logger.info(
            "%s : Start handling in manage_run_in_sample_sent_processing_state function",
            experiment_name,
        )
        platform = run_process_obj.get_run_platform()
        if platform == "None":
            string_message = (
                experiment_name + " : Used sequencer or the platform is not defined"
            )
            wetlab.utils.common.logging_errors(string_message, False, False)
            wetlab.utils.crontab_process.handling_errors_in_run(experiment_name, 24)
            logger.info(
                "%s ERROR in manage_run_in_sample_sent_processing_state function",
                experiment_name,
            )
            logger.debug(
                "%s End manage_run_in_sample_sent_processing_state function",
                experiment_name,
            )
            continue
        run_folder = (
            wetlab.models.RunningParameters.objects.filter(run_name_id=run_process_obj)
            .last()
            .get_run_folder()
        )
        number_of_cycles = (
            wetlab.models.RunningParameters.objects.filter(run_name_id=run_process_obj)
            .last()
            .get_number_of_cycles()
        )
        (
            run_status,
            run_completion_date,
        ) = wetlab.utils.crontab_process.check_sequencer_run_is_completed(
            conn, run_folder, platform, number_of_cycles, experiment_name
        )

        if run_status == "completed":
            run_process_obj.set_run_state("Processed Run")
            run_process_obj.set_run_completion_date(run_completion_date)
            logger.info("%s changed to Processed Run state", experiment_name)
            logger.debug(
                "%s End manage_run_in_sample_sent_processing_state function",
                experiment_name,
            )
        elif run_status == "cancelled":
            run_process_obj.set_run_state("Processed Run")
            logger.info("%s changed to Processed Run state", experiment_name)
            string_message = experiment_name + "was cancelled on the sequencer"
            wetlab.utils.common.logging_warnings(string_message, True)
            logger.debug(
                "%s End manage_run_in_sample_sent_processing_state function",
                experiment_name,
            )
        elif "ERROR" in run_status:
            if run_status["ERROR"] == 18:
                string_message = (
                    experiment_name
                    + " : Unable to fetch logs files for checking run completion status "
                    + run_folder
                )
            else:
                string_message = (
                    experiment_name
                    + " : platform "
                    + platform
                    + " is not defined in wetlab.config.py file (on PLATFORM_WAY_TO_CHECK_RUN_COMPLETION variable) "
                )
            wetlab.utils.common.logging_errors(string_message, False, False)
            wetlab.utils.crontab_process.handling_errors_in_run(
                experiment_name, run_status["ERROR"]
            )
            logger.debug(
                "%s : End manage_run_in_sample_sent_processing_state function",
                experiment_name,
            )
        else:
            maximun_time = (
                wetlab.models.ConfigSetting.objects.filter(
                    configuration_name__exact="MAXIMUM_TIME_WAIT_RUN_COMPLETION"
                )
                .last()
                .get_configuration_value()
            )
            time_to_check = run_process_obj.get_run_generated_date_no_format().date()
            if not wetlab.utils.crontab_process.waiting_time_expired(
                time_to_check, maximun_time, experiment_name
            ):
                logger.info(
                    "%s : Waiting more time to get Sequencer completion",
                    experiment_name,
                )
                run_process_obj.set_run_state("Processing Run")
                logger.info("%s : changed to Processing Run state", experiment_name)
                logger.debug(
                    "%s  : End manage_run_in_sample_sent_processing_state function",
                    experiment_name,
                )
            else:
                string_message = (
                    experiment_name
                    + " : Expired time for waiting for Sequencer completion file on folder "
                    + run_folder
                )
                wetlab.utils.common.logging_errors(string_message, False, False)
                wetlab.utils.crontab_process.handling_errors_in_run(experiment_name, 9)
                logger.debug(
                    "%s  : End manage_run_in_sample_sent_processing_state function",
                    experiment_name,
                )

    logger.debug("End function manage_run_in_sample_sent_processing_state")
    return


def manage_run_in_processed_run_state(conn, run_process_objs):
    """
    Description:
        The funtion get the runs in processed run state. In this state it will collect run Metrics
        files and store in StatsRunSummary table
    Input:
        conn                # samba connection instance
        run_process_objs    # list of runProcess objects that are in processed_run
    Constants:
        RUN_TEMP_DIRECTORY_PROCESSING
    Functions:
        delete_existing_run_metrics_table_processed
        get_run_metric_files
        handling_errors_in_run
        delete_run_metric_files
        parsing_run_metrics_files
        create_run_metric_graphics
    Return:
        None
    """
    logger = logging.getLogger(__name__)
    logger.debug(" Starting function manage_run_in_processed_run_state")
    for run_process_obj in run_process_objs:
        experiment_name = run_process_obj.get_run_name()
        logger.info(
            "%s : Start handling in manage_run_in_processed_run_state function",
            experiment_name,
        )
        run_folder = (
            wetlab.models.RunningParameters.objects.filter(run_name_id=run_process_obj)
            .last()
            .get_run_folder()
        )
        # delete existing information to avoid having duplicated tables
        wetlab.utils.crontab_process.delete_existing_run_metrics_table_processed(
            run_process_obj, experiment_name
        )
        # Check run_folder time creation
        f_created_date = wetlab.utils.common.get_samba_atribute_data(
                    conn,
                    wetlab.utils.crontab_process.get_samba_shared_folder(),
                    run_folder,
                    "created_time",
                )
        time_to_check = datetime.datetime.utcfromtimestamp(
            f_created_date
        ).date()
        # Check maximum time for waiting run metric files
        max_time_for_run_parameters = (
            wetlab.models.ConfigSetting.objects.filter(
                configuration_name__exact="MAXIMUM_TIME_WAIT_RUN_PARAMETERS"
            )
            .last()
            .get_configuration_value()
        )

        # Get run metric files
        run_metric_files = wetlab.utils.crontab_process.get_run_metric_files(
            conn, run_folder, experiment_name
        )

        if "ERROR" in run_metric_files:
            if wetlab.utils.crontab_process.waiting_time_expired(
                time_to_check, max_time_for_run_parameters, experiment_name
            ):
                string_message = (
                    experiment_name + " : Unable to collect all files for run metrics"
                )
                wetlab.utils.common.logging_errors(string_message, True, False)
                wetlab.utils.crontab_process.handling_errors_in_run(
                    experiment_name, run_metric_files["ERROR"]
                )
                wetlab.utils.crontab_process.delete_run_metric_files(experiment_name)
                logger.debug(
                    "%s : Deleting run metric files", experiment_name
                )
            else:
                logger.debug(
                        "%s : RunParameter not in folder run. Allowing more time",
                        experiment_name,
                    )
            logger.debug(
                    "%s : End manage_run_in_processed_run_state function", experiment_name
                )
            continue
        
        # Parse run metric files
        (
            parsed_run_stats_summary,
            parsed_run_stats_read,
        ) = wetlab.utils.crontab_process.parsing_run_metrics_files(
            wetlab.config.RUN_TEMP_DIRECTORY_PROCESSING,
            run_process_obj,
            experiment_name,
        )

        for run_stat_summary in parsed_run_stats_summary:
            wetlab.models.StatsRunSummary.objects.create_stats_run_summary(
                run_stat_summary, run_process_obj
            )
        logger.info("%s : run metrics summary data saved to database", experiment_name)
        for run_stat_read in parsed_run_stats_read:
            wetlab.models.StatsRunRead.objects.create_stats_run_read(
                run_stat_read, run_process_obj
            )
        logger.info("%s :run metrics read data saved to database", experiment_name)

        # create run graphics
        run_graphics = wetlab.utils.crontab_process.create_run_metric_graphics(
            wetlab.config.RUN_TEMP_DIRECTORY_PROCESSING,
            run_process_obj,
            run_folder,
            experiment_name,
        )
        if "ERROR" in run_graphics:
            string_message = (
                experiment_name + " : Unable to save graphics for run metrics"
            )
            wetlab.utils.common.logging_errors(string_message, True, False)
            wetlab.utils.crontab_process.handling_errors_in_run(
                experiment_name, run_graphics["ERROR"]
            )
            wetlab.utils.crontab_process.delete_existing_run_metrics_table_processed(
                run_process_obj, experiment_name
            )
            wetlab.utils.crontab_process.delete_run_metric_files(experiment_name)
            logger.debug(
                "%s : End manage_run_in_processed_run_state function", experiment_name
            )
            continue
        logger.info(
            "%s : run metrics graphics processed and copied to plot image folder",
            experiment_name,
        )
        # deleting temporary run metrics files
        wetlab.utils.crontab_process.delete_run_metric_files(experiment_name)
        # return the state to Processed Run
        run_process_obj.set_run_state("Processing Bcl2fastq")

    logger.debug(" End function manage_run_in_processed_run_state")
    return


def manage_run_in_processing_bcl2fastq_state(conn, run_process_objs):
    """
    Description:
        The funtion get the runs in bcl2fastq processing run state. In this state it will check if
        bcl2fastq process is completed by checking if report folder exists
    Input:
        conn                # samba connection instance
        run_process_objs    # list of runProcess objects that are in processing_bcl2fastq
    Constants:
        RUN_TEMP_DIRECTORY_PROCESSING
    Functions:
        check_run_metrics_processed
        waiting_time_expired
        check_demultiplexing_folder_exists
    Return:
        None
    """
    logger = logging.getLogger(__name__)
    logger.debug("Starting function manage_run_in_processing_bcl2fastq_state")
    for run_process_obj in run_process_objs:
        experiment_name = run_process_obj.get_run_name()
        logger.info(
            "%s : Start handling in manage_run_in_processing_bcl2fastq_state function",
            experiment_name,
        )
        run_folder = wetlab.models.RunningParameters.objects.get(
            run_name_id=run_process_obj
        ).get_run_folder()

        bcl2fastq_finish_date = (
            wetlab.utils.crontab_process.check_demultiplexing_folder_exists(
                conn, run_folder, experiment_name
            )
        )
        if "ERROR" in bcl2fastq_finish_date:
            if (
                bcl2fastq_finish_date["ERROR"] == 29
                or bcl2fastq_finish_date["ERROR"] == 31
            ):
                maximun_time = (
                    wetlab.models.ConfigSetting.objects.filter(
                        configuration_name__exact="MAXIMUM_TIME_WAIT_TO_RUN_BCL2FASTQ"
                    )
                    .last()
                    .get_configuration_value()
                )
                try:
                    time_to_check = (
                        run_process_obj.get_run_completion_date_no_format().date()
                    )
                except Exception:
                    string_message = (
                        experiment_name
                        + " :  Aborting the process. No Run completion date was defined."
                    )
                    wetlab.utils.common.logging_errors(string_message, True, False)
                    wetlab.utils.crontab_process.handling_errors_in_run(
                        experiment_name, 30
                    )
                    logger.debug(
                        "%s : Aborting the process. No Run completion date was defined",
                        experiment_name,
                    )
                    continue
                if not wetlab.utils.crontab_process.waiting_time_expired(
                    time_to_check, maximun_time, experiment_name
                ):
                    logger.debug(
                        "%s : End the process. Waiting more time to get bcf2fastq stats.",
                        experiment_name,
                    )
                    continue
                else:
                    string_message = (
                        experiment_name
                        + " :  Aborting the process. Exceeded the waiting time for fetching demultiplexing files "
                        + run_folder
                    )
                    wetlab.utils.common.logging_errors(string_message, True, False)
                    wetlab.utils.crontab_process.handling_errors_in_run(
                        experiment_name, bcl2fastq_finish_date["ERROR"]
                    )
                    continue
            else:
                string_message = (
                    experiment_name
                    + " : Aborting the process. Error in bcf2fastq process."
                    + run_folder
                )
                wetlab.utils.common.logging_errors(string_message, True, False)
                wetlab.utils.crontab_process.handling_errors_in_run(
                    experiment_name, bcl2fastq_finish_date["ERROR"]
                )
                continue
        run_process_obj.set_run_bcl2fastq_finished_date(bcl2fastq_finish_date)
        run_process_obj.set_run_state("Processed Bcl2fastq")
        logger.info("%s : Updated to Processed Bcl2Fastq state", experiment_name)
        logger.info(
            "%s : End handling in manage_run_in_processing_bcl2fastq_state function",
            experiment_name,
        )
    logger.debug(" End function manage_run_in_processing_bcl2fastq_state")
    return


def manage_run_in_processed_bcl2fastq_state(conn, run_process_objs):
    """
    Description:
        The funtion get the runs in bcl2fastq processed run state. In this state it will get the
        demultiplexing data generated by bcl2fastq.
    Input:
        conn                # samba connection instance
        run_process_objs    # list of runProcess objects that are in processed_bcl2fastq
    Constants:
        RUN_TEMP_DIRECTORY_PROCESSING
    Functions:
        delete_existing_bcl2fastq_table_processed
        get_demultiplexing_files
        parsing_demux_and_conversion_files
        parsing_demux_sample_project
        process_and_store_raw_demux_project_data
        process_and_store_fl_summary_data
        process_and_store_lane_summary_data
        process_and_store_unknown_barcode_data
        process_and_store_samples_projects_data
        get_run_disk_utilization

    Return:
        None
    """
    logger = logging.getLogger(__name__)
    logger.debug(" Starting function manage_run_in_processed_bcl2fastq_state")
    for run_process_obj in run_process_objs:
        experiment_name = run_process_obj.get_run_name()
        logger.info(
            "%s : Start handling in manage_run_in_processed_bcl2fastq_state function",
            experiment_name,
        )
        run_param_obj = wetlab.models.RunningParameters.objects.get(
            run_name_id=run_process_obj
        )
        run_folder = run_param_obj.get_run_folder()
        # delete existing information to avoid having duplicated tables
        wetlab.utils.crontab_process.delete_existing_bcl2fastq_table_processed(
            run_process_obj, experiment_name
        )
        demux_files = wetlab.utils.crontab_process.get_demultiplexing_files(
            conn, run_folder, experiment_name
        )
        if "ERROR" in demux_files:
            if demux_files["ERROR"] == 29 or demux_files["ERROR"] == 31:
                string_message = (
                    experiment_name
                    + " : Aborting the process. Unable to reach demultiplexing files "
                )
                wetlab.utils.common.logging_errors(string_message, True, False)
                wetlab.utils.crontab_process.handling_errors_in_run(
                    experiment_name, demux_files["ERROR"]
                )
                continue

        number_of_lanes = run_param_obj.get_number_of_lanes()
        try:
            number_of_lanes = int(number_of_lanes)
        except Exception:
            string_message = (
                experiment_name
                + " : Sequencer used in the run has not defined the number of lanes"
            )
            wetlab.utils.common.logging_errors(string_message, True, False)
            wetlab.utils.crontab_process.handling_errors_in_run(experiment_name, 32)
            logger.debug(
                "%s : Aborting the process. Number of lanes not defined on the sequencer",
                experiment_name,
            )
            continue

        # parsing the files to get the xml Stats
        logger.info("%s : Start parsing  demultiplexing files", experiment_name)
        project_parsed_data = (
            wetlab.utils.crontab_process.parsing_demux_and_conversion_files(
                demux_files, number_of_lanes, experiment_name
            )
        )

        logger.info("%s : Start parsing  samples demultiplexing", experiment_name)
        sample_project_parsed_data = (
            wetlab.utils.crontab_process.parsing_demux_sample_project(
                demux_files, number_of_lanes, experiment_name
            )
        )

        try:
            # clean up the fetched files in the local temporary folder
            os.remove(demux_files["demux_stats"])
            os.remove(demux_files["conversion_stats"])
            logger.info(
                "%s : Deleted temporary demultiplexing and conversion files",
                experiment_name,
            )
        except Exception:
            string_message = (
                experiment_name
                + " : Unable to delete the Conversion Stats file"
                + run_folder
            )
            wetlab.utils.common.logging_errors(string_message, True, False)
            logger.debug(
                "%s : Allowing to contine the parssing process", experiment_name
            )

        wetlab.utils.crontab_process.process_and_store_raw_demux_project_data(
            project_parsed_data, run_process_obj, experiment_name
        )

        wetlab.utils.crontab_process.process_and_store_fl_summary_data(
            project_parsed_data, run_process_obj, number_of_lanes, experiment_name
        )

        wetlab.utils.crontab_process.process_and_store_lane_summary_data(
            project_parsed_data, run_process_obj, number_of_lanes, experiment_name
        )

        wetlab.utils.crontab_process.process_and_store_unknown_barcode_data(
            project_parsed_data, run_process_obj, number_of_lanes, experiment_name
        )

        try:
            wetlab.utils.crontab_process.process_and_store_samples_projects_data(
                sample_project_parsed_data, run_process_obj, experiment_name
            )
        except KeyError as key_error:
            string_message = (
                experiment_name
                + " : Error when processing and storing samples in projects."
            )
            wetlab.utils.common.logging_errors(string_message, True, False)
            if key_error.args[0] == 33:
                wetlab.utils.crontab_process.handling_errors_in_run(
                    experiment_name, "33"
                )
            else:
                string_message = (
                    experiment_name
                    + " : Unknown error when processing and storing samples in projects."
                )
                wetlab.utils.common.logging_errors(string_message, True, False)
            logger.debug(
                "%s : End function manage_run_in_processed_bcl2fast2_run with error",
                experiment_name,
            )
            continue
        except Exception:
            string_message = (
                experiment_name
                + " : Unknown error when processing and storing samples in projects."
            )
            wetlab.utils.common.logging_errors(string_message, True, False)
            continue

        # Get the disk space utilization for this run
        try:
            disk_utilization = wetlab.utils.crontab_process.get_run_disk_utilization(
                conn, run_folder, experiment_name
            )
        except Exception:
            string_message = (
                experiment_name + " : Error when fetching the disk utilization"
            )
            wetlab.utils.common.logging_errors(string_message, True, False)
            wetlab.utils.crontab_process.handling_errors_in_run(experiment_name, "17")
            logger.debug(
                "%s : End function manage_run_in_processed_bcl2fast2_run with error",
                experiment_name,
            )
            continue

        run_process_obj.set_used_space(disk_utilization)
        finish_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        run_process_obj.set_run_finish_date(finish_date)
        # Update the run state to completed
        run_process_obj.set_run_state("Completed")
        logger.info("%s : is Completed", experiment_name)
    logger.debug(" End function manage_run_in_processed_bcl2fastq_state")
