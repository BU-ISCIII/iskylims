#!/usr/bin/env python3

import sys, os, re

from iSkyLIMS_wetlab.models import RunProcess, RunStates, RunningParameters

import logging

from iSkyLIMS_wetlab import wetlab_config

from .handling_crontab_common_functions import *
from .handling_crontab_manage_run_states import *
from .generic_functions import *
from django.conf import settings

from .common_run_functions import manage_run_in_processed_run, manage_run_in_processing_bcl2fastq, manage_run_in_processed_bcl2fastq
from .generic_functions import get_userid_list

from .sample_sheet_utils import validate_userid_in_user_iem_file


from django.conf import settings
from django_utils.models import Center


def get_list_processed_runs():
    '''
    Description:
        The function get the run folder id from the Running Parameter
        table. This list will be used to compare agains the folder on
        remote server.
    Return:
        processed_runs
    '''
    logger = logging.getLogger(__name__)
    logger.debug('Starting function get_list_processed_runs')
    processed_runs = []
    r_parameters_objects = RunningParameters.objects.all()

    for r_parameter in r_parameters_objects:
        processed_runs.append(r_parameter.get_run_folder())

    logger.info('run processed list is filled')
    logger.debug('End function get_list_processed_runs')
    return processed_runs


def search_update_new_runs(request_reason):
    '''
    Description:
        The function will check if there are new run folders in the remote
        server.
        Get the runParameter file to identify if run is NextSeq or miSeq
        to execute its dedicate handler process.
    Input:
        request_reason      # define if crontab process request the search or was for testing
    Functions:
        assign_projects_to_run               # located at utils.handling_crontab_common_functions
        assign_used_library_in_run           # located at utils.handling_crontab_common_functions
        copy_sample_sheet_to_remote_folder   # located at utils.handling_crontab_common_functions
        fetch_remote_file                    # located at utils.handling_crontab_common_functions
        open_samba_connection                # located in utils.generic_functions.py
        get_list_processed_runs              # located at this file
        get_new_runs_on_remote_server        # located at utils.generic_functions.py
        get_new_runs_from_remote_server      # located at utils.handling_crontab_common_functions
        get_experiment_name_from_file        # located at utils.generic_functions.py
        get_remote_sample_sheet              # located at utils.handling_crontab_common_functions
        get_samba_application_shared_folder  # located at utils.handling_crontab_common_functions.py
        get_samba_shared_folder              # located at utils.handling_crontab_common_functions.py
        fetch_remote_file                    # located at utils.handling_crontab_common_functions
        parsing_run_info_and_parameter_information  #  located at utils.handling_crontab_common_functions
        store_sample_sheet_if_not_defined_in_run     #  located at utils.handling_crontab_common_functions
    Constants:
        PROCESSED_RUN_FILE
        RUN_TEMP_DIRECTORY
        SAMBA_SHARED_FOLDER_NAME
        SAMPLE_SHEET
    Return:
        new_processed_runs # List with all run successfully processed
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function for searching new runs')
    #processed_run_file = os.path.join( wetlab_config.RUN_TEMP_DIRECTORY, wetlab_config.PROCESSED_RUN_FILE)
    #processed_runs = read_processed_runs_file (processed_run_file)
    processed_runs = get_list_processed_runs()
    #process_run_file_update = False
    new_processed_runs = []
    run_with_error = []
    try:
        conn = open_samba_connection()
        logger.info('Sucessfully  SAMBA connection for search_update_new_runs')
    except Exception:
        string_message = 'Unable to open SAMBA connection for the process search update runs'
        # raising the exception to stop crontab
        logging_errors(string_message, True, True)
        raise

    new_runs = get_new_runs_from_remote_server(processed_runs, conn, get_samba_shared_folder())

    if len(new_runs) > 0:
        for new_run in new_runs :
            l_run_parameter_path = os.path.join(wetlab_config.RUN_TEMP_DIRECTORY, wetlab_config.RUN_PARAMETER_FILE)
            s_run_parameter_path = os.path.join(get_samba_application_shared_folder(), new_run, wetlab_config.RUN_PARAMETER_FILE)
            try:
               l_run_parameter = fetch_remote_file (conn, new_run, s_run_parameter_path, l_run_parameter_path)
               logger.info('%s : Sucessfully fetch of RunParameter file', new_run)
               experiment_name = get_experiment_name_from_file (l_run_parameter)
            except Exception:
                error_message = 'Unable to fetch RunParameter file for folder :' + new_run
                logging_errors(error_message, True, False)
                running_parameters["run_date"] = ""
                run_process_obj = get_run_process_obj_or_create_if_not_exists(running_parameters, experiment_name)
                handling_errors_in_run(experiment_name, '21')
                continue

            if request_reason == 'crontab_request':
                if experiment_name == '' or experiment_name == 'NOT FOUND' or 'test' in experiment_name.lower():
                    if experiment_name == '':
                        string_message = new_run + ' : Experiment name is empty'
                        logging_errors(string_message, False, False)
                    elif experiment_name == 'NOT FOUND':
                        string_message = new_run + ' : Experiment name field was not found in file'
                        logging_errors(string_message, False, False)
                    else:
                        string_message = new_run + ' : Ignoring test folder ' + experiment_name
                        logger.info(string_message)
                    os.remove(l_run_parameter)
                    logger.info(' %s  : Deleted temporary run parameter file', new_run)
                    continue
            else:
                if experiment_name != request_reason:
                    logger.info('ignoring test folder %s', experiment_name)
                    os.remove(l_run_parameter)
                    logger.info(' %s  : Deleted temporary run parameter file', new_run)
                    continue

            logger.debug('%s : Found the experiment name called : %s', new_run, experiment_name)
            exclude_states = ["Error", "Recorded"]
            if RunProcess.objects.filter(runName__exact=experiment_name).exclude(state__runStateName__in=exclude_states).exists():
                # This situation should not occurr. The run_processed file should  have this value. To avoid new iterations with this run
                # we update the run process file with this run and continue  with the next item
                run_state = RunProcess.objects.get(runName__exact = experiment_name).get_state()
                string_message = new_run  + ' :  experiment name  state ' + experiment_name + 'in incorrect state. Run state is ' + run_state
                logging_errors( string_message, False, True)
                logger.info('%s : Deleting temporary runParameter file' , experiment_name)
                os.remove(l_run_parameter)
                logger.info('%s : RunParameter file. Local copy deleted',experiment_name)
                continue

            # Fetch run info
            l_run_info_path = os.path.join(wetlab_config.RUN_TEMP_DIRECTORY, wetlab_config.RUN_INFO)
            s_run_info_path = os.path.join(get_samba_application_shared_folder(), new_run, wetlab_config.RUN_INFO)
            try:
                l_run_info = fetch_remote_file (conn, new_run, s_run_info_path, l_run_info_path)
                logger.info('%s : Sucessfully fetch of RunInfo file', experiment_name)
            except Exception:
                string_message = experiment_name + ' : Unable to fetch the RunInfo file on folder ' + new_run
                logging_errors(string_message, True, False)
                running_parameters["run_date"] = ""
                run_process_obj = get_run_process_obj_or_create_if_not_exists(running_parameters, experiment_name)
                handling_errors_in_run(experiment_name, '20')
                # cleaning up the RunParameter in local temporaty file
                logger.debug ('%s : Deleting RunParameter file', experiment_name)
                os.remove(l_run_parameter)
                logger.debug ('%s : Aborting the process for this run. Continue with the next.', experiment_name)
                continue

            running_parameters = parsing_run_info_and_parameter_information(l_run_info, l_run_parameter, experiment_name)
            logger.info('%s  : Deleting runParameter file', experiment_name)
            os.remove(l_run_parameter)
            logger.info('%s  : Deleting runInfo file', experiment_name)
            os.remove(l_run_info)

            run_process_obj = get_run_process_obj_or_create_if_not_exists(running_parameters, experiment_name)
            sequencer_obj = get_sequencer_obj_or_create_if_no_exists(running_parameters, experiment_name)
            run_process_obj = run_process_obj.set_used_sequencer(sequencer_obj)
            run_process_obj.set_run_date(running_parameters["run_date"])
            logger.info('%s : Sequencer  stored on database', experiment_name)

            if run_process_obj.get_sample_file() == '':
                # Fetch sample Sheet from remote server
                l_sample_sheet_path = get_remote_sample_sheet(conn, new_run, experiment_name)
                if not l_sample_sheet_path:
                    logger.debug('%s : End the process. Waiting more time to get Sample Sheet file', experiment_name)
                    # Delete the previous collected data to create again when file is availabñe
                    run_process_obj.delete()
                    continue
                # check if sampleSheet contains userID in description
                if get_configuration_value("DESCRIPTION_IN_SAMPLE_SHEET_MUST_HAVE_USERNAME") == "TRUE":
                    user_id_list = get_userid_list()
                    users = validate_userid_in_user_iem_file(l_sample_sheet_path, user_id_list)
                    if "user_ids" in users and len(users["user_ids"]) == 0:
                        string_message = experiment_name + ' : Description field does not contains userid'
                        logging_errors(string_message, True, True)
                        handling_errors_in_run(experiment_name, '1')
                        continue
                assign_projects_to_run(run_process_obj, l_sample_sheet_path, experiment_name)
                assign_used_library_in_run (run_process_obj,l_sample_sheet_path, experiment_name)
                store_sample_sheet_if_not_defined_in_run (run_process_obj,l_sample_sheet_path, experiment_name)
            else :
                if wetlab_config.COPY_SAMPLE_SHEET_TO_REMOTE and  'NextSeq' in running_parameters['running_data'][wetlab_config.APPLICATION_NAME_TAG]:
                    sample_sheet = run_process_obj.get_sample_file()
                    sample_sheet_path = os.path.join(settings.MEDIA_ROOT, sample_sheet)
                    run_folder = RunningParameters.objects.get( runName_id__exact = run_process_obj).get_run_folder()
                    try:
                        copy_sample_sheet_to_remote_folder(conn, sample_sheet_path, run_folder,experiment_name)
                    except Exception as e:
                        string_message = experiment_name + ' : Unable to copy Sample Sheet to Remote folder' + new_run
                        logging_errors(string_message, True, True)
                        handling_errors_in_run(experiment_name, '23')
                        logger.debug ('%s : Aborting the process. Exiting with exception', experiment_name)
                        continue   # returning to handle next run folder
            save_run_parameters_data_to_database(running_parameters['running_data'], run_process_obj, experiment_name)
            logger.info('%s : RunParameters information  stored on database', experiment_name)
            run_process_obj.set_run_state('Sample Sent')

    logger.info ('Clossing SAMBA connection')
    conn.close()
    logger.debug ('End function searching new runs. Returning handle runs ' )
    return


def handle_not_completed_run ():
    '''
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
        open_samba_connection   # located in utils.generic_functions.py

        save_new_miseq_run      # located at this file
    Constants:
        PROCESSED_RUN_FILE
        RUN_TEMP_DIRECTORY
        SAMBA_SHARED_FOLDER_NAME
        SAMPLE_SHEET
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function for search_not_completed_run')
    try:
        conn = open_samba_connection()
        logger.info('Sucessfully  SAMBA connection for the process_run_in_recorded_state')
    except Exception as e:
        string_message = 'Unable to open SAMBA connection for the process search update runs'
        # raising the exception to stop crontab
        logging_errors(string_message, True, False)
        logger.debug ('End function for search_not_completed_run')
        raise Exception

    runs_to_handle = {}

    #state_list_be_processed = ['Sample Sent','Processing Run','Processed Run', 'Processing Bcl2fastq',
    #                                'Processed Bcl2fastq', 'Recorded']
    state_list_be_processed = ['Recorded','Sample Sent', 'Processing Run','Processed Run', 'Processing Bcl2fastq', 'Processed Bcl2fastq']
    # get the list for all runs that are not completed
    for state in state_list_be_processed:
        #run_state_obj = RunStates.objects.filter(runStateName__exact = state).last()

        if RunProcess.objects.filter(state__runStateName__exact = state).exists():
            runs_to_handle[state] = []
            runs_in_state_objs = RunProcess.objects.filter(state__runStateName__exact = state)
            for run_in_state_obj in runs_in_state_objs :
                runs_to_handle[state].append(run_in_state_obj)

    for state in runs_to_handle.keys():
        logger.info ('Start processing the run found for state %s', state)
        if state == 'Recorded':
            manage_run_in_recorded_state(conn, runs_to_handle[state])

        elif state == 'Sample Sent':
            manage_run_in_sample_sent_processing_state(conn, runs_to_handle[state])
        elif state == 'Processing Run':
            manage_run_in_sample_sent_processing_state(conn, runs_to_handle[state])
        elif state == 'Processed Run':
            manage_run_in_processed_run_state(conn, runs_to_handle[state])
        elif state == 'Processing Bcl2fastq':
            manage_run_in_processing_bcl2fastq_state(conn, runs_to_handle[state])
        elif state == 'Processed Bcl2fastq':
            manage_run_in_processed_bcl2fastq_state(conn, runs_to_handle[state])
        else:
            for run_obj in runs_to_handle[state]:
                experiment_name = run_obj.get_run_name()
                string_message = experiment_name + ' : Is in state not supported by crontab process'
                logging_errors(string_message, False, True)

    logger.debug ('End function for search_not_completed_run')
    return
