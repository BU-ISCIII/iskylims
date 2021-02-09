import os
import logging
from logging.config import fileConfig
from logging.handlers import RotatingFileHandler

from iSkyLIMS_wetlab.models import RunProcess, RunningParameters, Projects, ConfigSetting
from iSkyLIMS_wetlab.wetlab_config import *

from .sample_sheet_utils import get_projects_in_run, get_index_library_name
from .handling_crontab_common_functions import *


def manage_run_in_recorded_state(conn, run_process_objs):
    '''
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
        get_remote_sample_sheet                 # located in utils.handling_crontab_common_functions.py
        waiting_time_expired                    # located in utils.handling_crontab_common_functions.py
        logging_errors                          # located in utils.handling_crontab_common_functions.py
        handling_errors_in_run                  # located in utils.handling_crontab_common_functions.py
        assign_projects_to_run                  # located in utils.handling_crontab_common_functions.py
        assign_used_library_in_run              # located in utils.handling_crontab_common_functions.py
        store_sample_sheet_if_not_defined_in_run    # located in utils.handling_crontab_common_functions.py
        copy_sample_sheet_to_remote_folder      # located in utils.handling_crontab_common_functions.py
    Return:
        None
    '''
    logger = logging.getLogger(__name__)
    logger.debug (' Starting function manage_run_in_recorded_state')
    for run_process_obj in run_process_objs:
        experiment_name = run_process_obj.get_run_name()
        if not  RunningParameters.objects.filter(runName_id = run_process_obj).exists():
            logger.info('%s : Ignore this run. Waiting for existing folder on remote server', experiment_name)
            continue
        run_folder = RunningParameters.objects.filter(runName_id = run_process_obj).last().get_run_folder()
        if run_process_obj.get_sample_file() == '':
            # sample sheet does not included yet in run process. Try to get it from remote server
            l_sample_sheet_path = get_remote_sample_sheet(conn, run_folder ,experiment_name)
            if not l_sample_sheet_path :
                maximun_time = ConfigSetting.objects.filter(configurationName__exact = 'MAXIMUM_TIME_WAIT_SAMPLE_SHEET').last().get_configuration_value()
                if not waiting_time_expired(run_process_obj, maximun_time ,experiment_name):
                    logger.debug ('%s : End the process. Waiting more time to get Sample Sheet file', experiment_name)
                    continue
                else:
                    string_message = experiment_name + ' : Expired time for waiting for Sample Sheet file on folder ' + run_folder
                    logging_errors(string_message, False, True)
                    handling_errors_in_run (experiment_name, 19)
                    logger.debug ('%s : Aborting the process. Exceeded the waiting time for fetching ths sample sheet', experiment_name)
                    continue
            assign_projects_to_run(run_process_obj, l_sample_sheet_path, experiment_name)
            assign_used_library_in_run (run_process_obj,l_sample_sheet_path, experiment_name )
            store_sample_sheet_if_not_defined_in_run (run_process_obj,l_sample_sheet_path, experiment_name)

        if COPY_SAMPLE_SHEET_TO_REMOTE and 'NextSeq' in run_process_obj.get_run_platform() :
            sample_sheet_path = run_process_obj.get_sample_file()

            try:
                copy_sample_sheet_to_remote_folder(conn, sample_sheet_path, run_folder ,experiment_name)
            except:
                logger.info('%s : Aborting process, Unable to copy sammple sheet to remote server', experiment_name)

                continue

        run_process_obj.set_run_state('Sample Sent')
        logger.info('%s  : is now on Sample Sent state' , experiment_name)
    logger.debug (' End function manage_run_in_recorded_state')
    return

def manage_run_in_sample_sent_processing_state(conn, run_process_objs):
    '''
    Description:
        The funtion get the runs in sample sent state.
    Input:
        conn                # samba connection instance
        run_process_objs    # list of runProcess objects that are in recorded
    Functions:
        check_log_for_run_completions           # located in utils.handling_crontab_common_functions.py
        waiting_time_expired                    # located in utils.handling_crontab_common_functions.py
    Return:
        None
    '''
    logger = logging.getLogger(__name__)
    logger.debug (' Starting function manage_run_in_sample_sent_processing_state')
    for run_process_obj in run_process_objs:
        experiment_name = run_process_obj.get_run_name()
        logger.info('%s : Start handling in manage_run_in_sample_sent_processing_state function', experiment_name)
        platform = run_process_obj.get_run_platform()
        if platform == 'None':
            string_message = experiment_name + ' : Used sequencer or the platform is not defined'
            logging_errors(string_message, False, True)
            handling_errors_in_run (experiment_name, 24)
            logger.info('%s ERROR in manage_run_in_sample_sent_processing_state function', experiment_name)
            logger.debug('%s End manage_run_in_sample_sent_processing_state function', experiment_name)
            continue
        run_folder = RunningParameters.objects.filter(runName_id = run_process_obj).last().get_run_folder()
        number_of_cycles =  RunningParameters.objects.filter(runName_id = run_process_obj).last().get_number_of_cycles()
        run_status, run_completion_date = check_sequencer_run_is_completed(conn, run_folder , platform , number_of_cycles, experiment_name)

        if run_status == 'completed':
            run_process_obj.set_run_state('Processed Run')
            run_process_obj.set_run_completion_date(run_completion_date)
            logger.info('%s changed to Processed Run state', experiment_name)
            logger.debug('%s End manage_run_in_sample_sent_processing_state function', experiment_name)
        elif run_status == 'cancelled':
            run_process_obj.set_run_state('Processed Run')
            logger.info('%s changed to Processed Run state', experiment_name)
            string_message = experiment_name + 'was cancelled on the sequencer'
            logging_warnings(string_message, True)
            logger.debug('%s End manage_run_in_sample_sent_processing_state function', experiment_name)
        elif 'ERROR' in run_status:
            if run_status['ERROR'] == 18:
                string_message = experiment_name + ' : Unable to fetch logs files for checking run completion status ' + run_folder
            else :
                string_message = experiment_name + ' : platform ' + platform + ' is not defined in wetlab_config.py file (on PLATFORM_WAY_TO_CHECK_RUN_COMPLETION variable) '
            import pdb; pdb.set_trace()
            logging_errors(string_message, False, True)
            handling_errors_in_run (experiment_name, run_status['ERROR'])
            logger.debug('%s : End manage_run_in_sample_sent_processing_state function', experiment_name)
        else:
            maximun_time = ConfigSetting.objects.filter(configurationName__exact = 'MAXIMUM_TIME_WAIT_RUN_COMPLETION').last().get_configuration_value()
            if not waiting_time_expired(run_process_obj, maximun_time ,experiment_name):
                logger.info ('%s : Waiting more time to get Sequencer completion', experiment_name)
                run_process_obj.set_run_state('Processing Run')
                logger.info('%s : changed to Processing Run state', experiment_name)
                logger.debug('%s  : End manage_run_in_sample_sent_processing_state function', experiment_name)
            else:
                string_message = experiment_name + ' : Expired time for waiting for Sequencer completion file on folder ' + run_folder
                logging_errors(string_message, False, True)
                handling_errors_in_run (experiment_name, 9)
                logger.debug('%s  : End manage_run_in_sample_sent_processing_state function', experiment_name)

    logger.debug (' End function manage_run_in_sample_sent_processing_state')
    return
