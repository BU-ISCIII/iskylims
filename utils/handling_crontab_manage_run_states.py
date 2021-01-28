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
    '''
    logger = logging.getLogger(__name__)
    logger.debug (' Starting function manage_run_in_recorded_state')
    for run_process_obj in run_process_objs:
        experiment_name = run_process_obj.get_run_name()
        run_folder = RunningParameters.objects.filter(runName_id = run_process_obj).last().get_run_folder()
        l_sample_sheet_path = get_remote_sample_sheet(conn, run_folder ,experiment_name)
        if not l_sample_sheet_path :
            maximun_time = ConfigSetting.objects.filter(configurationName__exact = 'MAXIMUM_TIME_WAIT_SAMPLE_SHEET').last().get_configuration_value()
            if not waiting_time_expired(run_process_obj, maximun_time ,experiment_name):
                logger.debug ('%s : End the process. Waiting more time to get Sample Sheet file', experiment_name)
                continue
            else:
                string_message = experiment_name + ' : Expired time for waiting for Sample Sheet file on folder ' + new_run
                logging_errors(string_message, False, True)
                handling_errors_in_run (experiment_name, 19)
                logger.debug ('%s : Aborting the process. Exceeded the waiting time for fetching ths sample sheet', experiment_name)
                continue
        assign_projects_to_run(run_process_obj, l_sample_sheet_path, experiment_name)
        assign_used_library_in_run (run_process_obj,l_sample_sheet_path, experiment_name )
        store_sample_sheet_if_not_defined_in_run (run_process_obj,l_sample_sheet_path, experiment_name)
        import pdb; pdb.set_trace()
        if COPY_SAMPLE_SHEET_TO_REMOTE :
            sample_sheet_file = run_process_obj.get_sample_file()
            sample_sheet_path = os.path.join(wetlab_config.RUN_SAMPLE_SHEET_DIRECTORY, sample_sheet_file)
            try:
                copy_sample_sheet_to_remote_folder(conn, sample_sheet_path, run_folder ,experiment_name)
            except:
                continue

        run_process_obj.set_run_state('Sample Sent')
        logger.info('%s  : is now on Sample Sent state' , experiment_name)
    logger.debug (' End function manage_run_in_recorded_state')
    return
