import os
import logging
from logging.config import fileConfig
from logging.handlers import RotatingFileHandler

from iSkyLIMS_wetlab.models import RunProcess, RunningParameters, Projects, ConfigSetting
from iSkyLIMS_wetlab.wetlab_config import *

from .sample_sheet_utils import get_projects_in_run, get_index_library_name
from .handling_crontab_common_functions import *
from .handling_crontab_run_metric import *
from .handling_crontab_bcl2fastq import *

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
                time_to_check = run_process_obj.get_run_generated_date_no_format().date()
                if not waiting_time_expired(run_process_obj,time_to_check, maximun_time ,experiment_name):
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
            logging_errors(string_message, False, True)
            handling_errors_in_run (experiment_name, run_status['ERROR'])
            logger.debug('%s : End manage_run_in_sample_sent_processing_state function', experiment_name)
        else:
            maximun_time = ConfigSetting.objects.filter(configurationName__exact = 'MAXIMUM_TIME_WAIT_RUN_COMPLETION').last().get_configuration_value()
            time_to_check = run_process_obj.get_run_generated_date_no_format().date()
            if not waiting_time_expired(run_process_obj, time_to_check, maximun_time ,experiment_name):
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


def manage_run_in_processed_run_state(conn, run_process_objs):
    '''
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
        get_run_metric_files          # located in utils.handling_crontab_run_metric.py
        handling_errors_in_run        # located in utils.handling_crontab_common_functions.py
        delete_run_metric_files       # located in utils.handling_crontab_run_metric.py
        parsing_run_metrics_files     # located in utils.handling_crontab_run_metric.py
        create_run_metric_graphics    # located in utils.handling_crontab_run_metric.py
    Return:
        None
    '''
    logger = logging.getLogger(__name__)
    logger.debug (' Starting function manage_run_in_processed_run_state')
    for run_process_obj in run_process_objs:
        experiment_name = run_process_obj.get_run_name()
        logger.info('%s : Start handling in manage_run_in_processed_run_state function', experiment_name)
        run_folder = RunningParameters.objects.filter(runName_id = run_process_obj).last().get_run_folder()
        # delete existing information to avoid having duplicated tables
        delete_existing_run_metrics_table_processed(run_process_obj, experiment_name)
        run_metric_files = get_run_metric_files (conn, run_folder, experiment_name)

        if 'ERROR'in run_metric_files :
            string_message = experiment_name + ' : Unable to collect all files for run metrics'
            logging_errors(string_message, False, True)
            handling_errors_in_run (experiment_name, run_metric_files['ERROR'])
            delete_run_metric_files (experiment_name)
            logger.debug('%s : End manage_run_in_processed_run_state function', experiment_name)
            continue
        parsed_run_stats_summary, parsed_run_stats_read = parsing_run_metrics_files(RUN_TEMP_DIRECTORY_PROCESSING, run_process_obj, experiment_name)

        for run_stat_summary in parsed_run_stats_summary :
            saved_run_stat_summary = StatsRunSummary.objects.create_stats_run_summary(run_stat_summary, run_process_obj)
        logger.info('%s : run metrics summary data saved to database', experiment_name)
        for run_stat_read in parsed_run_stats_read :
            saved_run_stat_read = StatsRunRead.objects.create_stats_run_read(run_stat_read, run_process_obj)
        logger.info('%s :run metrics read data saved to database',experiment_name)

        # create run graphics
        run_graphics = create_run_metric_graphics (RUN_TEMP_DIRECTORY_PROCESSING, run_process_obj, run_folder, experiment_name)
        if 'ERROR' in run_graphics:
            string_message = experiment_name + ' : Unable to save graphics for run metrics'
            logging_errors(string_message, False, True)
            handling_errors_in_run (experiment_name, run_graphics['ERROR'])
            delete_existing_run_metrics_table_processed(run_process_obj, experiment_name)
            delete_run_metric_files (experiment_name)
            logger.debug('%s : End manage_run_in_processed_run_state function', experiment_name)
            continue
        logger.info('%s : run metrics graphics processed and copied to plot image folder',experiment_name)
        # deleting temporary run metrics files
        delete_run_metric_files ( experiment_name)
        # return the state to Processed Run
        run_state = run_process_obj.set_run_state('Processing Bcl2fastq')

    logger.debug (' End function manage_run_in_processed_run_state')
    return


def manage_run_in_processing_bcl2fastq_state(conn, run_process_objs):
    '''
    Description:
        The funtion get the runs in bcl2fastq processing run state. In this state it will check if
        bcl2fastq process is completed by checking if report folder exists
    Input:
        conn                # samba connection instance
        run_process_objs    # list of runProcess objects that are in processing_bcl2fastq
    Constants:
        RUN_TEMP_DIRECTORY_PROCESSING
    Functions:
        check_run_metrics_processed             # located in utils.handling_crontab_run_metric.py
        waiting_time_expired                    # located in utils.handling_crontab_common_functions.py
    Return:
        None
    '''
    logger = logging.getLogger(__name__)
    logger.debug (' Starting function manage_run_in_processing_bcl2fastq_state')
    for run_process_obj in run_process_objs:
        experiment_name = run_process_obj.get_run_name()
        logger.info('%s : Start handling in manage_run_in_processing_bcl2fastq_state function', experiment_name)
        run_folder = RunningParameters.objects.get(runName_id = run_process_obj).get_run_folder()

        bcl2fastq_finish_date = check_demultiplexing_folder_exists(conn, run_folder, experiment_name)
        if 'ERROR' in bcl2fastq_finish_date:
            if bcl2fastq_finish_date['ERROR'] == 29:
                maximun_time = ConfigSetting.objects.filter(configurationName__exact = 'MAXIMUM_TIME_WAIT_TO_RUN_BCL2FASTQ').last().get_configuration_value()
                try:
                    time_to_check = run_process_obj.get_run_completion_date_no_format().date()
                except:
                    string_message = experiment_name + ' :  No Run completion date defined '
                    logging_errors(string_message, False, True)
                    handling_errors_in_run (experiment_name, 30)
                    logger.debug ('%s : Aborting the process. No Run completion date was defined', experiment_name)
                    continue
                if not waiting_time_expired(run_process_obj,time_to_check, maximun_time ,experiment_name):
                    logger.debug ('%s : End the process. Waiting more time to get Sample Sheet file', experiment_name)
                    continue
                else:
                    string_message = experiment_name + ' : Expired time for waiting for demultiplexing folder ' + run_folder
                    logging_errors(string_message, False, True)
                    handling_errors_in_run (experiment_name, 29)
                    logger.debug ('%s : Aborting the process. Exceeded the waiting time for fetching ths Demultiplexing files', experiment_name)
                    continue
            else:
                string_message = experiment_name + ' : Unable to fetch the Conversion Stats file' + run_folder
                logging_errors(string_message, False, True)
                handling_errors_in_run (experiment_name, 31)
                logger.debug ('%s : Aborting the process. Unable to fetch Confersion Stats file', experiment_name)
                continue
        run_process_obj.set_run_bcl2fastq_finished_date(bcl2fastq_finish_date)
        run_process_obj.set_run_state('Processed Bcl2fastq')
        logger.info('%s : Updated to Processed Bcl2Fastq state', experiment_name)
        logger.info('%s : End handling in manage_run_in_processing_bcl2fastq_state function', experiment_name)
    logger.debug (' End function manage_run_in_processing_bcl2fastq_state')
    return

def manage_run_in_processed_bcl2fastq_state(conn, run_process_objs):
    '''
    Description:
        The funtion get the runs in bcl2fastq processed run state. In this state it will get the
        demultiplexing data generated by bcl2fastq.
    Input:
        conn                # samba connection instance
        run_process_objs    # list of runProcess objects that are in processed_bcl2fastq
    Constants:
        RUN_TEMP_DIRECTORY_PROCESSING
    Functions:
        delete_existing_bcl2fastq_table_processed    # located at utils.handling_crontab_bcl2fastq.py
        get_demultiplexing_files                     # located at utils.handling_crontab_bcl2fastq.py
        parsing_demux_and_conversion_files           # located at utils.handling_crontab_bcl2fastq.py
        parsing_demux_sample_project                 # located at utils.handling_crontab_bcl2fastq.py
        process_and_store_raw_demux_project_data     # located at utils.handling_crontab_bcl2fastq.py
        process_and_store_fl_summary_data            # located at utils.handling_crontab_bcl2fastq.py
        process_and_store_lane_summary_data          # located at utils.handling_crontab_bcl2fastq.py
        process_and_store_unknown_barcode_data       # located at utils.handling_crontab_bcl2fastq.py
        process_and_store_samples_projects_data      # located at utils.handling_crontab_bcl2fastq.py
        get_run_disk_utilization                     # located at utils.handling_crontab_bcl2fastq.py

    Return:
        None
    '''
    logger = logging.getLogger(__name__)
    logger.debug (' Starting function manage_run_in_processed_bcl2fastq_state')
    for run_process_obj in run_process_objs:
        experiment_name = run_process_obj.get_run_name()
        logger.info('%s : Start handling in manage_run_in_processed_bcl2fastq_state function', experiment_name)
        run_param_obj = RunningParameters.objects.get(runName_id = run_process_obj)
        run_folder = run_param_obj.get_run_folder()
        # delete existing information to avoid having duplicated tables
        delete_existing_bcl2fastq_table_processed(run_process_obj, experiment_name)
        demux_files = get_demultiplexing_files(conn, run_folder, experiment_name)
        if 'ERROR' in demux_files:
            if demux_files['ERROR'] == 29:
                string_message = experiment_name + ' :  Unable to reach demultiplexing folder '
                logging_errors(string_message, False, True)
                handling_errors_in_run (experiment_name, 30)
                logger.debug ('%s : Aborting the process. No demultiplexing folder exists', experiment_name)
                continue
            else:
                string_message = experiment_name + ' : Unable to fetch the Conversion Stats file' + run_folder
                logging_errors(string_message, False, True)
                handling_errors_in_run (experiment_name, 31)
                logger.debug ('%s : Aborting the process. Unable to fetch Confersion Stats files', experiment_name)
                continue
        number_of_lanes = run_param_obj.get_number_of_lanes()
        try:
            number_of_lanes = int(number_of_lanes)
        except:
            string_message = experiment_name + ' : Sequencer used in the run has not defined the number of lanes'
            logging_errors(string_message, False, True)
            handling_errors_in_run (experiment_name, 32)
            logger.debug ('%s : Aborting the process. Number of lanes not defined on the sequencer', experiment_name)
            continue
        # parsing the files to get the xml Stats
        logger.info('%s : Start parsing  demultiplexing files',experiment_name)
        project_parsed_data = parsing_demux_and_conversion_files(demux_files, number_of_lanes, experiment_name)

        logger.info('%s : Start parsing  samples demultiplexing',experiment_name)
        sample_project_parsed_data = parsing_demux_sample_project (demux_files, number_of_lanes, experiment_name)

        try:
            # clean up the fetched files in the local temporary folder
            os.remove(demux_files['demux_stats'])
            os.remove(demux_files['conversion_stats'])
            logger.info ('%s : Deleted temporary demultiplexing and conversion files', experiment_name)
        except:
            string_message = experiment_name + ' : Unable to delete the Conversion Stats file' + run_folder
            logging_errors(string_message, False, True)
            logger.debug ('%s : Allowing to contine the parssing process', experiment_name)

        process_and_store_raw_demux_project_data(project_parsed_data, run_process_obj, experiment_name)

        process_and_store_fl_summary_data(project_parsed_data, run_process_obj,number_of_lanes, experiment_name )

        process_and_store_lane_summary_data(project_parsed_data, run_process_obj,number_of_lanes, experiment_name )

        process_and_store_unknown_barcode_data(project_parsed_data, run_process_obj,number_of_lanes, experiment_name )

        process_and_store_samples_projects_data(sample_project_parsed_data, run_process_obj, experiment_name)

        ## Get the disk space utilization for this run
        try:
            disk_utilization = get_run_disk_utilization (conn, run_folder, experiment_name)
        except:
            string_message = experiment_name + ' : Error when fetching the disk utilization'
            logging_errors (string_message, True, True)
            handling_errors_in_run (experiment_name, '17' )
            logger.debug('%s : End function manage_run_in_processed_bcl2fast2_run with error', experiment_name)
            continue

        result_store_usage = run_process_obj.set_used_space (disk_utilization)
        finish_date = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        result_set_finish_date = run_process_obj.set_run_finish_date(finish_date)
        # Update the run state to completed
        run_state = run_process_obj.set_run_state('Completed')
        logger.info('%s : is Completed',experiment_name)
    logger.debug (' End function manage_run_in_processed_bcl2fastq_state')
