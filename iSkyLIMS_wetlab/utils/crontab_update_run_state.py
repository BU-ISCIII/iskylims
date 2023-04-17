import os
import logging
from logging.config import fileConfig
from logging.handlers import RotatingFileHandler

from iSkyLIMS_wetlab.models import *
from iSkyLIMS_wetlab.wetlab_config import *

from .sample_sheet_utils import *
from .common import *
from .crontab_process import *

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
        open_samba_connection                # located in utils.common.py
        get_list_processed_runs              # located at this file
        get_new_runs_on_remote_server        # located at utils.common.py
        get_new_runs_from_remote_server      # located at utils.handling_crontab_common_functions
        get_experiment_name_from_file        # located at utils.common.py
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
            except Exception as e:
                error_message = 'Unable to fetch RunParameter file for folder :' + new_run
                logging_errors(error_message, True, False)
                # we don't have experiment name when there is no run_parameters file. We used the run folder name instead.
                run_process_obj = get_run_process_obj_or_create_if_not_exists(new_run)
                handling_errors_in_run(new_run, '21')
                logger.debug ('%s : Deleting RunParameter file', experiment_name)
                logger.debug ('%s : Aborting the process for this run. Continue with the next.', experiment_name)
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

            # Fetch run info
            l_run_info_path = os.path.join(wetlab_config.RUN_TEMP_DIRECTORY, wetlab_config.RUN_INFO)
            s_run_info_path = os.path.join(get_samba_application_shared_folder(), new_run, wetlab_config.RUN_INFO)
            try:
                l_run_info = fetch_remote_file (conn, new_run, s_run_info_path, l_run_info_path)
                logger.info('%s : Sucessfully fetch of RunInfo file', experiment_name)
            except Exception:
                string_message = experiment_name + ' : Unable to fetch the RunInfo file on folder ' + new_run
                logging_errors(string_message, True, False)
                run_process_obj = get_run_process_obj_or_create_if_not_exists(experiment_name)
                handling_errors_in_run(experiment_name, '20')
                # cleaning up the RunParameter in local temporaty file
                logger.debug ('%s : Deleting RunParameter file', experiment_name)
                os.remove(l_run_parameter)
                logger.debug ('%s : Aborting the process for this run. Continue with the next.', experiment_name)
                continue
            
            logger.debug('%s : Found the experiment name called : %s', new_run, experiment_name)
            exclude_states = ["Error", "Recorded"]
            if RunProcess.objects.filter(runName__exact=experiment_name).exclude(state__runStateName__in=exclude_states).exists():
                # This situation should not occurr. The run_processed file should  have this value. To avoid new iterations with this run
                # we update the run process file with this run and continue  with the next item
                run_state = RunProcess.objects.get(run_name__exact = experiment_name).get_state()
                string_message = new_run  + ' :  experiment name  state ' + experiment_name + ' in incorrect state. Run state is ' + run_state
                logging_errors( string_message, False, False)
                logger.info('%s : Deleting temporary runParameter file' , experiment_name)
                os.remove(l_run_parameter)
                logger.info('%s : RunParameter file. Local copy deleted',experiment_name)
                continue

            running_parameters = parsing_run_info_and_parameter_information(l_run_info, l_run_parameter, experiment_name)
            logger.info('%s  : Deleting runParameter file', experiment_name)
            os.remove(l_run_parameter)
            logger.info('%s  : Deleting runInfo file', experiment_name)
            os.remove(l_run_info)

            run_process_obj = get_run_process_obj_or_create_if_not_exists(experiment_name)
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
                    file_read = read_user_iem_file(l_sample_sheet_path)
                    users = validate_userid_in_user_iem_file(file_read, user_id_list)
                
                    if 'ERROR' in users:
                        string_message = experiment_name + ' : Description field does not contains userid.'
                        logging_errors(string_message, True, False)
                        handling_errors_in_run(experiment_name, '1')
                        continue
                
                assign_projects_to_run(run_process_obj, l_sample_sheet_path, experiment_name)
                assign_used_library_in_run (run_process_obj,l_sample_sheet_path, experiment_name)
                store_sample_sheet_if_not_defined_in_run (run_process_obj,l_sample_sheet_path, experiment_name)
            else :
                if wetlab_config.COPY_SAMPLE_SHEET_TO_REMOTE and  'NextSeq' in running_parameters['running_data'][wetlab_config.APPLICATION_NAME_TAG]:
                    sample_sheet = run_process_obj.get_sample_file()
                    sample_sheet_path = os.path.join(settings.MEDIA_ROOT, sample_sheet)
                    run_folder = RunningParameters.objects.get( run_name_id__exact = run_process_obj).get_run_folder()
                    try:
                        copy_sample_sheet_to_remote_folder(conn, sample_sheet_path, run_folder,experiment_name)
                    except Exception as e:
                        string_message = experiment_name + ' : Unable to copy Sample Sheet to Remote folder' + new_run
                        logging_errors(string_message, True, False)
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
        open_samba_connection   # located in utils.common.py

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
                logging_errors(string_message, False, False)

    logger.debug ('End function for search_not_completed_run')
    return


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
                    logging_errors(string_message, False, False)
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
        check_log_for_run_completions           
        waiting_time_expired                    
        handling_errors_in_run                 
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
            logging_errors(string_message, False, False)
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
            logging_errors(string_message, False, False)
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
                logging_errors(string_message, False, False)
                handling_errors_in_run (experiment_name, 9)
                logger.debug('%s  : End manage_run_in_sample_sent_processing_state function', experiment_name)

    logger.debug ('End function manage_run_in_sample_sent_processing_state')
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
            logging_errors(string_message, True, False)
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
            logging_errors(string_message, True, False)
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
        check_demultiplexing_folder_exists      # located in utils.handling_crontab_bcl2fastq.py
    Return:
        None
    '''
    logger = logging.getLogger(__name__)
    logger.debug('Starting function manage_run_in_processing_bcl2fastq_state')
    for run_process_obj in run_process_objs:
        experiment_name = run_process_obj.get_run_name()
        logger.info('%s : Start handling in manage_run_in_processing_bcl2fastq_state function', experiment_name)
        run_folder = RunningParameters.objects.get(run_name_id = run_process_obj).get_run_folder()

        bcl2fastq_finish_date = check_demultiplexing_folder_exists(conn, run_folder, experiment_name)
        if 'ERROR' in bcl2fastq_finish_date:
            if bcl2fastq_finish_date['ERROR'] == 29 or bcl2fastq_finish_date['ERROR'] == 31:
                maximun_time = ConfigSetting.objects.filter(configurationName__exact = 'MAXIMUM_TIME_WAIT_TO_RUN_BCL2FASTQ').last().get_configuration_value()
                try:
                    time_to_check = run_process_obj.get_run_completion_date_no_format().date()
                except:
                    string_message = experiment_name + ' :  Aborting the process. No Run completion date was defined.'
                    logging_errors(string_message, True, False)
                    handling_errors_in_run (experiment_name, 30)
                    logger.debug ('%s : Aborting the process. No Run completion date was defined', experiment_name)
                    continue
                if not waiting_time_expired(run_process_obj,time_to_check, maximun_time ,experiment_name):
                    logger.debug ('%s : End the process. Waiting more time to get bcf2fastq stats.', experiment_name)
                    continue
                else:
                    string_message = experiment_name + ' :  Aborting the process. Exceeded the waiting time for fetching demultiplexing files ' + run_folder
                    logging_errors(string_message, True, False)
                    handling_errors_in_run (experiment_name, bcl2fastq_finish_date['ERROR'])
                    continue
            else:
                string_message = experiment_name + ' : Aborting the process. Error in bcf2fastq process.' + run_folder
                logging_errors(string_message, True, False)
                handling_errors_in_run (experiment_name, bcl2fastq_finish_date['ERROR'])
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
    '''
    logger = logging.getLogger(__name__)
    logger.debug (' Starting function manage_run_in_processed_bcl2fastq_state')
    for run_process_obj in run_process_objs:
        experiment_name = run_process_obj.get_run_name()
        logger.info('%s : Start handling in manage_run_in_processed_bcl2fastq_state function', experiment_name)
        run_param_obj = RunningParameters.objects.get(run_name_id = run_process_obj)
        run_folder = run_param_obj.get_run_folder()
        # delete existing information to avoid having duplicated tables
        delete_existing_bcl2fastq_table_processed(run_process_obj, experiment_name)
        demux_files = get_demultiplexing_files(conn, run_folder, experiment_name)
        if 'ERROR' in demux_files:
            if demux_files['ERROR'] == 29 or demux_files['ERROR'] == 31:
                string_message = experiment_name + ' : Aborting the process. Unable to reach demultiplexing files '
                logging_errors(string_message, True, False)
                handling_errors_in_run (experiment_name, demux_files['ERROR'])
                continue
           
        number_of_lanes = run_param_obj.get_number_of_lanes()
        try:
            number_of_lanes = int(number_of_lanes)
        except:
            string_message = experiment_name + ' : Sequencer used in the run has not defined the number of lanes'
            logging_errors(string_message, True, False)
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
            logging_errors(string_message, True, False)
            logger.debug ('%s : Allowing to contine the parssing process', experiment_name)

        process_and_store_raw_demux_project_data(project_parsed_data, run_process_obj, experiment_name)

        process_and_store_fl_summary_data(project_parsed_data, run_process_obj,number_of_lanes, experiment_name )

        process_and_store_lane_summary_data(project_parsed_data, run_process_obj,number_of_lanes, experiment_name )

        process_and_store_unknown_barcode_data(project_parsed_data, run_process_obj,number_of_lanes, experiment_name )

        try:
            process_and_store_samples_projects_data(sample_project_parsed_data, run_process_obj, experiment_name)
        except KeyError as key_error:
            string_message = experiment_name + ' : Error when processing and storing samples in projects.'
            logging_errors (string_message, True, False)
            if key_error.args[0] == 33:
                handling_errors_in_run (experiment_name, '33' )
            else:
                string_message = experiment_name + ' : Unknown error when processing and storing samples in projects.'
                logging_errors (string_message, True, False) 
            logger.debug('%s : End function manage_run_in_processed_bcl2fast2_run with error', experiment_name)
            continue
        except Exception as e:
            string_message = experiment_name + ' : Unknown error when processing and storing samples in projects.'
            logging_errors (string_message, True, False)
            continue

        ## Get the disk space utilization for this run
        try:
            disk_utilization = get_run_disk_utilization (conn, run_folder, experiment_name)
        except:
            string_message = experiment_name + ' : Error when fetching the disk utilization'
            logging_errors (string_message, True, False)
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
