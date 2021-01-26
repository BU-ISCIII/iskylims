#!/usr/bin/env python3

import sys, os, re
#import xml.etree.ElementTree as ET
#import shutil
#import locale
#import datetime, time
from iSkyLIMS_wetlab.models import RunProcess, RunStates, RunningParameters

import logging

from iSkyLIMS_wetlab import wetlab_config


from .handling_crontab_common_functions import *
from .generic_functions import *

from .miseq_run_functions import  handle_miseq_run , manage_miseq_in_samplesent,  manage_miseq_in_processing_run
from .nextseq_run_functions import handle_nextseq_recorded_run, manage_nextseq_in_samplesent, manage_nextseq_in_processing_run
from .common_run_functions import manage_run_in_processed_run, manage_run_in_processing_bcl2fastq, manage_run_in_processed_bcl2fastq


from django.conf import settings
from django_utils.models import Center



def read_processed_runs_file (processed_run_file) :
    '''
    Description:
        The function reads the file that contains all the processed runs
        and return a list with the run folder names
    Input:
        processed_run_file # full path and name of the file
    Variable:
        processed_runs  # list of the folder names read from file
    Return:
        Error when file can not be readed
        processed_runs variable for successful file read
    '''
    logger = logging.getLogger(__name__)
    logger.debug('Starting for reading processed run file' )
    w_dir = os.getcwd()
    logger.debug('Folder for fetching the proccess run file is %s', w_dir)
    logger.debug('processed_run_file is %s', processed_run_file)
    processed_runs = []
    if os.path.exists(processed_run_file):
        try:
            fh = open (processed_run_file,'r')
            for line in fh:
                line=line.rstrip()
                processed_runs.append(line)
        except Exception as e:
            string_message = 'Unable to open the processed run file. '
            logging_errors(string_message, True, True)
            raise
        fh.close()
        logger.info('run processed file have been read')
        logger.debug('Exiting sucessfully the function read_processed_runs_file' )
        return processed_runs
    else:
        logger.debug('No found processed run file. Exiting the function' )
        return 'Error'


def get_list_processed_runs () :
    '''
    Description:
        The function get the run folder id from the Running Parameter
        table. This list will be used to compare agains the folder on
        remote server.
    Return:
        raise exception when not access to database
        processed_runs
    '''
    logger = logging.getLogger(__name__)
    logger.debug('Starting function get_list_processed_runs' )
    processed_runs = []
    r_parameters_objects = RunningParameters.objects.all().exclude(runName_id__state__runStateName = "Recorded")

    for r_parameter in r_parameters_objects :
        processed_runs.append(r_parameter.get_run_folder())

    logger.info('run processed list is filled')
    logger.debug('End function get_list_processed_runs' )
    return processed_runs


def update_processed_run_file (processed_run_file, processed_runs) :
    '''
    Description:
        The function write the file that contains all the processed runs
        return True if the file was sucessfully write
    Input:
        processed_run_file # full path and name of the file
    Variable:
        processed_runs  # list of the folder names to write in the file
    Return:
        Error when file can not be create
        True if sucessfully
    '''
    logger = logging.getLogger(__name__)
    logger.debug('Starting update_processed_run_file' )
    #import pdb; pdb.set_trace()
    fh =open (processed_run_file,'w')
    # update the process_run_file with new runs
    for processed in processed_runs:
        fh.write(processed)
        fh.write('\n')
    fh.close()


    logger.debug('Exit update_processed_run_file' )
    return True

def search_update_new_runs ():
    '''
    Description:
        The function will check if there are new run folder in the remote
        server.
        Get the runParameter file to identify if run is NextSeq or miSeq
        to execute its dedicate handler process.
    Functions:
        fetch_remote_file                    # located at utils.handling_crontab_common_functions
        open_samba_connection                # located in utils.generic_functions.py
        get_list_processed_runs              # located at this file
        get_new_runs_on_remote_server        # located at utils.generic_functions.py
        get_experiment_name_from_file        # located at utils.generic_functions.py
        get_samba_application_shared_folder  # located at utils.handling_crontab_common_functions.py
        handle_nextseq_recorded_run          # located at utils.nextseq_run_functions.py
        validate_sample_sheet                # located at this file
        save_new_miseq_run                   # located at this file
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

        conn=open_samba_connection()
        logger.info('Sucessfully  SAMBA connection for search_update_new_runs')
    except Exception as e:
        string_message = 'Unable to open SAMBA connection for the process search update runs'
        # raising the exception to stop crontab
        logging_errors(string_message, True, True)
        raise

    new_runs = get_new_runs_from_remote_server (processed_runs, conn, get_samba_shared_folder())

    if len (new_runs) > 0 :
        for new_run in new_runs :
            l_run_parameter_path = os.path.join(wetlab_config.RUN_TEMP_DIRECTORY, wetlab_config.RUN_PARAMETER_NEXTSEQ)
            s_run_parameter_path = os.path.join(get_samba_application_shared_folder(), new_run, wetlab_config.RUN_PARAMETER_NEXTSEQ)
            try:
               l_run_parameter = fetch_remote_file (conn, new_run, s_run_parameter_path, l_run_parameter_path)
               logger.info('%s : Sucessfully fetch of RunParameter file', new_run)
            except Exception as e:
                error_message = 'Unable to fetch RunParameter file for folder :' +  new_run
                logging_errors(error_message, True, True)
                continue

            experiment_name = get_experiment_name_from_file (l_run_parameter)

            if experiment_name == ''  or experiment_name == 'NOT FOUND':
                if experiment_name == '':
                    string_message = new_run + ' : Experiment name is empty'
                else:
                    string_message = new_run + ' : Experiment name field was not found in file'
                logging_errors(string_message, False, True)
                os.remove(l_run_parameter)
                logger.info(' %s  : Deleted temporary run parameter file', new_run)
                continue
            logger.debug('%s : Found the experiment name called : %s', new_run , experiment_name)
            if RunProcess.objects.filter(runName__exact = experiment_name).exclude(state__runStateName ='Recorded').exists():
                # This situation should not occurr. The run_processed file should
                # have this value. To avoid new iterations with this run
                # we update the run process file with this run and continue
                # with the next item
                run_state = RunProcess.objects.get(runName__exact = experiment_name).get_state()
                string_message = new_run  + ' :  experiment name  state ' + experiment_name + 'in incorrect state. Run state is ' + run_state
                logging_errors( string_message, False, True)
                logger.info('%s : Deleting temporary runParameter file' , experiment_name)
                os.remove(l_run_parameter)
                logger.info('%s : RunParameter file. Local copy deleted',experiment_name)
                continue


            # Fetch run info
            l_run_info_path = os.path.join(wetlab_config.RUN_TEMP_DIRECTORY, wetlab_config.RUN_INFO)
            s_run_info_path = os.path.join(get_samba_application_shared_folder() , new_run, wetlab_config.RUN_INFO)
            try:
                l_run_info = fetch_remote_file (conn, new_run, s_run_info_path, l_run_info_path)
                logger.info('%s : Sucessfully fetch of RunInfo file', experiment_name)
            except Exception as e:
                string_message = experiment_name + ' : Unable to fetch the RunInfo file on folder ' + new_run
                logging_errors(string_message, True, True)
                handling_errors_in_run(experiment_name, '20')
                # cleaning up the RunParameter in local temporaty file
                logger.debug ('%s : Deleting RunParameter file', experiment_name)
                os.remove(l_run_parameter)
                logger.debug ('%s : End function for handling NextSeq run with exception', experiment_name)
                raise Exception   # returning to handle next run folder

            # Finding out the platform to continue the run processing
            run_platform =  get_run_platform_from_file (l_run_parameter)
            # branch according platform to continue the run processing

            logger.debug('%s : Found platform name  , %s', experiment_name, run_platform)
            if run_platform == 'NOT FOUND':
                string_message = new_run + ': Exting this run becuase Platform tag  was not found RunParameter file'
                logging_errors (string_message, False, True)
                continue
            if 'MiSeq' in run_platform :
                logger.debug('%s  : MiSeq run found. Executing miseq handler ', experiment_name)
                try:
                    update_miseq_process_run =  handle_miseq_run (conn, new_run, l_run_parameter, experiment_name)
                    if update_miseq_process_run != '' :
                        new_processed_runs.append(experiment_name)
                        logger.info('%s : was successfully processed ', experiment_name)
                        logger.debug('%s : Finished miSeq handling process', experiment_name)
                    continue
                except ValueError as e :
                    logger.warning(' %s : Error found when processing miSeq run %s ', experiment_name, e)
                    run_with_error.append(experiment_name)
                    logger.debug('%s : Finished miSeq handling process with error', experiment_name)
                    continue
                except Exception as e:
                    string_message = experiment_name + " : Unexpected Error " + str (e)
                    logging_errors(string_message, True, True)
                    continue
                '''
                except :
                    logger.warning('miSeq run  %s does not have all required files. Giving more time for the sequencer to write them.', experiment_name )
                    logger.info('Continue processing next item ')
                    logger.debug('Finished miSeq handling process with error')
                    continue
                '''
            elif 'NextSeq' in run_platform :
                logger.debug('%s  : Executing NextSeq handler ', experiment_name)
                try:
                    update_nextseq_process_run =  handle_nextseq_recorded_run (conn, new_run, l_run_parameter, l_run_info, experiment_name)
                except Exception as e:
                    string_message = experiment_name +  ' : Error when processing the handle_nextseq_recorded_run function'
                    logging_errors(string_message, True, True)
                    logger.debug('%s :  Finished NextSeq handling process', experiment_name)
                    continue
                else:
                    import pdb; pdb.set_trace()
                    if update_nextseq_process_run != '' :
                        new_processed_runs.append(experiment_name)
                        logger.info('%s : was successfully processed ', experiment_name)
                        logger.debug('%s : Finished miSeq handling process', experiment_name)
                    continue
                logger.debug('%s :  Finished NextSeq handling process', experiment_name)
            else:
                string_message = 'Platform for this run is not supported'
                logging_errors(string_message, False , True)
                # Set run to error state
                os.remove(l_run_parameter)

    else:
        logger.info('No found new run folders on the remote server')

    '''
    if process_run_file_update :
        logger.info('Updating the process_run_file with the new runs')
        try:
            update_processed_run_file (processed_run_file, processed_runs)
            logger.info('File for processed runs was updted ')
        except:
            string_message = 'Unable to write the processed runs file'
            logging_errors(string_message, True, True)
    '''
    logger.info ('Clossing SAMBA connection')
    conn.close()
    logger.debug ('End function searching new runs. Returning handle runs ' )
    return new_processed_runs , run_with_error



def search_not_completed_run ():
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
        open_samba_connection # located in utils.generic_functions.py

        save_new_miseq_run # located at this file
    Constants:
        PROCESSED_RUN_FILE
        RUN_TEMP_DIRECTORY
        SAMBA_SHARED_FOLDER_NAME
        SAMPLE_SHEET
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function for search_not_completed_run')
    try:
        conn=open_samba_connection()
        logger.info('Sucessfully  SAMBA connection for the process_run_in_recorded_state')
    except Exception as e:
        string_message = 'Unable to open SAMBA connection for the process search update runs'
        # raising the exception to stop crontab
        raise logging_errors(string_message, True, False)

    runs_to_handle = {}
    updated_run={}
    runs_with_error = {}
    state_list_be_processed = ['Sample Sent','Processing Run','Processed Run', 'Processing Bcl2fastq',
                                'Processed Bcl2fastq']
    # get the list for all runs that are not completed
    for state in state_list_be_processed:
        run_state = RunStates.objects.get(runStateName__exact = state)
        #import pdb; pdb.set_trace()
        if RunProcess.objects.filter(state__exact = run_state).exists():
            runs_found_in_state = RunProcess.objects.filter(state__exact = run_state)
            runs_to_handle[state] = []
            for run_found_in_state in runs_found_in_state :
                runs_to_handle[state].append(run_found_in_state)

    for state in runs_to_handle:
        logger.info ('Start processing the run found for state %s', state)

        if state == 'Sample Sent':
            updated_run['Sample Sent'] = []
            runs_with_error['Sample Sent'] = []
            logger.debug('--------Start handling the runs in Sample Sent state--------')
            for run_in_sample_sent in runs_to_handle[state] :
                experiment_name = run_in_sample_sent.get_run_name()
                # get platform
                try:
                    run_platform = run_in_sample_sent.get_run_platform()
                except:
                    string_message = experiment_name + ' : Platform not defined'
                    logging_errors (string_message , False, True)
                    continue
                experiment_name = run_in_sample_sent.runName
                logger.info(' %s : Start handling in Sample Sent state', experiment_name)
                try:
                    if 'NextSeq' in run_platform :
                        logger.info(' %s : Handle on NextSeq branch', experiment_name)
                        updated_run['Sample Sent'].append(manage_nextseq_in_samplesent(conn, run_in_sample_sent))
                    elif 'MiSeq' in run_platform :
                        logger.info(' %s : Handle on MiSeq branch', experiment_name)
                        updated_run['Sample Sent'].append(manage_miseq_in_samplesent(conn, run_in_sample_sent))
                    else:
                        string_message = experiment_name + ' : Platform ' + run_platform +' is not supported '
                        logging_errors (string_message , False, True)
                        continue
                except :
                    runs_with_error[state].append(run_in_sample_sent.get_run_name())
                    logger.info('%s : Handling the exception to continue with the next item', experiment_name)
                    continue
            logger.debug('--------End runs in Sample Sent state--------')

        elif state == 'Processing Run':
            updated_run[state] = []
            runs_with_error[state] = []
            logger.debug('--------Start handling the runs in  Processing Run state------')
            for run_in_processing_run in runs_to_handle[state] :
                experiment_name = run_in_processing_run.get_run_name()
                try:
                    run_platform = run_in_processing_run.get_run_platform()
                except:
                    string_message = experiment_name + ' : Platform not defined'
                    logging_errors (string_message , False, True)
                    continue
                experiment_name = run_in_processing_run.get_run_name()
                logger.info(' %s : Start handling in Processing run state', experiment_name)

                try:
                    if 'NextSeq' in run_platform :
                        updated_run[state].append(manage_nextseq_in_processing_run(conn, run_in_processing_run))
                    elif 'MiSeq' in run_platform :
                        updated_run[state].append(manage_miseq_in_processing_run(conn, run_in_processing_run))
                    else:
                        string_message = 'Platform ' + run_platform +' is not supported '
                        logging_errors (string_message , False, False)
                        continue
                except :
                    runs_with_error[state].append(run_in_processing_run.get_run_name())
                    logger.info('%s : Handling the exception to continue with the next item', experiment_name)
                    continue
            logger.debug('--------End runs in Processing Run state--------')

        elif state == 'Processed Run':
            updated_run[state] = []
            runs_with_error[state] = []
            logger.debug('--------Start handling the runs in  Processed Run state--------')
            for run_in_processed_run in runs_to_handle[state]:
                try:
                    updated_run[state].append(manage_run_in_processed_run(conn, run_in_processed_run))
                except :
                    runs_with_error[state].append(run_in_processed_run.get_run_name())
                    logger.info('Handling the exception to continue with the next item')
                    continue
            logger.debug('--------End runs in Processed Run state--------')

        elif state == 'Processing Bcl2fastq':
            updated_run[state] = []
            runs_with_error[state] = []
            logger.debug('--------Start handling the runs in  Processing Bcl2fastq--------')
            for run_in_processing_bcl2fastq_run in runs_to_handle[state] :
                try:
                    updated_run[state].append( manage_run_in_processing_bcl2fastq (conn, run_in_processing_bcl2fastq_run))
                except :
                    runs_with_error[state].append(run_in_processing_bcl2fastq_run.get_run_name())
                    logger.info('Handling the exception on Processing Bcl2fastq.  Continue with the next item')
                    continue
            logger.debug('--------End handling the runs in  Processing Bcl2fastq--------')
        elif state == 'Processed Bcl2fastq':
            updated_run[state] = []
            runs_with_error[state] = []
            logger.debug('--------Start handling the runs in  Processed Bcl2fastq--------')
            for run_in_processed_bcl2fastq_run in runs_to_handle[state] :
                try:
                    updated_run[state].append( manage_run_in_processed_bcl2fastq (conn, run_in_processed_bcl2fastq_run))
                except :
                    runs_with_error[state].append(run_in_processed_bcl2fastq_run.get_run_name())
                    logger.info('Handling the exception on Processed Bcl2fastq.  Continue with the next item')
                    continue
            logger.debug('--------End handling the runs in  Processed Bcl2fastq--------')
        else:
            string_message = 'Run in unexpected state. ' + state
            logging_errors (string_message , False, False)
            continue
    logger.debug ('End function for search_not_completed_run')
    return updated_run, runs_with_error
