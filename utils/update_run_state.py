#!/usr/bin/env python3

import sys, os, re
#import xml.etree.ElementTree as ET
#import shutil
#import locale
#import datetime, time
from iSkyLIMS_wetlab.models import RunProcess
#from .interop_statistics import *
import logging


from iSkyLIMS_wetlab import wetlab_config
from iSkyLIMS_drylab.models import Machines, Platform

from .run_common_functions import *
from .miseq_run_functions import handle_miseq_run , manage_miseq_in_samplesent,  manage_miseq_in_processing_run
from .nextseq_run_functions import handle_nextseq_recorded_run, manage_nextseq_in_samplesent, manage_nextseq_in_processing_run

#from .sample_sheet_utils import get_experiment_library_name, get_projects_in_run

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
            logging_errors(logger, string_message, True, True)
            
            raise
        fh.close()
        logger.info('run processed file have been read')
        logger.debug('Exiting sucessfully the function read_processed_runs_file' )
        return processed_runs
    else:
        logger.debug('No found processed run file. Exiting the function' )
        return 'Error'
    

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
    import pdb; pdb.set_trace()
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
        open_samba_connection # located in utils.wetlab_misc_utilities.py 
        get_new_runs_on_remote_server # located at this file
        validate_sample_sheet   # located at this file
        save_new_miseq_run # located at this file
    Constants:
        PROCESSED_RUN_FILE
        RUN_TEMP_DIRECTORY
        SAMBA_SHARED_FOLDER_NAME
        SAMPLE_SHEET
    Variables:
        handle_new_miseq_run # list with all new miseq runs
        l_sample_sheet  # full path for storing sample sheet file on 
                        tempary local folder
        new_processed_runs # list with the new folder run processed
        s_sample_sheet  # full path for remote sample sheet file
        processed_run_file # path
        processed_runs  # list of run that are moved to Sample Sent state
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function for searching new runs')
    processed_run_file = os.path.join( wetlab_config.RUN_TEMP_DIRECTORY, wetlab_config.PROCESSED_RUN_FILE)
    processed_runs = read_processed_runs_file (processed_run_file)
    process_run_file_update = False
    new_processed_runs = []
    try:
        conn=open_samba_connection()
        logger.info('Sucessfully  SAMBA connection for the process_run_in_recorded_state')
    except Exception as e:
        string_message = 'Unable to open SAMBA connection for the process search update runs'
        # raising the exception to stop crontab
        raise logging_errors(logger, string_message, True, False)
    new_runs = get_new_runs_from_remote_server (processed_runs, conn,
                            wetlab_config.SAMBA_SHARED_FOLDER_NAME)
    
    if len (new_runs) > 0 :
        for new_run in new_runs :
            l_run_parameter = os.path.join(wetlab_config.RUN_TEMP_DIRECTORY, wetlab_config.RUN_PARAMETER_NEXTSEQ)
            s_run_parameter = os.path.join(new_run,wetlab_config.RUN_PARAMETER_NEXTSEQ)
            try:
               l_run_parameter = fetch_remote_file (conn, new_run, s_run_parameter, l_run_parameter)
               logger.info('Sucessfully fetch of RunParameter file')
            except:
                logger.info('Exception fetched to continue with the next run')
                continue

            experiment_name = get_experiment_name_from_file (l_run_parameter)
            logger.debug('found the experiment name  , %s', experiment_name)
            if experiment_name == '' :
                string_message = 'Experiment name for ' + new_run + ' is empty'
                logging_errors(logger, string_message, False, True)
                logger.info('Exceptional condition reported on log. Continue with the next run')
                continue

            if RunProcess.objects.filter(runName__exact = experiment_name).exclude(runState__exact ='Recorded').exists():
                # This situation should not occurr. The run_processed file should
                # have this value. To avoid new iterations with this run 
                # we update the run process file with this run and continue
                # with the next item
                run_state = RunProcess.objects.get(runName__exact = experiment_name).get_state()
                string_message = 'The run ' + experiment_name + 'is in state ' + run_state + '. Should be in Recorded'
                logging_errors(logger, string_message, False, False)
                process_run_file_update = True
                processed_runs.append(new_run)
                new_processed_runs.append(new_run)
                continue
            # Run is new or it is in Recorded state. 
            # Finding out the platform to continue the run processing
            run_platform =  get_run_platform_from_file (l_run_parameter)
            logger.debug('found the platform name  , %s', run_platform)
            if 'MiSeq' in run_platform :
                logger.debug('Executing miseq handler ')
                try:
                    update_miseq_process_run =  handle_miseq_run (conn, new_run, l_run_parameter, experiment_name)
                    if update_miseq_process_run != '' :
                        process_run_file_update = True
                        processed_runs.append(new_run)
                        new_processed_runs.append(new_run)
                    continue
                except ValueError as e :
                    # Include the run in the run processed file 
                    logger.warning('Error found when processing miSeq run %s ', e)
                    logger.info('Including this run in the run processed file ')
                    process_run_file_update = True
                    processed_runs.append(new_run)
                    new_processed_runs.append(new_run)
                    continue
                except :
                    logger.warning('miSeq run is waiting for sequencer run to have all files')
                    logger.info('Continue processing next item ')
                    continue
                logger.debug('Finished miSeq handling process')

            elif 'NextSeq' in run_platform :

                logger.debug('Executing NextSeq handler ')

                try:
                    update_nextseq_process_run =  handle_nextseq_recorded_run (conn, new_run, l_run_parameter, experiment_name)
                    if update_nextseq_process_run != '' :
                        process_run_file_update = True
                        processed_runs.append(new_run)
                        new_processed_runs.append(new_run)
                    continue
                except ValueError as e :
                    # Include the run in the run processed file 
                    logger.warning('Error found when processing miSeq run %s ', e)
                    logger.info('Including this run in the run processed file ')
                    process_run_file_update = True
                    processed_runs.append(new_run)
                    new_processed_runs.append(new_run)
                    continue
                except :
                    logger.warning('NextSeq run is waiting for sequencer to have all files')
                    logger.info('Continue processing next item ')
                    continue
                logger.debug('Finished miSeq handling process')
            else:
                string_message = 'Platform for this run is not supported'
                logging_errors(logger, string_message)
                # Set run to error state
                os.remove(l_run_parameter)

    else:
        logger.info('No found new run folders on the remote server')
        return ''

    if process_run_file_update :
        logger.info('Updating the process_run_file wit the new runs')
        try:
            update_processed_run_file (processed_run_file, processed_runs)
            logger.info('File for processed runs was updted ')
        except:
            string_message = 'Unable to write the processed runs file'
            logging_errors(logger, string_message, True, True)
    logger.info ('Clossing SAMBA connection')
    conn.close()
    logger.debug ('End function searching new runs. Returning handle runs ' )
    return new_processed_runs
    


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
        open_samba_connection # located in utils.wetlab_misc_utilities.py 
        get_new_runs_on_remote_server # located at this file
        validate_sample_sheet   # located at this file
        save_new_miseq_run # located at this file
    Constants:
        PROCESSED_RUN_FILE
        RUN_TEMP_DIRECTORY
        SAMBA_SHARED_FOLDER_NAME
        SAMPLE_SHEET
        
    Variables:
    
    
    
        runs_to_handle # dictionary contains the run state as key and
                        run objects in a value list
        run_platform  # platform used on the run in sample sent state
                        
        state_list_be_processed # list contains the string state that have
                                to be handle to move to complete state
        
        handle_new_miseq_run # list with all new miseq runs
        l_sample_sheet  # full path for storing sample sheet file on 
                        tempary local folder
        s_sample_sheet  # full path for remote sample sheet file
        processed_run # dictionary having state as key and the runs 
                        processed as value list
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function for search_not_completed_run')
    try:
        conn=open_samba_connection()
        logger.info('Sucessfully  SAMBA connection for the process_run_in_recorded_state')
    except Exception as e:
        string_message = 'Unable to open SAMBA connection for the process search update runs'
        # raising the exception to stop crontab
        raise logging_errors(logger, string_message, True, False)
        
    runs_to_handle = {}
    state_list_be_processed = ['Sample Sent','Processing Run','Processed Run', 'Processing Bcl2Fastq',
                                'Processed Bcl2Fastq']
    # get the list for all runs that are not completed
    for state in state_list_be_processed:
        if RunProcess.objects.filter(runState__exact = state).exists():
            runs_to_handle[state]=RunProcess.objects.filter(runState__exact = state)
    
    updated_run={}
    for state in runs_to_handle:
        logger.info ('Start processing the run found for state %s', state)
        
        if state == 'Sample Sent':
            updated_run['Sample Sent'] = []
            logger.debug('Start handling the runs in Sample Sent state')
            for run_in_sample_sent in runs_to_handle[state] :
                # get platform
                run_platform = run_in_sample_sent.get_run_platform()
                experiment_name = run_in_sample_sent.runName
                logger.info('Start handling the runs  %s in Sample Sent state', experiment_name)
                try:
                    if 'Next-Seq' in run_platform :
                        updated_run['Sample Sent'].append(manage_nextseq_in_samplesent(conn, run_in_sample_sent))
                    elif 'Mi-Seq' in run_platform :
                        updated_run['Sample Sent'].append(manage_miseq_in_samplesent(conn, run_in_sample_sent))
                    else:
                        string_message = 'Platform ' + run_platform +' is not supported '
                        logging_errors (logger, string_message , False, False)
                        continue
                except :
                    logger.info('Handling the exception to continue with the next item')
                    continue
            logger.debug('End runs in Sample Sent state')

        elif state == 'Processing Run':
            updated_run['Processing run'] = []
            logger.debug('Start handling the runs in  Processing Run state')
            for run_in_processing_run in runs_to_handle[state] :
                run_platform = run_in_processing_run.get_run_platform()
                try:
                    if 'Next-Seq' in run_platform :
                        updated_run['Processing run'].append(manage_nextseq_in_processing_run(conn, run_in_processing_run))
                    elif 'Mi-Seq' in run_platform :
                        updated_run['Processing run'].append(manage_miseq_in_processing_run(conn, run_in_processing_run))
                    else:
                        string_message = 'Platform ' + run_platform +' is not supported '
                        logging_errors (logger, string_message , False, False)
                        continue
                except :
                    logger.info('Handling the exception to continue with the next item')
                    continue
            logger.debug('End runs in Processing Run state')

        elif state == 'Processed Run':
            updated_run['Processed run'] = []
            logger.debug('Start handling the runs in  Processed Run state') 
            for run_in_processed_run in runs_to_handle[state]:
                run_platform = run_in_processed_run.get_run_platform()
                try:
                    if 'Next-Seq' in run_platform :
                        updated_run['Processed run'].append(manage_nextseq_in_processing_run(conn, run_in_processed_run))
                    elif 'Mi-Seq' in run_platform :
                        updated_run['Processed run'].append(manage_miseq_in_processing_run(conn, run_in_processed_run))
                    else:
                        string_message = 'Platform ' + run_platform +' is not supported '
                        logging_errors (logger, string_message , False, False)
                        continue
                except :
                    logger.info('Handling the exception to continue with the next item')
                    continue
            logger.debug('End runs in Processed Run state')

        elif state == 'Processing Bcl2fastq':
            for run_in_processing_bcl2fastq_run in runs_to_handle[state] :
                pass
                #updated_run[state].append( manage_run_in_processing_bcl2fast2_run(conn, run_in_processing_bcl2fastq_run))
        elif state == 'Processed Bcl2fastq':
            for run_in_processed_bcl2fastq_run in runs_to_handle[state] :
                pass
                #updated_run[state].append( manage_run_in_processed_bcl2fastq_run(conn, run_in_bcl2fastq_processed_run))
        else:
            string_message = 'Run in unexpected state. ' + state
            logging_errors (logger, string_message , False, False)
            logger.debug ('Run name is : $s ', run_to_handle[state])
            continue

    return (processed_run)

'''

def perform_xml_stats (xml_statistics, run_name_value):
    for project in xml_statistics:
        print (project)
        ### Flowcell Summary
        fl_pf_yield_sum=0
        fl_raw_yield_sum=0
        fl_mbases=0
        for values in xml_statistics[project]['PF_Yield']:
            fl_pf_yield_sum+= int(values)
        for values in xml_statistics[project]['RAW_Yield']:
            fl_raw_yield_sum+= int(values)
        for values in xml_statistics[project]['']:
            print()
'''
