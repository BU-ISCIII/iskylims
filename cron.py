
from datetime import datetime
from .utils.stats_calculation import *
from .utils.parsing_run_info import *

import os , sys
import logging
from logging.handlers import RotatingFileHandler

def open_log(log_name):
    
    log_name=os.path.join('/srv/iSkyLIMS/wetlab/log/', log_name)
    #def create_log ():
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    #create the file handler
    handler = logging.handlers.RotatingFileHandler(log_name, maxBytes=40000, backupCount=5)
    handler.setLevel(logging.DEBUG)
    
    #create a Logging format
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    #add the handlers to the logger
    logger.addHandler(handler)
    
    return logger
    

def check_recorded_folder ():
    time_start= datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(time_start )
    working_path = '/srv/iSkyLIMS'
    print('Starting the process for recorded_folder ')
    logger=open_log('check_recorded_folder.log')
    os.chdir(working_path)
    path=os.path.join(working_path,'documents/wetlab/tmp/recorded/')
    logger.info('Looking for new runs in directory %s', path)
    
    dir_wetlab=os.getcwd()
    logger.debug('check_recorder_folder function is running on directory  %s', dir_wetlab)
    # true if there are folders under the recorded directory
    if os.listdir(path):
        # There are sample sheet files that need to be processed
        updated_run=process_run_in_recorded_state(logger)
        if updated_run == 'Error':
            logger.error('No connection is available to Flavia')
            logger.error('Exiting the process for searching run in recorded state ')
            time_stop= datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            print(time_stop)
            print('******* Exiting the check_recorder_folder module due to error when connecting to NGS_Data_test')
        else:
            for run_changed in updated_run:
                logger.info('The run  %s is now on Sample Sent state', run_changed)
            time_stop= datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            print(time_stop)
            logger.info('Exiting the check_recorded_folder')
    else:
        time_stop= datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        print(time_stop)
        logger.info( 'Exiting the crontab for record_folder. No directories under recorded folder have been found')

def check_not_finish_run():
    time_start= datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(time_start )
    print('Starting the process for searching not completed runs ')
    logger=open_log('checking_uncompleted_run.log')
    logger.info('starting execute the crontab for not finish run')
    working_path = '/srv/iSkyLIMS'
    os.chdir(working_path)
    dir_wetlab=os.getcwd()
    logger.info('Running the check_not_finish_run  module in directory %s', dir_wetlab)
    updated_run=find_not_completed_run(logger)
    for run in updated_run :
        logger.debug('Display the list of the updated_run %s', run)
    
    count=0
    for state in updated_run:
        if (updated_run[state] == "" ):
            logger.debug('found runs on %s but not found the conditions to upgrade the state', state)
        elif (updated_run[state]== 'Error'):
            logger.error('Not connection was available for state %s', state)
        else:
            for run_changed in updated_run[state]:
                logger.info('the run  %s was changed from %s', run_changed, state)
                count +=1
         
    if count == 0:
        logger.info('***** Exiting the crontab without performing any changes')
    time_stop= datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(time_stop)
    print ('****** Exiting the process for searching not completed runs')

def delete_unregister_run ():
    from datetime import datetime, timedelta
    
    time_start= datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(time_start )
    print('Starting the process for deleting runs in register state older than ', today_date)
    date_for_removing = datetime.today() - timedelta(days=days_to_subtract)
    run_found_for_deleting = RunProcess.objects.filter(runState__exact ='Pre-Recorded', generatedat__lte = date_for_removing)
    for run_found in run_found_for_deleting:
        run_id = run_found.id
        if Projects.objects.filter(runprocess_id__exact = run_id).exists():
            projects_to_be_deleted = Projects.objects.filter(runprocess_id__exact = run_id)
            for projects in projects_to_be_deleted:
                projects.delete()
        sample_sheet_file = os.paht.join(settings.MEDIA_ROOT, run_found.sampleSheet)
        os.remove(sample_sheet_file)
        print('deleting run ' , run_found.runName,'\n')
        run_found.delete()
    end_start= datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(end_time)
    print('End of deleting process ')
