from datetime import datetime
from django.conf import settings
from django.contrib.auth.models import User
from django.conf import settings
from .wetlab_config import *
from .utils.update_run_state import search_update_new_runs, hanlding_not_completed_run


from .utils.run_common_functions import  open_log

import os , sys, traceback,errno



def looking_for_new_runs ():
    '''
    Description:
        The function is called from crontab to find and update new runs
        It is split in 2 main functions. 

        The first one  "search_update_new_runs" will look for new miSeq
        runs and move to Recorded state. For NextSeq runs will copy the
        sample sheet to remote folder.
        For both types of run information is collected from the run folder
        files to store in database.
        
        The second one "hanlding_not_completed_run" will check the files
        on the run remote folder to fetch the information and moving 
        into run steps from Sample Sent towards Completed
    Functions:
        search_new_miseq_runs # located at utils.update_run_state
        hanlding_not_completed_run # located at utils.update_run_state
        open_log    # located in utils.wetlab_misc_utilities
    Constants:
        LOG_NAME_MISEQ_FETCH_SAMPLE_SHEET
        MEDIA_ROOT
    Variables:
        logger          # contain the log object 
        new_miseq_runs  # will contains the run names for the runs that 
                        were in Recorded state and they have been updated
                        to Sample Sent state
    Return:
        None
    '''
    working_path = settings.MEDIA_ROOT
    os.chdir(working_path)
    logger=open_log(LOG_NAME_MISEQ_FETCH_SAMPLE_SHEET)
    logger.info('###########---Start Crontab-----############')
    logger.info('Start searching for new/updating runs')

    new_runs_updated = search_update_new_runs ()
    if len(runs_updated) > 0:
        for new_run in new_runs_updated :
            logger.info('%s has been updated in database', new_run)
    else:
        logger.info ('No new runs have been found ')
    logger.info('Exiting the proccess for  new/updating runs')
    
    # looking in database for the runs that are not completed
    logger.info('----------------------------------')
    logger.info('Start looking for uncompleted runs')
    working_path = settings.MEDIA_ROOT
    os.chdir(working_path)
    
    updated_runs = hanlding_not_completed_run()


    for state in updated_runs:
        if (updated_runs[state] == "" ):
            logger.debug('found runs on %s but not found the conditions to upgrade the state', state)
        elif (updated_run[state]== 'Error'):
            logger.error('Not connection was available for state %s', state)
        else:
            for run_changed in updated_run[state]:
                logger.info('the run  %s was changed from %s', run_changed, state)
                count +=1


    time_stop= datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(time_stop)
    print ('****** Exiting the process for searching not completed runs')
    logger.info('###########-----End Crontab--######################')
    
def update_run_in_recorded_state ():
    '''
    Description:
        The function is called from crontab to check if there are runs 
        in recorded state.
    Functions:
        handle_run_in_recorded_state # located in utils.update_run_state
        open_log    # located in utils.wetlab_misc_utilities
    Variables:
        logger          # contain the log object 
        updated_runs    # will contains the run names for the runs that 
                        were in Recorded state and they have been updated
                        to Sample Sent state
    Return:
        None
    '''
    working_path = settings.MEDIA_ROOT
    os.chdir(working_path)
    logger=open_log(LOG_NAME_RUN_IN_RECORDED_STATE)
    logger.info('Starting to check runs in recorded state')
    # check if there are runs in recorded state
    if RunProcess.objects.filter(runState__exact = 'Recorded').exists():
        logger.info('Processing the run in recorded state.')
        updated_runs = handle_runs_in_recorded_state(logger)
        if updated_runs == 'Error':
            print('******* ERROR ********')
            print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
            print('When processing run in recorded state. Check log for detail information')
        else:
            for run_changed in updated_runs:
                logger.info('The run  %s is now on Sample Sent state', run_changed)
            print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
            logger.info('Exiting update_run_in_recorded_state')
    else:
        logger.info( 'Exiting the crontab for record_folder. No runs in recorded state')


def check_not_finish_run():
    time_start= datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(time_start )
    print('Starting the process for searching not completed runs ')
    logger=open_log('checking_uncompleted_run.log')
    logger.info('starting execute the crontab for not finish run')
    working_path = settings.MEDIA_ROOT
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
    print('Starting the process for deleting runs in register state older than ', datetime.today())
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
