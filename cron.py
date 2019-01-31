from datetime import datetime
from django.conf import settings
#from django.contrib.auth.models import User

from .wetlab_config import *
from .utils.update_run_state import search_update_new_runs, search_not_completed_run

from .utils.generic_functions import  open_log

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

        The second one "search_not_completed_run" will check different files
        (depending on the state of the run ) on the run remote folder,
        to fetch the required information and moving the run into steps from
        Sample Sent towards Completed
    Functions:
        search_update_new_runs # located at utils.update_run_state
        search_not_completed_run # located at utils.update_run_state
        open_log    # located in utils.run_common_functions
    Constants:
        LOG_NAME_MISEQ_FETCH_SAMPLE_SHEET
        MEDIA_ROOT
    Variables:
        logger          # contain the log object
        new_runs_updated  # will have the run names for the runs that
                        were in Recorded state and they have been updated
                        to Sample Sent state
        updated_runs  # will have the run names for the runs that were
                        processed. It will use to have a summary of the
                        updated runs.
        working_path    # contains the path folder defined on MEDIA_ROOT
    Return:
        None
    '''
    working_path = settings.MEDIA_ROOT
    os.chdir(working_path)
    logger=open_log()
    logger.info('###########---Start Crontab-----############')
    logger.info('Start searching for new/updating runs')
    
    new_runs_updated, run_with_error = search_update_new_runs ()
    logger.info('------- Printing summary for search_update_new_runs -----')
    if len (new_runs_updated) > 0:
        for new_run in new_runs_updated :
            logger.info('%s has been updated in database', new_run)
    else:
        logger.info('No Run has been updated ')

    if len (run_with_error) > 0:
        for error_run in run_with_error :
            logger.info('%s has been set to error state', error_run)
    logger.info('------- End summary for search_update_new_runs -----')
    logger.info('Exiting the proccess for  new/updating runs')

    # looking in database for the runs that are not completed
    logger.info('----------------------------------')
    logger.info('----------------------------------')
    logger.info('Start looking for uncompleted runs')
    working_path = settings.MEDIA_ROOT
    os.chdir(working_path)

    updated_runs, run_with_error = search_not_completed_run()
    logger.info('------ Printing the summary result for the manage runs -----')
    if len (updated_runs) > 0 :
        for state in updated_runs:
            if len (updated_runs[state]) > 0 :
                for run_changed in updated_runs[state]:
                    logger.info('Run  %s was  processed on  %s  state', run_changed, state)
            else:
                logger.info('There is no updated run for %s ', state)
    else:
        logger.info('There are no updated runs ')

    if len (run_with_error) > 0 :
        for state in run_with_error:
            if len (run_with_error[state]) > 0 :
                for run_error in run_with_error[state]:
                    logger.info('Run  %s was set to Error when processing run on %s state', run_error, state)

    logger.info('------- End summary for search_update_new_runs -----')
    time_stop= datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(time_stop)
    print ('****** Exiting the process for searching not completed runs')
    logger.info('###########-----End Crontab--######################')
    return



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
