from datetime import datetime
from django.conf import settings

# from django.contrib.auth.models import User
import os
import traceback
from wetlab import config
from .utils.crontab_update_run import (
    search_update_new_runs,
    handle_not_completed_run,
)

from .utils.common import open_log


def looking_for_new_runs():
    """
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
    """
    time_start = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(time_start)
    working_path = settings.MEDIA_ROOT
    os.chdir(working_path)
    config_file = os.path.join(
        settings.BASE_DIR, "iSkyLIMS_wetlab", config.LOGGING_CONFIG_FILE
    )
    logger = open_log(config_file)
    logger.info("----------------------------------")
    logger.info("###########---Start Crontab-----############")
    logger.info("----------------------------------")
    logger.info("Start searching for new/updating runs")
    try:
        search_update_new_runs("crontab_request")
    except Exception as e:
        print("****** Exiting abnormal the crontab  process with errors")
        print(e)
        logger.debug("*************** ERROR *********************")
        logger.debug(e)
        logger.debug("################################")
        logger.debug(traceback.format_exc())
        logger.debug("################################")
        logger.debug("*************** ERROR *********************")
        logger.info("###########-----End Crontab--######################")
        print(traceback.print_exc())
    logger.info("Exiting the proccess for  new/updating runs")

    # looking in database for the runs that are not completed
    logger.info("----------------------------------")
    logger.info(" ***** LOOKING FOR RUNS ALREADY DEFINED ON DATABASE *****")
    logger.info("----------------------------------")
    logger.info("Start looking for uncompleted runs")
    working_path = settings.MEDIA_ROOT
    os.chdir(working_path)
    try:
        handle_not_completed_run()
    except Exception as e:
        print("****** Exiting abnormal the crontab  process with errors")
        print(e)
        logger.debug("*************** ERROR *********************")
        logger.debug(e)
        logger.debug("################################")
        logger.debug(traceback.format_exc())
        logger.debug("################################")
        logger.debug("*************** ERROR *********************")
        print(traceback.print_exc())

    time_stop = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(time_stop)
    print("****** Exiting the process for searching not completed runs")
    logger.info("###########-----End Crontab--######################")

    return


def delete_invalid_run():
    from datetime import datetime, timedelta
    import os
    from django.conf import settings
    from wetlab import config
    from wetlab.models import RunProcess, Projects

    time_start = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(time_start)
    print(
        "Starting the process for deleting runs in register state older than ",
        datetime.today(),
    )
    date_for_removing = datetime.today() - timedelta(
        days=int(config.RETENTION_TIME)
    )
    run_found_for_deleting = RunProcess.objects.filter(
        state__run_state_name="Pre-Recorded", generate_dat__lte=date_for_removing
    )

    for run_found in run_found_for_deleting:
        if Projects.objects.filter(runprocess_id=run_found).exists():
            projects_to_be_deleted = Projects.objects.filter(runprocess_id=run_found)
            for projects in projects_to_be_deleted:
                projects.delete()
                print("deleted project", projects)
        sample_sheet_file = os.path.join(
            settings.MEDIA_ROOT, run_found.get_sample_file()
        )
        if os.path.isfile(sample_sheet_file):
            os.remove(sample_sheet_file)
            print("deleted sample sheet file ")
        run_found.delete()
        print("deleted run", run_found)
    all_runs = RunProcess.objects.all()
    sample_sheet_valid_files = []

    for run in all_runs:
        sample_sheet_valid_files.append(os.path.basename(run.get_sample_file()))
    sample_sheet_folder = os.path.join(
        settings.MEDIA_ROOT, os.path.dirname(run.get_sample_file())
    )

    file_list_ss_folder = os.listdir(sample_sheet_folder)
    print("Deleting invalid sample sheets")

    for file_ss in file_list_ss_folder:
        if file_ss not in sample_sheet_valid_files:
            file_to_delete = os.path.join(sample_sheet_folder, file_ss)
            os.remove(file_to_delete)

    print("Completed deletion of invalid sample sheets")
    time_end = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(time_end)
    print("End of deleting process ")
