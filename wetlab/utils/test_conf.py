# Generic imports
import grp
import os
import pwd
# import shutil

from django.conf import settings

# Local imports
import wetlab.config
import wetlab.models
import wetlab.utils.common
import wetlab.utils.crontab_process
import wetlab.utils.crontab_update_run


def get_config_file(config_file):
    c_file = []
    try:
        with open(config_file, "r") as fh:
            for line in fh:
                line = line.replace("\n", "")
                c_file.append(line)
    except Exception:
        return
    return c_file


def get_files_attribute(directory):
    attr_files = []
    try:
        for dirpath, dirnames, filenames in os.walk(directory, topdown=True):
            for d in dirnames:
                for entry in os.scandir(directory):
                    attr_files.append(
                        [
                            oct(os.stat(entry).st_mode)[-3:],
                            pwd.getpwuid(entry.stat().st_uid).pw_name,
                            grp.getgrgid(entry.stat().st_gid).gr_name,
                            os.path.join(dirpath, entry.name),
                        ]
                    )
    except Exception:
        return
    return attr_files


def get_iSkyLIMS_settings():
    s_file = []
    settings_file = os.path.join(settings.BASE_DIR, "iSkyLIMS", "settings.py")
    try:
        with open(settings_file, "r") as fh:
            for line in fh:
                line = line.replace("\n", "")
                if "PASSWORD" in line or "SECRET_KEY" in line:
                    if "VALIDATORS" not in line:
                        split_line = line.split("=")
                        if len(split_line) > 1:
                            split_line[1] = "XXXXXXXXXXXXXXXXXX"
                            line = " = ".join(split_line)
                        else:
                            split_line = line.split(":")
                            if len(split_line) > 1:
                                split_line[1] = "XXXXXXXXXXXXXXXXXX"
                            line = " : ".join(split_line)
                s_file.append(line)
    except Exception:
        return

    return s_file


def check_access_database():
    try:
        wetlab.models.RunProcess.objects.all()
        return "OK"
    except Exception:
        return "NOK"


def check_samba_connection():
    try:
        wetlab.utils.common.open_samba_connection()
        return "OK"
    except Exception:
        return "NOK"


def folder_test_exists(folder_run_name):
    """
    Description:
        The funtion check if remote folder exists
    Input:
        folder_run_name     # folder test name in remote server

    Functions:
        open_samba_connection   # located at utils.common.py
        get_samba_application_shared_folder
    Return:
        True if folder exists
    """
    result = False
    conn = wetlab.utils.common.open_samba_connection()
    run_data_root_folder = os.path.join("/", wetlab.utils.crontab_process.get_samba_application_shared_folder())
    shared_folder = wetlab.utils.crontab_process.get_samba_shared_folder()
    run_folder_list = conn.listPath(shared_folder, run_data_root_folder)
    for sfh in run_folder_list:
        if sfh.isDirectory:
            folder_run = sfh.filename
            if folder_run == folder_run_name:
                result = True
                break
    return result


def delete_test_run(run_obj):
    """
    Description:
        The funtion delete the database information of the run test
    Input:
        run_obj
    """
    run_obj.delete()
    return


def execute_test_for_testing_run(run_test_name):
    """
    Description:
        The funtion call the functions used in the crontab to collect run information
    Input:
        run_test_name     # folder test name in remote server

    Functions:
        open_samba_connection
        get_samba_application_shared_folder
        get_samba_shared_folder


    Return:
        True if folder exists
    """
    run_result = {}
    if wetlab.models.RunProcess.objects.filter(run_name__exact=run_test_name).exists():
        run_obj = wetlab.models.RunProcess.objects.filter(run_name__exact=run_test_name).last()
        run_obj.delete()
    working_path = settings.MEDIA_ROOT
    os.chdir(working_path)
    config_file = os.path.join(
        settings.BASE_DIR, "wetlab", wetlab.config.LOGGING_CONFIG_FILE
    )
    logger = wetlab.utils.common.open_log(config_file)
    logger.info("----------------------------------")
    logger.info("###########---Start RUN Testing  -----############")
    logger.info("----------------------------------")
    wetlab.utils.crontab_update_run.search_update_new_runs(run_test_name)
    conn = wetlab.utils.common.open_samba_connection()

    # Execute 6 times to be sure it has completed all steps
    state_run_test = ["Sample Sent", "Processed Run", "Processed Bcl2fastq"]
    for state_run in state_run_test:
        run_result[state_run] = "NOK"
    if not wetlab.models.RunProcess.objects.filter(run_name__exact=run_test_name).exists():
        run_result["ERROR"] = wetlab.config.ERROR_NO_RUN_TEST_WAS_CREATED
        return run_result
    run_obj = wetlab.models.RunProcess.objects.filter(run_name__exact=run_test_name).last()
    for step in range(6):
        state = run_obj.get_state()
        if state == "ERROR":
            run_result["ERROR"] = "error"
            return run_result
        elif state == "Sample Sent":
            wetlab.utils.crontab_update_run.manage_run_in_sample_sent_processing_state(conn, [run_obj])
            if run_obj.get_state() == "ERROR":
                run_result["ERROR"] = "Error when processing run in Sample Sent state"
                break
            else:
                run_result["Sample Sent"] = "OK"
        elif state == "Processing Run":
            wetlab.utils.crontab_update_run.manage_run_in_sample_sent_processing_state(conn, [run_obj])
            if run_obj.get_state() == "ERROR":
                run_result[
                    "ERROR"
                ] = "Error when processing run in Processing Run state"
                break
            else:
                run_result["Processing Run"] = "OK"
        elif state == "Processed Run":
            wetlab.utils.crontab_update_run.manage_run_in_processed_run_state(conn, [run_obj])
            if run_obj.get_state() == "ERROR":
                run_result["ERROR"] = "Error when processing run in Processed Run state"
                break
            else:
                run_result["Processed Run"] = "OK"
        elif state == "Processing Bcl2fastq":
            wetlab.utils.crontab_update_run.manage_run_in_processing_bcl2fastq_state(conn, [run_obj])
            if run_obj.get_state() == "ERROR":
                run_result[
                    "ERROR"
                ] = "Error when processing run in Processing Bcl2fastq state"
                break
            else:
                run_result["Processing Bcl2fastq"] = "OK"
        elif state == "Processed Bcl2fastq":
            wetlab.utils.crontab_update_run.manage_run_in_processed_bcl2fastq_state(conn, [run_obj])
            if run_obj.get_state() == "ERROR":
                run_result[
                    "ERROR"
                ] = "Error when processing run in Processed Bcl2fastq state"
                break
            else:
                run_result["Processed Bcl2fastq"] = "OK"
        elif state == "Completed":
            run_result["Completed"] = "OK"
            break
        logger.info("----------------------------------")
        logger.info("###########---End RUN Testing  -----############")
        logger.info("----------------------------------")
    return run_result
