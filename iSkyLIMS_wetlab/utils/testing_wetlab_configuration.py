import os
from django.contrib.auth.models import User
from django.conf import settings
import pwd, stat, os, grp, shutil
import logging

#from pathlib import Path
from iSkyLIMS_wetlab.models import *
from iSkyLIMS_wetlab.wetlab_config import *
from iSkyLIMS_wetlab.utils.generic_functions import *
from iSkyLIMS_wetlab.utils.handling_crontab_manage_run_states import *
from iSkyLIMS_wetlab.utils.handling_crontab_common_functions import *
from iSkyLIMS_wetlab.utils.update_run_state import search_update_new_runs
#from .common_run_functions import manage_run_in_processed_run, manage_run_in_processing_bcl2fastq, manage_run_in_processed_bcl2fastq

def get_config_file (config_file):
    c_file = []
    try:
        with open (config_file ,'r') as fh:
            for line in fh:
                line = line.replace('\n', '')
                c_file.append(line)
    except:
        return
    return c_file

def get_files_attribute (directory):
    attr_files = []
    try:
        for (dirpath, dirnames, filenames) in os.walk(directory, topdown=True):
            for d in dirnames:
               for entry in os.scandir(directory):
                    attr_files.append([oct(os.stat(entry).st_mode)[-3:],
                        pwd.getpwuid(entry.stat().st_uid).pw_name,
                        grp.getgrgid(entry.stat().st_gid).gr_name,
                        os.path.join(dirpath, entry.name)])
    except:
        return
    return attr_files

def get_iSkyLIMS_settings():
    s_file = []
    settings_file = os.path.join(settings.BASE_DIR, 'iSkyLIMS' , 'settings.py')
    try:
        with open (settings_file ,'r') as fh:
            for line in fh:
                line = line.replace('\n', '')
                s_file.append(line)
    except:
        return

    return s_file


def check_access_database ():
    try:
        runs = RunProcess.objects.all()
        return 'OK'
    except:
        return 'NOK'

def check_samba_connection() :

    try:
        open_samba_connection()
        return 'OK'
    except  :
        return 'NOK'


def run_exists_in_db (run_name):

    if RunProcess.objects.filter(runName__exact = run_name).exists() :
        return True
    else:
        return False


def folder_test_exists (folder_run_name):
    '''
    Description:
        The funtion check if remote folder exists
    Input:
        folder_run_name     # folder test name in remote server

    Functions:
        open_samba_connection   # located at utils.generic_functions.py
        get_samba_application_shared_folder
    Return:
        True if folder exists
    '''
    result = False
    conn = open_samba_connection()
    run_data_root_folder = os.path.join('/', get_samba_application_shared_folder() )
    shared_folder = get_samba_shared_folder()
    run_folder_list = conn.listPath( shared_folder, run_data_root_folder)
    for sfh in run_folder_list:
        if sfh.isDirectory:
            folder_run = sfh.filename
            if (folder_run == folder_run_name):
                result = True
                break
    return result


def delete_graphic_folder_if_exists(run_name):
    run_object = RunProcess.objects.get(runName__exact = run_name)
    if RunningParameters.objects.filter(runName_id = RunProcess.objects.get(runName__exact = run_name)).exists() :
        folder_graphic = RunningParameters.objects.get(runName_id = RunProcess.objects.get(runName__exact = run_name)).get_run_folder()
        if os.path.isdir(os.path.join(settings.MEDIA_ROOT,'wetlab', 'images_plot', folder_graphic)) :
            shutil.rmtree(os.path.join(settings.MEDIA_ROOT,'wetlab', 'images_plot', folder_graphic))
    return True


def execute_test_for_testing_run(run_test_name, run_test_folder):
    '''
    Description:
        The funtion call the functions used in the crontab to collect run information
    Input:
        run_test_name     # folder test name in remote server

    Functions:
        open_samba_connection               # located at utils.generic_functions.py
        get_samba_application_shared_folder # located at utils.handling_crontab_common_functions.py
        get_samba_shared_folder             # located at utils.handling_crontab_common_functions.py


    Return:
        True if folder exists
    '''
    run_result = {}

    if RunProcess.objects.filter(runName__exact = run_test_name).exists():
        run_obj = RunProcess.objects.filter(runName__exact = run_test_name).last()
        run_obj.delete()
    working_path = settings.MEDIA_ROOT
    os.chdir(working_path)
    config_file = os.path.join(settings.BASE_DIR,'iSkyLIMS_wetlab',  wetlab_config.LOGGING_CONFIG_FILE )
    logger=open_log(config_file)
    logger.info('----------------------------------')
    logger.info('###########---Start RUN Testing  -----############')
    logger.info('----------------------------------')
    search_update_new_runs()
    conn = open_samba_connection()

    # Execute 6 times to be sure it has completed all steps
    state_run_test = ['Sample Sent', 'Processed Run', 'Processed Bcl2fastq']
    for state_run in state_run_test:
        run_result[state_run] = 'NOK'
    if not RunProcess.objects.filter(runName__exact = run_test_name).exists():
        run_result ['ERROR'] = wetlab_config.ERROR_NO_RUN_TEST_WAS_CREATED
        return run_result
    run_obj = RunProcess.objects.filter(runName__exact = run_test_name).last()
    for step in range(6):
        state = run_obj.get_state()
        if state == 'ERROR' :
            run_result['ERROR'] = 'error'
            return run_result
        elif state == 'Sample Sent':
            manage_run_in_sample_sent_processing_state(conn, [run_obj])
            if run_obj.get_state() == 'ERROR':
                run_result ['ERROR'] = 'Error when processing run in Sample Sent state'
                break
            else:
                run_result['Sample Sent'] = 'OK'
        elif state == 'Processing Run':
            manage_run_in_sample_sent_processing_state(conn, [run_obj])
            if run_obj.get_state() == 'ERROR':

                run_result ['ERROR'] = 'Error when processing run in Processing Run state'
                break
            else:
                run_result['Processing Run'] = 'OK'
        elif state == 'Processed Run':
            manage_run_in_processed_run_state(conn, [run_obj])
            if run_obj.get_state() == 'ERROR':
                run_result ['ERROR'] = 'Error when processing run in Processed Run state'
                break
            else:
                run_result['Processed Run'] = 'OK'
        elif state == 'Processing Bcl2fastq':
            manage_run_in_processing_bcl2fastq_state(conn, [run_obj])
            if run_obj.get_state() == 'ERROR':
                run_result ['ERROR'] = 'Error when processing run in Processing Bcl2fastq state'
                break
            else:
                run_result['Processing Bcl2fastq'] = 'OK'
        elif state == 'Processed Bcl2fastq':
            manage_run_in_processed_bcl2fastq_state(conn, [run_obj])
            if run_obj.get_state() == 'ERROR':
                run_result ['ERROR'] = 'Error when processing run in Processed Bcl2fastq state'
                break
            else:
                run_result['Processed Bcl2fastq'] = 'OK'
        elif state == 'Completed':
            run_result['Completed'] = 'OK'
            break
        logger.info('----------------------------------')
        logger.info('###########---End RUN Testing  -----############')
        logger.info('----------------------------------')
    return run_result
