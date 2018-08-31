from django.conf import settings
from iSkyLIMS_wetlab import wetlab_config

from datetime import datetime
from .utils.stats_calculation import *
from .utils.parsing_run_info import *

import os , sys
import logging
from logging.handlers import RotatingFileHandler


def open_log(log_name):

    log_name=os.path.join(settings.BASE_DIR, wetlab_config.LOG_DIRECTORY, log_name)
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



def get_miseqruns_samplesheets_list():
    ## Function to identify stored MiSeq runs whose samplesheets must be treated

    ## 1st: scan of wetlab_config.SAMBA_SHARED_FOLDER_NAME to build a list
    ## with runs containing 'M0' (MiSeq)

    ## 2nd: construction of sublist with runs which fullfil:
    ## a) have a samplesheet b) are not terminated sequences (i.e. already
    ## finished primary analysis) c) are not runs featuring faulty samplesheets





def getSampleSheetFromSequencer():
    ## This function is used for sequencers (as of today, MiSeq) for which
    ## the system do not interact with the wetlab manager via web forms
    ## when dealing with the samplesheet:it fetches it straight from the
    ## sequencer storage directory. The process is periodically kicked off by 'crontab'
    ## So far, we just consider the case of just one "library index name"




    ## Search for new runs which have finished primary analysis
    ##    primary_analysis_run_list=
    ##      list of runs in SAMBA remote dir for our type(s) of sequencer(s):MiSeq
    ##      and still not present in the DB (otherwise, present in state
    ##          CANCELLED, RECORDED or the following ones until COMPLETED}

    ## For each run in primary_analysis_run_list:
    ##   -local copy of samplesheet
    ##    run_name, index_library_name = get_experiment_library_name(stored_file)
    ##   -sanity checks
    ##    if ko, delete sampleSheet (cómo avisar al usuario (mail)? vs CANCELLED)
    ##
    ##   -Ensure unique ids for sampleSheetç
    ##   (-No BaseSpace formatting needed for samplesheet)
    ##   -update DB tables  (run state= "Recorded"  (run state= "Recorded"))
    ##      Since the protocol of the preparation of the library (library kit) is
    ##      not provided via a form, it will be stored as "Unknown"




            ## Calculate project_list (--> profileUserIDs) -- take center of the 1st user by default
            center_requested_id = Profile.objects.get(profileUserID = request.user).profileCenter.id
            center_requested_by = Center.objects.get(pk = center_requested_id)
            run_proc_data = RunProcess(
                runName=run_name,sampleSheet= file_name, runState='Pre-Recorded',
                centerRequestedBy = center_requested_by)
            run_proc_data.save()
            experiment_name = '' if run_name == timestr else run_name

            ## create new project record based on the project involved in the run and
            ## include the project information in projects variable to build the new FORM
            run_info_values ={}
            run_info_values['experiment_name'] = experiment_name
            run_info_values['index_library_name'] = index_library_name
            for key, val  in project_list.items():
                userid=User.objects.get(username__exact = val)
                p_data=Projects(runprocess_id=RunProcess.objects.get(runName =run_name), projectName=key, user_i
                p_data.save()
                projects.append([key, val])
            run_info_values['projects_user'] = projects
            run_info_values['runname']= run_name
            ## Get the list of the library kit used (libraryKit)
            #import pdb; pdb.set_trace()
            used_libraries = []
            list_libraries = LibraryKit.objects.order_by().valuejs_list('libraryName', flat=True)
            run_info_values['used_libraryKit'] =  list_libraries



           for p in range(len( projects)):
               my_project = projects [p]
               my_name = user_name[p]
               my_libkit = library_kit[p]
               library_kit_id = LibraryKit.objects.get(libraryName__exact = library_kit[p])
               update_info_proj=Projects.objects.get(projectName = my_project)
               update_info_proj.libraryKit=project_index_kit[p]
               update_info_proj.baseSpaceFile=bs_file[project_index_kit[p]]
               update_info_proj.LibraryKit_id = library_kit_id
               update_info_proj.proState='Recorded'
               update_info_proj.save()
           results.append(['runname', experiment_name])
           ## save the sample sheet file under tmp/recorded to be processed when run folder was created
           subfolder_name=str(run_p.id)
           #base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


           # update the sample sheet with the experiment name
           if run_name != experiment_name :
               update_sample_sheet (in_file, experiment_name)
           ## update the Experiment name and the state of the run to 'Recorded'
           run_p.runName = experiment_name
           run_p.index_library = run_index_library_name
           run_p.runState='Recorded'
           run_p.save()
           ## update the project state to Recorded
           project_to_be_updated = Projects.objects.filter(runprocess_id__exact = run_p.id)
           for project in project_to_be_updated :
               project.procState='Recorded'
               project.save()




            ##
            ##       + state= "SampleSent"
            ##

    ##  2.2.if ko, log problem on file and show message on console  ???? TODO .. EndTODO
    ##  2.3 if ok, update database tables: RunProcess, Projects. Run state="Recorded??"
    ##
    ## TODO
    ## ...
    ## endTODO




def check_recorded_folder ():
    ## This function will intend to take runs to the "Sample Sent" state
    ## A run in state="SampleSent" implies that :
    ##  the run primary analysis has been succesfully executed  (secondary
    ##      analysis will be performed subsequently via 'Bcl2Fastq')
    ##  and the file 'PROCESSED_RUN_FILE' has been updated with the run directory name
    ##
    ## TODO
    ## if case 'MiSeq':
    ##  -checks, database update and update state (getSampleSheetFromSequencer())
    ##



    ## case "NextSeq" (and all???):
    ## end TODO
    time_start= datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(time_start )
    working_path = settings.MEDIA_ROOT
    print('Starting the process for recorded_folder ')
    logger=open_log('check_recorded_folder.log')
    os.chdir(working_path)
    path=os.path.join(working_path,wetlab_config.RUN_TEMP_DIRECTORY_RECORDED )
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
