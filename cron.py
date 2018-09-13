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



def get_miseqruns_samplesheets():
    ## Identification of runs whose samplesheets must be treated:
    ## 1st: scan of wetlab_config.SAMBA_SHARED_FOLDER_NAME to build a list
    ## with MiSeq runs

    ## 2nd: construction of sublist with runs which fullfill:
    ## a) have a samplesheet b) have not been already processed
    ## c) are not runs featuring faulty samplesheets


    ## Reading wetlab_config.SAMBA_SHARED_FOLDER_NAME (NGS_Data in production):
    try:
        conn=open_samba_connection()
        if (conn != True):  ##https://pysmb.readthedocs.io/en/latest/api/smb_SMBConnection.html
            logger.error('Error when trying to set up SMB (samba) connection. Potential authentication error')
            raise Exception ('Error when trying to set up SMB (samba) connection. Potential authentication error')

        logger.info('Succesfully SAMBA connection for get_miseqruns_samplesheetsr'
        file_list = conn.listPath(wetlab_config.SAMBA_SHARED_FOLDER_NAME, '/')
        if len(file_list) < 1:
            logger.error('Unexpected empty folder: nº of elements= ',len(file_list))
            raise Exception ('Unexpected empty folder: nº of elements= ',len(file_list))

    except: ##
        print('Exception when trying to set up SMB (samba) connection')
        raise
    finally: #always: either if execution of try was OK or KO
        conn.close()

    ## total available miSeq runs format
    ## {run_dir=... ,{samplesheet_filename=..., sequencer_family=..., sequencer_model= ....}}
    temp_run_folders= {}
    target_run_folders={}: ##subset of runs of temp_run_folders with MiSeq runs retained (Same format)
    unexpected_samplesheet_miseqruns = os.path.join(
        wetlab_config.RUN_TEMP_DIRECTORY, wetlab_config.UNEXPECTED_SAMPLESHEET_MISEQRUNS_FILE)


    for sfh in file_list:
        if sfh.isDirectory:
            run_dir=(sfh.filename)
            if ('.' == run_dir or '..'== run_dir):
                continue
            sequencer_info=re.search('_M0\d+_', sfh.filename) ## MiSeq run_dir_path

            if None != sequencer_info
                run_dir_path = os.path.join(wetlab_config.SAMBA_SHARED_FOLDER_NAME ,run_dir)
                file_list = conn.listPath( run_dir_path, '/')
                samplesheet_found=false
                for file in file_list:
                    ##  E samplesheet.csv   ?
                    ## usually present as "SampleSheet.csv".
                    ## Once as "samplesheet.csv" ( 180725_M03352_0112_000000000-D38LV).
                    if file.filename.lower() == "samplesheet.csv":
                        samplesheet_found = True
                        temp_run_folders.append(run_dir)
                        temp_run_folders[run_dir]={}
                        temp_run_folder[run_dir]['samplesheet_filename']=file.filename
                        temp_run_folder[run_dir]['sequencer_......]= ## MiSeq
                        break

                    else: ## no sample sheet found:
                        logger.error('Run',run_dir,'without samplesheet').
                        ####Treatment of samplesheet non existent... TODO
                        try:
                            with open (unexpected_samplesheet_miseqruns, "a") as unexpected_samplesheet_file:
                                unexpected_samplesheet_file.write(run_dir)
                        except:
                            print('Exception when trying to record run ', run_dir,
                                ' without sample sheet in file ',
                                wetlab_config.UNEXPECTED_SAMPLESHEET_MISEQRUNS_FILE)
                            raise

                        logger.error('Recorded in: ', wetlab_config.PROCESSED_RUN_FILE)
                        print('Run',run_dir,'without samplesheet'. Recorded in: ', wetlab_config.PROCESSED_RUN_FILE)
                        ### End TODO

    if len(temp_run_folders) < 1:  ## There are not MiSeq runs
        print("No MiSeq runs at this moment")
        logger.info("No MiSeq runs at this moment")
        return miseqruns_update_dir_list


    process_run_file_miseqelements=[]
    process_run_file = os.path.join(wetlab_config.RUN_TEMP_DIRECTORY, wetlab_config.PROCESSED_RUN_FILE)

    try:
        with open (process_run_file,'r') as fh:
            for line in fh:
                line=line.rstrip()
                if 'M0' in line:
                    process_run_file_miseqelements.append(line)
    except:
        time_stop= datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        logger.error('Exception when opening (and scanning) the file containing the processed runs. Time stop=  ', time_stop)
        print('Exception when opening (and scanning) the file containing the processed runs. Time stop= ', time_stop)
        raise

    ##TODO
    ## Getting info about runs with unexpected samplesheets
    try:
        with open (unexpected_samplesheet_miseqruns,'r') as fh:
            for line in fh:
                line=line.rstrip()
                if 'M0' in line: ## TODO endTODO
                    process_run_file_miseqelements.append(line)
    except:
        time_stop= datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        logger.error('Exception when opening (and scanning) the file containing the processed runs. Time stop=  ', time_stop)
        print('Exception when opening (and scanning) the file containing the processed runs. Time stop= ', time_stop)
        raise
###End TODO

    for run,val in temp_run_folders:
        if (run in process_run_file_miseqelements) or (run in unexpected_samplesheet_ miseqruns):
            continue;
        else:
            run_dir_name = os.path.join(wetlab_config.SAMBA_SHARED_FOLDER_NAME ,run)
            file_list = conn.listPath( run_dir_name, '/')
                for file in file_list:
                    if file == "SampleSheet.csv"exactly:
                        target_run_folders.append(run)=val
                        break
                    else:
                        continue

    return target_run_folders



def fetch_remote_samplesheets(run_dir_dict):
    ## run_dir_dict= {run_dir:{samplesheet_filename:..., sequencer_family:..., sequencer_model:...}}
    ##Open samba connection and copy the remote samplesheet locally. If exception happens, delete the (locally) copied files

    ##Treatment of the selected runs:
    local_samplesheet_filepaths=[] ##filepaths of transfered (locally copied) samplesheets. Used to clean up in case things go wrong
    try:
        conn=open_samba_connection()
        for run_index, run_info_dict in run_dir_dict:
            ##Storage of the original samplesheet:
            samba_samplesheet_filepath = os.path.join(run_index,runinfodict[samplesheet_filename])
            ##Construction of local_samplesheet_filename:
            split_filename=re.search('(.*)(\.\w+$)',run_info_dict[samplesheet_filename])
            sequencer=  re.search('_M0\d+_', run_index)
            if None==sequencer: ##this should be impossible...
                print ('Sequencer == None in run name= ',run_index)
                logger.error
                raise exception RunTimeError?
            ext_filename=split_filename.group(2)
            local_samplesheet_filename = str(split_filename.group(1)+ timestr +ext_file)
            local_samplesheet_filepath= os.path.join(
                settings.MEDIA_ROOT, wetlab_config.RUN_SAMPLE_SHEET_DIRECTORY,
                local_samplesheet_filename)

            run_info_dict.append['samplesheet_filepath']=local_samplesheet_filepath ## used later in the function
            ## now run_dir_dict= {run_dir:{samplesheet_filename:..., sequencer_family:...,
            ##              sequencer_model:..., local_samplesheet_filepath:...}}
            try:
                with open (local_samplesheet_filepath, 'wb') as file_handler:
                    conn.retrieveFile(wetlab_config.SAMBA_SHARED_FOLDER_NAME,
                        samba_samplesheet_filepath, file_handler)
                    transfered_samplesheet_filepaths.append(local_samplesheet_filepath)
            except:
                exception treatment:
                raise

    except:
        exception treatment:
            for transfered_file in transfered_samplesheet_filepaths:
                os.delete(transfered_file)
                logger
            raise

    finally: #always
        conn.close()

    database_info={} ## Information to be returned
    for run_index, run_info_dict in run_dir_list:

        run_name,index_library_name=get_experiment_library_name(
            run_info_dict['local_samplesheet_filepath'])
        if run_name='' :
            run_name= timestr #unique value
            print('empty run_name in samplesheet: ',local_samplesheet_filepath,
                '. run_name allocated timestamp value= ',timestr)
            logger.error

        ##samplesheet checks:
        ## if KO: tratamiento:
        ##  registrar en 'unexpected_samplesheet_ miseq_runs'
        ##  error behaviour: log? something else? TODO
        ## mail angel y nosotros
        ## filtro wetlab: "Get Incompleted" + causa error
        check_run_name_free_to_use(run_name)
        check_run_projects_in_samplesheet(run_info_dict['local_samplesheet_filepath'])
        check_run_users_definition(run_info_dict['local_samplesheet_filepath'])

        project_dict=get_projects_in_run(stored_file) ##project_dict is a dict
        check_run_projects_definition(project_dict)

        database_info['run_name']={}
        database_info['run_name']['samplesheet']=wetlab_config.RUN_SAMPLE_SHEET_DIRECTORY+samplesheet_filename
        database_info['run_name']['run_projects']=project_dict

        ## For MiSeq runs we take the 1st researcher as user for the 'center requested'
        key= next.iter(project_dict) ## 1st (and maybe only) project in run
        researcher=project_dict[key]
        database_info['run_name']'userId']=User.objects.get(username_exact=researcher)
        database_info['run_name']['index_library']=index_library_name
        database.info['sequencer_model']=run_info_dict['sequencer_model']


    return database_info







def getSampleSheetFromSequencer():
    ## This function is used for sequencers (as of today, MiSeq) for which
    ## the system do not interact with the wetlab manager via web forms
    ## when dealing with the run samplesheet:it fetches it straight from the
    ## sequencer storage directory. The process is periodically kicked off by 'cron'
    ## So far, we just consider the case of just one "library index name"

    ## For each run in primary_analysis_run_list:
    ##   -local copy of samplesheet
    ##    run_name, index_library_name = get_experiment_library_name(stored_file)
    ##   -sanity checks
    ##    if ko, delete sampleSheet (cómo avisar al usuario (mail)? vs CANCELLED)
    ##
    ##   -At this stage, no need to ensure unique ids for sampleSheet
    ##   -No BaseSpace formatting needed for samplesheet)
    ##   -Since the protocol of the preparation of the library (LibraryKit_id) is
    ##      not provided via a form, it will be stored as "Unknown"

    try:
        ## Launch elaboration of the list of the MiSeq samplesheets to study:
        target_run_folders= get_miseqruns_samplesheets()


        if len(target_run_folders) < 1:
            time_stop= datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            logger.info('No new MiSeq runs available. Time stop= ',time_stop)
            logger.info('****** Leaving crontab')
            print('No new MiSeq runs available. Time stop= ',time_stop)
            print('****** Leaving crontab')
        else:
            ##Launch treatment of selected runs
            database_info=fetch_remote_samplesheets(temp_miseqrun_folders)
            if database_info: ## information fetched :)
                ##Store in DB the information corresponding to the fetched set of samplesheets / runs
                for key, val in database_info:

                    ##1.- table 'RunProcess' data:
                    runName= key
                    sampleSheet= val[file_name] # file_name = wetlab_config.RUN_SAMPLE_SHEET_DIRECTORY
                                                # + split_filename.group(1)+timestr +

                    center_requested_id = Profile.objects.get(
                        profileUserID = val[userId]).profileCenter.id
                    centerRequestedBy=Center.objects.get(pk=center_requested_id)
                    index_library = val['index_library']
                    sequencerPlatformModel= val['sequencerPlatformModel']
                    runState='Recorded'
                    new_run_info=RunProcess(runName,sampleSheet,centerRequestedBy,index_library,
                        sequencerPlatformModel,runState)
                    new_run_info.save()

                    ##2.- table 'Projects' data:
                    project_dict=val['run_projects']
                    for key2, val2 in project_dict.items()
                        runprocess_id=RunProcess.objects.get(runName=key)
                        projectName= key2
                        user_id=User.objects.get(username_exact=val2['Description'])
                        ##for the moment, no info about the library prep protocol
                        LibraryKit_id=LibraryKit.objects.get(libraryName__exact = 'Unknown')
                        ## as of today, only one set of indexes is being considered
                        libraryKit=val['index_library']
                        procState='Recorded'
                        new_project_info=Projects(
                            runprocess_id,projectName,user_id,LibraryKit_id,procState)
                        new_project_info.save()

            else: ## no information fetched
                time_stop= datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                logger.info('No MiSeq runs to introduce in database. Time stop= ',time_stop)
                logger.info('****** Leaving crontab')
                print('No MiSeq runs to introduce in database. Time stop= ',time_stop)
                print('****** Leaving crontab')

    except:
        ## generic exception handling (exception info).
        var =traceback.format_exc()  " "
        time_stop= datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        print('getSampleSheetFromSequencer: exception handling. Time stop= ', time_stop,'.  Traceback info= ',var)
        logger.error ('getSampleSheetFromSequencer: exception handling. Time stop= ', time_stop,'.  Traceback info= ',var)








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
            print('******* Exiting the check_recorder_folder module due to error when connecting to ',wetlab_config.SAMBA_SHARED_FOLDER_NAME)
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
