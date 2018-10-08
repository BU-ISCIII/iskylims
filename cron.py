from django.conf import settings
from iSkyLIMS_wetlab import wetlab_config

from datetime import datetime
from time import strftime
from .utils.stats_calculation import *
from .utils.parsing_run_info import *
from .utils.sample_convertion import get_experiment_library_name


import os , sys, traceback
import logging
from .utils.samplesheet_checks import *
from django.conf import settings
from logging.handlers import RotatingFileHandler


def timestamp_print(message):
    starting_time= datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(starting_time+' '+message)

### TBD
def open_samba_connection():
    ## needed for testing in cuadrix
    ## to be commented-out when delivery (will use instead the function
    ## located in utils/parsing_run_info.py)

    timestamp_print('Starting the process for open_samba_connection() (cron.py- cuadrix testing)')

    ###logger.info('user ID= '+wetlab_config.SAMBA_USER_ID+'. domain= '+wetlab_config.SAMBA_DOMAIN)
    conn=SMBConnection(wetlab_config.SAMBA_USER_ID, wetlab_config.SAMBA_USER_PASSWORD,
        wetlab_config.SAMBA_SHARED_FOLDER_NAME,wetlab_config.SAMBA_REMOTE_SERVER_NAME,
        use_ntlm_v2=wetlab_config.SAMBA_NTLM_USED,domain=wetlab_config.SAMBA_DOMAIN)
    if True != conn.connect(wetlab_config.SAMBA_IP_SERVER, int(wetlab_config.SAMBA_PORT_SERVER)):
        logger=open_log('open_samba_connection_testing.log')
        logger.error('Cannot set up SMB connection with '+wetlab_config.SAMBA_REMOTE_SERVER_NAME)
        assert 1==2, 'Cannot set up SMB connection with '+wetlab_config.SAMBA_REMOTE_SERVER_NAME
    timestamp_print('Leaving open_samba_connection() (cron.py- cuadrix testing)')
    return conn

### End TBD


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



def determine_target_miseqruns(logger):
    ## Identification of runs whose samplesheets must be treated:
    ## 1st: scan of wetlab_config.SAMBA_SHARED_FOLDER_NAME to build a list
    ## with MiSeq runs

    ## 2nd: construction of sublist with runs which fullfill:
    ## a) have a samplesheet b) have not been already processed
    ## c) are not runs featuring faulty samplesheets

    timestamp_print('Starting the process to determine_target_miseqruns()')
    logger.info('Starting the process to determine_target_miseqruns()')

    ## Reading wetlab_config.SAMBA_SHARED_FOLDER_NAME (NGS_Data in production):
    try:
        conn=open_samba_connection()
        logger.info('Succesfully SAMBA connection for determine_target_miseqruns')
        file_list= conn.listPath(wetlab_config.SAMBA_SHARED_FOLDER_NAME, '/')

        file_list_filenames_debug=[x.filename for x in file_list] ##debug
        logger.debug('length of file list= '+str(len(file_list))
            +'. file list=\n'+'\n'.join(file_list_filenames_debug))
        if len(file_list) < 1:
            logger.error('Unexpected empty folder: nº of elements = len(file_list)= ',len(file_list))
            assert len(file_list) >= 1, 'Unexpected empty folder: nº of elements=len(file_list)= '+len(file_list)

    except: ##
        logger.exception('Exception when trying to set up SMB (samba) connection')
        conn.close()
        raise


    ## total available miSeq runs format
    ## {'run_dir#1':{samplesheet_filename:'samplesheet_filename#1', sequencer_family:
    ##  'sequencer_family#1', sequencer_model:'sequencer_model#1'},
    ##  'run_dir#2':{samplesheet_filename:'samplesheet_filename#2', sequencer_family:
    ##  'sequencer_family#2',sequencer_model:'sequencer_model#2'},
    ##  ...}

    temp_run_folders= {}
    target_run_folders={} ##subset of runs of temp_run_folders with MiSeq runs retained (Same format)
    faulty_samplesheet_miseqruns_file = os.path.join(
        settings.MEDIA_ROOT,wetlab_config.RUN_TEMP_DIRECTORY, wetlab_config.FAULTY_SAMPLESHEET_MISEQRUNS_FILE)


    for sfh in file_list:
        if sfh.isDirectory:
            run_dir=(sfh.filename)
            logger.debug('run_dir= '+run_dir)
            if ('.' == run_dir or '..'== run_dir):
                continue
            sequencer=re.search('_M0\d+_', run_dir) ## MiSeq run_dir_path

            if None != sequencer: ##found a MiSeq run dir
                logger.debug('sequencer information= '+sequencer.group())
                samplesheet_found=False
                sequencer_info=sequencer.group() ##sequencer string
                run_dir_file_list = conn.listPath(wetlab_config.SAMBA_SHARED_FOLDER_NAME, run_dir)
                run_dir_file_list_filenames_debug=[x.filename for x in run_dir_file_list] ##debug
                logger.debug('\tlength of run dir file list= '+str(len(run_dir_file_list))
                    +'. file list=\n\t'+'\n\t'.join(run_dir_file_list_filenames_debug)+'\n')

                for file in run_dir_file_list:
                    ## SampleSheet usually present as "SampleSheet.csv".
                    ## Once as "samplesheet.csv" ( 180725_M03352_0112_000000000-D38LV).
                    if file.filename.lower() == "samplesheet.csv":
                        samplesheet_found=True
                        temp_run_folders[run_dir]={}
                        temp_run_folders[run_dir]['samplesheet_filename']=file.filename
                        temp_run_folders[run_dir]['sequencer_model']= sequencer_info[1:-1] # MiSeq
                        logger.debug('temp_run_folders['+run_dir+']: '+str(temp_run_folders[run_dir]))
                        break
                    else:
                        continue

                if False==samplesheet_found: ## no sample sheet found in run_dir:
                    ## TODO
                    assert 1==2,'Debugging...'
                    ## endTODO
                    samplesheet_check_error_dict={'run_name':run_dir,'error':'Run without samplesheet'}
                    logger.error('Run: '+samplesheet_check_error_dict['run_name']
                        + '   '+'Error: '+samplesheet_check_error_dict['error'])
                    ####Treatment of samplesheet non existent... TODO
                    try:
                        with open (faulty_samplesheet_miseqruns_file, 'a+') as faulty_samplesheet_file:
                            faulty_samplesheet_file.write(
                                'Run: '+samplesheet_check_error_dict['run_name']
                                + '   '+'Error: ',samplesheet_check_error_dict['error'])
                            logger.error('Run: '+samplesheet_check_error_dict['run_name']
                                + 'recorded in: '+wetlab_config.FAULTY_SAMPLESHEET_MISEQRUNS_FILE)

                    except:
                        logger.exception('Exception when trying to record run '+ run_dir+
                            ' without sample sheet in file '+
                            wetlab_config.FAULTY_SAMPLESHEET_MISEQRUNS_FILE)
                        raise

            else: ##unexpected case: no sequencer information in the run directory
                logger.info(run_dir+' does not carry MiSeq sequencer info in the run name. Different sequencer? ')

        else: ##unexpected case
            logger.warning(run_dir+' is not a directory....')
        ### End TODO

    conn.close()
    timestamp_print('Samba connection closed')

    if len(temp_run_folders) <1 :  ## There are not MiSeq runs
        logger.info("No MiSeq runs at this moment")
        timestamp_print('Leaving the process to determine_target_miseqruns()')

    else:  ## analysis of the MiSeq runs list built to select the final target ones
        logger.debug('temp_run_folders: '+str(temp_run_folders))
        process_run_file_miseqelements=[]
        faulty_samplesheet_miseqruns=[]
        process_run_file = os.path.join(settings.MEDIA_ROOT,wetlab_config.RUN_TEMP_DIRECTORY, wetlab_config.PROCESSED_RUN_FILE)

        try:
            with open (process_run_file,'r') as fh:
                for line in fh:
                    line=line.rstrip()
                    sequencer=  re.search('_M0\d+_', line)
                    if None != sequencer: ##MiSeq found
                        process_run_file_miseqelements.append(line)
        except FileNotFoundError:
            logger.warning(process_run_file+ ' does not exist. Creating it...')
            try:
                with open(process_run_file,'x') as fh:
                    os.chmod(process_run_file,0o664)
                    logger.debug('Creation of empty '+process_run_file)
                    timestamp_print('Creation of empty '+process_run_file)
            except:
                logger.exception('Exception when creating empty '+process_run_file)
                raise
        except:
            logger.exception('Exception when reading file containing'
                ' processed runs. Time stop=  '+ time_stop)
            raise

        ## Getting info about runs with unexpected samplesheets
        try:
            with open (faulty_samplesheet_miseqruns_file,'r') as fh:
                for line in fh:
                    line=line.rstrip()
                    faulty_samplesheet_miseqruns.append(line)
        except FileNotFoundError:
            logger.warning(faulty_samplesheet_miseqruns_file+ ' does not exist. Creating it...')
            try:
                with open (faulty_samplesheet_miseqruns_file,'x') as fh:
                    os.chmod(faulty_samplesheet_miseqruns_file,0o664)
                    logger.info('Creation of empty '+faulty_samplesheet_miseqruns_file)
                    timestamp_print('Creation of empty '+faulty_samplesheet_miseqruns_file)
            except:
                logger.exception('Exception when creating empty '+faulty_samplesheet_miseqruns_file)
                raise

        except:
            logger.exception(
                'Exception when reading (or creating) the file containing MiSeq runs',
                ' with faulty samplesheets . Time stop=  ')
            raise

        for run,val in temp_run_folders.items():
            if (run in process_run_file_miseqelements) or (run in faulty_samplesheet_miseqruns):
                continue;
            else:
                target_run_folders[run]=val

    timestamp_print('Leaving the process to determine_target_miseqruns()')
    logger.info('Leaving the process to determine_target_miseqruns()')
    return target_run_folders



def fetch_remote_samplesheets(run_dir_dict,logger):
    ## run_dir_dict= {run_dir:{samplesheet_filename:..., sequencer_family:..., sequencer_model:...}}
    ##Open samba connection and copy the remote samplesheet locally. If exception happens, delete the (locally) copied files

    timestamp_print('Starting the process for fetch_remote_samplesheets()')
    logger.info('Starting the process for fetch_remote_samplesheets()')
    ##Treatment of the selected runs:
    transfered_samplesheet_filepaths=[] ##filepaths of transfered (locally copied) samplesheets. Used to clean up in case things go wrong
    try:
        conn=open_samba_connection()
        counter=0
        for run_index, run_info_dict in run_dir_dict.items():
            counter+=1
            ##Storage of the original samplesheet:
            samba_samplesheet_filepath = os.path.join(run_index,run_info_dict['samplesheet_filename'])
            ##Construction of local_samplesheet_filename:
            split_filename=re.search('(.*)(\.\w+$)',run_info_dict['samplesheet_filename'])
            ext_filename=split_filename.group(2)
            timestr=time.strftime('%Y%m%d-%H%M%S')
            local_samplesheet_filename = str(split_filename.group(1)+ timestr+'-#'
                +str(counter) +ext_filename)
            local_samplesheet_filepath= os.path.join(
                settings.MEDIA_ROOT, wetlab_config.RUN_SAMPLE_SHEET_DIRECTORY,
                local_samplesheet_filename)

            ## now run_dir_dict= {run_dir:{'samplesheet_filename':...,'sequencer_model':...,
            ##              ,'local_samplesheet_filepath':...}
            ##                    run_dir2:{...}...}
            run_info_dict['local_samplesheet_filepath']=local_samplesheet_filepath ## used later in the function
            logger.debug('local_samplesheet_filepath= '+local_samplesheet_filepath)
            try:
                with open (local_samplesheet_filepath, 'wb') as file_handler:
                    conn.retrieveFile(wetlab_config.SAMBA_SHARED_FOLDER_NAME,
                        samba_samplesheet_filepath, file_handler)
                    transfered_samplesheet_filepaths.append(local_samplesheet_filepath)
                    logger.info('run: '+run_index+'. Local copy of samplesheet: '+local_samplesheet_filepath)
            except:
                logger.exception('Problem when retrieving MiSeq samplesheets from '
                    + wetlab_config.SAMBA_SHARED_FOLDER_NAME+samba_samplesheet_filepath
                    + ' to local storage')
                raise

    except:
        logger.exception('Problem when opening samba connection to retrieve samplesheets')
        for transfered_file in transfered_samplesheet_filepaths:
            os.delete(transfered_file)
            logger.info('Deleted from local storage: ',transfered_file)
        raise

    finally: #always
        conn.close()

    database_info={} ## Information to be returned
    for run_index, run_info_dict in run_dir_dict.items():
        counter=0
        experiment_run_name,index_library_name=get_experiment_library_name(
            run_info_dict['local_samplesheet_filepath'])
        logger.debug('For local_samplesheet_filepath: '+local_samplesheet_filepath
            + '(samplesheet) experiment_name: ' + experiment_run_name +'. index_library: '+index_library_name)
        if experiment_run_name=='' :
            timestr=time.strftime('%Y%m%d-%H%M%S')
            experiment_run_name= timestr+'-#'+str(counter) #unique value
            logger.info('empty run_name in samplesheet: ',local_samplesheet_filepath,
                '. run_name allocated value= ',experiment_run_name)

        ##samplesheet checks:
        ## if KO: tratamiento:
        ##  registrar en 'faulty_samplesheet_ miseq_runs'
        ##  error behaviour: log? something else? TODO
        ## mail angel y nosotros
        ## filtro wetlab: "Get Incompleted" + causa error

        project_dict=get_projects_in_run(
            run_info_dict['local_samplesheet_filepath'])
        logger.debug('project_dict= '+str(project_dict))
        samplesheet_check_error_dict={} ## to register the faulty run and the error cause
        message_output=''
        if 'Free' != check_run_name_free_to_use(experiment_run_name):
            message_output= check_run_name_free_to_use(experiment_run_name)
        elif 'OK_projects_samplesheet' != check_run_projects_in_samplesheet(
                run_info_dict['local_samplesheet_filepath']):
            message_output=check_run_projects_in_samplesheet(
                run_info_dict['local_samplesheet_filepath'])
        elif 'OK_users' != check_run_users_definition(
                run_info_dict['local_samplesheet_filepath']):
            message_output=check_run_users_definition(
                run_info_dict['local_samplesheet_filepath'])
        elif 'OK_projects_db' != check_run_projects_definition(project_dict):
            message_output=check_run_projects_definition(project_dict)

        if '' != message_output:
            samplesheet_check_error_dict={'run_name':run_index,
                'experiment_run_name':experiment_run_name, 'error':message_output}
            logger.error('Run: '+samplesheet_check_error_dict['run_name']
                +'experiment run name'+experiment_run_name
                +'Error: ',samplesheet_check_error_dict['error'])
            try:
                with open (faulty_samplesheet_miseqruns_file, 'a+') as faulty_samplesheet_file:
                    faulty_samplesheet_file.write(samplesheet_check_error_dict['run_name'])
                    ### TODO To store the cause
                    ### faulty_samplesheet_file.write('Run: ',samplesheet_check_error_dict['run_name'],
                    ###    '   ','Error: ',samplesheet_check_error_dict['error'])
                    ### EndTODO
                    logger.debug('Run: '+samplesheet_check_error_dict['run_name']
                        + 'recorded in: ',faulty_samplesheet_miseqruns_file)
            except:
                logger.exception('Exception when trying to record run ', run_dir,
                    ' with faulty samplesheet in file ',
                    wetlab_config.FAULTY_SAMPLESHEET_MISEQRUNS_FILE)
                raise

        ##TODO
        ## mail angel y nosotros
        ## filtro wetlab: "Get Incompleted" + causa error
        #EndTODO

        ## TODO
        assert 1==2,'Debugging...'
        ## endTODO


        database_info[run_name]={}
        database_info[run_name]['samplesheet']=wetlab_config.RUN_SAMPLE_SHEET_DIRECTORY+samplesheet_filename
        database_info[run_name]['run_projects']=project_dict

        ## For MiSeq runs we take the 1st researcher as user for the 'center requested'
        key= next.iter(project_dict) ## 1st (and maybe only) project in run
        researcher=project_dict[key]
        database_info[run_name]['userId']=User.objects.get(username_exact=researcher)
        database_info[run_name]['index_library']=index_library_name
        database.info[run_name]['sequencer_model']=run_info_dict['sequencer_model']


    timestamp_print('Leaving the process for fetch_remote_samplesheets()')
    logger.info('Leaving the process for fetch_remote_samplesheets()')
    return database_info







def getSampleSheetFromSequencer():
    ## TODO
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
    ## EndTODO

    ##TBD
    assert 1==2, 'quitar para empezar :)'
    #End TBD
    logger=open_log('getSampleSheetFromSequencer.log')
    timestamp_print('Starting the process for getSampleSheetFromSequencer()')
    logger.info('Starting the process for getSampleSheetFromSequencer()')
    try:
        ## Launch elaboration of the list of the MiSeq samplesheets to study:
        target_run_folders= determine_target_miseqruns(logger)


        if len(target_run_folders) < 1:
            time_stop= datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            print('No new MiSeq runs available. Time stop= ',time_stop)
        else:
            ##Launch treatment of selected runs
            logger.debug('target_run_folders: '+str(target_run_folders))
            database_info=fetch_remote_samplesheets(target_run_folders,logger)
            if database_info: ## information fetched :)
                ##Store in DB the information corresponding to the fetched set of samplesheets / runs
                for key, val in database_info:

                    ##1.- table 'RunProcess' data:
                    runName= key
                    sampleSheet= val['file_name'] # file_name = wetlab_config.RUN_SAMPLE_SHEET_DIRECTORY
                                                # + split_filename.group(1)+timestr +

                    center_requested_id = Profile.objects.get(
                        profileUserID = val['userId']).profileCenter.id
                    centerRequestedBy=Center.objects.get(pk=center_requested_id)
                    index_library = val['index_library']
                    sequencerPlatformModel= val['sequencerPlatformModel']
                    runState='Recorded'
                    new_run_info=RunProcess(runName,sampleSheet,centerRequestedBy,index_library,
                        sequencerPlatformModel,runState)
                    new_run_info.save()

                    ##2.- table 'Projects' data:
                    project_dict=val['run_projects']
                    for key2, val2 in project_dict.items():
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



        timestamp_print('Leaving the process for getSampleSheetFromSequencer()')
        logger.info('Leaving the process for getSampleSheetFromSequencer()')
        timestamp_print('****** Leaving crontab')

    except:
        ## generic exception handling (exception info).
        var =traceback.format_exc()
        time_stop= datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        logger.error ('getSampleSheetFromSequencer: exception handling. Time stop= '+
            time_stop+'.  Traceback info= '+var)








def check_recorded_folder ():
    ## TODO
    ## This function will intend to take runs to the "Sample Sent" state
    ## A run in state="SampleSent" implies that :
    ##  the run primary analysis has been succesfully executed  (secondary
    ##      analysis will be performed subsequently via 'Bcl2Fastq')
    ##  and the file 'PROCESSED_RUN_FILE' has been updated with the run directory name
    ##
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
