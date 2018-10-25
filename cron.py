from django.conf import settings
from django.contrib.auth.models import User
from django.conf import settings
from .wetlab_config import *
from .models import RunProcess
from django_utils.models import Profile

from datetime import datetime
from .utils.stats_calculation import *
from .utils.parsing_run_info import *
from .utils.sample_convertion import get_experiment_library_name
from .utils.samplesheet_checks import *
from .utils.wetlab_misc_utilities import timestamp_print

import os , sys, traceback
import logging
from logging.handlers import RotatingFileHandler


### TBD ### moved to utils/wetlab_misc_utilities.py To be deleted here
'''
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
        timestamp_print('Cannot set up SMB connection with '+wetlab_config.SAMBA_REMOTE_SERVER_NAME)

    timestamp_print('Leaving open_samba_connection() (cron.py- cuadrix testing)')
    return conn
'''
### End TBD


def managed_open_file(logger, file_path,mode='r'):
    ## Function to open text files managing exception logging if something happen
    file_text_lines_list=[]
    try:
        with open (file_path,mode) as fh:
            file_text_lines_list=[line.rstrip() for line in fh]
    except FileNotFoundError:
        logger.warning(file_path+ ' does not exist. Creating it...')
        try:
            with open(file_path,'x') as fh:
                os.chmod(file_path,0o664)
                logger.debug('Creation of empty '+file_path)
                timestamp_print('Creation of empty '+file_path)
        except:
            logger.exception('Exception when creating empty '+file_path)
            raise
    except:
        logger.exception('Exception when accessing file: '+file_path)
        raise

    return file_text_lines_list



def open_log(log_name):

    log_name=os.path.join(settings.BASE_DIR, wetlab_config.LOG_DIRECTORY, log_name)
    #def create_log ():
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    #create the file handler
    ##TBD ## maxBytes=40000, backupCount=5
    handler = logging.handlers.RotatingFileHandler(log_name, maxBytes=400000, backupCount=2)
    ##EndTBD
    handler.setLevel(logging.DEBUG)

    #create a Logging format
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    #add the handlers to the logger
    logger.addHandler(handler)

    return logger



def determine_target_miseqruns(logger):
    #Determination of potential new MISEQ runs in the repository, thus excluding :
    # a) runs already BCL-finished or cancelled;
    # b) runs which presented a faulty' samplesheet in a previous iteration
    # c) runs which presented a samplesheet lacking the 'experiment name' in a previous iteration
    timestamp_print('Starting the process to determine_target_miseqruns()')
    logger.info('Starting the process to determine_target_miseqruns()')

    temp_run_folders= {}
    target_run_folders={}

    ##subset of runs of temp_run_folders with MiSeq runs retained (Same format)
    ## Reading wetlab_config.SAMBA_SHARED_FOLDER_NAME (NGS_Data in production):
    try:
        conn=open_samba_connection()
        logger.info('Succesfully SAMBA connection for determine_target_miseqruns')
        file_list= conn.listPath(wetlab_config.SAMBA_SHARED_FOLDER_NAME, '/')
        file_list=file_list[0:7] ##TBDDebugEndDebug
        file_list_filenames_debug=[x.filename for x in file_list] ##debug
        logger.debug(
            'number of existing directory runs of any kind= '+str(len(file_list)-2))## -2 ->"." and ".."
        logger.debug('run dir list=\n'+'\n'.join(file_list_filenames_debug))

    except: ##
        logger.exception('==>Exception when trying to set up SMB (samba) connection')
        timestamp_print('==>Exception when trying to set up SMB (samba) connection')
        conn.close()
        return target_run_folders##to avoid exception: we just log it

    ## total available miSeq runs format
    ## {'run_dir#1':{samplesheet_filename:'samplesheet_filename#1', sequencer_family:
    ##  'sequencer_family#1', sequencer_model:'sequencer_model#1'},
    ##  'run_dir#2':{samplesheet_filename:'samplesheet_filename#2', sequencer_family:
    ##  'sequencer_family#2',sequencer_model:'sequencer_model#2'},
    ##  ...}

    for sfh in file_list:
        if sfh.isDirectory:
            run_dir=(sfh.filename)
            #logger.debug('run_dir= '+run_dir)
            if ('.' == run_dir or '..'== run_dir):
                continue
            sequencer=re.search('_M0\d+_', run_dir) ## MiSeq run_dir_path

            if None != sequencer: ##found a MiSeq run dir
                #logger.debug('sequencer information= '+sequencer.group())
                samplesheet_found=False
                sequencer_info=sequencer.group() ##sequencer string
                run_dir_file_list = conn.listPath(wetlab_config.SAMBA_SHARED_FOLDER_NAME, run_dir)
                run_dir_file_list_filenames_debug=[x.filename for x in run_dir_file_list] ##debug
                #logger.debug('\tlength of run dir file list= '+str(len(run_dir_file_list)))
                #logger.debug('. file list=\n\t'+'\n\t'.join(run_dir_file_list_filenames_debug)+'\n')

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

                if False==samplesheet_found:
                    samplesheet_check_error_dict={'run_name':run_dir,'error':'Run without samplesheet'}
                    logger.error('Run: '+samplesheet_check_error_dict['run_name']
                        + '   '+'Error: '+samplesheet_check_error_dict['error'])

                    try:
                        registered_faulty_runs=managed_open_file(logger, wetlab_config.FAULTY_SAMPLESHEET_MISEQRUNS_FILEPATH,'r')
                        logger.debug('Existing faulty runs:\n'+'\n'.join(registered_faulty_runs))

                        if samplesheet_check_error_dict['run_name'] not in registered_faulty_runs:
                            with open (wetlab_config.FAULTY_SAMPLESHEET_MISEQRUNS_FILEPATH, 'a+') as faulty_samplesheet_file:
                                faulty_samplesheet_file.write(samplesheet_check_error_dict['run_name']+'\n')
                                logger.error('Run: '+samplesheet_check_error_dict['run_name']
                                + ' recorded in: '+wetlab_config.FAULTY_SAMPLESHEET_MISEQRUNS_FILEPATH)

                    except:
                        raise

            else:##No MiSeq
                continue

        else: #No directory
            continue

    conn.close()
    logger.debug('SMB connection closed')

    if len(temp_run_folders) <1 :  ## There are not MiSeq runs
        logger.info("No MiSeq runs at this moment")

    else:  ## analysis of the MiSeq runs list built to select the final target ones
        logger.debug('temp_run_folders: '+str(temp_run_folders))
        process_run_file_miseqelements=[]
        faulty_samplesheet_miseqruns=[]
        miseqruns_sequencing_in_progress=[]

        try:##runs which already finished the sequencing (primary analysis)
            ##in the past (can be in later states...)
            ##or have been CANCELLED
            process_run_file_miseqelements=managed_open_file(logger, wetlab_config.MISEQ_PROCESSED_RUN_FILEPATH,'r')
            logger.debug('Existing processed miseq runs:\n'+'\n'.join(process_run_file_miseqelements))

            ##runs with 'problematic'/absent samplesheets registered in previous iterations
            faulty_samplesheet_miseqruns=managed_open_file(logger,wetlab_config.FAULTY_SAMPLESHEET_MISEQRUNS_FILEPATH,'r')
            logger.debug('Existing faulty runs:\n'+'\n'.join(faulty_samplesheet_miseqruns))

            ##runs with no 'experiment name' in samplesheet registered iterations
            noexpname_samplesheet_miseqruns=managed_open_file(
                logger,wetlab_config.SAMPLESHEET_NOEXPNAME_MISEQRUNS_FILEPATH,'r')
            logger.debug('Existing no exp name runs:\n'+'\n'.join(noexpname_samplesheet_miseqruns))
        except:
            raise


        for run,val in temp_run_folders.items():
            if run in process_run_file_miseqelements:
                logger.debug('run: '+run+' is in '+wetlab_config.MISEQ_PROCESSED_RUN_FILEPATH+': ignoring...')
                continue

            elif run in faulty_samplesheet_miseqruns:
                logger.debug('run: '+run+' is in '+wetlab_config.FAULTY_SAMPLESHEET_MISEQRUNS_FILEPATH+': ignoring...')
                continue

            elif run in noexpname_samplesheet_miseqruns:
                logger.debug('run: '+run+' is in '
                    +wetlab_config.SAMPLESHEET_NOEXPNAME_MISEQRUNS_FILEPATH+': ignoring...')
                continue


            else:
                target_run_folders[run]=val

    timestamp_print('Leaving the process to determine_target_miseqruns()')
    logger.info('Leaving the process to determine_target_miseqruns()')
    return target_run_folders



def fetch_remote_samplesheets(run_dir_dict,logger):
    ## run_dir_dict= {run_dir:{samplesheet_filename:..., sequencer_model:...}}
    ##Open samba connection and copy the remote samplesheet locally. If exception happens, delete the (locally) copied files

    timestamp_print('Starting the process for fetch_remote_samplesheets()')
    logger.info('Starting the process for fetch_remote_samplesheets()')
    ##Treatment of the selected runs:
    transfered_samplesheet_filepaths=[] ##filepaths of transfered (locally copied) samplesheets. Used to clean up in case things go wrong
    try:
        conn=open_samba_connection()
        for run_index, run_info_dict in run_dir_dict.items():
            ##Storage of the original samplesheet:
            samba_samplesheet_filepath = os.path.join(run_index,run_info_dict['samplesheet_filename'])
            ##Construction of local_samplesheet_filename:
            split_filename=re.search('(.*)(\.\w+$)',run_info_dict['samplesheet_filename'])
            ext_filename=split_filename.group(2)
            timestr=datetime.datetime.now().strftime('%Y%m%d-%H%M%S.%f') ## till milliseconds :)
            local_samplesheet_filename = str(split_filename.group(1)+ timestr +ext_filename)
            local_samplesheet_filepath= os.path.join(
                settings.MEDIA_ROOT, wetlab_config.RUN_SAMPLE_SHEET_DIRECTORY,
                local_samplesheet_filename)

            run_info_dict['local_samplesheet_filepath']=local_samplesheet_filepath ## used later in the function
            ##  run_dir_dict= {run_dir:{'samplesheet_filename':...,'sequencer_model':...,
            ##              ,'local_samplesheet_filepath':...}
            ##                    run_dir2:{...}...}
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
            os.remove(transfered_file)
            logger.info('Deleted from local storage: ',transfered_file)
        raise

    finally: #always
        conn.close()
        logger.debug('SMB connection closed')

    database_info={} ## Information to be returned
    for run_index, run_info_dict in run_dir_dict.items():
        experiment_run_name,index_library_name=get_experiment_library_name(
            run_info_dict['local_samplesheet_filepath'])
        logger.debug('============================\n')
        logger.debug('run_index: '+run_index+'. run_info_dict: '+str(run_info_dict))
        logger.debug('For local_samplesheet_filepath( '+run_info_dict['local_samplesheet_filepath']
            + '):\n(samplesheet) experiment_name: ' + experiment_run_name +'\nindex_library: '
            +index_library_name)
        if ''==experiment_run_name :
            timestr=datetime.datetime.now().strftime('%Y%m%d-%H%M%S.%f') ## till milliseconds :)
            experiment_run_name= timestr

            logger.info('==> empty experiment_run_name in samplesheet: '
                +run_info_dict['local_samplesheet_filepath']
                +'.\nexperiment_run_name allocated value= '+experiment_run_name)
            try:
                with open (
                    wetlab_config.SAMPLESHEET_NOEXPNAME_MISEQRUNS_FILEPATH, 'a+') as noexpname_file:
                    noexpname_file.write(run_index+'\n')
                    logger.info('Run: '+run_index
                        + ' recorded in: '+wetlab_config.SAMPLESHEET_NOEXPNAME_MISEQRUNS_FILEPATH)

            except:
                logger.exception('Exception when trying to record run '+ run_index
                    +' with NO EXPERIMENT NAME in '
                    + wetlab_config.SAMPLESHEET_NOEXPNAME_MISEQRUNS_FILEPATH)
                raise

        ##samplesheet checks:
        logger.debug('\nBeginning Samplesheet checks:\n')
        project_dict=get_projects_in_run(
            run_info_dict['local_samplesheet_filepath'])
        logger.debug('project_dict previous to check:=\n' +str(project_dict))
        samplesheet_check_error_dict={} ## to register the faulty run and the error cause
        message_output='OK_check'

        ##pre-generation of check results so to have them available
        message_output_run_name_free_to_use=''.join(check_run_name_free_to_use(experiment_run_name))
        #logger.debug('message_output_run_name_free_to_use= '+message_output_run_name_free_to_use)

        message_output_run_projects_in_samplesheet=''.join(check_run_projects_in_samplesheet(
                run_info_dict['local_samplesheet_filepath']))
        #logger.debug('message_output_run_projects_in_samplesheet= '+''.join(message_output_run_projects_in_samplesheet))

        message_output_run_users_definition=''.join(check_run_users_definition(
                run_info_dict['local_samplesheet_filepath']))
        #logger.debug('message_output_run_users_definition= '+message_output_run_users_definition)

        message_output_run_projects_definition=''.join(check_run_projects_definition(project_dict))
        #logger.debug('message_output_run_projects_definition= '+message_output_run_projects_definition)


        if 'Free' != message_output_run_name_free_to_use:
            state=RunProcess.objects.get(runName=experiment_run_name).runState
            if 'Recorded'==state:
                logger.info('Sequencing of Run: '+run_index  + '(experiment_run_name: '
                    +experiment_run_name+' is still in progress...')
                timestamp_print('Sequencing of Run: '+run_index  + '(experiment_run_name: '
                    +experiment_run_name+') is still in progress...')

            elif 'CANCELLED'==state:
                logger.info(' Run: '+run_index  + '(experiment_run_name: '
                    +experiment_run_name+') was stored as CANCELLED')

            if 'Recorded'==state or 'CANCELLED'==state:
                message_output='OK_check' #this is a normal behaviour
                os.remove(run_info_dict['local_samplesheet_filepath'])
                logger.info('Deleted file: '+ run_info_dict['local_samplesheet_filepath'])
                continue
            else: ##having checked before if run was among sequenced ones, this should not happen
                message_output= message_output_run_name_free_to_use
                logger.info ('Problem detected in check_run_name_free_to_use')
                raise ValueError('Unexpected status value for: '+experiment_run_name)


        elif 'OK_projects_samplesheet' != message_output_run_projects_in_samplesheet:
            message_output= message_output_run_projects_in_samplesheet
            logger.info ('Problem detected in check_run_projects_in_samplesheet')

        elif 'OK_users' != message_output_run_users_definition:
            message_output=message_output_run_users_definition
            logger.info ('Problem detected in check_run_users_definition')

        elif 'OK_projects_db' !=message_output_run_projects_definition:

            message_output=message_output_run_projects_definition
            logger.info ('Problem detected in check_run_projects_definition')

        else: ##everything alright
            logger.debug('Samplesheet checks message: '+message_output)

        if 'OK_check'!= message_output: ##there has been a samplesheet check error..
            samplesheet_check_error_dict={'run_name':run_index,
                'experiment_run_name':experiment_run_name, 'error':message_output}
            timestamp_print('Run: '+samplesheet_check_error_dict['run_name']
                +'. (experiment) run name: '+experiment_run_name
                +'. Error: '+str(samplesheet_check_error_dict['error']))
            logger.error('Run: '+samplesheet_check_error_dict['run_name']
                +'. (experiment) run name: '+experiment_run_name
                +'. Error: '+str(samplesheet_check_error_dict['error']))
            try:
                with open (wetlab_config.FAULTY_SAMPLESHEET_MISEQRUNS_FILEPATH, 'a+') as faulty_samplesheet_file:
                    faulty_samplesheet_file.write(samplesheet_check_error_dict['run_name']+'\n')
                logger.debug('Run: '+samplesheet_check_error_dict['run_name']
                        + 'recorded in: ',wetlab_config.FAULTY_SAMPLESHEET_MISEQRUNS_FILEPATH)
                os.remove(run_info_dict['local_samplesheet_filepath'])
                logger.info('Deleted file: '+ run_info_dict['local_samplesheet_filepath'])

            except:
                logger.exception('Exception when trying to record run '+ run_index
                    +' with faulty samplesheet in file '
                    + wetlab_config.FAULTY_SAMPLESHEET_MISEQRUNS_FILEPATH)
                raise

        else:#samplesheet checks ok
            ## database_info= {experiment_run_name:{'run_dir':..., 'relative_samplesheet_filepath':...,
            ##      'run_projects':{experiment_run_name:researcher,...},
            ##      'userId':...,'index_library':..., 'sequencer_model':...},
            ##
            ##      experiment_run_name2:{...},
            ##      ...}
            database_info[experiment_run_name]={}
            database_info[experiment_run_name]['run_dir']=run_index
            ## RunProcess keeps samplesheet paths below '.../documents/'
            database_info[experiment_run_name]['relative_samplesheet_filepath']= run_info_dict['local_samplesheet_filepath'][len(
                os.path.join(settings.MEDIA_ROOT, wetlab_config.RUN_SAMPLE_SHEET_DIRECTORY)):]
            database_info[experiment_run_name]['run_projects']={}
            database_info[experiment_run_name]['run_projects']=project_dict

            ## For MiSeq runs we take the "1st" (and potentially only) researcher as user
            ## for later calculation of the center from which the request came

            key= next(iter(project_dict))
            researcher=project_dict[key]
            database_info[experiment_run_name]['userId']=User.objects.get(username__exact = researcher)
            database_info[experiment_run_name]['index_library']=index_library_name
            database_info[experiment_run_name]['sequencerPlatformModel']=run_info_dict['sequencer_model']

    #logger.debug('\ndatabase_info= '+str(database_info)+'\n')
    timestamp_print('Leaving the process for fetch_remote_samplesheets()')
    logger.info('Leaving the process for fetch_remote_samplesheets()')
    return database_info




def getSampleSheetFromSequencer():
    ## This function is used for sequencers (as of today, MiSeq) for which
    ## the system do not interact with the wetlab manager via web forms
    ## when dealing with the run samplesheet:it fetches it straight from the
    ## sequencer storage directory. The process is periodically kicked off by 'cron'

    ## For each new run:
    ##   -local copy of samplesheet
    ##   -samplesheet sanity checks
    ##
    ##   -At this stage, no need to ensure unique ids for sampleSheet
    ##   -No BaseSpace formatting needed for samplesheet)
    ##   -Since the protocol of the preparation of the library (LibraryKit_id) is
    ##      not provided via a form, it will be stored as "Unknown"

    logger=open_log('getSampleSheetFromSequencer.log')
    timestamp_print('Starting the process for getSampleSheetFromSequencer()')
    logger.info('Starting the process for getSampleSheetFromSequencer()')
    #assert 0==2, 'comentar para arrancar'
    try:
        ## Launch elaboration of the list of the MiSeq samplesheets to study:
        target_run_folders= determine_target_miseqruns(logger)
        logger.debug('target_run_folders: '+str(target_run_folders))

        if len(target_run_folders) < 1:
            time_stop= datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            timestamp_print('No new MiSeq runs available. Time stop= '+time_stop)
            logger.info('No new MiSeq runs available. Time stop= '+time_stop)
        else:
            ##Launch treatment of selected runs
            database_info=fetch_remote_samplesheets(target_run_folders,logger)
            logger.debug('\nnew run information to be stored in the database:\n'
                +str(database_info))
            project_dict={}
            if database_info: ## not empty
                ##Store in DB the information corresponding to the fetched set of samplesheets / runs
                for key, val in database_info.items():

                    ##1.- table 'RunProcess' data:
                    relative_samplesheet_filepath= val['relative_samplesheet_filepath']
                    center_requested_id = Profile.objects.get(
                        profileUserID = val['userId']).profileCenter.id

                    new_run_info=RunProcess(
                        runName=key,
                        sampleSheet=relative_samplesheet_filepath,
                        centerRequestedBy=Center.objects.get(pk=center_requested_id),
                        index_library = val['index_library'],
                        sequencerPlatformModel= val['sequencerPlatformModel'],
                        runState='Recorded')
                    new_run_info.save()
                    logger.info('------------------------------')
                    logger.info('new record saved in RunProcess: '
                        +RunProcess.get_runprocess_info_debug(new_run_info))
                    logger.info('------------------------------')

                    ##2.- table 'Projects' data:
                    project_dict=val['run_projects']
                    for key2, val2 in project_dict.items():
                        logger.debug('key2: '+key2)
                        logger.debug('val2: '+val2)
                        run_id=RunProcess.objects.get(runName=key)

                        new_project_info=Projects(
                            runprocess_id=run_id,
                            projectName=key2,
                            user_id=User.objects.get(username__exact=val2),
                            ##for the moment, no info about the library prep protocol
                            LibraryKit_id=LibraryKit.objects.get(libraryName__exact = 'Unknown'),
                            ## as of today, only one set of indexes is being considered
                            libraryKit=val['index_library'],
                            procState='Recorded',
                            baseSpaceFile="")


                        new_project_info.save()
                        logger.info('------------------------------')
                        logger.info('new record saved in Projects: '
                            +Projects.get_project_info_debug(new_project_info))
                        logger.info('------------------------------')

                    ##3.- Update of MISEQ runs in RECORDED state
                        try:
                            with open (wetlab_config.RECORDED_MISEQRUNS_FILEPATH, 'a+') as recorded_miseqruns_file:
                                recorded_miseqruns_file.write(val['run_dir']+'\n')
                            os.chmod(wetlab_config.RECORDED_MISEQRUNS_FILEPATH, 0o664)
                            logger.debug('Run: '+val['run_dir']
                                    + 'recorded in: ',wetlab_config.FAULTY_SAMPLESHEET_MISEQRUNS_FILEPATH)
                        except:
                            logger.exception('Exception when trying to record run '+val['run_dir']
                                +' in file '
                                + wetlab_config.RECORDED_MISEQRUNS_FILEPATH)
                            raise



            else: ## no information fetched
                time_stop= datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                logger.info('No MiSeq runs to introduce in database. Time stop= '+time_stop)
                print('No MiSeq runs to introduce in database. Time stop= '+time_stop)


        timestamp_print('Leaving the process for getSampleSheetFromSequencer()')
        logger.info('Leaving the process for getSampleSheetFromSequencer()')
        timestamp_print('****** Leaving crontab')
        logger.info('****** Leaving crontab')
    except:
        ## generic exception handling (exception info).
        var =traceback.format_exc()
        time_stop= datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        logger.error ('getSampleSheetFromSequencer: exception handling. Time stop= '+
            time_stop+'.  Traceback info= '+var)



def test_parsing_xml_files():
    logger=open_log('test_parsing_xml_files.log')
    timestamp_print('Starting the process for test_parsing_xml_files()')
    logger.info('Starting the process for test_parsing_xml_files()')
    #assert 0==3, 'comentar para arrancar'
    #open connection for a rundir
    try:
        conn=open_samba_connection()
        ##Path of original files:
        run_index='180713_M03352_0109_000000000-B29DK' ##TBDDebugEndDebug
        samba_runinfoxml_filepath = os.path.join(run_index,'RunInfo.xml')
        samba_runparametersxml_filepath = os.path.join(run_index,'runParameters.xml')
        ##Construction of local_filepaths:
        timestr=datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
        local_runinfoxml_filepath= os.path.join(
            settings.MEDIA_ROOT, wetlab_config.RUN_TEMP_DIRECTORY,'RunInfo.xml'+timestr)
        local_runparametersxml_filepath= os.path.join(
            settings.MEDIA_ROOT, wetlab_config.RUN_TEMP_DIRECTORY,'runParameters.xml'+timestr)
        try:
            with open (local_runinfoxml_filepath, 'wb') as file_handler:
                conn.retrieveFile(wetlab_config.SAMBA_SHARED_FOLDER_NAME,
                    samba_runinfoxml_filepath, file_handler)
                logger.info('run: '+run_index+'. Local copy of RunInfo.xml: '
                    +local_runinfoxml_filepath)

            with open (local_runparametersxml_filepath, 'wb') as file_handler:
                conn.retrieveFile(wetlab_config.SAMBA_SHARED_FOLDER_NAME,
                    samba_runparametersxml_filepath, file_handler)
                logger.info('run: '+run_index+'. Local copy of RunParameters.xml: '
                    +local_runparametersxml_filepath)

        except:
            logger.exception('Problem when retrieving RunInfo.xml/RunParameters.xml from:'
                + wetlab_config.SAMBA_SHARED_FOLDER_NAME+samba_runinfoxml_filepath+'\t'
                + wetlab_config.SAMBA_SHARED_FOLDER_NAME+samba_runparametersxml_filepath
                + ' to local storage')
            raise

    except:
        logger.exception('Problem when opening samba connection to retrieve samplesheets')
        for transfered_file in transfered_samplesheet_filepaths:
            os.remove(transfered_file)
            logger.info('Deleted from local storage: ',transfered_file)
        raise

    finally: #always
        conn.close()
        logger.debug('SMB connection closed')    #fetch RunInfo y RunParameters

    #parse of files: check trazas
    ### TODO see 'id' next line EndTODO
    save_miseq_run_info(local_runinfoxml_filepath,local_runparametersxml_filepath,248,logger)
    return



##TODO

def miseq_check_recorded(): ## to be integrated in common flow...
    timestamp_print('Starting the process for miseq_check_recorded')
    logger=open_log('miseq_check_recorded')

    #Build list of runs in RECORDED state
    # a) check run BCL finished?
    # b) check run CANCELLED?
    # c) ...BCL still in progress

    try:
        ##runs with no 'experiment name' in samplesheet registered iterations
        recorded_miseqruns=managed_open_file(
            logger,wetlab_config.RECORDED_MISEQRUNS_FILEPATH,'r')
        logger.debug('Existing MISEQ run dirs (Recorded state):\n'+'\n'.join(
            recorded_miseqruns))
    except:
        raise

    for run_dir in recorded_miseqruns:
        # a) check run BCL finished?

        if checkTODO:
            continue
        # b) check run CANCELLED?
        if checkTODO:
            continue

        # c) ...BCL still in progress
        #update file!!! :-)

    return




'''
def check_recorded_folder ():
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
            print('******* Exiting the check_recorder_folder module due to error when connecting to '+wetlab_config.SAMBA_SHARED_FOLDER_NAME)
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
'''
##EndTODO



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
