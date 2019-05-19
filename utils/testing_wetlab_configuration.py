from django.contrib.auth.models import User
from django.conf import settings
import pwd, stat, os, grp, shutil
import logging

#from pathlib import Path
from iSkyLIMS_wetlab.models import *
from iSkyLIMS_wetlab.wetlab_config import *
from iSkyLIMS_wetlab.utils.generic_functions import *
from iSkyLIMS_wetlab.utils.update_run_state import handle_miseq_run
from iSkyLIMS_wetlab.utils.miseq_run_functions import *
from iSkyLIMS_wetlab.utils.nextseq_run_functions import *
from iSkyLIMS_wetlab.utils.common_run_functions import manage_run_in_processed_run, manage_run_in_processing_bcl2fastq, manage_run_in_processed_bcl2fastq
from iSkyLIMS_wetlab.utils.sample_sheet_utils import *

def get_config_file (config_file):
    c_file = []
    try:
        with open (config_file ,'r') as fh:
            for line in fh:
                if 'PASSWORD' in line:
                    hide_passwd = line.split('=')
                    hide_passwd[1] = 'XXXXXXXXXXXXXXXXX'
                    line = ' = '.join(hide_passwd)
                line = line.replace('\n', '')
                c_file.append(line)
    except:
        return
    return c_file

def get_files_attribute (directory):
    attr_files = []
    '''
    #try:
    for entry in os.scandir(directory):
        attr_files.append([oct(os.stat(entry).st_mode)[-3:], 
                    pwd.getpwuid(entry.stat().st_uid).pw_name, 
                    grp.getgrgid(entry.stat().st_gid).gr_name, 
                    entry.name])
    import pdb; pdb.set_trace()
    '''
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

def check_access_database ():
    try:
        runs = RunProcess.objects.all()
        return True
    except:
        return

def check_samba_connection() :
    try:
        open_samba_connection()
        return True
    except :
        return 



def run_exists_in_db (run_name):
    
    if RunProcess.objects.filter(runName__exact = run_name).exists() :
        return True
    else:
        return False

def delete_graphic_folder_if_exists(run_name):
    run_object = RunProcess.objects.get(runName__exact = run_name)
    if RunningParameters.objects.filter(runName_id = RunProcess.objects.get(runName__exact = run_name)).exists() :
        #import pdb; pdb.set_trace()
        folder_graphic = RunningParameters.objects.get(runName_id = RunProcess.objects.get(runName__exact = run_name)).get_run_folder()
        if os.path.isdir(os.path.join(settings.MEDIA_ROOT,'wetlab', 'images_plot', folder_graphic)) :
            shutil.rmtree(os.path.join(settings.MEDIA_ROOT,'wetlab', 'images_plot', folder_graphic))
    return True

def delete_run_in_db (run_name):
    
    if RunProcess.objects.filter(runName__exact = run_name).exists() :
        delete_run = RunProcess.objects.get(runName = run_name)
        sample_sheet_file = delete_run.get_sample_file()
        if sample_sheet_file != '' :
            full_path_sample_sheet_file = os.path.join(settings.MEDIA_ROOT, sample_sheet_file)
            os.remove(full_path_sample_sheet_file)
        delete_run.delete()
    return True

def delete_test_run (run_name) :
    if run_exists_in_db(run_name):
        delete_graphic_folder_if_exists(run_name)
        delete_run_in_db(run_name)
    return True

def folder_run_exists (conn, folder_run_name):
    
    try:
        run_folder_list = conn.listPath( wetlab_config.SAMBA_SHARED_FOLDER_NAME, run_data_root_folder)
    except :
         raise
    
    return conn

def create_project (p_name, run_name, bs_file):
    #import pdb; pdb.set_trace()
    n_project = Projects( runprocess_id = RunProcess.objects.get(runName__exact = run_name),
                    LibraryKit_id = LibraryKit.objects.get(libraryName__exact = 'Nextera XT Index Kit v2 Set B'),
                    projectName = p_name , libraryKit='Nextera XT V2 ', baseSpaceFile = bs_file,
                    user_id = User.objects.get(username__exact = 'Test_user1'))
    n_project.save()
    
    
    return True

def create_run_test_nextseq_in_recorded (run_folder, experiment_name ) :
    recorded_results = []
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function create_run_test_nextseq')
    # if run exist, then delete if before created it again
    if RunProcess.objects.filter(runName__exact = experiment_name).exists() :
        run_name = RunProcess.objects.get(runName__exact = experiment_name)
        sample_sheet = run_name.get_sample_file()
        temp_folder_sample_sheet = run_name.get_run_id()
        try:
            os.remove(os.path.join(wetlab_config.RUN_SAMPLE_SHEET_DIRECTORY, sample_sheet))
        except:
            logger.info('unable to delete sample sheet')
        try:
            shutil.rmtree(os.path.join(settings.MEDIA_ROOT, wetlab_config.RUN_TEMP_DIRECTORY_RECORDED, temp_folder_sample_sheet))
        except:
            logger.info('unable to delete tmp forlder with the smample sheet')
            
        delete_test_run(experiment_name)
    
    conn = open_samba_connection()
    
    l_sample = os.path.join(wetlab_config.RUN_TEMP_DIRECTORY,wetlab_config.SAMPLE_SHEET)
    s_sample= os.path.join(run_folder, wetlab_config.SAMPLE_SHEET)

    try:
        l_sample = fetch_remote_file (conn, run_folder, s_sample, l_sample)
        logger.info('Sucessfully fetch of Sample Sheet file')
    except Exception as ex:
        string_message = 'Unable to fetch the Sample Sheet file'
        logging_errors(string_message, False, True)
        recorded_results.append((str(ex) , 'NOK'))
        return recorded_results, 'NOK'
    projects = get_projects_in_run(l_sample)
    recorded_results.append(('Successful creation Run in database', 'OK'))
    
    center_requested_by = Center.objects.get(pk = 2)
    run_process = RunProcess(runName=experiment_name,sampleSheet= '',
                                state = RunStates.objects.get(runStateName__exact = 'Recorded'),
                                index_library = 'Nextera XT', centerRequestedBy = center_requested_by)
    run_process.save()
    
    run_process_id = run_process.get_run_id()
    dir_to_store_sample_sheet = os.path.join(settings.MEDIA_ROOT, wetlab_config.RUN_TEMP_DIRECTORY_RECORDED, run_process_id)
    os.makedirs(dir_to_store_sample_sheet)
    shutil.copyfile(l_sample, os.path.join(dir_to_store_sample_sheet, 'samplesheet.csv'))
    
    recorded_results.append(('Successful copy of Sample Sheet to temporary folder', 'OK'))

    new_sample_sheet_name = store_sample_sheet_in_run (l_sample, experiment_name )
    recorded_results.append(('Successful update of Sample Sheet in database', 'OK'))

    new_sample_sheet_file = os.path.join (settings.MEDIA_ROOT, wetlab_config.RUN_SAMPLE_SHEET_DIRECTORY, new_sample_sheet_name)
    bs_file = run_process.get_sample_file()
    # create projects in DDBB
    for project in projects :
        create_project(project, experiment_name, bs_file)
    recorded_results.append(('Successful projects creation in database', 'OK'))
    
    return recorded_results, 'OK'
    

def  run_nextseq_test_rec_to_sample_sent (run_folder, experiment_name) :
    recorded_results = []
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function run_nextseq_test_rec_to_sample_sent')
    
    conn = open_samba_connection()
    # check if test folder exists
    
    run_data_root_folder = os.path.join('/', wetlab_config.SAMBA_APPLICATION_FOLDER_NAME , run_folder)
    try:   
        run_folder_list = conn.listPath( wetlab_config.SAMBA_SHARED_FOLDER_NAME, run_data_root_folder)
        recorded_results.append(('Check if test Run folder exists ', 'OK'))
        logger.info('Test Run folder exists on remote server')
    except Exception as ex :
        string_message = 'Test run folder does not exists '
        logging_errors(string_message, False, True)
        recorded_results.append((str(ex) , 'NOK'))
        logger.debug ('End function run_nextseq_test_rec_to_sample_sent with error')
        return recorded_results, 'NOK'
    
    # Fetch RunParameter.xml file from remote server 
    l_run_parameter = os.path.join(wetlab_config.RUN_TEMP_DIRECTORY, wetlab_config.RUN_PARAMETER_NEXTSEQ)
    s_run_parameter = os.path.join(wetlab_config.SAMBA_APPLICATION_FOLDER_NAME, run_folder, wetlab_config.RUN_PARAMETER_NEXTSEQ)
    try:
       l_run_parameter = fetch_remote_file (conn, run_folder, s_run_parameter, l_run_parameter)
       logger.info('Sucessfully fetch of RunParameter file')
       
       recorded_results.append(('Successful RunParameter.xml file ', 'OK'))
    except:
        string_message = 'Unable to fetch the RunParameter.xml file '
        logging_errors(string_message, False, True)
        recorded_results.append(('Successful RunParameter.xml file ', 'NOK'))
        logger.debug ('End function create_run_test_nextseq_in_recorded with error')
        return recorded_results, 'NOK'
    # Test MiSeq run on recorded state
    try:  
        update_nextseq_in_recorded_state = handle_nextseq_recorded_run (conn, run_folder, l_run_parameter, experiment_name)
    
    except ValueError as e :
        recorded_results.append((str(e) , 'NOK'))
        logger.debug ('End function create_run_test_nextseq_in_recorded with error')
        return recorded_results, 'NOK'
    except Exception as ex :
        recorded_results.append((str(ex) , 'NOK'))
        logger.debug ('End function create_run_test_nextseq_in_recorded with error')
        return recorded_results, 'NOK'

    if   experiment_name == update_nextseq_in_recorded_state.get_run_name():
        if 'Sample Sent' == RunProcess.objects.get(runName__exact = experiment_name).get_state():
            recorded_results.append(('NextSeq Run is in Sample Sent state', 'OK'))
            logger.debug ('End function create_run_test_nextseq_in_recorded')
            return recorded_results, 'OK'
        else :
            recorded_results.append(('NextSeq Run is wrong state', 'NOK'))
            logger.debug ('End function create_run_test_nextseq_in_recorded with error')
            return recorded_results, 'NOK'
    else:
        recorded_results.append(('NextSeq function returns an invalid run ', 'NOK'))
        logger.debug ('End function create_run_test_nextseq_in_recorded with error')
        return recorded_results, 'NOK'
    
    return recorded_results, 'OK'



def  run_nextseq_test_sample_sent_to_Processing_Run (experiment_name) :
    sample_sent_results = []
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function run_nextseq_test_sample_sent_to_Processing_Run')
    run_name = RunProcess.objects.get(runName__exact = experiment_name)
    conn = open_samba_connection()
    
    try:
        update_in_sample_sent = manage_nextseq_in_samplesent (conn, run_name)
    
    except Exception as ex :
        sample_sent_results.append(('MiSeq Run is in Processing Run state', 'NOK'))
        logger.debug ('End function run_nextseq_test_sample_sent_to_Processing_Run with error')
        return sample_sent_results, 'NOK'
    
    import pdb; pdb.set_trace()
    if update_in_sample_sent == experiment_name:
        if 'Processing Run' == RunProcess.objects.get(runName__exact = experiment_name).get_state():
            sample_sent_results.append(('Check log files from the sequencer', 'OK'))
            sample_sent_results.append(('MiSeq Run is in Processing Run state', 'OK'))
            logger.debug ('End function run_nextseq_test_sample_sent_to_Processing_Run')
            return sample_sent_results, 'OK'
        else :
            sample_sent_results.append(('MiSeq Run is wrong state', 'NOK'))
            logger.debug ('End function run_nextseq_test_sample_sent_to_Processing_Run with error')
            return sample_sent_results, 'NOK'
    else:
        sample_sent_results.append(('MiSeq function returns an invalid run', 'NOK'))
        logger.debug ('End function run_nextseq_test_sample_sent_to_Processing_Run with error')
        return sample_sent_results, 'NOK'


def run_nextseq_test_Processing_Run_to_Processed_Run (experiment_name):
    processing_run_results = []
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function run_nextseq_test_Processing_Run_to_Processed_Run')
    run_name = RunProcess.objects.get(runName__exact = experiment_name)
    conn = open_samba_connection()
    
    try:
        update_in_procesing_run = manage_nextseq_in_processing_run (conn, run_name)
    
    except Exception as ex :
        processing_run_results.append(('NextSeq Run is in Processing Run state', 'NOK'))
        logger.debug ('End function run_nextseq_test_Processing_Run_to_Processed_Run with error')
        return processing_run_results, 'NOK'
    import pdb; pdb.set_trace()
    if update_in_procesing_run == experiment_name:
        if 'Processed Run' == RunProcess.objects.get(runName__exact = experiment_name).get_state():
            processing_run_results.append(('Successful completion on Sequencer for MiSeq Run', 'OK'))
            processing_run_results.append(('NextSeq Run is in Processed Run state', 'OK'))
            logger.debug ('End function run_nextseq_test_Processing_Run_to_Processed_Run')
            return processing_run_results, 'OK'
        else :
            processing_run_results.append(('NextSeq Run is wrong state', 'NOK'))
            logger.debug ('End function run_nextseq_test_Processing_Run_to_Processed_Run with error')
            return processing_run_results, 'NOK'
    else:
        processing_run_results.append(('NextSeq function returns an invalid run', 'NOK'))
        logger.debug ('End function run_nextseq_test_Processing_Run_to_Processed_Run with error')
        return processing_run_results, 'NOK'




















def run_miseq_test_rec_to_sample_sent (new_run, experiment_name) :
    recorded_results = []
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function run_miseq_test_rec_to_sample_sent')
    working_path = settings.MEDIA_ROOT
    os.chdir(working_path)
    #import pdb; pdb.set_trace()
    
    try:        
        conn = open_samba_connection()
        recorded_results.append(('Connection to remote server ', 'OK'))
        logger.info('Success connection to remote server')
    except Exception as ex :
        string_message = 'Unable to connect to remote server '
        logging_errors(string_message, False, True)
        recorded_results.append((str(ex) , 'NOK'))
        return recorded_results, 'NOK'
    
    # check if test folder exists
    run_data_root_folder = os.path.join('/', wetlab_config.SAMBA_APPLICATION_FOLDER_NAME , new_run)
    try:   
        run_folder_list = conn.listPath( wetlab_config.SAMBA_SHARED_FOLDER_NAME, run_data_root_folder)
        recorded_results.append(('Check if test Run folder exists ', 'OK'))
        logger.info('Test Run folder exists on remote server')
    except Exception as ex :
        string_message = 'Test run folder does not exists '
        logging_errors(string_message, False, True)
        recorded_results.append((str(ex) , 'NOK'))
        return recorded_results, 'NOK'
    
    # Check user is defined in database
    if not User.objects.filter(username__exact = 'test_user1').exists():
        user = User.objects.create_user(username='test_user1',
                                 email='test_user1@iSkyLIMS.com',
                                 password='test_user1')

    # Fetch RunParameter.xml file from remote server 
    l_run_parameter = os.path.join(wetlab_config.RUN_TEMP_DIRECTORY, wetlab_config.RUN_PARAMETER_NEXTSEQ)
    s_run_parameter = os.path.join(wetlab_config.SAMBA_APPLICATION_FOLDER_NAME, new_run,wetlab_config.RUN_PARAMETER_NEXTSEQ)
    try:
       l_run_parameter = fetch_remote_file (conn, new_run, s_run_parameter, l_run_parameter)
       logger.info('Sucessfully fetch of RunParameter file')
       
       recorded_results.append(('Successful RunParameter.xml file ', 'OK'))
    except:
        string_message = 'Unable to fetch the RunParameter.xml file '
        logging_errors(string_message, False, True)
        recorded_results.append(('Successful RunParameter.xml file ', 'NOK'))
        return recorded_results, 'NOK'
    # Test MiSeq run on recorded state
    try:
        update_miseq_in_recorded_state =  handle_miseq_run (conn, new_run, l_run_parameter, experiment_name)
        
        recorded_results.append(('Valid Sample Sheet ', 'OK'))
        recorded_results.append(('Stored Sample Sheet ', 'OK'))
        recorded_results.append(('Updated library name ', 'OK'))
        recorded_results.append(('Stored Projects ', 'OK'))
    except ValueError as e :
        recorded_results.append((str(e) , 'NOK'))
        logger.debug ('End function run_miseq_test_rec_to_sample_sent with error')
        return recorded_results, 'NOK'
    except Exception as ex :
        recorded_results.append((str(ex) , 'NOK'))
        logger.debug ('End function run_miseq_test_rec_to_sample_sent with error')
        return recorded_results, 'NOK'
    
    if update_miseq_in_recorded_state == new_run :
        if 'Sample Sent' == RunProcess.objects.get(runName__exact = experiment_name).get_state():
            recorded_results.append(('MiSeq Run is in Sample Sent state', 'OK'))
            logger.debug ('End function run_miseq_test_rec_to_sample_sent')
            return recorded_results, 'OK'
        else :
            recorded_results.append(('MiSeq Run is wrong state', 'NOK'))
            logger.debug ('End function run_miseq_test_rec_to_sample_sent with error')
            return recorded_results, 'NOK'
    else:
        recorded_results.append(('MiSeq function returns an invalid run ', 'NOK'))
        logger.debug ('End function run_miseq_test_rec_to_sample_sent with error')
        return recorded_results, 'NOK'


def run_miseq_test_sample_sent_to_Processing_Run (experiment_name):
    sample_sent_results = []
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function run_miseq_test_sample_sent_to_Processing_Run')
    run_name = RunProcess.objects.get(runName__exact = experiment_name)
    conn = open_samba_connection()
    
    try:
        update_in_sample_sent = manage_miseq_in_samplesent (conn, run_name)
    
    except Exception as ex :
        sample_sent_results.append(('MiSeq Run is in Processing Run state', 'NOK'))
        logger.debug ('End function run_miseq_test_sample_sent_to_Processing_Run with error')
        return sample_sent_results, 'NOK'

    if update_in_sample_sent == experiment_name:
        if 'Processing Run' == RunProcess.objects.get(runName__exact = experiment_name).get_state():
            sample_sent_results.append(('Check log files from the sequencer', 'OK'))
            sample_sent_results.append(('MiSeq Run is in Processing Run state', 'OK'))
            logger.debug ('End function run_miseq_test_sample_sent_to_Processing_Run')
            return sample_sent_results, 'OK'
        else :
            sample_sent_results.append(('MiSeq Run is wrong state', 'NOK'))
            logger.debug ('End function run_miseq_test_sample_sent_to_Processing_Run with error')
            return sample_sent_results, 'NOK'
    else:
        sample_sent_results.append(('MiSeq function returns an invalid run', 'NOK'))
        logger.debug ('End function run_miseq_test_sample_sent_to_Processing_Run with error')
        return sample_sent_results, 'NOK'

def run_miseq_test_Processing_Run_to_Processed_Run (experiment_name):
    processing_run_results = []
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function run_miseq_test_Processing_Run_to_Processed_Run')
    run_name = RunProcess.objects.get(runName__exact = experiment_name)
    conn = open_samba_connection()
    
    try:
        update_in_procesing_run = manage_miseq_in_processing_run (conn, run_name)
    
    except Exception as ex :
        processing_run_results.append(('MiSeq Run is in Processing Run state', 'NOK'))
        logger.debug ('End function run_miseq_test_Processing_Run_to_Processed_Run with error')
        return processing_run_results, 'NOK'

    if update_in_procesing_run == experiment_name:
        if 'Processed Run' == RunProcess.objects.get(runName__exact = experiment_name).get_state():
            processing_run_results.append(('Successful completion on Sequencer for MiSeq Run', 'OK'))
            processing_run_results.append(('MiSeq Run is in Processed Run state', 'OK'))
            logger.debug ('End function run_miseq_test_Processing_Run_to_Processed_Run')
            return processing_run_results, 'OK'
        else :
            processing_run_results.append(('MiSeq Run is wrong state', 'NOK'))
            logger.debug ('End function run_miseq_test_Processing_Run_to_Processed_Run with error')
            return processing_run_results, 'NOK'
    else:
        processing_run_results.append(('MiSeq function returns an invalid run', 'NOK'))
        logger.debug ('End function run_miseq_test_Processing_Run_to_Processed_Run with error')
        return processing_run_results, 'NOK'

def run_test_Processed_Run_to_Processing_Bcl2fastq (experiment_name):
    processed_run_results = []
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function run_test_Processed_Run_to_Processing_Bcl2fastq')
    run_name = RunProcess.objects.get(runName__exact = experiment_name)
    conn = open_samba_connection()
    
    try:
        update_in_processing_run = manage_run_in_processed_run (conn, run_name)
    
    except Exception as ex :
        processed_run_results.append(('Run in Processed Run state', 'NOK'))
        logger.debug ('End function run_test_Processed_Run_to_Processing_Bcl2fastq with error')
        return processed_run_results, 'NOK'

    if update_in_processing_run == experiment_name:
        if 'Processing Bcl2fastq' == RunProcess.objects.get(runName__exact = experiment_name).get_state():
            processed_run_results.append(('Successful completion on Sequencer for MiSeq Run', 'OK'))
            processed_run_results.append(('Run is in Processing Bcl2fastq state', 'OK'))
            logger.debug ('End function run_test_Processed_Run_to_Processing_Bcl2fastq')
            return processed_run_results, 'OK'
        else :
            processed_run_results.append(('Run is wrong state', 'NOK'))
            logger.debug ('End function run_test_Processed_Run_to_Processing_Bcl2fastq with error')
            return processed_run_results, 'NOK'
    else:
        processed_run_results.append(('Function returns an invalid run', 'NOK'))
        logger.debug ('End function run_test_Processed_Run_to_Processing_Bcl2fastq with error')
        return processed_run_results, 'NOK'

def run_test_Processing_Bcl2fastq_to_Processed_Bcl2fastq (experiment_name):
    processing_bcl2fastq_results = []
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function run_test_Processing_Bcl2fastq_to_Processed_Bcl2fastq')
    run_name = RunProcess.objects.get(runName__exact = experiment_name)
    conn = open_samba_connection()
    
    try:
        update_in_procesing_bclfastq = manage_run_in_processing_bcl2fastq (conn, run_name)
    
    except Exception as ex :
        processing_bcl2fastq_results.append(('Run is in Processing Run state', 'NOK'))
        logger.debug ('End function run_test_Processing_Bcl2fastq_to_Processed_Bcl2fastq with error')
        return processing_bcl2fastq_results, 'NOK'

    if update_in_procesing_bclfastq == experiment_name:
        if 'Processed Bcl2fastq' == RunProcess.objects.get(runName__exact = experiment_name).get_state():
            processing_bcl2fastq_results.append(('Successful completion of Bcl2fastq conversion', 'OK'))
            processing_bcl2fastq_results.append(('Bcl2fastq conversion date updated', 'OK'))
            processing_bcl2fastq_results.append(('Run is in Processed Bcl2fastq state', 'OK'))
            logger.debug ('End function run_test_Processing_Bcl2fastq_to_Processed_Bcl2fastq')
            return processing_bcl2fastq_results, 'OK'
        else :
            processing_bcl2fastq_results.append(('Run is wrong state', 'NOK'))
            logger.debug ('End function run_test_Processing_Bcl2fastq_to_Processed_Bcl2fastq with error')
            return processing_bcl2fastq_results, 'NOK'
    else:
        processing_bcl2fastq_results.append(('Function returns an invalid run', 'NOK'))
        logger.debug ('End function run_test_Processing_Bcl2fastq_to_Processed_Bcl2fastq with error')
        return processing_bcl2fastq_results, 'NOK'


def run_test_Processed_Bcl2fastq_to_Completed (experiment_name):
    processed_bcl2fastq_results = []
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function run_test_Processed_Bcl2fastq_to_Completed')
    run_name = RunProcess.objects.get(runName__exact = experiment_name)
    conn = open_samba_connection()
    
    try:
        update_in_processed_bclfastq = manage_run_in_processed_bcl2fastq (conn, run_name)
    
    except Exception as ex :
        processed_bcl2fastq_results.append(('Run in Processed Bcl2fastq state', 'NOK'))
        logger.debug ('End function run_test_Processed_Bcl2fastq_to_Completed with error')
        return processed_bcl2fastq_results, 'NOK'

    if update_in_processed_bclfastq == experiment_name:
        if 'Completed' == RunProcess.objects.get(runName__exact = experiment_name).get_state():
            processed_bcl2fastq_results.append(('Successful demultiplexing files', 'OK'))
            processed_bcl2fastq_results.append(('Updated database with demultiplexion data', 'OK'))
            processed_bcl2fastq_results.append(('Projects in the run successfuly updated', 'OK'))
            processed_bcl2fastq_results.append(('Stored graphics ', 'OK'))
            processed_bcl2fastq_results.append(('Run is in Processing Completed state', 'OK'))
            logger.debug ('End function run_test_Processed_Bcl2fastq_to_Completed')
            return processed_bcl2fastq_results, 'OK'
        else :
            processed_bcl2fastq_results.append(('Run is wrong state', 'NOK'))
            logger.debug ('End function run_test_Processed_Bcl2fastq_to_Completed with error')
            return processed_bcl2fastq_results, 'NOK'
    else:
        processed_bcl2fastq_results.append(('Function returns an invalid run', 'NOK'))
        logger.debug ('End function run_test_Processed_Bcl2fastq_to_Completed with error')
        return processed_bcl2fastq_results, 'NOK'



