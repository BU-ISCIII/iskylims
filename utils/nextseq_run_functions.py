import xml.etree.ElementTree as ET
from datetime import datetime
import os
import shutil
from .generic_functions import *
from .handling_crontab_common_functions import *
from iSkyLIMS_wetlab.models import *
#from iSkyLIMS_drylab.models import Machines, Platform
from iSkyLIMS_wetlab.wetlab_config import *

from django_utils.models import Center

def check_completion_success (l_run_completion, experiment_name):
    '''
    Description:
        The function will check if the run was success
    Input:
        l_run_completion  # local path for the run completion file
        experiment_name   # Contains the experiment name
    Functions:
        find_xml_tag_text # located at utils.generic_functions
    Constant:
        COMPLETION_SUCCESS
    Variables
        status_run # contain the text string of tag completion
    Return
        True if successfuly completed
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting function for check_completion_success', experiment_name)
    # check if NextSEq run have been successful completed
    status_run = find_xml_tag_text (l_run_completion, wetlab_config.COMPLETION_TAG )
    if  status_run != wetlab_config.COMPLETION_SUCCESS:
        string_message = 'Run status was ' + status_run
        logging_errors (string_message, False, False)
        return False
    else:
        logger.info ('%s : Run successfuly completed ', experiment_name)
    logger.debug ('%s : End function for check_completion_success', experiment_name)
    return True




def manage_nextseq_in_samplesent(conn, run_object_name) :
    '''
    Description:
        The function will look for the run completion file to check if run is
        completed.
        If file is not found, the sequencer is  still processing the run
    Input:
        run_in_sample_sent  #  object name of the run
        experiment_name     # name of the run
        run_folder          # run folder name on the remote server
    Functions:
        check_completion_success # located as this file
        copy_to_remote_file   # located at utils.generic_functions
        handling_errors_in_run # located at utils.generic_functions
        nextseq_parsing_run_information # located as this file
        get_samba_application_shared_folder     # located at utils.handling_crontab_common_functions.py
    Variable:
        l_run_completion    # completion file name in the local temporary folder
        run_object_name     # RunProcess object
        s_run_completion    # completion file name at the remote server
        updated_library # return value after updating the library
    Return
        Exception if file does not exist to process next run.
        Experiment name if run was completed
    '''
    logger = logging.getLogger(__name__)
    experiment_name = run_object_name.get_run_name()
    logger.debug ('% : Starting function manage_nextseq_in_samplesent', experiment_name)

    logger.info('%s : Start handling in Sample Sent state ', experiment_name)
    run_folder = RunningParameters.objects.get(runName_id = run_object_name).get_run_folder()
    l_run_completion = os.path.join(wetlab_config.RUN_TEMP_DIRECTORY, wetlab_config.RUN_COMPLETION)
    s_run_completion = os.path.join(get_samba_application_shared_folder() , run_folder, wetlab_config.RUN_COMPLETION)

    try:
        l_run_completion = fetch_remote_file (conn, run_folder, s_run_completion, l_run_completion)
        logger.info('%s : Sucessfully fetch of Completion status file',experiment_name)
        # Get status value for run completion
        completion_status = check_completion_success (l_run_completion, experiment_name)

    except Exception as e:
        logger.info ('%s : Completion status file is not present on the run folder', experiment_name)
        logger.info ('%s : Moving run to Processing run state', experiment_name)
        run_updated = run_object_name.set_run_state('Processing Run')
        logger.debug ('%s : End function for handling NextSeq run with exception', experiment_name)
        raise

    logger.info ('%s : Deleting local copy of completion status', experiment_name)
    # cleaning up the completion  in local temporary file
    os.remove(l_run_completion)
    if not completion_status :
        string_message = experiment_name + ' : Run status was ' + status_run
        logging_errors (string_message, False, False)
        # Set run to cancelled  state
        run_object_name.set_run_state('Cancelled')
        logger.debug ('%s : End function for handling NextSeq run Cancelled',experiment_name)
        raise ValueError ('Run was CANCELLED')
    else:
        run_updated = run_object_name.set_run_state('Processing Run')
        logger.info('%s : is now on Processing Run state', experiment_name)
        return experiment_name


def manage_nextseq_in_processing_run(conn, run_object_name) :
    '''
    Description:
        The function will look for the run completion file to check if run is
        completed.
        If file is not found, the run will keep in processing the run.
        If after the waiting time defined on MAXIMUM_TIME_WAIT_FOR_RUN_COMPLETION
        the runn will move to error state
    Input:
        run_in_sample_sent  #  object name of the run
        experiment_name     # name of the run
        run_folder          # run folder name on the remote server
    Functions:
        check_completion_success    # located in this file
        get_attributes_remote_file  # Located at utils.generic_functions
    Variable:
        l_run_completion    # completion file name in the local temporary folder
        run_object_name     # RunProcess object
        s_run_completion    # completion file name at the remote server
        updated_library # return value after updating the library
    Return
        Exception if file does not exist to process next run.
        Experiment name if run was completed
    '''
    logger = logging.getLogger(__name__)
    experiment_name = run_object_name.get_run_name()
    logger.debug ('%s : Starting function manage_nextseq_in_processing_run', experiment_name)

    logger.info('%s :  Handling in Processing Run state ', experiment_name)
    run_folder = RunningParameters.objects.get(runName_id = run_object_name).get_run_folder()
    l_run_completion = os.path.join(wetlab_config.RUN_TEMP_DIRECTORY, wetlab_config.RUN_COMPLETION)
    s_run_completion = os.path.join(wetlab_config.SAMBA_APPLICATION_FOLDER_NAME , run_folder, wetlab_config.RUN_COMPLETION)

    try:
        l_run_completion = fetch_remote_file (conn, run_folder, s_run_completion, l_run_completion)
        logger.info('%s : Sucessfully fetch of Completion status file', experiment_name)
        # Get status value for run completion
        completion_status = check_completion_success (l_run_completion, experiment_name)
        # Get the date and time when the RunCompletionStatus is created
        run_completion_attributes = get_attributes_remote_file (conn, run_folder, s_run_completion)
        run_completion_date = datetime.datetime.fromtimestamp(int(run_completion_attributes.create_time)).strftime('%Y-%m-%d %H:%M:%S')
        run_completion_date = run_object_name.set_run_completion_date(run_completion_date)
        logger.info('%s : Updated with completion date', experiment_name)

    except Exception as e:
        logger.info ('%s : Completion status file is still not present waiting more time', experiment_name)
        logger.debug ('%s : End function for handling manage_nextseq_in_processing_run', experiment_name)
        raise Exception ( 'Completion file not present')
        #return
    logger.info ('%s : Deleting local copy of completion status' , experiment_name)
    # cleaning up the completion  in local temporary file
    os.remove(l_run_completion)

    if not completion_status :
        string_message = experiment_name +  ' : Run status was ' + completion_status
        logging_errors (string_message, False, False)
        # Set tun to Cancelled state
        logger.info('Set run %s  to cancelled state', experiment_name)
        logger.debug ('End function for manage_nextseq_in_processing_run cancelled')
        raise ValueError ('Run was CANCELLED')
    else:
        run_updated = run_object_name.set_run_state('Processed Run')
        logger.info('%s : is now on Processed Run state', experiment_name)
        logger.debug ('%s : End function for handling manage_nextseq_in_processing_run', experiment_name)
        return experiment_name



def handle_nextseq_recorded_run (conn, new_run, l_run_parameter, l_run_info, experiment_name):
    '''
    Description:
        The function will copy the sample sheet to the remote run folder
    Input:
        conn        # samba connection object
        new_run     # folder remote directory for miseq run
        l_run_parameter  # local path for the run parameter file
        l_run_info      # local path for run info file
        experiment_name  # name used on miseq run
    Functions:
        copy_to_remote_file   # located at utils.generic_functions
        fetch_remote_file               # loacted at handling_crontab_common_functions.py
        handling_errors_in_run          # located at utils.generic_functions
        nextseq_parsing_run_information # located as this file
        create_new_sequencer_lab_not_defined    # located at utils/generic_functions.py

    Constant:
        RUN_TEMP_DIRECTORY
        RUN_INFO
        SAMPLE_SHEET
    Return:
        number_of_cycles
        file_content
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting function handle_nextseq_recorded_run',experiment_name)
    if RunProcess.objects.filter(runName__exact = experiment_name , state__runStateName ='Recorded').exists():
        logger.info('%s : Processing on Recorded state', experiment_name)
        run_process_obj = RunProcess.objects.filter(runName__exact = experiment_name ).last()

        running_parameters = nextseq_parsing_run_info_and_parameter_information(l_run_info, l_run_parameter, experiment_name)
        logger.info('%s  : Deleting runParameter file', experiment_name)
        os.remove(l_run_parameter)
        logger.info('%s  : Deleting runInfo file', experiment_name)
        os.remove(l_run_info)

        if RunningParameters.objects.filter(runName_id = run_process_obj).exists():
            run_parameter_objs = RunningParameters.objects.filter(runName_id = run_process_obj)
            for run_parameter_obj in run_parameter_objs:
                logger.info('%s  : Deleting RunParameters object from database', experiment_name)
                run_parameter_obj.delete()
        run_parameters = RunningParameters.objects.create_running_parameters(running_parameters['running_data'], run_process_obj)
        logger.info( '%s  : Created RunParameters object on database', experiment_name)

        if SequencerInLab.objects.filter(sequencerName__exact = running_parameters['instrument']).exists():
            sequencer_obj = SequencerInLab.objects.filter(sequencerName__exact = running_parameters['instrument']).last()

        else:
            string_message = experiment_name + ' : ' + instrument + ' no sequencer defined '
            logging_errors(string_message, False, True)
            sequencer_obj = create_new_sequencer_lab_not_defined ( running_parameters['instrument'],  running_parameters['NumLanes'], experiment_name)
            logger.info('%s : Continue the proccess after creating the new sequencer' ,experiment_name)
        run_process_obj.set_used_sequencer(sequencer_obj)
        logger.info('%s : Sequencer  stored on database', experiment_name)

        sample_sheet_tmp_dir = os.path.join(wetlab_config.RUN_TEMP_DIRECTORY_RECORDED,run_process_obj.get_run_id())
        l_sample = os.path.join(sample_sheet_tmp_dir, wetlab_config.SAMPLE_SHEET)

        if os.path.exists(l_sample):
            if wetlab_config.COPY_SAMPLE_SHEET_TO_REMOTE :
                # copy Sample heet file to remote directory
                logger.info('%s : Copy sample sheet to remote folder %s', experiment_name, new_run)
                s_sample= os.path.join(new_run, wetlab_config.SAMPLE_SHEET)
                try:
                    sample_sheet_copied = copy_to_remote_file  (conn, run_folder, s_sample, l_sample)
                    logger.info('%s : Sucessfully copy Sample sheet to remote folder', experiment_name)
                except Exception as e:
                    string_message = experiment_name + ': Unable to copy the Sample Sheet to remote folder ' + new_run
                    logging_errors (string_message, True, False)
                    handling_errors_in_run (experiment_name, '23')
                    logger.debug ('%s : End function for handling NextSeq run with exception', experiment_name)
                    raise Exception
            # Update run to Sample Sent

            run_process_obj.set_run_state('Sample Sent')
            # Cleaning up the temporary folder with the sample sheet
            try:
                shutil.rmtree(sample_sheet_tmp_dir)
                logger.debug('%s : Deleted temporary folder containing the samplesheet', experiment_name)
            except Exception as e:
                string_message = experiment_name +' : Unable to delete temporary folder with the sample sheet ' + sample_sheet_tmp_dir
                logging_errors (string_message, True, True)
                logger.debug ('%s : End function for handling NextSeq run with exception', experiment_name)
                raise Exception

        else:
            string_message = experiment_name + ' : sample sheet not found on local Directory ' + sample_sheet_tmp_dir
            logging_errors (string_message, False, True)
            handling_errors_in_run (experiment_name, '22')
            logger.debug ('%s : End function for handling NextSeq run with exception', experiment_name)
            raise Exception

        logger.debug ('%s : End function for handling NextSeq run ', experiment_name)
        return run_process_obj
    else:
        if not RunProcess.objects.filter(runName__exact = experiment_name).exists():
            string_message =  experiment_name + ' : is not yet defined on database'
            logging_errors (string_message, False, False)
            logger.debug ('%s : End function handle_nextseq_recorded_run with exception', experiment_name)
            raise ValueError ('Run not defined yet')
        else:
            run_state = RunProcess.objects.filter(runName__exact = experiment_name).get_state()
            string_message = 'Run ' + experiment_name +' is in state ' + run_state + '.  Should be in Recorded'
            logging_errors (string_message, False, False)
            logger.debug ('%s : End function fhandle_nextseq_recorded_run with exception', experiment_name)
            raise ValueError ('Run on wrong state')
    return
