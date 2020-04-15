import xml.etree.ElementTree as ET
from datetime import datetime
import os
import shutil
from .generic_functions import *
from iSkyLIMS_wetlab.models import *
from iSkyLIMS_drylab.models import Machines, Platform
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

def nextseq_parsing_run_information(l_run_info, l_run_parameter, experiment_name) :
    '''
    Description:
        The function is called for parsing the RunInfo and RunParameter
        files.
        After parsing the RunningParameters database table will be
        updated with a new row having the parsed data
        Empty values will be set for MiSeq runs that exist on NextSeq
        but not in MiSeq runs
    Input:
        run_info    # contains the path for RunInfo.xml file
        run_parameter # contains the path for RunParameter.xml file
        experiment_name      # contains the experiment name
    CONSTANTS:
        FIELDS_TO_COLLECT_FROM_RUN_INFO_FILE
    Variables:
        image_channel   # list containing the image channel values
        logger          # contains the logger object to write information
                        on log file
        p_parameter     # element tree object for run parameter
        p_run           # element tree for run info
        projects_to_update # list of projects that belong to the same run
        running_data    # dictionary with  the  information parsed
        running_parameters # new RunningParameter object to store parsed data
        sequencer       # sequencer index of the machine
     Return:
        running_data
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting function nextseq_parsing_run_information', experiment_name)
    running_data={}
    image_channel=[]

    #################################################
    ## parsing RunInfo.xml file
    #################################################
    run_data=ET.parse(l_run_info)
    run_root=run_data.getroot()
    logger.info('%s : parsing the runInfo.xml file ', experiment_name)
    p_run=run_root[0]
    ## getting the common values NextSeq and MiSeq


    running_data['Flowcell']=p_run.find('Flowcell').text
    running_data['FlowcellLayout']=p_run.find('FlowcellLayout').attrib
    #################################################
    ## parsing RunParameter.xml file
    #################################################
    logger.info('%s : Parsing the runParameter.xml file  ',experiment_name )
    parameter_data=ET.parse(l_run_parameter)
    parameter_data_root=parameter_data.getroot()
    p_parameter=parameter_data_root[1]
    ## getting the common values NextSeq and MiSeq
    for field in FIELDS_TO_COLLECT_FROM_RUN_INFO_FILE:
        try:
            running_data[field] = parameter_data_root.find(field).text
        except:
            running_data[field] = ''
            string_message = experiment_name + ' : Parameter ' + field  + ' not found in RunParameter.xml'
            logging_warnings(string_message, False)


    '''
    running_data['RunID']=parameter_data_root.find('RunID').text
    running_data['ExperimentName']=parameter_data_root.find('ExperimentName').text
    running_data['RTAVersion']=parameter_data_root.find('RTAVersion').text
    running_data['Chemistry']=parameter_data_root.find('Chemistry').text
    running_data['RunStartDate']=parameter_data_root.find('RunStartDate').text
    ###
    try:
        running_data['RunManagementType']=parameter_data_root.find('RunManagementType').text
    except:
        running_data['RunManagementType']=''
    running_data['ApplicationVersion']=p_parameter.find('ApplicationVersion').text
    running_data['NumTilesPerSwath']=p_parameter.find('NumTilesPerSwath').text
    running_data['SystemSuiteVersion']=parameter_data_root.find('SystemSuiteVersion').text
    ##
    running_data['LibraryID']=parameter_data_root.find('LibraryID').text
    running_data['AnalysisWorkflowType']=parameter_data_root.find('AnalysisWorkflowType').text
    running_data['PlannedRead1Cycles']=parameter_data_root.find('PlannedRead1Cycles').text
    running_data['PlannedRead2Cycles']=parameter_data_root.find('PlannedRead2Cycles').text
    running_data['PlannedIndex1ReadCycles']=parameter_data_root.find('PlannedIndex1ReadCycles').text
    running_data['PlannedIndex2ReadCycles']=parameter_data_root.find('PlannedIndex2ReadCycles').text
    '''

    for i in run_root.iter('Name'):
        image_channel.append(i.text)
    running_data['ImageChannel']=image_channel
    running_data['ImageDimensions']=p_run.find('ImageDimensions').attrib

    ## get the instrument for NextSeq run
    instrument = parameter_data_root.find('InstrumentID').text

    ##############################################
    ## updating the date fetched from the Date tag for run and project
    ##############################################
    date = p_run.find('Date').text
    logger.debug('%s : Found date that was recorded the Run %s', experiment_name , date)
    run_date = datetime.datetime.strptime(date, '%y%m%d')

    logger.debug('%s : End function nextseq_parsing_run_information', experiment_name)
    return running_data, run_date, instrument


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
    s_run_completion = os.path.join(wetlab_config.SAMBA_APPLICATION_FOLDER_NAME , run_folder, wetlab_config.RUN_COMPLETION)

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



def handle_nextseq_recorded_run (conn, new_run, l_run_parameter, experiment_name):
    '''
    Description:
        The function will copy the sample sheet to the remote run folder
    Input:
        conn        # samba connection object
        new_run     # folder remote directory for miseq run
        l_run_parameter  # local path for the run parameter file
        experiment_name  # name used on miseq run
    Functions:
        copy_to_remote_file   # located at utils.generic_functions
        handling_errors_in_run # located at utils.generic_functions
        nextseq_parsing_run_information # located as this file
    Import:
        os
        shutil
        RunningParameters
        RunProcess
    Constant:
        RUN_TEMP_DIRECTORY
        RUN_INFO
        SAMPLE_SHEET
    Variables:
        instrument  # name of the sequencer used in the run
        l_run_info  # local temporary copy of RunInfo
        l_sample_sheet  # local temporary copy of sample sheet
        new_run_parameters # new RunningParameters object created for this run
        run_process  # RunProcess object related to this run

        run_date        # date of starting run
        run_folder      # run foloder on remote server
        running_parameters  # dictionary with parsing information from
                            RunParameter and RunInfo
        s_run_info  # path of RunInfo on the remote server
        s_sample_sheet  # path of SampleSheet on the remote server
    Functions:
        create_new_sequencer_lab_not_defined    # located at utils/generic_functions.py
    Return:
        number_of_cycles
        file_content
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting function handle_nextseq_recorded_run',experiment_name)
    if RunProcess.objects.filter(runName__exact = experiment_name , state__runStateName ='Recorded').exists():
        logger.info('%s : Processing on Recorded state', experiment_name)
        run_process = RunProcess.objects.get(runName__exact = experiment_name )
        run_folder = os.path.join(wetlab_config.SAMBA_APPLICATION_FOLDER_NAME , new_run)
        # Fetch run info
        l_run_info = os.path.join(wetlab_config.RUN_TEMP_DIRECTORY, wetlab_config.RUN_INFO)
        s_run_info = os.path.join(run_folder, wetlab_config.RUN_INFO)
        try:
            l_run_info = fetch_remote_file (conn, new_run, s_run_info, l_run_info)
            logger.info('%s : Sucessfully fetch of RunInfo file', experiment_name)
        except Exception as e:
            string_message = experiment_name + ' : Unable to fetch the RunInfo file'
            logging_errors(string_message, True, True)
            handling_errors_in_run(experiment_name, '20')

            # cleaning up the RunParameter in local temporaty file
            os.remove(l_run_parameter)
            logger.debug ('%s : End function for handling NextSeq run with exception', experiment_name)
            raise ValueError ('Unable to fetch RunInfo')  # returning to handle next run folder

        # Parsing RunParameter and Run Info

        running_parameters, run_date, instrument = nextseq_parsing_run_information(l_run_info, l_run_parameter, experiment_name)

        run_process.set_run_date(run_date)

        if  Projects.objects.filter(runprocess_id = run_process).exists():
            project_objs = Projects.objects.filter(runprocess_id = run_process)
            for project_obj in project_objs :
                project_obj.set_project_run_date(run_date)

        if not RunningParameters.objects.filter(runName_id = run_process).exists():
            run_parameters = RunningParameters.objects.create_running_parameters(running_parameters, run_process)
        else:
            string_message = experiment_name + " : RunParameters already in database for " + run_process.get_run_name()
            logging_warnings(string_message,False)

        if SequencerInLab.objects.filter(sequencerName__exact = instrument).exists():
            sequencer = SequencerInLab.objects.filter(sequencerName__exact = instrument)

        else:
            string_message = experiment_name + ' : ' + instrument + ' no sequencer defined '
            logging_errors(string_message, False, True)
            sequencer = create_new_sequencer_lab_not_defined (instrument,l_run_info, experiment_name)
            logger.info('%s : Continue the proccess after creating the new sequencer' ,experiment_name)
        run_process.set_used_sequencer(sequencer)
        logger.info('%s : Sequencer  stored on database', experiment_name)

        # deleting RunParameter and RunInfo
        logger.info('%s : Deleting RunParameter and RunInfo',experiment_name)
        os.remove(l_run_parameter)
        os.remove(l_run_info)
        # Handling Sample Sheet file
        #exp_name_id = RunProcess.objects.get(runName__exact = experiment_name).id
        sample_sheet_tmp_dir = os.path.join(wetlab_config.RUN_TEMP_DIRECTORY_RECORDED,run_process.get_run_id())
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
                    raise ValueError ('Unable to copy Sample Sheet on remote server')
            # Update run to Sample Sent

            run_process.set_run_state('Sample Sent')
            # Cleaning up the temporary folder with the sample sheet
            try:
                shutil.rmtree(sample_sheet_tmp_dir)
                logger.debug('%s : Deleted temporary folder containing the samplesheet', experiment_name)
            except Exception as e:
                string_message = experiment_name +' : Unable to delete temporary folder with the sample sheet ' + sample_sheet_tmp_dir
                logging_errors (string_message, True, True)
                logger.debug ('%s : End function for handling NextSeq run with exception', experiment_name)
                raise ValueError ('Unable to delete local Sample Sheet')

        else:
            string_message = experiment_name + ' : sample sheet not found on local Directory ' + sample_sheet_tmp_dir
            logging_errors (string_message, False, True)
            handling_errors_in_run (experiment_name, '22')
            logger.debug ('%s : End function for handling NextSeq run with exception', experiment_name)
            raise ValueError ('Sample Sheet not found in local folder')

        logger.debug ('%s : End function for handling NextSeq run ', experiment_name)
        return run_process
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
