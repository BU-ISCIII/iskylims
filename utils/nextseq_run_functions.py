import xml.etree.ElementTree as ET
from datetime import datetime
import os
import shutil
from .generic_functions import *
from iSkyLIMS_wetlab.models import *
from iSkyLIMS_drylab.models import Machines, Platform

from django_utils.models import Center
from .sample_sheet_utils import get_experiment_library_name, get_projects_in_run

def check_completion_success (l_run_completion):
    '''
    Description:
        The function will check if the run was success 
    Input:
        l_run_completion  # local path for the run completion file
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
    logger.debug ('Starting function for check_completion_success')
    # check if NextSEq run have been successful completed
    status_run = find_xml_tag_text (l_run_completion, wetlab_config.COMPLETION_TAG )
    if  status_run != wetlab_config.COMPLETION_SUCCESS:
        string_message = 'Run status was ' + status_run  
        logging_errors (string_message, False, False)
        return False
    else:
        logger.info ('Run successfuly completed ')
    logger.debug ('End function for check_completion_success')
    return True

def nextseq_parsing_run_information(l_run_info, l_run_parameter) :
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
        run_id      # contains the id value of the run
    Import:
        xml.etree.ElementTree 
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
    logger.debug ('Starting function nextseq_parsing_run_information')
    running_data={}
    image_channel=[]
    
    #################################################
    ## parsing RunInfo.xml file
    #################################################
    run_data=ET.parse(l_run_info)
    run_root=run_data.getroot()
    logger.info('parsing the runInfo.xml file ')
    p_run=run_root[0]
    ## getting the common values NextSeq and MiSeq
    running_data['Flowcell']=p_run.find('Flowcell').text
    running_data['FlowcellLayout']=p_run.find('FlowcellLayout').attrib
    #################################################
    ## parsing RunParameter.xml file
    #################################################
    logger.info('Parsing the runParameter.xml file  ' )
    parameter_data=ET.parse(l_run_parameter)
    parameter_data_root=parameter_data.getroot()
    p_parameter=parameter_data_root[1]
    ## getting the common values NextSeq and MiSeq
    running_data['RunID']=parameter_data_root.find('RunID').text
    running_data['ExperimentName']=parameter_data_root.find('ExperimentName').text
    running_data['RTAVersion']=parameter_data_root.find('RTAVersion').text
    running_data['Chemistry']=parameter_data_root.find('Chemistry').text
    running_data['RunStartDate']=parameter_data_root.find('RunStartDate').text
    running_data['RunManagementType']=parameter_data_root.find('RunManagementType').text
    running_data['ApplicationVersion']=p_parameter.find('ApplicationVersion').text
    running_data['NumTilesPerSwath']=p_parameter.find('NumTilesPerSwath').text
    running_data['SystemSuiteVersion']=parameter_data_root.find('SystemSuiteVersion').text 
    running_data['LibraryID']=parameter_data_root.find('LibraryID').text 
    running_data['AnalysisWorkflowType']=parameter_data_root.find('AnalysisWorkflowType').text  
    running_data['PlannedRead1Cycles']=parameter_data_root.find('PlannedRead1Cycles').text  
    running_data['PlannedRead2Cycles']=parameter_data_root.find('PlannedRead2Cycles').text  
    running_data['PlannedIndex1ReadCycles']=parameter_data_root.find('PlannedIndex1ReadCycles').text  
    running_data['PlannedIndex2ReadCycles']=parameter_data_root.find('PlannedIndex2ReadCycles').text    

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
    logger.debug('Found the de date that was recorded the Run %s', date)
    run_date = datetime.datetime.strptime(date, '%y%m%d')
    
    logger.debug('End function nextseq_parsing_run_information')
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
    logger.debug ('Starting function manage_nextseq_in_samplesent')
    experiment_name = run_object_name.get_run_name()
    logger.info('Manage %s in Sample Sent state ', experiment_name)
    run_folder = RunningParameters.objects.get(runName_id = run_object_name).get_run_folder()
    l_run_completion = os.path.join(wetlab_config.RUN_TEMP_DIRECTORY, wetlab_config.RUN_COMPLETION)
    s_run_completion = os.path.join(run_folder, wetlab_config.RUN_COMPLETION)
    
    try:
        l_run_completion = fetch_remote_file (conn, run_folder, s_run_completion, l_run_completion)
        logger.info('Sucessfully fetch of Completion status file')
        # Get status value for run completion
        completion_status = check_completion_success (l_run_completion)
        
    except Exception as e:
        logger.info ('Completion status file is not present on the run folder')
        logger.info ('Moving run to Processing run state')
        run_updated = run_object_name.set_run_state('Processing Run')
        logger.debug ('End function for handling NextSeq run with exception')
        raise
    finally :
        logger.info ('Deleting local copy of completion status')
        # cleaning up the completion  in local temporary file
        os.remove(l_run_completion)
    
    if not completion_status :
        string_message = 'Run status was ' + status_run  
        logging_errors (string_message, False, False)
        # Set tun to error state
        run_updated = handling_errors_in_run (experiment_name)
        logger.debug ('End function for handling NextSeq run with exception')
        raise ValueError ('Run was CANCELLED')
    else:
        run_updated = run_object_name.set_run_state('Processing Run')
        logger.info('Run %s is now on Processing Run state', experiment_name)
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
    logger.debug ('Starting function manage_nextseq_in_processing_run')
    experiment_name = run_object_name.get_run_name()
    logger.info('Manage %s in Processing Run state ', experiment_name)
    run_folder = RunningParameters.objects.get(runName_id = run_object_name).get_run_folder()
    l_run_completion = os.path.join(wetlab_config.RUN_TEMP_DIRECTORY, wetlab_config.RUN_COMPLETION)
    s_run_completion = os.path.join(run_folder, wetlab_config.RUN_COMPLETION)
    
    try:
        l_run_completion = fetch_remote_file (conn, run_folder, s_run_completion, l_run_completion)
        logger.info('Sucessfully fetch of Completion status file')
        # Get status value for run completion
        completion_status = check_completion_success (l_run_completion)
        # Get the date and time when the RunCompletionStatus is created
        run_completion_attributes = get_attributes_remote_file (conn, run_folder, s_run_completion)
        run_completion_date = datetime.datetime.fromtimestamp(int(run_completion_attributes.create_time)).strftime('%Y-%m-%d %H:%M:%S')
        run_completion_date = run_object_name.set_run_completion_date(run_completion_date)
        
    except Exception as e:
        logger.info ('Completion status file is not present on the run folder')
        logger.info ('Moving run to Processing run state')
        logger.debug ('End function for handling manage_nextseq_in_processing_run with exception')
        raise
    finally :
        logger.info ('Deleting local copy of completion status')
        # cleaning up the completion  in local temporary file
        os.remove(l_run_completion)

    if not completion_status :
        string_message = 'Run status was ' + status_run  
        logging_errors (string_message, False, False)
        # Set tun to error state
        run_updated = handling_errors_in_run (experiment_name)
        logger.debug ('End function for manage_nextseq_in_processing_run with exception')
        raise ValueError ('Run was CANCELLED')
    else:
        run_updated = run_object_name.set_run_state('Processed Run')
        logger.info('Run %s is now on Processed Run state', experiment_name)
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
        get_projects_in_run # located at utils.sample_sheet_utils
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
        running_parameters  # dictionary with parsing information from 
                            RunParameter and RunInfo
        s_run_info  # path of RunInfo on the remote server
        s_sample_sheet  # path of SampleSheet on the remote server
    Return:
        number_of_cycles 
        file_content
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function handle_nextseq_recorded_run')
    if RunProcess.objects.filter(runName__exact = experiment_name , runState__exact ='Recorded').exists():
        run_process = RunProcess.objects.get(runName__exact = experiment_name )
        
        # Fetch run info
        l_run_info = os.path.join(wetlab_config.RUN_TEMP_DIRECTORY, wetlab_config.RUN_INFO)
        s_run_info = os.path.join(new_run,wetlab_config.RUN_INFO)
        try:
            l_run_info = fetch_remote_file (conn, new_run, s_run_info, l_run_info)
            logger.info('Sucessfully fetch of RunInfo file')
        except Exception as e:
            logger.info('Exception fetched $s ' , e)
            logger.info('Skiping the run $s , due to the error', new_run)
            # cleaning up the RunParameter in local temporaty file
            os.remove(l_run_parameter)
            logger.debug ('End function for handling NextSeq run with exception')
            raise   # returning to handle next run folder
        
        # Parsing RunParameter and Run Info

        running_parameters, run_date, instrument = nextseq_parsing_run_information(l_run_info, l_run_parameter)
        
        run_date = run_process.set_run_date(run_date)
        run_process_id = run_process.id
        run_parameters = RunningParameters.objects.create_running_parameters(running_parameters, run_process_id)
        if  Machines.objects.filter(machineName__exact = instrument).exists() :
            sequencer = Machines.objects.get(machineName__exact = instrument)
            instrument = run_process.set_sequencer(sequencer)
            logger.info('Sequencer  stored on database')
        else:
            string_message = instrument + ' has been not defined on machines '
            logging_errors(string_message, False, True)
        #update the project state 
        projects = Projects.objects.filter(runprocess_id = run_process)
        for project in projects:
                p_state = project.set_project_state ( 'Sample Sent')
        # deleting RunParameter and RunInfo
        os.remove(l_run_parameter)
        os.remove(l_run_info)

        # Handling Sample Sheet file
        #exp_name_id = RunProcess.objects.get(runName__exact = experiment_name).id
        sample_sheet_tmp_dir = os.path.join(wetlab_config.RUN_TEMP_DIRECTORY_RECORDED,str(run_process_id))
        l_sample = os.path.join(sample_sheet_tmp_dir, wetlab_config.SAMPLE_SHEET)
        
        if os.path.exists(l_sample):
            
            if wetlab_config.COPY_SAMPLE_SHEET_TO_REMOTE :
                # copy Sample heet file to remote directory
                logger.info('Copy sample sheet to remote folder %s', new_run)
                s_sample= os.path.join(new_run, wetlab_config.SAMPLE_SHEET)
                try:
                    sample_sheet_copied = copy_to_remote_file  (conn, new_run, s_sample, l_sample)
                    logger.info('Sucessfully copy Sample sheet to remote folder')
                except Exception as e:
                    string_message = 'Unable to copy the Sample Sheet to remote folder ' + new_run
                    logging_errors (string_message, True, False)
                    logger.info ('Deleting local copy of completion status ')
                    logger.debug ('End function for handling NextSeq run with exception')
                    raise 
            # Update run to Sample Sent
            
            run_process.set_run_state('Sample Sent')
            # Cleaning up the temporary folder with the sample sheet
            try:
                shutil.rmtree(sample_sheet_tmp_dir)
                logger.debug('Deleted temporary folder containing the samplesheet')
            except Exception as e:
                string_message = 'Unable to delete temporary folder with the sample sheet ' + sample_sheet_tmp_dir
                logging_errors (string_message, True, True)
                logger.debug ('End function for handling NextSeq run with exception')
                raise ValueError ('Unable to delete local Sample Sheet')

        else:
            string_message = 'sample sheet not found on local Directory ' + sample_sheet_tmp_dir
            logging_errors (string_message, False, True)
            run_updated = handling_errors_in_run (experiment_name)
            logger.debug ('End function for handling NextSeq run with exception')
            raise ValueError ('Sample Sheet not found in local folder')

        logger.debug ('End function for handling NextSeq run ')
        return run_process
    else:
        if not RunProcess.objects.filter(runName__exact = experiment_name).exists():
            string_message = 'Run ' + experiment_name +' is not yet defined on database'
            logging_errors (string_message, False, False)
            logger.debug ('End function handle_nextseq_recorded_run with exception')
            raise ValueError ('Run not defined yet')
        else:
            run_state = RunProcess.objects.filter(runName__exact = experiment_name).get_state()
            string_message = 'Run ' + experiment_name +' is in state ' + run_state + '.  Should be in Recorded'
            logging_errors (string_message, False, False)
            logger.debug ('End function fhandle_nextseq_recorded_run with exception')
            raise ValueError ('Run on wrong state')
    return




