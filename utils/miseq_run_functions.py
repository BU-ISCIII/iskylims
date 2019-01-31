import xml.etree.ElementTree as ET
from datetime import datetime

from .generic_functions import *
from iSkyLIMS_wetlab.models import *
from iSkyLIMS_drylab.models import Machines, Platform

from django_utils.models import Center
from .sample_sheet_utils import get_experiment_library_name, get_projects_in_run

def check_miseq_completion_run (conn, experiment_name, log_folder):
    '''
    Description:
        The function will check the run log files to check if all cycles
        were running on the sequencer and if the last log does not
        mention that run was canceled
    Input:
        conn    # contains the samba connection object
        log_folder  # remote path where are located the logs for miseq
        experiment_name # name for this run
    Functions:
        get_miseq_run_cycles # located at utils.wetlab_misc_utilities
        get_latest_miseq_log #  located at utils.wetlab_misc_utilities
    Constants:
        COMPLETION_SUCCESS
    Variables:
        latest_log # string containing the latest log information
        run_cycles # number of cycles to be completed in the run
        log_cycles # number of maximum cycles logged on the Log folder
        status_run # status of the run
        run_completion_date # date and time of the completion
    Return:
        status_run and run_completion_date
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function check_miseq_completion_run')
    run_completion_date =''
    run_object = RunProcess.objects.get(runName__exact = experiment_name)
    number_of_cycles =  RunningParameters.objects.get(runName_id__exact = run_object).get_number_of_cycles()

    try:
        log_cycles, log_file_content = get_latest_miseq_log(conn, log_folder)
    except :
        string_message = 'Unable to fetch the log files for run ' + experiment_name
        logging_errors( string_message, False, False)
        handling_errors_in_run (experiment_name, '18' )
        logger.debug('End function check_miseq_completion_run with IOError exception')
        raise
        #raise IOError ('Unable to fetch log file')
    if 'Cancel' in log_file_content :
        status_run = 'Canceled'
    elif log_cycles != number_of_cycles :
        status_run = 'still_running'
    else:
        status_run = wetlab_config.COMPLETION_SUCCESS
        last_line_in_file = log_file_content.split('\n')[-2]
        last_log_time = last_line_in_file.split(' ')[0:2]
        last_log_time[0] = str('20'+ last_log_time[0])
        string_completion_date = ' '.join(last_log_time)
        run_completion_date = datetime.datetime.strptime(string_completion_date, '%Y-%m-%d %H:%M:%S.%f')
    logger.info('Status of the run is : %s ', status_run)
    logger.debug('End function check_miseq_completion_run')
    return status_run, run_completion_date



def get_latest_miseq_log(conn, log_folder) :
    '''
    Description:
        The function will find the latest log file for the input folder
    Input:
        conn        # samba connection object
        log_folder  # folder name
    Constant:
        RUN_TEMP_DIRECTORY
    Variables:
        file_content # text information in the latest log
        file_remote # file name on the remote server
        remote_file_list # samba object with the list of files
        max_cycle   # Counter with the maximum cycle found
        number_of_cycles # will add the partial number of cycles
    Return:
        number_of_cycles
        file_content
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function for fetching remote file')
    remote_file_list = conn.listPath( wetlab_config.SAMBA_SHARED_FOLDER_NAME, log_folder)
    max_cycle = -1
    for sfh in remote_file_list:
        if sfh.isDirectory:
            continue
        file_remote = sfh.filename
        if file_remote.endswith('.log') :
            log_file = re.search('.*_Cycle(\d+)_.*',file_remote)

            cycle_number = int(log_file.group(1))
            if cycle_number > max_cycle :
                max_cycle = cycle_number
                latest_log = file_remote

    temporary_log = os.path.join(wetlab_config.RUN_TEMP_DIRECTORY,'miseq_cycle.log')
    s_latest_log = os.path.join(log_folder,latest_log)
    with open(temporary_log ,'wb') as log_fp :
        temporary_log = fetch_remote_file (conn, log_folder, s_latest_log, temporary_log)

    with open (temporary_log, 'r') as fh :
        log_file_content = fh.read()
    os.remove(temporary_log)
    return max_cycle, log_file_content


def miseq_parsing_run_information(run_info, run_parameter):
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
    logger.debug ('Starting function for parsing xml file for miSeq run')
    running_data={}
    running_data_fields = ['Flowcell','FlowcellLayout','RunID','ExperimentName','RTAVersion',
            'Chemistry','RunStartDate','RunManagementType','ApplicationVersion','NumTilesPerSwath',
            'PlannedRead1Cycles','PlannedIndex1ReadCycles','PlannedIndex2ReadCycles','PlannedRead2Cycles',
            'SystemSuiteVersion','LibraryID','AnalysisWorkflowType','ImageChannel','ImageDimensions']
    # Initialize all fields to default -> empty
    for fields in running_data_fields:
        running_data[fields] = ""

    #################################################
    ## parsing RunInfo.xml file
    #################################################
    run_data=ET.parse(run_info)
    run_root=run_data.getroot()
    logger.info('parsing the runInfo.xml file ' )
    p_run=run_root[0]
    ## getting the common values NextSeq and MiSeq
    running_data['Flowcell']=p_run.find('Flowcell').text
    running_data['FlowcellLayout']=p_run.find('FlowcellLayout').attrib

    #################################################
    ## parsing RunParameter.xml file
    #################################################
    logger.info('Parsing the runParameter.xml file ')
    parameter_data=ET.parse(run_parameter)
    parameter_data_root=parameter_data.getroot()
    p_parameter=parameter_data_root[1]
    ## getting the common values NextSeq and MiSeq
    running_data['RunID']=parameter_data_root.find('RunID').text
    running_data['ExperimentName']=parameter_data_root.find('ExperimentName').text
    running_data['RTAVersion']=parameter_data_root.find('RTAVersion').text
    running_data['Chemistry']=parameter_data_root.find('Chemistry').text
    running_data['RunStartDate']=parameter_data_root.find('RunStartDate').text
    running_data['RunManagementType']=parameter_data_root.find('RunManagementType').text
    running_data['ApplicationVersion']=parameter_data_root.find('Setup').find('ApplicationVersion').text
    running_data['NumTilesPerSwath']=parameter_data_root.find('Setup').find('NumTilesPerSwath').text

    ## get the length index number for reads and indexes
    for run_info_read in parameter_data_root.iter('RunInfoRead'):
        if run_info_read.attrib['Number'] == '1' :
            running_data['PlannedRead1Cycles'] = run_info_read.attrib['NumCycles']
        elif run_info_read.attrib['Number'] == '2' :
            running_data['PlannedIndex1ReadCycles'] = run_info_read.attrib['NumCycles']
        elif run_info_read.attrib['Number'] == '3' :
            running_data['PlannedIndex2ReadCycles'] = run_info_read.attrib['NumCycles']
        elif run_info_read.attrib['Number'] == '4' :
            running_data['PlannedRead2Cycles'] = run_info_read.attrib['NumCycles']

    ## setting empty values for MiSeq
    ## get the instrument for MiSeq run located on RunInfo.xml
    instrument =p_run.find('Instrument').text


    ##############################################
    ## updating the date fetched from the Date tag for run and project
    ##############################################
    date = p_run.find('Date').text
    logger.debug('Found the de date that was recorded the Run %s', date)
    run_date = datetime.datetime.strptime(date, '%y%m%d')

    logger.debug ('End function for parsing xml file for miSeq run')
    return running_data, run_date, instrument


def run_waiting_for_sample_sheet (experiment_name):
    '''
    Description:
        The function get the value of the sample sheet field in database
        If value is empty , run is waiting for smaple sheet.
    Input:
        experiment_name  # experiment name to be checked
    Import:
        RunProccess     # from iSkyLIMS_wetlab.models
    Variable:
        sample_sheet_value  # contains the value of the sample sheet
                            stored in database
    Return:
        True if run is waiting to fetch the sample sheet
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function run_waiting_for_sample_sheet')
    sample_sheet_value = RunProcess.objects.get(runName__exact = experiment_name).get_sample_file()
    if sample_sheet_value == '':
        logger.debug ('End function. Run is waiting for sample sheet')
        return True
    else:
        logger.debug ('End function. Sample sheet already fetched')
        return False


def save_new_miseq_run ( experiment_name, run_date, instrument) :
    '''
    Description:
        The function will create a new instance for RunProcess
        The information stored is :  experiment_name, run_date, sequencerModel
        run_state, centerRequestedBy, index_library . Sample Sheet will
        be empty . This field will be updated when
        fetching the sample sheet from remote server
    Input:
        experiment_name # name of the experiment
        run_date        # date of run creation
        instrument      # sequencer name
    Functions:
        logging_errors   # locate at utils.generic_functions
    Imports:
        Center      # from django_utils.models
        Machines    # from iSkyLIMS_drylab.models
    Constants:
        DEFAULT_CENTER
    Variables:
        center_requested_by # Center object for run requested center
        sequencer       # sequencer object that contains the reference to
                        machines class on iSkyLIMS_wetlab.models
    '''
    logger = logging.getLogger(__name__)
    logger.debug('Executing the function save_new_miseq_run' )
    ## create a new entry on runProcess database
    if Center.objects.filter(centerAbbr__exact = wetlab_config.DEFAULT_CENTER).exists():
        center_requested_by = Center.objects.get(centerAbbr__exact = wetlab_config.DEFAULT_CENTER)
    else:
        string_message = 'The requested center ' +  wetlab_config.DEFAULT_CENTER + ' defined in config wetlab file does not exist'
        logging_errors(string_message, False, False)
        logger.info('Using the first center defined in database')
        center_requested_by = Center.objects.all().first()

    # Get the machines reference to add to sequencer field in the run
    if  Machines.objects.filter(machineName__exact = instrument).exists() :
        sequencer = Machines.objects.get(machineName__exact = instrument)
        logger.info('Updated the run date and sequencer used for the runProcess table ')
    else:
        string_message = instrument + ' has been not defined on machines '
        logging_errors(string_message, False, False)

    run_state = RunStates.objects.get(runStateName__exact = 'Recorded')
    run_process = RunProcess(runName=experiment_name,sampleSheet= '',
                            run_date = run_date, index_library = '',
                            state = run_state, centerRequestedBy = center_requested_by,
                            sequencerModel = sequencer)
    run_process.save()
    logger.info('Create the new runProcess into database')
    logger.debug('Exiting the function save_new_miseq_run' )
    return run_process


def save_miseq_projects_found (projects_users , experiment_name, library_name):
    '''
    Description:
        The function will save the project information fetched from
        sample sheet in Projects database table
    Input:
        projects_users # dictionary contains projects and users
        library_name # library name fetched from sample sheet
    Imports:
        LibraryKit  # from iSkyLIMS_wetlab.models
        RunProcess  # from iSkyLIMS_wetlab.models
    Constants:
        DEFAULT_LIBRARY_KIT
    Variable:
        base_space_file # path modification from sample_sheet_on_database
                        to store in the Project
        library_kit # library kit object used in the project
        p_data      # New Project object to save into database
        run         # RunProcess object for miseq run
        sample_sheet_on_database  # sample sheet value stored in RunProcess
        string_message # message to be showed in the log file

    Return:
        Returns an empty value
    '''
    logger = logging.getLogger(__name__)
    logger.debug('Executing the function save_miseq_projects_found' )

    if LibraryKit.objects.filter(libraryName__exact = wetlab_config.DEFAULT_LIBRARY_KIT).exists():
        library_kit = LibraryKit.objects.get(libraryName__exact = wetlab_config.DEFAULT_LIBRARY_KIT)
    else:
        string_message = 'The default library ' +  wetlab_config.DEFAULT_LIBRARY_KIT + ' defined in config wetlab file does not exist'
        logger_errors(logger, string_message)
        logger.info('Using the first library kit defined in database')
        library_kit = LibraryKit.objects.all().first()

    run_process = RunProcess.objects.get(runName = experiment_name)
    sample_sheet_on_database = run_process.get_sample_file()
    base_space_file = os.path.join(settings.MEDIA_URL.replace('/',''), sample_sheet_on_database)

    for project, user  in projects_users.items():
        userid=User.objects.get(username__exact = user)
        p_data=Projects(runprocess_id = run_process, projectName = project,
                        user_id = userid, procState = 'Sample Sent',
                        baseSpaceFile = base_space_file,
                        LibraryKit_id = library_kit, libraryKit = library_name)
        p_data.save()
    logger.info('Updated Projects table with the new projects found')
    logger.debug('End the function save_miseq_projects_found' )
    return ''


def store_sample_sheet_in_run (l_sample_sheet, experiment_name ) :
    '''
    Description:
        The function will move the sample sheet from the local temporary
        folder to the folder destination defined in RUN_SAMPLE_SHEET_DIRECTORY
        It will update the field sampleSheet in database
    Input:
        l_sample_sheet  # local copy of sample sheet file
        experiment_name  # name of the run
    Import:
        datetime
        os
    Variable:
        run_updated     # RunProcess object
        new_sample_sheet_name  # renamed sample sheet including time to have
                                a unique file name
        now     # Present time of the server
        run_updated # RunProcess object for miseq run
        sample_sheet_on_database # sample sheet path used in database
        timestr     # Present time including 3 digits for miliseconds
    Return:
        sample_sheet_on_database
    '''
    logger = logging.getLogger(__name__)
    logger.debug('Starting the function store_sample_sheet_in_run')
    # Get the present time in miliseconds to add to have an unique file name
    now = datetime.datetime.now()
    timestr = now.strftime("%Y%m%d-%H%M%S.%f")[:-3]
    new_sample_sheet_name = 'SampleSheet' + timestr + '.csv'

    new_sample_sheet_file = os.path.join (settings.MEDIA_ROOT, wetlab_config.RUN_SAMPLE_SHEET_DIRECTORY, new_sample_sheet_name)
    logger.debug('new sample sheet name %s', new_sample_sheet_file)
    # Path to be included in database
    sample_sheet_on_database = os.path.join(wetlab_config.RUN_SAMPLE_SHEET_DIRECTORY, new_sample_sheet_name)
    ## Move sample sheet to final folder
    os.rename(l_sample_sheet, new_sample_sheet_file)
    # Update the run with the sample sheet information
    run_updated = RunProcess.objects.get(runName__exact = experiment_name)
    run_updated.sampleSheet = sample_sheet_on_database
    run_updated.save()

    logger.info('Updated runProccess table with the sample sheet')
    logger.debug('End function store_sample_sheet_in_run')
    return sample_sheet_on_database


def update_library_name_in_run (experiment_name, library_name) :
    '''
    Description:
        The function will update the library name in database
    Input:
        experiment_name  # name of the run
        library_name    # libary used in the run
    Variable:
        run         # RunProcess object
        updated_library # return value after updating the library
    Return
        updated_library
    '''
    logger = logging.getLogger(__name__)
    logger.debug('Starting the function update_library_name_in_run')
    run = RunProcess.objects.get(runName = experiment_name)
    updated_library = run.update_library(library_name)

    logger.debug('Endfunction update_library_name_in_run')
    return updated_library


def validate_sample_sheet (sample_sheet):
    '''
    Description:
        The function get the sample sheet file and will do some checks
        to validate the file.
    Input:
        sample_sheet  # full path for smaple sheet file
    Functions:
        get_projects_in_run # located at this file
        logging_errors   # locate at utils.wetlab_misc_utilities
    Import:
        User        # from django_utils
        Projects    # from iSkyLIMS_wetlab.models
    Variable:
        projects_users  # dictionary containing projects and user owner
        experiment_name # contains the experiment name from the sample sheet
        library_name  # contains the library name from the sample sheet
    Return:
        True if all checking are successful False if any of the check fails
    '''
    logger = logging.getLogger(__name__)
    logger.debug('Starting the function validate_sample_sheet' )
    # get the projects and user owner form sample sheet
    projects_users = get_projects_in_run(sample_sheet)
    if len(projects_users) == 0 :
        string_message = 'No sample lines have been found '
        logging_errors(string_message, False, False)
        logger.debug('End the function validate sample_sheet with error')
        return False
    for project in projects_users.keys() :
        if project == "":
            string_message = 'Empty projects have been found '
            logging_errors(string_message, False, False)
            logger.debug('End the function validate sample_sheet with error')
            return False
        if Projects.objects.filter(projectName__exact = project).exists():
            string_message = 'project name %s , already been used ' + project
            logging_errors(string_message, False, False)
            logger.debug('Exiting the function validate sample_sheet with error')
            return False
    for user in projects_users.values():
        if user == "":
            string_message =  'Empty users have been found '
            logging_errors(string_message, False, False)
            logger.debug('End the function validate sample_sheet with error')
            return False
        if ( not User.objects.filter(username__icontains = user).exists()):
            string_message = 'user name ' +user  +' is not defined in the system'
            logging_errors(string_message, False, False)
            logger.debug('Exiting the function validate sample_sheet with error')
            return False

    logger.debug('End the function validate_sample_sheet' )
    return True

def manage_miseq_in_processing_run (conn, run_name):
    '''
    Description:
        The function will look the run log files to check if run is
        completed. Run will remain in processing run if the latest
        cycle log did not reach the cycle number in the run
    Input:
        conn            #  Samba connection object
        run_name        # run object
    Constant:
        MAXIMUM_TIME_WAIT_RUN_COMPLETION
    Functions:
        need_to_wait_more   # located in utils.generic_functions
        check_miseq_completion_run # located in this file
    Variable:
        run_completion_date # completion date of the run
        run_folder      # folder on the remote server
        run_updated     # Updated run object
        status_run         # RunProcess status (Cancelled, still_running,
                            Completed as planned)
    Return
        experiment_name
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function manage_miseq_in_processing_run')
    experiment_name = run_name.get_run_name()
    run_folder = RunningParameters.objects.get(runName_id = run_name).get_run_folder()
    log_folder = os.path.join(wetlab_config.SAMBA_APPLICATION_FOLDER_NAME, run_folder, wetlab_config.RUN_LOG_FOLDER)
    try: # waiting for sequencer run completion
        status_run, run_completion_date = check_miseq_completion_run (conn, experiment_name, log_folder)
        if status_run == 'Cancel' :
            run_updated = handling_errors_in_run (experiment_name)
            run_updated = run_name.set_run_state('CANCELLED')
            run_updated.save()
            logger.debug ('End function manage_miseq_in_processing_run was cancelled')
            raise ValueError ('Run was canceled')
        elif status_run == 'still_running':
            if need_to_wait_more (experiment_name, wetlab_config. MAXIMUM_TIME_WAIT_RUN_COMPLETION):
                logger.info('Run $s still waiting for sequencer to finish', experiment_name)
            raise # returning to handle next run folder
        else:
            run_updated = run_name.set_run_state('Processed Run')
            run_name.run_finish_date = run_completion_date
            run_name.save()
            logger.debug ('End function manage_miseq_in_processing_run %s' , experiment_name)
            return experiment_name
    except ValueError :
        logger.debug ('End function manage_miseq_in_processing_run with error')
        raise
    except:
        string_message = 'Error when fetching the log file for the run ' + new_run
        logging_errors (string_message, False, False)
        logger.debug ('End function manage_miseq_in_processing_run with error')
        raise




def manage_miseq_in_samplesent(conn, run_name) :
    '''
    Description:
        The function will look the run log files to check if run is
        completed. Run will remain in processing run if the latest
        cycle log did not reach the cycle number in the run
    Input:
        conn            #  Samba connection object
        run_name        # run object
    Functions:
        check_miseq_completion_run # located in this file
    Variable:
        run_completion_date # completion date of the run
        run_folder      # folder on the remote server
        run_updated     # Updated run object
        status_run         # RunProcess status (Cancelled, still_running,
                            Completed as planned)
    Return
        experiment_name
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function manage_miseq_in_samplesent')
    experiment_name = run_name.get_run_name()
    run_folder = RunningParameters.objects.get(runName_id = run_name).get_run_folder()
    try: # waiting for sequencer run completion
        log_folder = os.path.join(wetlab_config.SAMBA_APPLICATION_FOLDER_NAME, run_folder, wetlab_config.RUN_LOG_FOLDER)

        status_run, run_completion_date = check_miseq_completion_run (conn, experiment_name, log_folder)
        if status_run == 'Cancel' :
            run_updated = handling_errors_in_run (experiment_name)
            run_updated = run_name.set_run_state('CANCELLED')
            run_updated.save()
            logger.debug ('End function manage_miseq_in_samplesent with error')
            raise ValueError ('Run was canceled')
        else :
            run_updated = RunProcess.objects.get(runName__exact = experiment_name).set_run_state('Processing run')
            logger.debug ('End function manage_miseq_in_samplesent for run %s' , experiment_name)
            return experiment_name
    except ValueError :
        raise
    except:
        string_message = 'Error when fetching the log file for the run ' + new_run
        logging_errors (string_message, False, False)
        logger.debug ('End function manage_miseq_in_samplesent with error')
        raise


def handle_miseq_run (conn, new_run, l_run_parameter, experiment_name) :
    '''
    Description:
        The function will find the latest log file for the input folder
    Input:
        conn        # samba connection object
        new_run     # folder remote directory for miseq run
        l_run-parameter  # local path for the run parameter file
        experiment_name  # name used on miseq run
    Functions:
        get_projects_in_run # located at utils.sample_sheet_utils
        get_experiment_library_name # located at utils.sample_sheet_utils
        fetch_remote_file   # located at utils.generic_functions

        miseq_parsing_run_information # located as this file
        need_to_wait_more       # locate at utils.generic_functions
        run_waiting_for_sample_sheet  # located as this file
        save_new_miseq_run      # located as this file
        store_sample_sheet_in_run # located as this file
        validate_sample_sheet   # located as this file
    Import:
        os
        RunningParameters
        RunProcess
    Constant:
        RUN_TEMP_DIRECTORY
        RUN_INFO
        SAMPLE_SHEET
        MAXIMUM_TIME_WAIT_SAMPLE_SHEET
    Variables:
        instrument  # name of the sequencer used in the run
        l_run_info  # local temporary copy of RunInfo
        l_sample_sheet  # local temporary copy of sample sheet
        new_run_parameters # new RunningParameters object created for this run
        new_run_process  # new RunProcess object created for this run

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
    logger.debug ('Starting function for handling miSeq run')
    # Fetch run info from remote server
    l_run_info = os.path.join(wetlab_config.RUN_TEMP_DIRECTORY, wetlab_config.RUN_INFO)
    s_run_info = os.path.join(wetlab_config.SAMBA_APPLICATION_FOLDER_NAME, new_run,wetlab_config.RUN_INFO)
    try:
        l_run_info = fetch_remote_file (conn, new_run, s_run_info, l_run_info)
        logger.info('Sucessfully fetch of RunInfo file')
    except Exception as e:
        string_message = 'Unable to fetch the RunInfo file on folder : ' + new_run
        logging_errors(string_message, True, False)
        logger.info('Skiping the run %s , due to the error on fetching RunInfo.xml', new_run)
        # cleaning up the RunParameter in local temporaty file
        os.remove(l_run_parameter)
        raise  ValueError ('RunInfo.xml file not found on run folder') # returning to handle next run folder

    # Parsing RunParameter and RunInfo
    running_parameters, run_date, instrument = miseq_parsing_run_information(l_run_info, l_run_parameter)
    if not RunProcess.objects.filter(runName__exact = experiment_name).exists():
        # Save run data and set run to "Recorded" state
        new_run_process = save_new_miseq_run (experiment_name, run_date, instrument)
        logger.info ('New RunProcess %s was stored on database', experiment_name)
    else :
        new_run_process = RunProcess.objects.get(runName__exact = experiment_name)
    # Update the running parameter table with the information
    new_run_parameters = RunningParameters.objects.create_running_parameters(running_parameters, new_run_process)
    logger.info('Running parameters have been stored on database for %s', experiment_name )

    # deleting temporary copy of RunParameter and RunInfo files
    os.remove(l_run_parameter)
    os.remove(l_run_info)
    logger.info('Deleted RunParameter and RunInfo files')
    if run_waiting_for_sample_sheet (experiment_name) :
        # Fetch sample sheet from remote server
        l_sample_sheet = os.path.join(wetlab_config.RUN_TEMP_DIRECTORY, wetlab_config.SAMPLE_SHEET)
        s_sample_sheet = os.path.join(wetlab_config.SAMBA_APPLICATION_FOLDER_NAME, new_run, wetlab_config.SAMPLE_SHEET)
        try:
            l_sample_sheet = fetch_remote_file (conn, new_run, s_sample_sheet, l_sample_sheet)
            logger.info('Sucessfully fetch of Sample Sheet file')
        except Exception as e:
            logger.info('Error %s' , e)
            logger.info('Unable to get Sample sheet on folder  %s', new_run)
            if need_to_wait_more (experiment_name, wetlab_config.MAXIMUM_TIME_WAIT_SAMPLE_SHEET) :
                string_message = 'Sample Sheet for' + new_run + '. Extended more time to fetch it'
                logging_warnings(string_message, True)
                # Delete running paramters object for allowing  to look for,
                # in the romote folder, again next time the process is executed
                new_run_parameters.delete()
                logger.info("Deleted running parameters from database ")
                logger.debug ('End function for handling miSeq run with error')
                raise ValueError ('Sample sheet not found in run folder')
            else:
                # set run state to Error state
                string_message = 'Time expiration for getting sample Sheet for ' + experiment_name
                logging_errors( string_message, False, True)
                handling_errors_in_run(experiment_name, '19')
                logger.debug ('End function for handling miSeq run with error')
                raise ValueError ('Time expiration for getting sample sheet')

        if validate_sample_sheet (l_sample_sheet) :
            
            projects_users = get_projects_in_run (l_sample_sheet)
            logger.debug('Fetched projects from sample sheet')

            experiment_name, library_name = get_experiment_library_name(l_sample_sheet)
            logger.debug('Fetched experiment name and library name from sample sheet')

            sample_sheet_on_database = store_sample_sheet_in_run (l_sample_sheet, experiment_name )
            updated_libary_name = update_library_name_in_run (experiment_name, library_name)

            save_miseq_projects_found (projects_users, experiment_name, library_name)
            run_updated = RunProcess.objects.get(runName__exact = experiment_name).set_run_state('Sample Sent')
            
            logger.info('Run %s is now on sample sheet state', experiment_name)
            logger.debug ('End function for handling miSeq ')
            return new_run
        else:
            # set run state to ERROR
            string_message = "Invalid Sample Sheet for " + experiment_name + 'on folder ' + new_run
            logging_errors (string_message, False, True)
            run_updated = handling_errors_in_run (experiment_name,'1')
            logger.debug ('End function for handling miSeq run with error')
            raise ValueError ('Invalid sample sheet')
    else:
        run_state = RunProcess.objects.get(runName__exact = experiment_name).get_run_state()
        string_message = 'Wrong state: ' + run_state +' when calling to handle_miseq_run function '
        logging_errors(string_message, False, True)
        logger.debug ('End function for handling miSeq run with error')
        raise ValueError ('Wrong run state ')


