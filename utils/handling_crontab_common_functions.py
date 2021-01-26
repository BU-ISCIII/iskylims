import logging
from logging.config import fileConfig
from logging.handlers import RotatingFileHandler

import os
from smb.SMBConnection import SMBConnection
import socket
import traceback
import xml.etree.ElementTree as ET
from datetime import datetime

from iSkyLIMS_wetlab.models import RunProcess, RunStates, Projects, RunningParameters, SambaConnectionData, EmailData, ConfigSetting
from .generic_functions import get_samba_connection_data, get_email_data, send_error_email_to_user, find_xml_tag_text
from iSkyLIMS_wetlab.wetlab_config import *

def create_new_sequencer_lab_not_defined (sequencer_name,num_of_lanes, experiment_name):
    '''
    Description:

        creates a new entry in database wit only the sequencer name and the lane numbers
    Input:
        sequencer_name    # sequencer name
        num_of_lanes        # number of lanes
        experiment_name     # experiment name
    Functions:
        get_sequencer_lanes_number_from_file # located at this file
    Constants:
        EMPTY_FIELDS_IN_SEQUENCER
    Return:
        new_sequencer_obj
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting function create_new_sequencer_lab_not_defined', experiment_name)
    seq_data = {}
    for item in EMPTY_FIELDS_IN_SEQUENCER :
        seq_data[item] = None
    seq_data['sequencerNumberLanes'] = num_of_lanes
    seq_data['sequencerName'] = sequencer_name
    new_sequencer_obj = SequencerInLab.objects.create_sequencer_in_lab(seq_data)
    logger.info('%s : Created the new sequencer in database' , experiment_name )
    logger.debug ('%s : End function create_new_sequencer_lab_not_defined', experiment_name)
    return new_sequencer_obj


def fetch_remote_file (conn, run_dir, remote_file, local_file) :
    '''
    Description:
        Function will fetch the file from remote server and copy on local
        directory
    Input:
        conn    # Samba connection object
        run_dir # run folder to fetch the file
        remote_file # file name to fetch on remote server
        local_file # local copy of the file fetched
    Return:
        local_file if the file was successfuly copy.
        Exception if file could not be fetched
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting function for fetching remote file', run_dir)
    with open(local_file ,'wb') as r_par_fp :
        try:
            samba_folder = get_samba_shared_folder()
            conn.retrieveFile(samba_folder, remote_file, r_par_fp)
            logger.info('%s : Retrieving the remote %s file for %s ', run_dir, remote_file, local_file)
        except Exception as e:
            string_message = 'Unable to fetch the ' + local_file + ' file on folder : ' + run_dir
            logging_errors (string_message, False, True)
            os.remove(local_file)
            logger.debug ('%s : End function for fetching remote file', run_dir)
            raise Exception('File not found')
    logger.debug ('%s : End function for fetching remote file', run_dir)
    return local_file


def get_new_runs_from_remote_server (processed_runs, conn, shared_folder):
    '''
    Description:
        The function fetch the folder names from the remote server and
        returns a list containing the folder names that have not been
        processed yet.
    Input:
        processed_runs  # full path and name of the file
        conn # samba connection object
        shared_folder   # shared folder in the remote server
        logger          # log object
    Return:
        new runs
    '''
    logger = logging.getLogger(__name__)
    logger.debug('Starting function get_new_runs_on_remote_server' )
    new_runs = []
    run_data_root_folder = os.path.join('/', get_samba_application_shared_folder() )
    logger.debug('Shared folder  on remote server is : %s', run_data_root_folder)
    run_folder_list = conn.listPath( shared_folder, run_data_root_folder)
    for sfh in run_folder_list:
        if sfh.isDirectory:
            folder_run = sfh.filename
            if (folder_run == '.' or folder_run == '..'):
                continue
            # if the run folder has been already process continue searching
            if folder_run in processed_runs:
                continue
            else:
                logger.info (' %s  : Found new folder run ',folder_run)
                new_runs.append(folder_run)
    logger.debug('End function get_new_runs_on_remote_server' )
    return new_runs

def get_run_platform_from_file (l_run_parameter) :
    '''
    Description:
        The function will get the run platform for the xml element tag in the
        file and it will return the platform used
    Input:
        l_run_parameter  # file to find the tag
    Variables:
        platform # name of the experiment found in runParameter file
    Return:
        platform
    '''
    platform = find_xml_tag_text (l_run_parameter, APPLICATION_NAME_TAG)

    return platform



def get_samba_application_shared_folder():
    '''
    Description:
        The function get in database the shared application folder used for samba connections
        Application folder is a sub_directory of the shared folder. In many cases this value is empty
    Return:
        samba_folder_name
    '''
    return SambaConnectionData.objects.last().get_samba_application_folder_name()


def get_samba_shared_folder():
    '''
    Description:
        The function get in database the shared folder used for samba connections
    Return:
        samba_folder_name
    '''
    return SambaConnectionData.objects.last().get_samba_folder_name()


def handling_errors_in_run (experiment_name, error_code):
    '''
    Description:
        Function will manage the error situation where the run must be
        set to run state ERROR
    Input:
        experiment_name # name of the run to be updated
        error_code      # Error code
        state_when_error # Run state when error was found
    Constants:
        SAMBA_SHARED_FOLDER_NAME
    Return:
        True
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting function handling_errors_in_run' , experiment_name)
    logger.info('%s : Set run to ERROR state',  experiment_name)
    run_process_obj = RunProcess.objects.filter(runName__exact = experiment_name).last()
    run_process_obj.set_run_error_code(error_code)
    logger.debug ('%s : End function handling_errors_in_run' , experiment_name)
    return True


def logging_errors(string_text, showing_traceback , print_on_screen ):
    '''
    Description:
        The function will log the error information to file.
        Optional can send an email to inform about the issue
    Input:
        logger # contains the logger object
        string_text # information text to include in the log
    Functions:
        send_error_email_to_user # located on utils.generic_functions
    Constant:
        SENT_EMAIL_ON_ERROR
        EMAIL_USER_CONFIGURED
    Variables:
        subject # text to include in the subject email
    '''
    logger = logging.getLogger(__name__)
    logger.error('-----------------    ERROR   ------------------')
    logger.error(string_text )
    if showing_traceback :
        logger.error('################################')
        logger.error(traceback.format_exc())
        logger.error('################################')
    logger.error('-----------------    END ERROR   --------------')
    email_data = get_email_data()
    if email_data['SENT_EMAIL_ON_ERROR'] :
        subject = 'Error found on wetlab when running crontab'
        send_error_email_to_user (subject, string_text, email_data['USER_EMAIL'],
                            [email_data['USER_EMAIL']])
    if print_on_screen :
        from datetime import datetime
        print('********* ERROR **********')
        print(string_text)
        print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        print('When processing run . Check log for detail information')
        print('******* END ERROR ********')
    return ''

def logging_warnings(string_text, print_on_screen ):
    '''
    Description:
        The function will log the error information to file.
        Optional can send an email to inform about the issue
    Input:
        logger # contains the logger object
        string_text # information text to include in the log
    '''
    logger = logging.getLogger(__name__)
    logger.warning('-----------------    WARNING   ------------------')
    logger.warning(string_text )
    logger.warning('-----------------    END WARNING   --------------')
    if print_on_screen :
        from datetime import datetime
        print('******* WARNING ********')
        print(string_text)
        print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        print('When processing run . Check log for detail information')
        print('**** END WARNING *******')
    return ''




def open_log(config_file):
    '''
    Description:
        The function will create the log object to write all logging information
    Input:
        logger_name    # contains the logger name that will be included
                        in the log file
    Constant:
        LOGGING_CONFIG_FILE
    Return:
        logger object
    '''
    fileConfig(config_file)
    logger = logging.getLogger(__name__)
    return logger


def nextseq_parsing_run_info_and_parameter_information(l_run_info, l_run_parameter, experiment_name) :
    '''
    Description:
        The function is called for parsing the RunInfo and RunParameter  files.
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
        SETUP_TAG
        NUMBER_OF_LANES_TAG
        APPLICATION_NAME_TAG
     Return:
        parsing_data
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting function parsing_run_information', experiment_name)
    running_data={}
    parsing_data = {}
    image_channel=[]

    #################################################
    ## parsing RunInfo.xml file
    #################################################
    run_data=ET.parse(l_run_info)
    run_root=run_data.getroot()
    logger.info('%s : parsing the runInfo.xml file ', experiment_name)
    p_run=run_root[0]
    ## getting the common values NextSeq and MiSeq
    logger.info('%s  : Fetching Flowcell and FlowcellLayout data ', experiment_name)
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
            string_message = experiment_name + ' : Parameter ' + field  + ' unable to fetch in RunParameter.xml'
            logging_warnings(string_message, False)

    for i in run_root.iter('Name'):
        image_channel.append(i.text)
    running_data['ImageChannel']=image_channel
    running_data['ImageDimensions']=p_run.find('ImageDimensions').attrib

    ## get the instrument for NextSeq run
    parsing_data['instrument'] = parameter_data_root.find('InstrumentID').text
    ## get the nuber of lanes in case sequencer lab is not defined

    param_in_setup = ['ApplicationVersion', 'NumTilesPerSwath']
    for i in range(len(param_in_setup)):
        try:
            running_data[param_in_setup[i]] = parameter_data_root.find('Setup').find(param_in_setup[i]).text
        except:
            string_message = experiment_name + ' : Parameter ' + param_in_setup[i] + ' not found in RunParameter.xml'
            logging_warnings(string_message, False)
            continue
    for setup_field in FIELDS_TO_FETCH_FROM_SETUP_TAG:

        try:
            running_data[setup_field] = parameter_data_root.find(SETUP_TAG).find(setup_field).text

        except:
            running_data[setup_field] = ''
            string_message = experiment_name + ' : Parameter in Setup -- ' + setup_field + ' unable to fetch in RunParameter.xml'
            logging_errors(string_message, False, True)

    import pdb; pdb.set_trace()
    ##############################################
    ## updating the date fetched from the Date tag for run and project
    ##############################################
    date = p_run.find('Date').text
    logger.debug('%s : Found date that was recorded the Run %s', experiment_name , date)
    run_date = datetime.strptime(date, '%y%m%d')
    parsing_data['running_data'] = running_data
    parsing_data['run_date'] = run_date
    import pdb; pdb.set_trace()
    logger.debug('%s : End function nextseq_parsing_run_information', experiment_name)
    return parsing_data
