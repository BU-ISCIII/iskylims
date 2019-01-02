#!/usr/bin/env python3

import sys, os, re
import xml.etree.ElementTree as ET
import time
import shutil
import locale
import datetime, time
#from  ..models import *
from iSkyLIMS_wetlab.models import *
from .interop_statistics import *

from smb.SMBConnection import SMBConnection
from iSkyLIMS_wetlab import wetlab_config
from iSkyLIMS_drylab.models import Machines, Platform
from .wetlab_misc_utilities import open_samba_connection , logging_errors, get_miseq_run_cycles
from .sample_sheet_utils import get_experiment_library_name, get_projects_in_run

from django.conf import settings
from django_utils.models import Center


def get_size_dir (directory, conn, logger):
    '''
    Description:
        Recursive function to get the size of the run directory on the 
        remote server.
        Optional can send an email to inform about the issue
    Input:
        logger # contains the logger object 
        conn # Connectio samba object
        directory   # root folder to start the checking file size
    Variables:
        file_list # contains the list of file and subfolders
        count_file_size # partial size for the subfolder
    Return:
        count_file_size # in the last iteraction will return the total
                    size of the folder
    '''
    count_file_size = 0
    file_list = conn.listPath(wetlab_config.SAMBA_SHARED_FOLDER_NAME, directory)
    for sh_file in file_list:
        if sh_file.isDirectory:
            if (sh_file.filename == '.' or sh_file.filename == '..'):
                continue
            logger.debug('Checking space for directory %s', sh_file.filename)
            sub_directory = os.path.join (directory,sh_file.filename)
            count_file_size += get_size_dir (sub_directory, conn, logger)
        else:
            count_file_size += sh_file.file_size

    return count_file_size


def get_run_disk_utilization (conn, run_Id_used, run_processing_id, logger):
        '''
    Description:
        Recursive function to get the size of the run directory on the 
        remote server.
        Optional can send an email to inform about the issue
    Input:
        logger # contains the logger object 
        conn # Connectio samba object
        run_Id_used   # root folder to start the checking file size
        run_processing_id #
    Functions:
        get_size_dir    # Located on this file
    Variables:
        file_list # contains the list of file and subfolders
        count_file_size # partial size for the subfolder
    Return:
        count_file_size # in the last iteraction will return the total
                    size of the folder
    '''
    logger.debug('Executing the function get_run_disk_utilization')
    if RunProcess.objects.filter(pk = run_processing_id).exists():
        run_be_updated = RunProcess.objects.get(pk = run_processing_id)
        get_full_list = conn.listPath(wetlab_config.SAMBA_SHARED_FOLDER_NAME ,run_Id_used)
        rest_of_dir_size = 0
        data_dir_size = 0
        images_dir_size = 0
        in_mega_bytes = 1024*1024
        logger.info('Starting getting disk space utilization for runID  %s', run_Id_used)
        for item_list in get_full_list:
            if item_list.filename == '.' or item_list.filename == '..':
                continue
            if item_list.filename == 'Data':
                dir_data = os.path.join(run_Id_used,'Data')
                data_dir_size = get_size_dir(dir_data , conn,logger)
                continue
            elif item_list.filename == 'Images':
                dir_images = os.path.join(run_Id_used, 'Images')
                images_dir_size = get_size_dir(dir_images , conn,logger)
                continue
            if item_list.isDirectory:
                item_dir = os.path.join(run_Id_used, item_list.filename)
                rest_of_dir_size += get_size_dir(item_dir, conn,logger)
            else:
                rest_of_dir_size += item_list.file_size
        # format file space and save it into database
        data_dir_size_formated = '{0:,}'.format(round(data_dir_size/in_mega_bytes))
        images_dir_size_formated = '{0:,}'.format(round(images_dir_size/in_mega_bytes))
        rest_of_dir_size_formated = '{0:,}'.format(round(rest_of_dir_size/in_mega_bytes))
        run_be_updated.useSpaceImgMb= images_dir_size_formated
        run_be_updated.useSpaceFastaMb= data_dir_size_formated
        run_be_updated.useSpaceOtherMb= rest_of_dir_size_formated
        run_be_updated.save()
        logger.info('End  disk space utilization for runID  %s', run_Id_used)
    logger.debug('Exiting the function get_run_disk_utilization')


'''
def save_miseq_run_info(run_info,run_parameter,run_id,logger):
## Collecting information from MiSeq run to save it in our database
    running_data={}
    image_channel=[]
    pir=1 #PlannedIndexNread:1 or 2
    pr=1 #PlannedNRead: 1 or 2
    #################################################
    ## parsing RunInfo.xml file
    #################################################
    run_data=ET.parse(run_info)
    run_root=run_data.getroot()
    logger.info('Processing the MISEQ runInfo.xml file')
    p_run=run_root[0]
    running_data['Flowcell']=p_run.find('Flowcell').text
    running_data['FlowcellLayout']=p_run.find('FlowcellLayout').attrib
    running_data['ImageChannel']= None #Available in NextSeq but not in MiSeq
    running_data['ImageDimensions']=None #Available in NextSeq but not in MiSeq

    #################################################
    ## parsing RunParameter.xml file
    #################################################
    logger.info('Processing the MISEQ runParameter.xml file')
    parameter_data=ET.parse(run_parameter)
    parameter_data_root=parameter_data.getroot()
    running_data['RunID']=parameter_data_root.find('RunID').text
    running_data['ExperimentName']=parameter_data_root.find('ExperimentName').text
    running_data['RTAVersion']=parameter_data_root.find('RTAVersion').text
    running_data['SystemSuiteVersion']=None #Available in NextSeq but not in MiSeq
    running_data['LibraryID']=None #Available in NextSeq but not in MiSeq
    running_data['Chemistry']=parameter_data_root.find('Chemistry').text
    running_data['RunStartDate']=parameter_data_root.find('RunStartDate').text
    running_data['AnalysisWorkflowType']=(parameter_data_root.find('Workflow')).find('Analysis').text
    logger.debug('running_data information -intermediate!- only'+ str(running_data))

    reads=parameter_data_root.find('Reads')
    #initialization of expected structure
    running_data['PlannedIndex1ReadCycles']=0
    running_data['PlannedIndex2ReadCycles']=0
    running_data['PlannedRead1Cycles']=0
    running_data['PlannedRead2Cycles']=0

    for run_info_read in reads.iter('RunInfoRead'):
        if 'Y'== run_info_read.attrib['IsIndexedRead']:
            if 1==pir:
                running_data['PlannedIndex1ReadCycles']=run_info_read.attrib['NumCycles']
                pir=2
            else:
                running_data['PlannedIndex2ReadCycles']=run_info_read.attrib['NumCycles']
        else:
            if 1==pr:
                running_data['PlannedRead1Cycles']=run_info_read.attrib['NumCycles']
                pr=2
            else:
                running_data['PlannedRead2Cycles']=run_info_read.attrib['NumCycles']


    running_data['RunManagementType']=parameter_data_root.find('RunManagementType').text
    running_data['ApplicationVersion']=parameter_data_root.find('Setup').find('ApplicationVersion').text
    running_data['NumTilesPerSwath']=parameter_data_root.find('Setup').find('NumTilesPerSwath').text

    logger.debug('running_data information for table RunParameters'+ str(running_data))
    ###########################################
    ## saving data into database
    ###########################################
    logger.info ('Saving to database  the (MiSeq) run id %s', run_id)
    running_parameters= RunningParameters (runName_id=RunProcess.objects.get(pk=run_id),
                         RunID=running_data['RunID'],
                         ExperimentName=running_data['ExperimentName'],
                         RTAVersion=running_data['RTAVersion'],
                         SystemSuiteVersion= running_data['SystemSuiteVersion'],
                         LibraryID= running_data['LibraryID'],
                         Chemistry= running_data['Chemistry'],
                         RunStartDate= running_data['RunStartDate'],
                         AnalysisWorkflowType= running_data['AnalysisWorkflowType'],
                         RunManagementType= running_data['RunManagementType'],
                         PlannedRead1Cycles= running_data['PlannedRead1Cycles'],
                         PlannedRead2Cycles= running_data['PlannedRead2Cycles'],
                         PlannedIndex1ReadCycles= running_data['PlannedIndex1ReadCycles'],
                         PlannedIndex2ReadCycles= running_data['PlannedIndex2ReadCycles'],
                         ApplicationVersion= running_data['ApplicationVersion'],
                         NumTilesPerSwath= running_data['NumTilesPerSwath'],
                         ImageChannel= running_data['ImageChannel'],
                         Flowcell= running_data['Flowcell'],
                         ImageDimensions= running_data['ImageDimensions'],
                         FlowcellLayout= running_data['FlowcellLayout'])

    running_parameters.save()
    ##############################################
    ## updating the date fetched from the Date tag for run and project
    ##############################################
    date = p_run.find('Date').text
    logger.debug('Found the date that was recorded the Run %s', date)
    run_date = datetime.datetime.strptime(date, '%y%m%d')

    run_to_be_updated = RunProcess.objects.get(pk=run_id)
    run_to_be_updated.run_date = run_date
    logger.debug('run_to_be_updated: '+run_to_be_updated.runName)
    logger.debug('run_date: '+str(run_to_be_updated.run_date))
    run_to_be_updated.save()
    logger.info('Updated the run date for the runProcess table ')

    projects_to_update = Projects.objects.filter(runprocess_id__exact = run_id)
    for project in projects_to_update :
        project.project_run_date = run_date
        project.save()
        logger.debug('project_to_be_updated: '+project.projectName)
        logger.debug('project_run_date: '+str(project.project_run_date))
        logger.info('Updated the project date for the Project table ')

    return
'''



def parsing_and_save_run_info(run_info, run_parameter, run_id, logger):
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
        logger     # contains the logger object to write information
                    on the log file
    Variables:
        image_channel   # list containing the image channel values
        p_parameter     # element tree object for run parameter
        p_run           # element tree for run info
        projects_to_update # list of projects that belong to the same run
        running_data    # containg the parsed information to save into database
        running_parameters # new RunningParameter object to store parsed data
        sequencer       # sequencer index of the machine
        
     Return:
        No return values. the function store all the information in database
    '''
    logger.debug('Executing the parsing_and_save_run_info function')
    running_data={}
    image_channel=[]
    #################################################
    ## parsing RunInfo.xml file
    #################################################
    run_data=ET.parse(run_info)
    run_root=run_data.getroot()
    logger.info('parsing the runInfo.xml file for %s' , run_id)
    p_run=run_root[0]
    ## getting the common values NextSeq and MiSeq
    running_data['Flowcell']=p_run.find('Flowcell').text
    running_data['FlowcellLayout']=p_run.find('FlowcellLayout').attrib
    
    #################################################
    ## parsing RunParameter.xml file
    #################################################
    logger.info('Parsing the runParameter.xml file for %s' , run_id)
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
    running_data['ApplicationVersion']=p_parameter.find('ApplicationVersion').text
    running_data['NumTilesPerSwath']=p_parameter.find('NumTilesPerSwath').text
    ## fetch the application platform to get the machine used in the run
    if 'NextSeq' in p_parameter.find('ApplicationName').text :
        # getting the specific values for NextSeq
        running_data['SystemSuiteVersion']=parameter_data_root.find('SystemSuiteVersion').text ##NextSeq
        running_data['LibraryID']=parameter_data_root.find('LibraryID').text  ##NextSeq
        running_data['AnalysisWorkflowType']=parameter_data_root.find('AnalysisWorkflowType').text  ##NextSeq
        running_data['PlannedRead1Cycles']=parameter_data_root.find('PlannedRead1Cycles').text  ##NextSeq
        running_data['PlannedRead2Cycles']=parameter_data_root.find('PlannedRead2Cycles').text  ##NextSeq
        running_data['PlannedIndex1ReadCycles']=parameter_data_root.find('PlannedIndex1ReadCycles').text  ##NextSeq
        running_data['PlannedIndex2ReadCycles']=parameter_data_root.find('PlannedIndex2ReadCycles').text  ##NextSeq
        
        ## get the values on RunInfo that are only for NextSeq
        for i in run_root.iter('Name'):
            image_channel.append(i.text)
        running_data['ImageChannel']=image_channel  ##NextSeq
        running_data['ImageDimensions']=p_run.find('ImageDimensions').attrib  ##NextSeq
        
        ## get the instrument for NextSeq run
        instrument = parameter_data_root.find('InstrumentID').text  ##NextSeq
        
    else:
        ## get the length index number for reads and indexes 
        for run_info_read in parameter_data_root.iter('RunInfoRead'):
            if run_info_read.attrib['Number'] == 1 :
                running_data['PlannedRead1Cycles'] = run_info_read.attrib['NumCycles']
            elif run_info_read.attrib['Number'] == 2 :
                running_data['PlannedIndex1ReadCycles'] = run_info_read.attrib['NumCycles']
            elif run_info_read.attrib['Number'] == 3 :
                running_data['PlannedIndex2ReadCycles'] = run_info_read.attrib['NumCycles']
            elif run_info_read.attrib['Number'] == 4 :
                running_data['PlannedRead2Cycles'] = run_info_read.attrib['NumCycles']

        ## setting empty values for MiSeq
        running_data['SystemSuiteVersion'] = ''
        running_data['LibraryID'] = ''
        running_data['AnalysisWorkflowType'] = ''
        running_data['ImageChannel'] = ''
        running_data['ImageDimensions'] = ''
        ## get the instrument for MiSeq run located on RunInfo.xml
        instrument =p_run.find('Instrument').text
    '''
    Chemistry = parameter_data_root.find('Chemistry').text
    ## check if InstrumentID exists on database
    if not Machines.objects.filter(machineName__exact = InstrumentID).exists() :
        if not Platform.objects.filter(platformName__exact = Chemistry).exists():
            new_platform = Platform (platformName = Chemistry)
            new_platform.save()
        platform = Platform.objects.get(platformName__exact = Chemistry)
        new_machine = Machines(platformID= platform, machineName = InstrumentID )
        new_machine.save()
    '''
    

    logger.debug('running_data information', running_data)
    ###########################################
    ## saving data into database
    ###########################################
    logger.info ('Saving to database  the running_parameters for %s', run_id)

    running_parameters= RunningParameters (runName_id=RunProcess.objects.get(pk=run_id),
                         RunID=running_data['RunID'], ExperimentName=running_data['ExperimentName'],
                         RTAVersion=running_data['RTAVersion'], SystemSuiteVersion= running_data['SystemSuiteVersion'],
                         LibraryID= running_data['LibraryID'], Chemistry= running_data['Chemistry'],
                         RunStartDate= running_data['RunStartDate'], AnalysisWorkflowType= running_data['AnalysisWorkflowType'],
                         RunManagementType= running_data['RunManagementType'], PlannedRead1Cycles= running_data['PlannedRead1Cycles'],
                         PlannedRead2Cycles= running_data['PlannedRead2Cycles'], PlannedIndex1ReadCycles= running_data['PlannedIndex1ReadCycles'],
                         PlannedIndex2ReadCycles= running_data['PlannedIndex2ReadCycles'], ApplicationVersion= running_data['ApplicationVersion'],
                         NumTilesPerSwath= running_data['NumTilesPerSwath'], ImageChannel= running_data['ImageChannel'],
                         Flowcell= running_data['Flowcell'], ImageDimensions= running_data['ImageDimensions'],
                         FlowcellLayout= running_data['FlowcellLayout'])

    running_parameters.save()
    ##############################################
    ## updating the date fetched from the Date tag for run and project
    ##############################################
    date = p_run.find('Date').text
    logger.debug('Found the de date that was recorded the Run %s', date)
    run_date = datetime.datetime.strptime(date, '%y%m%d')

    run_to_be_updated = RunProcess.objects.get(pk=run_id)
    run_to_be_updated.run_date = run_date
    if  Machines.objects.filter(machineName__exact = instrument).exists() :
        ## Update the run information with the instrument used
        sequencer = Machines.objects.get(machineName__exact = instrument)
        run_to_be_updated.sequencerModel = sequencer
    else:
        logger.error('-----------------    ERROR   ------------------')
        logger.error('%s has been not defined on machines ', instrument)
        logger.error('-----------------    ERROR   ------------------')
    run_to_be_updated.save()
    logger.info('Updated the run date and sequencer used for the runProcess table ')

    projects_to_update = Projects.objects.filter(runprocess_id__exact = run_id)
    for project in projects_to_update :
        project.project_run_date = run_date
        project.save()
        logger.info('Updated the project date for the Project table ')
    
    logger.debug('Exiting the parsing_and_save_run_info function')
    return ''



def save_new_miseq_run (sample_sheet, logger) :
    '''
    Description:
        The function will get the sample sheet information and will 
        save the information in database by creating a new miseq run and
        new projects.
        The function will move the sample sheet from the temporary folder
        to the designed folder for all sample sheets.
    Input:
        sample_sheet  # full path for smaple sheet file
        logger # log object for logging 
    Functions:
        get_experiment_library_name # located at this file
        get_projects_in_run # located at this file
        logger_errors   # locate at utils.wetlab_misc_utilities
    Imports:
        Center # from django_utils.models
        datetime # to fetch present time in miliseconds
        os      # file utilities
    Constants:
        MEDIA_ROOT  
        RUN_SAMPLE_SHEET_DIRECTORY
        DEFAULT_CENTER
        DEFAULT_LIBRARY_KIT
    Variables:
        base_space_file # path formed with django settings.MEDIA_URL and
                    sample_sheet_on_database
        experiment_name # 
        center_requested_by # Center object for run requested center
        library_name # 
        projects_users # dictionary contains projects and users
        new_sample_sheet_file # sample sheet path on the final destination
        new_sample_sheet_name # sample sheet name including timestamp
        sample_sheet  # full path for storing sample sheet file on 
                        tempary local folder new_sample_sheet_name
        sample_sheet_on_database # contains the folder path where all
                        are stored and the 
        s_sample_sheet  # full path for remote sample sheet file
        timestr     # present time for adding to sample sheet name
    '''
    logger.debug('Executing the function save_new_miseq_run' )
    experiment_name, library_name = get_experiment_library_name (sample_sheet)
    projects_users = get_projects_in_run(sample_sheet)
    
    now = datetime.datetime.now()
    timestr = now.strftime("%Y%m%d-%H%M%S.%f")[:-3]
    new_sample_sheet_name = 'SampleSheet' + timestr + '.csv'
    new_sample_sheet_file = os.path.join (settings.MEDIA_ROOT, wetlab_config.RUN_SAMPLE_SHEET_DIRECTORY, new_sample_sheet_name)
    logger.debug('new sample sheet name %s', new_sample_sheet_file)
    sample_sheet_on_database = os.path.join(wetlab_config.RUN_SAMPLE_SHEET_DIRECTORY, new_sample_sheet_name)
    ## move sample sheet to final folder
    os.rename(sample_sheet, new_sample_sheet_file)
    #import pdb; pdb.set_trace()
    ## create a new entry on runProcess database
    if Center.objects.filter(centerAbbr__exact = wetlab_config.DEFAULT_CENTER).exists():
        center_requested_by = Center.objects.get(centerAbbr__exact = wetlab_config.DEFAULT_CENTER)
    else:
        string_message = 'The requested center ' +  wetlab_config.DEFAULT_CENTER + ' defined in config wetlab file does not exist'
        logger_errors(logger, string_message)
        logger.info('Using the first center defined in database')
        center_requested_by = Center.objects.all().first()
        
    run_proc_data = RunProcess(runName=experiment_name,sampleSheet= sample_sheet_on_database, 
                            runState='Recorded', centerRequestedBy = center_requested_by)
    run_proc_data.save()
    logger.info('Updated runProccess table with the new experiment name found')

    if LibraryKit.objects.filter(libraryName__exact = wetlab_config.DEFAULT_LIBRARY_KIT).exists():
        library_kit = LibraryKit.objects.get(libraryName__exact = wetlab_config.DEFAULT_LIBRARY_KIT)
    else:
        string_message = 'The default library ' +  wetlab_config.DEFAULT_LIBRARY_KIT + ' defined in config wetlab file does not exist'
        logger_errors(logger, string_message)
        logger.info('Using the first library kit defined in database')
        library_kit = LibraryKit.objects.all().first()
    base_space_file = os.path.join(settings.MEDIA_URL.replace('/',''), sample_sheet_on_database)

    for project, user  in projects_users.items():
        userid=User.objects.get(username__exact = user)
        p_data=Projects(runprocess_id=RunProcess.objects.get(runName =experiment_name), 
                        projectName=project, user_id=userid, procState ='Recorded',
                        baseSpaceFile = base_space_file, 
                        LibraryKit_id = library_kit, libraryKit = library_name)
        p_data.save()
    logger.info('Updated Projects table with the new projects found')
    
    logger.debug('Exiting the function save_new_miseq_run' )
    return False

def validate_sample_sheet (sample_sheet, logger):
    '''
    Description:
        The function get the sample sheet file and will make some chekings
        to validate the file.
    Input:
        sample_sheet  # full path for smaple sheet file
        logger          # log object 
    Functions:
        get_experiment_library_name # located at this file
        get_projects_in_run # located at this file
        logger_errors   # locate at utils.wetlab_misc_utilities
    Variable:
        projects_users  # dictionary containing projects and user owner
        experiment_name # contains the experiment name from the sample sheet
        library_name  # contains the library name from the sample sheet
    Return:
        True if all checking are successful False if a check fails
    '''
    logger.debug('Executing the function validate_sample_sheet' ) 
    # get experiment name
    experiment_name, library_name = get_experiment_library_name (sample_sheet)
    # check if experiment name is already used in the system
    if experiment_name != '' :
        if RunProcess.objects.filter(runName__exact = experiment_name).exists() :
            logger.warning('Experiment name %s , already been used', experiment_name)
            logger.debug('Exiting the function validate sample_sheet')
            return False
    else:
        string_message ='Experiment name is empty'
        logging_errors(logger, string_message, False)
        logger.debug('Exiting the function validate sample_sheet')
        return False
            
    # get the projects and user owner form sample sheet
    projects_users = get_projects_in_run(sample_sheet)
    if len(projects_users) == 0 :
        logging_errors(logger, 'No projects have been found ')
        logger.debug('Exiting the function validate sample_sheet')
        return False
    for project in projects_users.keys() :
        if Projects.objects.filter(projectName__exact = project).exists():
            string_message = 'project name %s , already been used ' + project
            logging_errors(logger, string_message)
            logger.debug('Exiting the function validate sample_sheet')
            return False
    for user in projects_users.values():
        if ( not User.objects.filter(username__icontains = user).exists()):
            string_message = 'user name ' +user  +' is not defined in the system'
            logging_errors(logger, string_message)
            logger.debug('Exiting the function validate sample_sheet')
            return False
        
    logger.debug('Exiting the function validate_sample_sheet' ) 
    return True

def get_new_runs_on_remote_server (processed_runs, conn, shared_folder, logger):
    '''
    Description:
        The function fetch the folder names from the remote server and 
        returns a list containing the folder names that have not been
        porcessed yet.
    Input:
        processed_runs  # full path and name of the file
        conn # samba connection object
        shared_folder   # shared folder in the remote server
        logger          # log object 
    Variable:
        new_runs        # list containing the new folder run names
        run_folder_list  # list of the folder names on remote server
    Return:
        new runs that have been found
    '''
    logger.debug('Executing the function get_new_runs_on_remote_server' ) 
    new_runs = []
    run_folder_list = conn.listPath( shared_folder, '/')
    for sfh in run_folder_list:
        if sfh.isDirectory:
            folder_run = sfh.filename
            if (folder_run == '.' or folder_run == '..'):
                continue
            # if the run folder has been already process continue searching
            if folder_run in processed_runs:
                logger.debug('folder run  %s already processed', folder_run)
                continue
            else:
                logger.info ('Found a new run  %s ',folder_run)
                new_runs.append(folder_run)
    logger.debug('Exiting the function get_new_runs_on_remote_server' )
    return new_runs


def read_processed_runs_file (processed_run_file, logger) :
    '''
    Description:
        The function reads the file that contains all the processed runs
        and return a list with the run folder names 
    Input:
        processed_run_file # full path and name of the file
        logger # log object 
    Variable:
        processed_runs  # list of the folder names read from file
    Return:
        Error when file can not be readed
        processed_runs variable for successful file read
    '''
    logger.debug('Executing the function read_processed_runs_file' )
    w_dir = os.getcwd()
    logger.debug('Folder for fetching the proccess run file is %s', w_dir)
    logger.debug('processed_run_file is %s', processed_run_file)
    processed_runs = []
    if os.path.exists(processed_run_file):
        try:
            fh = open (processed_run_file,'r')
            for line in fh:
                line=line.rstrip()
                processed_runs.append(line)
        except Exception as e:
            string_message = 'Unable to open the processed run file.  %s '
            logging_errors(logger, string_message, True)
            return 'Error'
        fh.close()
        logger.info('run processed file have been read')
        logger.debug('Exiting sucessfully the function read_processed_runs_file' )
        return processed_runs
    else:
        logger.debug('No found processed run file. Exiting the function' )
        return 'Error'
    

def search_new_miseq_runs (logger):
    '''
    Description:
        The function will check if there are new run folder in the remote
        server has a sample sheet file that was not processed. 
        When found it will perform checks on the file to validate it and
        store a new run in database
    Input:
        logger # log object for logging 
    Functions:
        open_samba_connection # located in utils.wetlab_misc_utilities.py 
        get_new_runs_on_remote_server # located at this file
        validate_sample_sheet   # located at this file
        save_new_miseq_run # located at this file
    Constants:
        PROCESSED_RUN_FILE
        RUN_TEMP_DIRECTORY
        SAMBA_SHARED_FOLDER_NAME
        SAMPLE_SHEET
        
    Variables:
        handle_new_miseq_run # list with all new miseq runs
        l_sample_sheet  # full path for storing sample sheet file on 
                        tempary local folder
        s_sample_sheet  # full path for remote sample sheet file
        processed_run_file # path
    '''
    handle_new_miseq_run = []
    processed_run_file = os.path.join( wetlab_config.RUN_TEMP_DIRECTORY, wetlab_config.PROCESSED_RUN_FILE)
    processed_runs = read_processed_runs_file (processed_run_file, logger)
    if processed_runs == 'Error' :
        return 'Error'
    try:
        conn=open_samba_connection()
        logger.info('Sucessfully  SAMBA connection for the process_run_in_recorded_state')
    except Exception as e:
        string_message = 'Unable to open SAMBA connection for the process_run_in_recorded_state %s'
        logging_errors(logger, string_message, True)
        return 'Error'

    new_runs = get_new_runs_on_remote_server (processed_runs, conn, 
                        wetlab_config.SAMBA_SHARED_FOLDER_NAME, logger)
    if len (new_runs) > 0 :
        for new_run in new_runs :
            l_sample_sheet = os.path.join(wetlab_config.RUN_TEMP_DIRECTORY, wetlab_config.SAMPLE_SHEET)
            s_sample_sheet = os.path.join(new_run, wetlab_config.SAMPLE_SHEET)
            with open(l_sample_sheet ,'wb') as s_sheet_fp :
                try:
                    conn.retrieveFile(wetlab_config.SAMBA_SHARED_FOLDER_NAME, s_sample_sheet, s_sheet_fp)
                    logger.info('Retrieving the remote %s/SampleSheet.csv', new_run)
                except Exception as e:
                    string_message = 'Unable to fetch ' + new_run + '/SampleSheet.csv'
                    logging_errors(logger, string_message, True)
                    continue
            
            if validate_sample_sheet (l_sample_sheet, logger) :
                logger.info('Successful sample sheet checking for folder %s ', new_run)
                save_new_miseq_run(l_sample_sheet, logger)
                handle_new_miseq_run.append(new_run)
            else:
                logger.info('Invalid sample sheet checking for folder %s ', new_run)
            
    else:
        logger.info('No found new run folders on the remote server')
        return ''
    
    return handle_new_miseq_run



def find_xml_tag_text (input_file, search_tag):
    '''
    Description:
        The function will look for the xml element tag in the 
        file and it will return the text value
    Input:
        input_file  # file to find the tag
        search_tag  # xml tag to be found in the input_file
    Variables:
        found_tag   # line containing the tag 
    '''
    fh = open (input_file, 'r')
    search_line = '<' + search_tag+ '>(.*)</' + search_tag+'>'
    for line in fh:
        found_tag = re.search('^\s+ %s' % search_line, line)
        if found_tag:
            fh.close()
            return found_tag.group(1)
    fh.close()
    return ''


def check_miseq_completion_run (conn, log_folder, l_run_parameter, l_logger):
    '''
    Description:
        The function will check the run log files to check if all cycles
        were running on the sequencer and if the last log does not
        mention that run was canceled
    Input:
        conn    # contains the samba connection object
        log_folder  # remote path where are located the logs for miseq
        logger     # contains the logger object to write information
                    on the log file     
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
    run_completion_date =''
    logger.debug('Executing the function check_miseq_completion_run')
    run_cycles = get_miseq_run_cycles (l_run_parameter)
    log_cycles, latest_log = get_latest_miseq_log(conn, log_folder)
    if 'Canceled' in latest_log :
        status_run = 'Canceled'
    elif log_cycles != run_cycles :
        status_run = 'still_processing'
    else:
        status_run = wetlab_config.COMPLETION_SUCCESS
        last_line_in_file = latest_log.split('/n')[-2]
        last_log_time = last_line_in_file.split(' ')[0:2]
        last_log_time[0] = str('20'+ last_log_time[0])
        run_completion_date = ''.join.(last_log_time)
    logger.debug('Exiting the function check_miseq_completion_run')
    return status_run, run_completion_date

def handle_runs_in_recorded_state(logger):
    '''
    Description:
        The function is called from crontab to check if NextSeq runs 
        are in recorded state.
        web, having 2 main parts:
            - User form with the information to add a new library
            - Result information as response of user submit
        
        Save a new library kit name in database if it is not already defined.
    Input:
        logger     # contains the logger object to write information
                    on the log file     
    Constants:
         
        RUN_INFO
        RUN_COMPLETION
        RUN_PARAMETER
        
        RUN_PARAMETER_NEXTSEQ
        RUN_PARAMETER_MISEQ
        RUN_TEMP_DIRECTORY
    Functions:
        find_xml_tag_text # located at this file
        read_processed_runs_file 
        parsing_and_save_run_info   # located at this file
        update_run_state # located at this file
        update_project_state # located at this file
    Variables:
        base_directory
        exp_name        # contains the text value for the experimentName
                        tag on runParameter.xml file
        new_runs # contains then new runs located on the remote server
        
        l_run_parameter # conatins the path for run parameter file
        l_run_info      # contains the path for run info file
        l_run_completion# contains the path for run completion status file
        processed_runs  # contain the list of the processed runs
        process_run_file_updated # boolean variable. False if there is not
                        need to update the process_run_file
        run_found       # run object that contains the run to be updated
        run_found_id    # Id number of the run object to be updated
        project_name_list # list of project objects that are in the run_found
        s_run_completion# contains the remote path for run completion status file 
        s_run_parameter # contains the remote path for run parameter file 

        s_sample        # contains the remote path for sample sheet 
        run_completion_attributes # contains the remote file attribute for 
                            completion_status
        run_completion_date   # contains the remote file date for 
                            completion_status
        run_found       # contains the object RunProccess for the found run
        run_found_id    # contains the id of the found run
        run_platform    # contains the platform name get from the application
                        name tag in the RunParameters.xml file
        recorded_dir
        status_run      # contains the text value for the competion status
                        tag on runCompletion.xml file
        sample_sheet_tmp_file # contains the path where is located the 
                        samplesheet file inside the temporary folder
    Return:
        return a list containing the experiment names that have been 
        sucessfully proccessed
    '''
    processed_run_file = []
    #run_list = [] 
    try:
        conn=open_samba_connection()
        logger.info('Sucessfully  SAMBA connection for the process_run_in_recorded_state')
    except:
        logger.error('Unable to open SAMBA connection for the process_run_in_recorded_state')
        return 'Error'

    l_run_parameter = os.path.join(wetlab_config.RUN_TEMP_DIRECTORY, wetlab_config.RUN_PARAMETER_NEXTSEQ)
    l_run_info = os.path.join(wetlab_config.RUN_TEMP_DIRECTORY, wetlab_config.RUN_INFO)
    l_run_completion = os.path.join(wetlab_config.RUN_TEMP_DIRECTORY, wetlab_config.RUN_COMPLETION)
    process_run_file = os.path.join(wetlab_config.RUN_TEMP_DIRECTORY, wetlab_config.PROCESSED_RUN_FILE)

    processed_runs=[]
    run_names_processed=[]
    ## get the list of the processed run

    processed_runs = read_processed_runs_file (processed_run_file, logger)
    if processed_runs == 'Error' :
        return 'Error'
    logger.info('processed_runs file was read')
    # Initialize the variable do not need to update the process run file.
    process_run_file_updated = False
    new_runs = get_new_runs_on_remote_server (processed_runs, conn, 
                        wetlab_config.SAMBA_SHARED_FOLDER_NAME, logger)
    if len (new_runs) > 0 :
        for new_run in new_runs :
            with open(l_run_parameter ,'wb') as r_par_fp :
                try: # get RunParameter.xml if NextSeq
                    s_run_parameter = os.path.join(run_dir,wetlab_config.RUN_PARAMETER_NEXTSEQ)
                    conn.retrieveFile(wetlab_config.SAMBA_SHARED_FOLDER_NAME, s_run_parameter, r_par_fp)
                    logger.info('Retrieving the remote RunParameter.xml file for %s', new_run)
                except:
                    try: # get runParameter.xml if miSeq
                        s_run_parameter = os.path.join(run_dir,wetlab_config.RUN_PARAMETER_MISEQ)
                        conn.retrieveFile(wetlab_config.SAMBA_SHARED_FOLDER_NAME, s_run_parameter, r_par_fp)
                        logger.info('Retrieving the remote runParameter.xml file for %s', new_run)

                    except Exception as e:
                        string_message = 'Unable to fetch the RunParameter.xml file on folder ' + new_run
                        logging_errors (logger, string_message, True)
                        os.remove(l_run_parameter)
                        continue
            # look for the experience name  for the new run folder. Then find the run_id valued for it
            exp_name = find_xml_tag_text (l_run_parameter, wetlab_config.EXPERIMENT_NAME_TAG)
            logger.debug('found the experiment name  , %s', exp_name)
            if exp_name == '':
                string_message = 'NO experiment name  was defined for run ' + new_run
                logging_errors (logger, string_message, False)
                continue
            run_platform = find_xml_tag_text (l_run_parameter, wetlab_config.APPLICATION_NAME_TAG)
            if 'MiSeq' in run_platform :
                logger.info('Found out that it is a miSeq run')
                log_folder = 
                status_run, run_completion_date = check_miseq_completion_run (conn, log_folder, l_run_parameter, logger)
                if 'Error' == status_run :
                    string_message = 'miseq run was canceled for the run ' + new_run
                    logging_errors (logger, string_message, False)
                    # lets continue to set the run and their projects to CANCELED state
                elif 'still_processing' == status_run :
                    logger.info('miSeq run process is still running on sequencer for %s', new_run)
                    logger.debug('Deleting local copy of the run parameter file')
                    os.remove(l_run_parameter)
                    continue
            else:
                # check if NextSEq run have been successful completed
                s_run_completion = os.path.join(new_run, wetlab_config.RUN_COMPLETION)
                try:
                    with open (l_run_completion, 'wb') as c_status_fp :
                        conn.retrieveFile(wetlab_config.SAMBA_SHARED_FOLDER_NAME, s_run_completion, c_status_fp )
                        # Get the date and time when the RunCompletionStatus is created
                        run_completion_attributes = conn.getAttributes(wetlab_config.SAMBA_SHARED_FOLDER_NAME , s_run_completion)
                        run_completion_date = datetime.datetime.fromtimestamp(int(run_completion_attributes.create_time)).strftime('%Y-%m-%d %H:%M:%S')
                except Exception as e:
                    string_message = 'Unable to fetch the completion status file on folder ' + new_run
                    logging_errors (logger, string_message, True)
                    logger.info ('Deleting local copy of completion status ')
                    os.remove(l_run_completion)
                    continue
                status_run = find_xml_tag_text (l_run_completion, wetlab_config.COMPLETION_TAG )
                if  status_run != wetlab_config.COMPLETION_SUCCESS:
                    string_message = 'Run status was ' + status_run  
                    logging_errors (logger, string_message, False)
                    # lets continue to set the run to CANCELED state
                else:
                    logger.info ('Run completed for Run ID %s ', new_run)

                logger.debug('Deleting RunCompletionStatus.xml file')
                os.remove(l_run_completion)

            # get the RunProcess object
            if  RunProcess.objects.filter(runName__icontains = exp_name, runState__exact = 'Recorded').exists():
                run_found_id=str(RunProcess.objects.get(runName__exact = exp_name).id)
                logger.debug('Matching the experimental name id %s with database ', exp_name_id)
                if status_run != wetlab_config.COMPLETION_SUCCESS:
                    # set the run in error state
                    run_found = RunProcess.objects.get(runName__exact = exp_name)
                    run_found.runState = 'CANCELLED'
                    run_found.run_finish_date = run_completion_date
                    run_found.save()
                    logger.warning('-----------------    WARNING   ------------------')
                    logger.warning('run was canceled for  %s', new_run)
                    logger.warning('-----------------    END WARNING   --------------')
                    project_name_list = Projects.objects.filter(runprocess_id__exact = run_found_id)
                    for project in project_name_list:
                        project.procState= 'CANCELLED'
                        project.save()
                    # delelete the runParameter file
                    os.remove(l_run_parameter)
                    # Updated the processed file  and continue with the next item in the list
                    process_run_file_updated = True
                    processed_runs.append(new_run)
                    continue

                if 'NextSeq' in run_platform :
                    sample_sheet_tmp_dir = os.path.join(recorded_dir,exp_name_id)
                    sample_sheet_tmp_file=os.path.join(sample_sheet_tmp_dir, wetlab_config.SAMPLE_SHEET)
                    if os.path.exists(sample_sheet_tmp_file):
                        if wetlab_config.COPY_SAMPLE_SHEET_TO_REMOTE :
                            # copy Sample heet file to remote directory
                            logger.info('Copy sample sheet to remote folder %s', new_run)
                            try:
                                with open(sample_sheet_tmp_file ,'rb') as  sample_samba_fp:
                                    s_sample= os.path.join(new_run, wetlab_config.SAMPLE_SHEET)
                                    logger.debug('Local dir for Shample Sheet %s', sample_sheet_tmp_file)
                                    logger.debug('Remote dir to copy Shample Sheet  is %s', new_run)
                                    conn.storeFile(wetlab_config.SAMBA_SHARED_FOLDER_NAME, s_sample, sample_samba_fp)
                                    logger.info('Samplesheet.csv file has been copied on the remote server')
                            except Exception as e:
                                string_message = 'Unable to copy Sample Sheet for run ' + new_run
                                logging_errors (logger, string_message, True)
                                print ('ERROR:: Unable to copy Sample Sheet for ', new_run)
                                print ('ERROR:: check log file for more information')
                                continue
                        try:
                            os.shutil(sample_sheet_tmp_dir)
                            logger.debug('Deleted temporary folder containing the samplesheet')
                        except Exception as e:
                            string_message = 'Unable to delete temporary folder with the sample sheet ' + sample_sheet_tmp_dir
                            logging_errors (logger, string_message, True)
                    else:
                        string_message = 'sample sheet not found on local Directory ' + sample_sheet_tmp_dir
                        logging_errors (logger, string_message, True)

                # get the runIfnfo.xml to collect the  information for this run
                try:
                    with open(l_run_info ,'wb') as r_info_fp :
                        samba_run_info_file=os.path.join(new_run,wetlab_config.RUN_INFO)
                        logger.debug('Local dir for RunInfo.xml %s', l_run_info)
                        logger.debug('Remote file of RunInfo.xml is located in %s', samba_run_info_file)
                        conn.retrieveFile(wetlab_config.SAMBA_SHARED_FOLDER_NAME, samba_run_info_file, r_info_fp)
                except Exception as e:
                    string_message = 'Unable to get the RunInfo.xml file for run ' + new_run
                    logging_errors (logger, string_message, True)
                    logger.info('Deleting RunParameter and RunInfo files ')
                    os.remove(l_run_parameter)
                    os.remove(l_run_info)
                    print ('ERROR:: Unable to get the RunInfo.xml filet for run ', new_run)
                    print ('Check log file for more information')
                    continue
                logger.info('copied to local the RunInfo.xml and start the parsing for RunInfo and RunParameter files')
                parsing_and_save_run_info (l_run_info, l_run_parameter, exp_name_id, logger)
                logger.info('Run information parsed and saved for %s' , exp_name)
                    # delete the copy of the run files
                os.remove(l_run_info)
                os.remove(l_run_parameter)
                logger.debug('Deleted runInfo and RunParameter files on local server')
                    # delete the file and folder for the sample sheet processed
                if 'NextSeq' in run_platform :
                    # for nextSeq run delete the temporary folder which contains the sample sheet file 
                    shutil.rmtree(os.path.join(recorded_dir, exp_name_id))
                    logger.debug('Deleted the recorded folder %s', exp_name_id)
                # change the run  to SampleSent state
                update_run_state(exp_name_id, 'Sample Sent', logger)
                update_project_state(exp_name_id, 'Sample Sent', logger)
                    # add the completion date in the run
                logger.info('Saving completion date for %s' , exp_name)
                run_update_date = RunProcess.objects.get(pk=exp_name_id)
                run_update_date.run_finish_date = run_completion_date
                run_update_date.save()
                    # add the run_dir inside the processed_runs file
                processed_runs.append(new_run)
                run_names_processed.append(exp_name)
                    # mark to update the processed_runs file with the new run folders
                process_run_file_updated = True
                logger.info('Finish the proccess for run folder %s', new_run)
            else:
                logger.warning('-----------------    WARNING   ------------------')
                logger.warning ('The run ID %s does not match any run in the RunProcess object.', new_run)
                logger.warning('-----------------   END  WARNING   ------------------')
                os.remove(l_run_parameter)
                continue
    conn.close()
    if process_run_file_update :
        fh =open (process_run_file,'w')
        # update the process_run_file with the new run
        for process in processed_runs:
            fh.write(process)
            fh.write('\n')
        fh.close()
    # check if all run process file are handled

    list_dir=os.listdir(recorded_dir)
    if list_dir:
        logger.warning('-----------------    WARNING   ------------------')
        logger.warning('There are runs in Recorded state that were not matching on the remote server')
        logger.debug('pending folder to be processed %s', list_dir)
        logger.warning('-----------------   END  WARNING   ------------------')
    return(run_names_processed)


def update_run_state(run_id, state, logger):
    run=RunProcess.objects.get(pk=run_id)
    logger.info('updating the run state for %s to %s ', run_id, state)
    run.runState= state
    run.save()

def update_project_state(run_id, state, logger):
    projects_to_update = Projects.objects.filter(runprocess_id__exact = run_id)
    for project in projects_to_update :
        logger.info('updating the project state for %s to %s ', project, state)
        project.procState = state
        project.save()

def parsing_statistics_xml(run_id, demux_file, conversion_file, logger):
    total_p_b_count=[0,0,0,0]
    stats_result={}

    demux_stat=ET.parse(demux_file)
    root=demux_stat.getroot()
    projects=[]
    logger.info('Starting conversion for demux file')
    for child in root.iter('Project'):
        projects.append(child.attrib['name'])
    total_samples = 0
    number_of_lanes = get_machine_lanes(run_id)
    for i in range(len(projects)):
        p_temp=root[0][i]
        samples=p_temp.findall('Sample')

        sample_all_index=len(samples)-1
        barcodeCount ,perfectBarcodeCount, b_count =[], [] ,[]
        p_b_count, one_mismatch_count =[], []

        dict_stats={}
        for c in p_temp[sample_all_index].iter('BarcodeCount'):
        #for c in p_temp[sample].iter('BarcodeCount'):
            #b_count.append(c.text)
            barcodeCount.append(c.text)
        for c in p_temp[sample_all_index].iter('PerfectBarcodeCount'):
            p_b_count.append(c.text)

        # look for One mismatch barcode

        if p_temp[sample_all_index].find('OneMismatchBarcodeCount') ==None:
             for  fill in range(number_of_lanes):
                one_mismatch_count.append('NaN')
        else:
            for c in p_temp[sample_all_index].iter('OneMismatchBarcodeCount'):
                one_mismatch_count.append(c.text)

        #one_mismatch_count.append(one_m_count)

        dict_stats['BarcodeCount']=barcodeCount
        dict_stats['PerfectBarcodeCount']=p_b_count
        dict_stats['sampleNumber']=len(samples) -1
        dict_stats['OneMismatchBarcodeCount']=one_mismatch_count
        stats_result[projects[i]]=dict_stats
        if projects[i] != 'default' and projects[i] != 'all':
            total_samples += len(samples) -1
        logger.info('Complete parsing from demux file for project %s', projects[i])
    # overwrite the value for total samples
    stats_result['all']['sampleNumber']=total_samples

    conversion_stat=ET.parse(conversion_file)
    root_conv=conversion_stat.getroot()
    projects=[]
    logger.info('Starting conversion for conversion file')
    for child in root_conv.iter('Project'):
        projects.append(child.attrib['name'])
    for i in range(len(projects)):
        p_temp=root_conv[0][i]
        samples=p_temp.findall('Sample')
        sample_all_index=len(samples)-1
        tiles=p_temp[sample_all_index][0][0].findall('Tile')
        tiles_index=len(tiles)-1
        list_raw_yield=[]
        list_raw_yield_q30=[]
        list_raw_qualityscore=[]
        list_pf_yield=[]
        list_pf_yield_q30=[]
        list_pf_qualityscore=[]

        for l_index in range(number_of_lanes):
            raw_yield_value = 0
            raw_yield_q30_value = 0
            raw_quality_value = 0
            pf_yield_value = 0
            pf_yield_q30_value = 0
            pf_quality_value = 0
            for t_index in range(tiles_index):

                     # get the yield value for RAW and for read 1 and 2
                for c in p_temp[sample_all_index][0][l_index][t_index][0].iter('Yield'):
                    raw_yield_value +=int(c.text)
                    # get the yield Q30 value for RAW  and for read 1 and 2
                for c in p_temp[sample_all_index][0][l_index][t_index][0].iter('YieldQ30'):
                    raw_yield_q30_value +=int(c.text)
                for c in p_temp[sample_all_index][0][l_index][t_index][0].iter('QualityScoreSum'):
                    raw_quality_value +=int(c.text)
                 # get the yield value for PF and for read 1 and 2
                for c in p_temp[sample_all_index][0][l_index][t_index][1].iter('Yield'):
                    pf_yield_value +=int(c.text)
                # get the yield Q30 value for PF and for read 1 and 2
                for c in p_temp[sample_all_index][0][l_index][t_index][1].iter('YieldQ30'):
                    pf_yield_q30_value +=int(c.text)
                for c in p_temp[sample_all_index][0][l_index][t_index][1].iter('QualityScoreSum'):
                    pf_quality_value +=int(c.text)
            list_raw_yield.append(str(raw_yield_value))
            list_raw_yield_q30.append(str(raw_yield_q30_value))
            list_raw_qualityscore.append(str(raw_quality_value))
            list_pf_yield.append(str(pf_yield_value))
            list_pf_yield_q30.append(str(pf_yield_q30_value))
            list_pf_qualityscore.append(str(pf_quality_value))

        stats_result[projects[i]]['RAW_Yield']=list_raw_yield
        stats_result[projects[i]]['RAW_YieldQ30']=list_raw_yield_q30
        stats_result[projects[i]]['RAW_QualityScore']=list_raw_qualityscore
        stats_result[projects[i]]['PF_Yield']=list_pf_yield
        stats_result[projects[i]]['PF_YieldQ30']=list_pf_yield_q30
        stats_result[projects[i]]['PF_QualityScore']=list_pf_qualityscore
        logger.info('completed parsing for xml stats for project %s', projects[i])

    unknow_lanes  = []
    unknow_barcode_start_index= len(projects)
    counter=0
    logger.info('Collecting the Top Unknow Barcodes')
    for un_child in root_conv.iter('TopUnknownBarcodes'):
        un_index= unknow_barcode_start_index + counter
        p_temp=root_conv[0][un_index][0]
        unknow_barcode_lines=p_temp.findall('Barcode')
        unknow_bc_count=[]
        for lanes in unknow_barcode_lines:
            unknow_bc_count.append(lanes.attrib)

        unknow_lanes.append(unknow_bc_count)
        counter +=1
    stats_result['TopUnknownBarcodes']= unknow_lanes
    logger.info('Complete XML parsing ')

    return stats_result


def store_raw_xml_stats(stats_projects, run_id,logger):
    error_found = False
    for project in stats_projects:
        if project == 'TopUnknownBarcodes':
            continue
        logger.info('processing project %s with rund_id = %s', project, run_id)
        if project == 'all' or project == 'default':
            logger.info('Found project %s setting the project_id to NULL', project)
            project_id= None
            default_all = project
        else:
            if Projects.objects.filter (projectName__exact = project).exists():
                p_name_id=Projects.objects.get(projectName__exact = project).id
            else:
                logger.error('ERROR:: Project name inside the report does not match with the one define in the sample sheet')
                print ('ERROR::: Project name ', project ,' is not include in the Sample Sheet ')
                error_found = True
                continue
            project_id= Projects.objects.get(pk=p_name_id)
            default_all = None
        # save the information when no error is found. This condition is set to avoid saving information on database
        # because the information stored before getting the error must be deleted
        if not error_found :
            raw_stats_xml = RawStatisticsXml (runprocess_id=RunProcess.objects.get(pk=run_id),
                                          project_id = project_id, defaultAll = default_all,
                                          rawYield= stats_projects[project]['RAW_Yield'], rawYieldQ30= stats_projects[project]['RAW_YieldQ30'],
                                          rawQuality= stats_projects[project]['RAW_QualityScore'], PF_Yield= stats_projects[project]['PF_Yield'],
                                          PF_YieldQ30= stats_projects[project]['PF_YieldQ30'], PF_QualityScore =stats_projects[project]['PF_QualityScore'])

            logger.info('saving raw stats for %s project', project)
            raw_stats_xml.save()
    if error_found :
        # delete all information stored on database
        if RawStatisticsXml.objects.filter (runprocess_id__exact = run_id).exists():
            projects_to_delete = RawStatisticsXml.objects.filter (runprocess_id__exact = run_id)
            logger.info('Deleting stored RawStatisticsXml information ')
            for project in projects_to_delete :
                project.delete()
                logger.debug('Deleted the RawStatisticsXml for project %s', project)
        return 'ERROR'

    logger.info('Raw XML data have been stored for all projects ')
    return ''


def process_xml_stats(stats_projects, run_id, logger):
    # get the total number of read per lane
    M_BASE=1.004361/1000000
    logger.debug('starting the process_xml_stats method')
    total_cluster_lane=(stats_projects['all']['PerfectBarcodeCount'])
    logger.info('processing flowcell stats for %s ', run_id)
    number_of_lanes=get_machine_lanes(run_id)
    for project in stats_projects:
        if project == 'TopUnknownBarcodes':
            continue
        flow_raw_cluster, flow_pf_cluster, flow_yield_mb = 0, 0, 0
        for fl_item in range(number_of_lanes):
             # make the calculation for Flowcell
            flow_raw_cluster +=int(stats_projects[project]['BarcodeCount'][fl_item])
            flow_pf_cluster +=int(stats_projects[project]['PerfectBarcodeCount'][fl_item])
            flow_yield_mb +=float(stats_projects[project]['PF_Yield'][fl_item])*M_BASE


        flow_raw_cluster='{0:,}'.format(flow_raw_cluster)
        flow_pf_cluster='{0:,}'.format(flow_pf_cluster)
        flow_yield_mb= '{0:,}'.format(round(flow_yield_mb))
        sample_number=stats_projects[project]['sampleNumber']

        if project == 'all' or project == 'default' :
            logger.info('Found project %s setting the project_id to NULL', project)
            project_id= None
            default_all = project
        else:
            p_name_id=Projects.objects.get(projectName__exact = project).id
            project_id= Projects.objects.get(pk=p_name_id)
            default_all = None

        #store in database
        logger.info('Processed information for flow Summary for project %s', project)
        ns_fl_summary = NextSeqStatsFlSummary(runprocess_id=RunProcess.objects.get(pk=run_id),
                                project_id=project_id, defaultAll=default_all, flowRawCluster=flow_raw_cluster,
                                flowPfCluster=flow_pf_cluster, flowYieldMb= flow_yield_mb,
                                sampleNumber= sample_number)


        ns_fl_summary.save()
        logger.info('saving processing flowcell xml data  for project %s', project)


    for project in stats_projects:
        if project == 'TopUnknownBarcodes':
            continue
        logger.info('processing lane stats for %s', project)

        for i in range (number_of_lanes):
            # get the lane information
            lane_number=str(i + 1)
            pf_cluster_int=(int(stats_projects[project]['PerfectBarcodeCount'][i]))
            pf_cluster='{0:,}'.format(pf_cluster_int)
            perfect_barcode=(format(int(stats_projects[project]['PerfectBarcodeCount'][i])*100/int(stats_projects[project]['BarcodeCount'][i]),'.3f'))
            percent_lane=  format(float(int(pf_cluster_int)/int(total_cluster_lane[i]))*100, '.3f')
            one_mismatch=stats_projects[project]['OneMismatchBarcodeCount'][i]
            yield_mb= '{0:,}'.format(round(float(stats_projects[project]['PF_Yield'][i])*M_BASE))

            bigger_q30=format(float(stats_projects[project]['PF_YieldQ30'][i])*100/float( stats_projects[project]['PF_Yield'][i]),'.3f')

            mean_quality=format(float(stats_projects[project]['PF_QualityScore'][i])/float(stats_projects[project]['PF_Yield'][i]),'.3f')

            # make the calculation for Flowcell
            flow_raw_cluster = stats_projects[project]['BarcodeCount'][i]
            flow_pf_cluster = stats_projects[project]['PerfectBarcodeCount'][i]
            flow_yield_mb ='{0:,}'.format(round(float(stats_projects[project]['PF_Yield'][i])*M_BASE))

            #store in database
            if project == 'all' or project == 'default':
                logger.info('Found project %s setting the project_id to NULL', project)
                project_id= None
                default_all =project
            else:
                p_name_id=Projects.objects.get(projectName__exact = project).id
                project_id= Projects.objects.get(pk=p_name_id)
                default_all = None

            #store in database
            logger.info('Processed information for Lane %s for project %s', lane_number, project)
            ns_lane_summary = NextSeqStatsLaneSummary(runprocess_id=RunProcess.objects.get(pk=run_id),
                                                 project_id=project_id, defaultAll = default_all, lane = lane_number,
                                                 pfCluster=pf_cluster, percentLane=percent_lane, perfectBarcode=perfect_barcode,
                                                 oneMismatch= one_mismatch, yieldMb=yield_mb,
                                                 biggerQ30=bigger_q30, meanQuality=mean_quality )

            ns_lane_summary.save()

    logger.info ('processing the TopUnknownBarcodes')
    for project in stats_projects:
        if project == 'TopUnknownBarcodes':
            for un_lane in range(number_of_lanes) :
                logger.info('Processing lane %s for TopUnknownBarcodes', un_lane)
                count_top=0
                lane_number=str(un_lane + 1)
                top_number =1
                for barcode_line in stats_projects[project][un_lane]:
                    barcode_count= '{0:,}'.format(int(barcode_line['count']))
                    barcode_sequence= barcode_line['sequence']

                    raw_unknow_barcode = RawTopUnknowBarcodes(runprocess_id=RunProcess.objects.get(pk=run_id),
                                                             lane_number = lane_number, top_number=str(top_number),
                                                             count=barcode_count, sequence=barcode_sequence)
                    raw_unknow_barcode.save()
                    top_number +=1


def parsing_sample_project_xml(run_id,demux_file, conversion_file, logger):
    total_p_b_count=[0,0,0,0]
    sample_result_dict={}
    #demux_file='example.xml'
    demux_stat=ET.parse(demux_file)
    root=demux_stat.getroot()
    projects=[]
    number_of_lanes=get_machine_lanes(run_id)
    logger.info('Starting parsing DemultiplexingStats.XML for getting Sample information')
    for child in root.iter('Project'):
        projects.append(child.attrib['name'])

    for i in range(len(projects)):
        if projects [i] == 'default' or projects [i] == 'all':
            continue
        p_temp=root[0][i]
        samples=p_temp.findall('Sample')
        sample_dict ={}
        for index in range (len(samples)):
            sample_name = samples[index].attrib['name']
            if sample_name == 'all':
                continue
            barcodeCount , perfectBarcodeCount = 0 , 0

            sample_stats={}
            sample_stats ['barcodeName'] = samples[index].find ('Barcode').attrib['name']

            for bar_count in p_temp[index][0].iter('BarcodeCount'):
                barcodeCount += int(bar_count.text)
            for p_bar_count in p_temp[index][0].iter('PerfectBarcodeCount'):
                perfectBarcodeCount += int(p_bar_count.text)
            sample_stats['BarcodeCount']=barcodeCount
            sample_stats['PerfectBarcodeCount']=perfectBarcodeCount
            sample_dict[sample_name] = sample_stats

            sample_result_dict[projects[i]]=sample_dict
    logger.info('Complete parsing from demux file for sample and for project %s', projects[i])


    conversion_stat=ET.parse(conversion_file)
    root_conv=conversion_stat.getroot()
    projects=[]
    logger.info('Starting conversion for conversion file')
    for child in root_conv.iter('Project'):
        projects.append(child.attrib['name'])
    for i in range(len(projects)):
        if projects [i] == 'default' or projects [i] == 'all':
            continue
        p_temp=root_conv[0][i]
        samples=p_temp.findall('Sample')

        for s_index in range (len (samples)):
            sample_name = samples[s_index].attrib['name']
            if sample_name == 'all':
                continue
            quality_per_sample = {}
            raw_yield_value = 0
            raw_yield_q30_value = 0
            raw_quality_value = 0
            pf_yield_value = 0
            pf_yield_q30_value = 0
            pf_quality_value = 0

            for l_index in range(number_of_lanes):
                tiles_index = len(p_temp[s_index][0][l_index].findall ('Tile'))
                for t_index in range(tiles_index):
                         # get the yield value for RAW and for read 1 and 2
                    for c in p_temp[s_index][0][l_index][t_index][0].iter('Yield'):
                        raw_yield_value +=int(c.text)
                        # get the yield Q30 value for RAW  and for read 1 and 2
                    for c in p_temp[s_index][0][l_index][t_index][0].iter('YieldQ30'):
                        raw_yield_q30_value +=int(c.text)
                    for c in p_temp[s_index][0][l_index][t_index][0].iter('QualityScoreSum'):
                        raw_quality_value +=int(c.text)
                     # get the yield value for PF and for read 1 and 2
                    for c in p_temp[s_index][0][l_index][t_index][1].iter('Yield'):
                        pf_yield_value +=int(c.text)
                    # get the yield Q30 value for PF and for read 1 and 2
                    for c in p_temp[s_index][0][l_index][t_index][1].iter('YieldQ30'):
                        pf_yield_q30_value +=int(c.text)
                    for c in p_temp[s_index][0][l_index][t_index][1].iter('QualityScoreSum'):
                        pf_quality_value +=int(c.text)

            sample_result_dict[projects[i]][sample_name]['RAW_Yield']=raw_yield_value
            sample_result_dict[projects[i]][sample_name]['RAW_YieldQ30']=raw_yield_q30_value
            sample_result_dict[projects[i]][sample_name]['RAW_QualityScore']=raw_quality_value
            sample_result_dict[projects[i]][sample_name]['PF_Yield']=pf_yield_value
            sample_result_dict[projects[i]][sample_name]['PF_YieldQ30']=pf_yield_q30_value
            sample_result_dict[projects[i]][sample_name]['PF_QualityScore']=pf_quality_value
        logger.info('completed parsing for xml stats for project %s', projects[i])


    logger.info('Complete XML parsing  for getting Samples')

    return sample_result_dict



def store_samples_projects(sample_project_stats, run_id, logger):
    # get the total number of read per lane
    M_BASE=1.004361/1000000
    logger.debug('starting store_sample_projects method')

    logger.info('processing flowcell stats for %s ', run_id)

    for project in sample_project_stats:
        # find the total number of PerfectBarcodeCount in the procjec to make percent calculations
        total_perfect_barcode_count = 0
        for sample in sample_project_stats[project]:
            total_perfect_barcode_count += sample_project_stats[project][sample] ['PerfectBarcodeCount']
        for sample in sample_project_stats[project]:
            sample_name = sample
            barcode_name = sample_project_stats[project][sample]['barcodeName']
            perfect_barcode = int(sample_project_stats[project][sample] ['PerfectBarcodeCount'])
            percent_in_project = format (float(perfect_barcode) *100 /total_perfect_barcode_count,'.2f')
            perfect_barcode = '{0:,}'.format(perfect_barcode)
            yield_mb = '{0:,}'.format(round(float(sample_project_stats[project][sample] ['PF_Yield'])*M_BASE))
            if sample_project_stats[project][sample] ['PF_Yield'] > 0:
                bigger_q30=format(float(sample_project_stats[project][sample]['PF_YieldQ30'])*100/float( sample_project_stats[project][sample]['PF_Yield']),'.3f')
                mean_quality=format(float(sample_project_stats[project][sample]['PF_QualityScore'])/float(sample_project_stats[project][sample]['PF_Yield']),'.3f')
            else:
                bigger_q30 = 0
                mean_quality =0
            p_name_id=Projects.objects.get(projectName__exact = project).id
            project_id= Projects.objects.get(pk=p_name_id)

            sample_to_store = SamplesInProject (project_id = project_id, sampleName = sample_name,
                            barcodeName = barcode_name, pfClusters = perfect_barcode,
                            percentInProject = percent_in_project , yieldMb = yield_mb,
                            qualityQ30 = bigger_q30, meanQuality = mean_quality)

            sample_to_store.save()
            logger.debug('Stored sample %s', sample_name)
        logger.info('Stored sample for the project %s', project)

        #store in database
        logger.info('Processed information sample %s', project)



def process_run_in_samplesent_state (process_list, logger):
     # prepare a dictionary with key as run_name and value the RunID
     processed_run=[]
     for run_item in process_list:
         logger.info ('Running the process sample sent state for %s', run_item)
         run_be_processed_id=RunProcess.objects.get(runName__exact=run_item).id
         logger.debug ('Run ID for the run process to be update is:  %s', run_be_processed_id)
         #run_Id_for_searching=RunningParameters.objects.get(runName_id= run_be_processed_id)
         update_run_state(run_be_processed_id, 'Process Running', logger)
         processed_run.append(run_be_processed_id)
     return processed_run

def process_run_in_processrunning_state (process_list, logger):
    processed_run=[]
    logger.debug('starting the process_run_in_processrunning_state method')
    try:
        conn=open_samba_connection()
        logger.info('check the Sucessful connection to NGS_Data before starting processing runing state method')

    except:
        return('Error')

    share_folder_name = wetlab_config.SAMBA_SHARED_FOLDER_NAME
    for run_item in process_list:
        logger.debug ('processing the run %s in process running state' , run_item)
        run_be_processed_id=RunProcess.objects.get(runName__exact=run_item).id
        run_Id_used=str(RunningParameters.objects.get(runName_id= run_be_processed_id))
        logger.debug ('found the run ID  %s' , run_Id_used )
        run_folder=os.path.join('/',run_Id_used,'Data/Intensities/BaseCalls')
        # check if runCompletion is avalilable
        logger.debug ('found the run ID  %s' , run_Id_used )
        try:  #####
            file_list = conn.listPath( share_folder_name, run_folder)
        except: #####
            logger.error('no folder found for run ID %s', run_Id_used)  #####
            continue  #####
        #import pdb; pdb.set_trace()
        found_report_directory = 0
        for sh in file_list:
            if sh.filename =='Reports' :
                logger.info('bcl2fastq has been completed for run %s', run_Id_used)
                processed_run.append(run_Id_used)
                update_run_state(run_be_processed_id, 'Bcl2Fastq Executed', logger)
                update_project_state(run_be_processed_id, 'B2FqExecuted', logger)
                # Get the time when  the Bcl2Fastq process is ending
                conversion_stats_file = os.path.join (run_Id_used,'Data/Intensities/BaseCalls/Stats/', 'ConversionStats.xml')
                conversion_attributes = conn.getAttributes(wetlab_config.SAMBA_SHARED_FOLDER_NAME ,conversion_stats_file)
                run_date = RunProcess.objects.get(pk=run_be_processed_id)
                run_date.bcl2fastq_finish_date = datetime.datetime.fromtimestamp(int(conversion_attributes.create_time)).strftime('%Y-%m-%d %H:%M:%S')
                run_date.save()
                logger.info ('Updated the Bcl2Fastq time in the run %s', run_item)
                break
            else:
                logger.debug('The directory %s has been found while looking for completion of the execution of bcl2fastq', sh.filename)
        if found_report_directory:
            logger.info('blc2fastq has been completed for the Run ID %s  it is now on Bcl2Fastq Executed state', run_Id_used)
        else:
            logger.info('blc2fastq was not finish for the Run ID %s  waiting for Bcl2Fastq to be completed', run_Id_used)


    # close samba connection
    conn.close()
    logger.info('Closing the remote connection ')
    return processed_run



def process_run_in_bcl2F_q_executed_state (process_list, logger):
    processed_run=[]
    # get the directory of samba to fetch the files
    share_folder_name = wetlab_config.SAMBA_SHARED_FOLDER_NAME
    local_dir_samba= wetlab_config.RUN_TEMP_DIRECTORY_PROCESSING
    remote_stats_dir= 'Data/Intensities/BaseCalls/Stats/'
    demux_file=os.path.join(local_dir_samba,'DemultiplexingStats.xml')
    conversion_file=os.path.join(local_dir_samba,'ConversionStats.xml')
    run_info_file=os.path.join(local_dir_samba, 'RunInfo.xml')

    logger.debug('Executing process_run_in_bcl2F_q_executed_state method')

    # check the connectivity to remote server
    try:
        conn=open_samba_connection()
        logger.info('Successful connection for updating run on bcl2F_q' )
    except:
        logger.error('ERROR:: Unable to connect to remote server')
        return 'Error'

    for run_item in process_list:

        logger.info ('Processing the process on bcl2F_q for the run %s', run_item)
        run_processing_id=RunProcess.objects.get(runName__exact=run_item).id
        run_Id_used=str(RunningParameters.objects.get(runName_id= run_processing_id))

        update_run_state(run_processing_id, 'Running Stats', logger)
        #copy the demultiplexingStats.xml file to wetlab/tmp/processing


        samba_demux_file=os.path.join('/',run_Id_used,remote_stats_dir, 'DemultiplexingStats.xml')
        logger.debug('path to fetch demultiplexingStats is %s',  samba_demux_file)
        try:
            with open(demux_file ,'wb') as demux_fp :
                conn.retrieveFile(share_folder_name, samba_demux_file, demux_fp)
        except:
            logger.error('Unable to fetch the DemultiplexingStats.xml file for RunID %s', run_Id_used)
            os.remove(demux_file)
            logger.debug('deleting DemultiplexingStats file  for RunID %s' , run_Id_used)
            continue
        logger.info('Fetched the DemultiplexingStats.xml')
        #copy the ConversionStats.xml file to wetlab/tmp/processing
        samba_conversion_file=os.path.join('/', run_Id_used,remote_stats_dir,'ConversionStats.xml')
        try:
            with open(conversion_file ,'wb') as conv_fp :
                conn.retrieveFile(share_folder_name, samba_conversion_file, conv_fp)
        except:
            logger.error('Unable to fetch the ConversionStats.xml file for RunID %s', run_Id_used)
            os.remove(conversion_file)
            os.remove(demux_file)
            logger.debug('deleting ConversionStats and DemultiplexingStats file  for RunID %s' , run_Id_used)
            continue
        logger.info('Fetched the conversionStats.xml')
        # copy RunInfo.xml  file to process the interop files
        try:
            with open(run_info_file ,'wb') as runinfo_fp :
                samba_conversion_file=os.path.join('/', run_Id_used,'RunInfo.xml')
                conn.retrieveFile(share_folder_name, samba_conversion_file, runinfo_fp)
        except:
            logger.error('Unable to fetch the RunInfo.xml file for RunID %s', run_Id_used)
            os.remove(run_info_file)
            os.remove(conversion_file)
            os.remove(demux_file)
            logger.debug('deleting RunInfo, ConversionStats and DemultiplexingStats file  for RunID %s' , run_Id_used)
            continue
        logger.info('Fetched the RunInfo.xml file')

        # copy all binary files in interop folder to local  documents/wetlab/tmp/processing/interop
        interop_local_dir_samba= os.path.join(local_dir_samba, 'InterOp')
        remote_interop_dir=os.path.join('/',run_Id_used,'InterOp')
        try:
            file_list = conn.listPath( share_folder_name, remote_interop_dir)
            logger.info('InterOp folder exists on the RunID %s', run_Id_used)
            run_parameters_file=os.path.join(local_dir_samba,'runParameters.xml')
            try:
                with open(run_parameters_file ,'wb') as runparam_fp :
                    samba_conversion_file=os.path.join('/', run_Id_used,'runParameters.xml')
                    conn.retrieveFile(share_folder_name, samba_conversion_file, runparam_fp)
                logger.info('Fetched the runParameters.xml file. Written as: '
                    +str(run_parameters_file))
            except:
                logger.error('Unable to fetch the runParameters.xml file for RunID %s', run_Id_used)
                os.remove(run_parameters_file)
                logger.debug('Deleting runParameters file  for RunID %s' , run_Id_used)

        except:
            logger.error('ERROR:: InterOP folder does not exist on RunID %s', run_Id_used)
            os.remove(run_info_file)
            os.remove(conversion_file)
            os.remove(demux_file)
            logger.debug('deleting RunInfo, ConversionStats and DemultiplexingStats file  for RunID %s' , run_Id_used)
            print('ERROR:: InterOP folder does not exist on RunID ', run_Id_used)
            continue
        error_in_interop = False
        for sh in file_list:
            if sh.isDirectory:
                continue
            else:
                interop_file_name=sh.filename
                remote_interop_file=os.path.join(remote_interop_dir, interop_file_name)
                copy_file=os.path.join(interop_local_dir_samba, interop_file_name)
                try:
                    with open(copy_file ,'wb') as cp_fp :
                        remote_file=os.path.join(remote_interop_dir,)
                        logger.debug('File %s to be copied on the local directory', interop_file_name)
                        conn.retrieveFile(share_folder_name, remote_interop_file, cp_fp)
                        logger.info('Copied %s to local Interop folder', interop_file_name)
                except:
                    logger.error("Not be able to fetch the file %s", interop_file_name)
                    os.remove(run_info_file)
                    os.remove(conversion_file)
                    os.remove(demux_file)
                    logger.debug('deleting files RunInfo, ConversionStats and DemultiplexingStats ')
                    logger.debug('because of error when fetching interop files for RunID  %s' , run_Id_used)
                    # deleting local Interop files
                    file_list_to_delete = os.list(interop_local_dir_samba)
                    logger.debug ('Deleting the local interop files ')
                    for file_name in file_list_to_delete:
                        if file_name == '.' or '..' :
                            continue
                        else:
                            os.remove(file_name)
                    error_in_interop= True
                    break
        if error_in_interop :
            continue
        else:
            # parsing the files to get the xml Stats
            logger.info('processing the XML files')
            xml_stats=parsing_statistics_xml(run_processing_id, demux_file, conversion_file, logger)
            result_of_raw_saving = store_raw_xml_stats(xml_stats,run_processing_id, logger)
            if result_of_raw_saving == 'ERROR':
                update_run_state(run_processing_id, 'ERROR-on-Raw-SavingStats', logger)
                update_project_state(run_processing_id, 'ERROR-on-Raw-SavingStats', logger)
                logger.error('Stopping process for this run an starting deleting the files')
            else:
                process_xml_stats(xml_stats,run_processing_id, logger)

                # parsing and processing the project samples
                sample_project_stats = parsing_sample_project_xml (run_processing_id,demux_file, conversion_file, logger)
                store_samples_projects (sample_project_stats, run_processing_id, logger)

                logger.info('processing interop files')
                # processing information for the interop files
                number_of_lanes = get_machine_lanes(run_processing_id)

                process_binStats(local_dir_samba, run_processing_id, logger, number_of_lanes)
                # Create graphics
                graphic_dir=os.path.join(settings.MEDIA_ROOT,wetlab_config.RUN_TEMP_DIRECTORY_PROCESSING)

                create_graphics(graphic_dir, run_processing_id, run_Id_used, logger)

                processed_run.append(run_Id_used)
                logger.info('run id %s is now on Completed state', run_Id_used)
                update_run_state(run_processing_id, 'Completed', logger)
                update_project_state(run_processing_id, 'Completed', logger)
            # clean up the used files and directories
            logger.info('starting the clean up for the copied files from remote server ')
            os.remove(demux_file)
            logger.debug('Demultiplexing file have been removed from %s', demux_file)
            os.remove(conversion_file)
            logger.debug('ConversionStats file have been removed from %s', conversion_file)
            os.remove(run_info_file)
            logger.debug('RunInfo file have been removed from %s', run_info_file)
            for file_object in os.listdir(interop_local_dir_samba):
                file_object_path = os.path.join(interop_local_dir_samba, file_object)
                if os.path.isfile(file_object_path):
                    logger.debug('Deleting file %s' , file_object_path)
                    os.remove(file_object_path)
            logger.info('xml files and binary files from InterOp folder have been removed')
            ## connect to server to get the disk space utilization of the folders
            get_run_disk_utilization (conn, run_Id_used, run_processing_id, logger)
            # Update the run with the date of the run completion
            completion_date = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            run_date_to_update = RunProcess.objects.get(pk = run_processing_id)
            run_date_to_update.process_completed_date = completion_date
            run_date_to_update.save()

    # close samba connection
    conn.close()
    logger.info('Samba connection close. Sucessful copy of files form remote server')
    return processed_run
'''
def find_state_and_save_data(run_name,run_folder):
    run_file='RunInfo.xml'
    run_parameter='RunParameters.xml'

    try:
        rn_found = RunProcess.objects.get(runName__exact=run_name)
    except:
        #os.chdir('wetlab/tmp/logs')
        with open('documents/wetlab/tmp/logs/wetlab.log', 'a') as log_file:
            time_log = str(datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y"),'\n')
            log_write(time_log)
            error_text= str('[ERROR]--  run name ',run_name, 'was not found in database  \n')
            log_file.write(error_text)
            return 'ERROR'
    rn_state = rn_found.get_state()
    if rn_state == 'Recorded':
        copy_sample_sheet(rn_found, run_folder)
    elif rn_state == 'Sample Sent':
        save_running_info(run_file, run_parameter, rn_found)
        rn_found.runState='Process Running'
    elif rn_state == 'Process Running':
        ## check if the run is completed by checking if RunCompletionStatus.xml exists
        rn_found.runState='Bcl2Fastq Executed'
    else:
        rn_found.runState='Completed'
'''

def find_not_completed_run (logger):
    working_list={}
    state_list = ['Sample Sent','Process Running','Bcl2Fastq Executed']
    # get the run that are not completed
    for state in state_list:
        logger.info('start looking for incompleted runs in state %s', state)

        if RunProcess.objects.filter(runState__exact = state).exists():
            working_list[state]=RunProcess.objects.filter(runState__exact = state)
            logger.debug ('found  %s not completed runs' , working_list  )

    processed_run={}
    for state in working_list:
        logger.debug('find_not_completed_run / working_list= '+str(working_list)) #Debug
        logger.info ('Start processing the run found for state %s', state)
        if state == 'Sample Sent':
            logger.debug ('found sample sent in state ')
            processed_run[state]=process_run_in_samplesent_state(working_list['Sample Sent'], logger)
        elif state == 'Process Running':
            logger.debug('Found runs for Process running %s', working_list['Process Running'])
            processed_run[state]=process_run_in_processrunning_state(working_list['Process Running'], logger)

        else:
            logger.debug('Found runs for Bcl2Fastq %s', working_list['Bcl2Fastq Executed'])
            processed_run[state]=process_run_in_bcl2F_q_executed_state(working_list['Bcl2Fastq Executed'], logger)

    return (processed_run)


'''

def copy_sample_sheet(run_name, run_folder):
    ## get the sample sheet file
    sample_file=rn_found.get_sample_file()
    ## send the sample sheet file to the run folder
    open_samba_connection()
    with open('wetlab/tmp/logs/wetlab.log', 'a') as log_file:
        time_log = str(datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y"),'\n')
        log_write(time_log)
        ## opening the samba connection
        info_text = str('[INFO]--  Openning the connection to samba server \n')
        log_file.write(info_text)
        #
        #
        #
        info_text = str('[INFO]--  Sending Sample Sheet to folder ',run_folder, ' for run ',run_name, '\n')
        log_file.write(info_text)
        #
        ## waiting for file copy completion
        info_text = str('[INFO]--  run name ',run_name, 'was sent to folder ',run_folder ,'\n')
        log_file.write(info_text)


        info_text = str('[INFO]--  run name ',run_name, 'state was changed to SampleSent \n')
        log_file.write(info_text)
        log_file.close()


def fetch_runID_parameter():
    runparameters_file='wetlab/tmp/tmp/processing/RunParameters.xml'
    data_from_runparameters=get_running_info(runparameters_file)
    run_name=data_from_runparameters['ExperimentName']
    ## to include the information on database we get the index first
    if RunProcess.object.filter(runName__exact = run_name).exists():
        r_name_id = RunProcess.object.filter(runName__exact = run_name).id
        r_name_id.RunID=data_from_runparameters['RunID']


def perform_xml_stats (xml_statistics, run_name_value):
    for project in xml_statistics:
        print (project)
        ### Flowcell Summary
        fl_pf_yield_sum=0
        fl_raw_yield_sum=0
        fl_mbases=0
        for values in xml_statistics[project]['PF_Yield']:
            fl_pf_yield_sum+= int(values)
        for values in xml_statistics[project]['RAW_Yield']:
            fl_raw_yield_sum+= int(values)
        for values in xml_statistics[project]['']:
            print()
'''
