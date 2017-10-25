import paramiko
from datetime import datetime
from .utils.stats_calculation import *
from .utils.parsing_run_info import *

import os , sys
import logging
from logging.handlers import RotatingFileHandler

def open_log():
    
    LOG_FILENAME = 'checking_uncompleted_run.log'
    log_name=os.path.join('iSkyLIMS/wetlab/log/', LOG_FILENAME)
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
    
'''
def check_crontab():
    time_start= datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(time_start ,'  checking the time for starting crontab')

def check_crontab2():
    time_start= datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(time_start ,'  checking the time for starting crontab2 has to be on min 28 ')


def createSSHClient(server, port, user, password):
    client = paramiko.SSHClient()
    client.load_system_host_keys()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(server, port, user, password)
    return client 




def fetching_stats_scheduled_job ():
    barbarroja = '10.15.60.54'
    remote_dir_stats_xml='161123_NS500454_0096_AHFGV5BGXY/fastq_files4/Stats/'
    remote_dir_run='161123_NS500454_0096_AHFGV5BGXY/'
    remote_interop = '161123_NS500454_0096_AHFGV5BGXY/Interop/'

    demux_file= remote_dir_stats_xml + 'DemultiplexingStats.xml'
    conversion_file=remote_dir_stats_xml + 'ConversionStats.xml'
    run_file=remote_dir_run + 'RunInfo.xml'
    run_parameter=remote_dir_run + 'RunParameters.xml'
    
    local_working_dir = '/home/bioinfo/web_carlosIII/polls/documents/uploadFromServer/'
    local_stats = '/home/bioinfo/web_carlosIII/polls/documents/uploadFromServer/Stats/'
    local_interop= '/home/bioinfo/web_carlosIII/polls/documents/uploadFromServer/Interop/'
    run_file = local_working_dir + 'RunInfo.xml'
    run_parameter=  local_working_dir + 'RunParameters.xml'
    copied_files = 0
    time_start= datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print('Transfer files started at:', time_start,'\n')
    #ssh = createSSHClient(barbarroja, 22, 'lchapado', 'chapadomaster')
    #ftp = ssh.open_sftp()
    #ftp.get('luis3_result.log', '/home/bioinfo/web_carlosIII/polls/uploadFromServer/barba-luis3-result.log')
    #ftp.put('localfile', 'remotefile')
    #for f in ftp.listdir(remote_interop):
    #    ftp.get(remote_interop+f, local_interop+f)
    #    copied_files +=1
    #ftp.close()
    print('Total files copied: ',copied_files)
    time_end= datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print('Transfer files ended at:', time_end,'\n')
    print('Working on getting statistics\n')

    running_data = get_running_data(run_file,run_parameter)
    store_in_db(running_data, 'running_table')
    #xml_statistics = get_statistics_xml(demux_file, conversion_file)
    #store_in_db(xml_statistics, 'nextSeq_xml_table')
    
    return True
'''
def check_recorded_folder ():
    time_start= datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(time_start )
    print('Starting the process for recorded_folder ')
    logger=open_log()
    
    logger.info('Looking for new run in directory on documents/wetlab/tmp/recorded ')
    path='iSkyLIMS/documents/wetlab/tmp/recorded/'
    dir_wetlab=os.getcwd()
    logger.debug('working directory is %s', dir_wetlab)
    if os.listdir(path):
        # There are sample sheet files that need to be processed
        updated_run=process_run_in_recorded_state (logger)
        if updated_run == 'Error':
            logger.error('No connection is available to Flavia')
            logger.error('Exiting the process for searching run in recorded state ')
            time_stop= datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            print(time_stop)
            print('Exiting the check_recorder_folder module due to error when connecting to Flavia')
        else:
            for run_changed in updated_run:
                logger.info('The run  %s is now on Sample Sent state', run_changed)
            time_stop= datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            print(time_stop)
    else:
        time_stop= datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        print(time_stop)
        logger.info( 'Exiting the crontab for record_folder. No directories have been found')

def check_not_finish_run():
    time_start= datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(time_start )
    print('Starting the process for searching not completed runs ')
    logger=open_log()
    logger.info('starting execute the crontab for not finish run')
    updated_run=find_not_completed_run(logger)
    logger.debug('Display the list of the updated_run %s', updated_run)
    
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
            
