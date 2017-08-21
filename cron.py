import paramiko
from datetime import datetime
from .utils.stats_calculation import *
from .utils.parsing_run_info import *

import os , sys

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

def check_recorded_folder ():
    time_start= datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(time_start ,'  looking for new run in directory on wetlab/tmp/recorded \n')
    path='iSkyLIMS/wetlab/tmp/recorded/'
    dir_wetlab=os.getcwd()
    print('woring dir = ', dir_wetlab)
    if os.listdir(path):
        # There are sample sheet files that need to be processed
        updated_run=process_run_in_recorded_state ()
        if updated_run == 'Error':
            print('No connection is available to Flavia \n')
            print ('Exiting the process for searching run in recorded state \n')
        else:
            for run_changed in updated_run:
                print('The run ', run_changed, 'is now on Sample Sent state\n')
    else:
        print( 'Exiting the crontab for record_folder. No directories have been found \n')

def check_not_finish_run():
    time_start= datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(time_start ,'  looking for not finish run  \n')
    updated_run=found_not_completed_run()
    if (not updated_run):
        print('Exiting the crontab for not finish run without performing any changes \n')
        return
    if updated_run != 'Error':
        print('Exiting the crontab for not finish run performing the changes in the following run \n')
        for state, value in updated_run.iteritems():
            print ('The following Run are now in the state ', state, '\n')
            for run_changed in value:
                print('--', run_changed, '\n')
    else:
            print('No connection is available to Flavia \n')

            
    