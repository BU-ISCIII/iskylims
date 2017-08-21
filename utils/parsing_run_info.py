#!/usr/bin/env python3

import sys, os, re
import xml.etree.ElementTree as ET
import time
import shutil
from  ..models import *

from smb.SMBConnection import SMBConnection

def open_samba_connection():
    ## open samba connection
    # There will be some mechanism to capture userID, password, client_machine_name, server_name and server_ip
    # client_machine_name can be an arbitary ASCII string
    # server_name should match the remote machine name, or else the connection will be rejected
    
    conn=SMBConnection('Luigi', 'Apple123', 'bioinfo', 'LUIGI-PC', use_ntlm_v2=True)
    conn.connect('192.168.1.3', 139)

    '''
    conn = SMBConnection(userid, password, client_machine_name, remote_machine_name, use_ntlm_v2 = True)
    conn.connect(server_ip, 139)
    '''
    return conn

def fetch_runID_parameter():
    runparameters_file='wetlab/tmp/tmp/processing/RunParameters.xml'
    data_from_runparameters=get_running_info(runparameters_file)
    run_name=data_from_runparameters['ExperimentName']
    ## to include the information on database we get the index first
    if RunProcess.object.filter(runName__exact = run_name).exists():
        r_name_id = RunProcess.object.filter(runName__exact = run_name).id
        r_name_id.RunID=data_from_runparameters['RunID']



def save_run_info(run_info, run_parameter, run_id):
    running_data={}
    image_channel=[]
    #################################################
    ## parsing RunInfo.xml file
    #################################################
    run_data=ET.parse(run_info)
    run_root=run_data.getroot()
    p_run=run_root[0]
    for i in run_root.iter('Name'):
        image_channel.append(i.text)
                 
    running_data['ImageChannel']=image_channel
    running_data['Flowcell']=p_run.find('Flowcell').text
    
    running_data['ImageDimensions']=p_run.find('ImageDimensions').attrib
    running_data['FlowcellLayout']=p_run.find('FlowcellLayout').attrib
    #################################################
    ## parsing RunParameter.xml file
    #################################################
    parameter_data=ET.parse(run_parameter)
    parameter_data_root=parameter_data.getroot()
    p_parameter=parameter_data_root[1]
    running_data['RunID']=parameter_data_root.find('RunID').text  
    running_data['ExperimentName']=parameter_data_root.find('ExperimentName').text
    running_data['RTAVersion']=parameter_data_root.find('RTAVersion').text
    running_data['SystemSuiteVersion']=parameter_data_root.find('SystemSuiteVersion').text
    running_data['LibraryID']=parameter_data_root.find('LibraryID').text
    running_data['Chemistry']=parameter_data_root.find('Chemistry').text
    running_data['RunStartDate']=parameter_data_root.find('RunStartDate').text
    running_data['AnalysisWorkflowType']=parameter_data_root.find('AnalysisWorkflowType').text
    running_data['RunManagementType']=parameter_data_root.find('RunManagementType').text
    running_data['PlannedRead1Cycles']=parameter_data_root.find('PlannedRead1Cycles').text
    running_data['PlannedRead2Cycles']=parameter_data_root.find('PlannedRead2Cycles').text
    running_data['PlannedIndex1ReadCycles']=parameter_data_root.find('PlannedIndex1ReadCycles').text
    running_data['PlannedIndex2ReadCycles']=parameter_data_root.find('PlannedIndex2ReadCycles').text
    running_data['ApplicationVersion']=p_parameter.find('ApplicationVersion').text
    running_data['NumTilesPerSwath']=p_parameter.find('NumTilesPerSwath').text
    ###########################################
    ## saving data into database
    ###########################################
    print ('Saving to database  the run id ', run_id)
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

    #running_parameters.save()
    
       
    return

def fetch_exp_name_from_run_info (local_run_parameter_file):

    ## look for   <ExperimentName>NextSeq_CNM_041</ExperimentName> in RunParameters.xml file
    fh=open(local_run_parameter_file ,'r')
    for line in fh:
        exp_name=re.search('^\s+<ExperimentName>(.*)</ExperimentName>',line)
        if exp_name:
            fh.close()
            return exp_name.group(1)


           
def process_run_in_recorded_state():
    try:
        conn=open_samba_connection()
        print('Sucessfully connection')
    except:
        return ('Error')
    processed_run_file, runlist = [] , []
    recorded_dir='iSkyLIMS/wetlab/tmp/recorded/'
    share_folder_name='Flavia'
    local_run_parameter_file='iSkyLIMS/wetlab/tmp/tmp/RunParameters.xml'
    local_run_info_file='iSkyLIMS/wetlab/tmp/tmp/RunInfo.xml'
    process_run_file='iSkyLIMS/wetlab/tmp/processed_run_file'
    processed_run=[]
    run_names_processed=[]
    ## get the list of the processed run
    if os.path.exists(process_run_file):
        fh = open (process_run_file,'r')
        for line in fh:
            line=line.rstrip()
            processed_run.append(line)
        fh.close()
    # Check if the directory from flavia has been processed
    file_list = conn.listPath( share_folder_name, '/')
    for sfh in file_list:
        if sfh.isDirectory:
            run_dir=(sfh.filename)
            if (run_dir == '.' or run_dir == '..'):
                continue
                # if the run folder has been already process continue searching
            if run_dir in processed_run_file:
                continue
            else:
                #copy the runParameter.xml file to wetlab/tmp/tmp
                with open(local_run_parameter_file ,'wb') as r_par_fp :
                    samba_run_parameters_file=os.path.join(run_dir,'RunParameters.xml')
                    conn.retrieveFile(share_folder_name, samba_run_parameters_file, r_par_fp)
                # look for the experience name  for the new run folder. Then find the run_id valued for it
                exp_name=fetch_exp_name_from_run_info(local_run_parameter_file)
                if  RunProcess.objects.filter(runName__icontains = exp_name).exists():
                #if True:
                    exp_name_id=str(RunProcess.objects.get(runName__exact = exp_name).id)
                    #exp_name_id ='2'
                    sample_sheet_tmp_dir=os.path.join('iSkyLIMS/wetlab/tmp/recorded',exp_name_id,'samplesheet.csv')
                    if os.path.exists(sample_sheet_tmp_dir):
                        # copy Sample heet file to samba directory
                        with open(sample_sheet_tmp_dir ,'rb') as  sample_samba_fp:
                            samba_sample_file= os.path.join(run_dir,'samplesheet.csv')
                            conn.storeFile(share_folder_name, samba_sample_file, sample_samba_fp)
                        # retrieve the runInfo.xml file from samba directory
                    # get the runIfnfo.xml to collect the  information for this run    
                    with open(local_run_info_file ,'wb') as r_info_fp :
                        samba_run_info_file=os.path.join(run_dir,'RunInfo.xml')
                        conn.retrieveFile(share_folder_name, samba_run_info_file, r_info_fp)
                    save_run_info (local_run_info_file, local_run_parameter_file, exp_name_id)
                        # delete the copy of the run files 
                    os.remove(local_run_info_file)
                    os.remove(local_run_parameter_file)
                        # delete the file and folder for the sample sheet processed
                    shutil.rmtree(os.path.join(recorded_dir, exp_name_id))
                        
                        # change the run  to SampleSent state
                    update_state(exp_name_id, 'Sample Sent')
                        # add the run_dir inside the processed_run file
                    processed_run.append(run_dir)
                    run_names_processed.append(exp_name)    
                else:
                    # error in log file
                    print ('The run ID ' , run_dir, 'does not match any run in the RunProcess object.\n') 
                    continue
    conn.close()
    fh =open (process_run_file,'w')
    for process in processed_run:
        fh.write(process)
        fh.write('\n')
    fh.close()
    # check if all run process file are handled
 
    list_dir=os.listdir(recorded_dir)
    if list_dir:
        print ('Warning: There are run in Recorded state that were not processed \n')
        for item in list_dir:
            print('directory for the run Id : ', item, '\n')
    return(run_names_processed)




def update_state(run_id, state):
    run=RunProcess.objects.get(pk=run_id)
    run.runState= state
    run.save()

def parsing_statistics_xml(demux_file, conversion_file):
 
    stats_result={}
    #demux_file='example.xml'
    demux_stat=ET.parse(demux_file)
    root=demux_stat.getroot()
    projects=[]
    for child in root.iter('Project'):
        projects.append(child.attrib['name'])
    
    for i in range(len(projects)):
        p_temp=root[0][i]
        samples=p_temp.findall('Sample')
        sample_all_index=len(samples)-1
        barcodeCount=[]
        perfectBarcodeCount=[]
        b_count=[]
        p_b_count=[]
        dict_stats={}
        for c in p_temp[sample_all_index].iter('BarcodeCount'):
            #b_count.append(c.text)
            barcodeCount.append(c.text)
        for c in p_temp[sample_all_index].iter('PerfectBarcodeCount'):
            p_b_count.append(c.text)
        perfectBarcodeCount.append(p_b_count)
        dict_stats['BarcodeCount']=barcodeCount
        dict_stats['PerfectBarcodeCount']=perfectBarcodeCount
        stats_result[projects[i]]=dict_stats
    
    #demux_file='example.xml'
    conversion_stat=ET.parse(conversion_file)
    root_conv=conversion_stat.getroot()
    projects=[]
    
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
        list_raw_quality=[]
        list_pf_yield=[]
        list_pf_yield_q30=[]
        list_pf_quality=[]
    
        for l_index in range(4):
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
            list_raw_quality.append(str(raw_quality_value))
            list_pf_yield.append(str(pf_yield_value))
            list_pf_yield_q30.append(str(pf_yield_q30_value))
            list_pf_quality.append(str(pf_quality_value))
                
        stats_result[projects[i]]['RAW_Yield']=list_raw_yield
        stats_result[projects[i]]['RAW_YieldQ30']=list_raw_yield_q30
        stats_result[projects[i]]['RAW_Quality']=list_raw_quality
        stats_result[projects[i]]['PF_Yield']=list_pf_yield
        stats_result[projects[i]]['PF_YieldQ30']=list_pf_yield_q30
        stats_result[projects[i]]['PF_Quality']=list_pf_quality
    return stats_result

def store_raw_xml_stats(stats_projects, run_id):
    for project in stats_projects:
        #project_id=Projects.objects.get(projectName__exact=project).id
        raw_stats_xml = RawStatisticsXml (runprocess_id=RunProgress.objects.get(pk=run_id),
                                          rawYield= project['RAW_Yield'], rawYieldQ30= project['RAW_YieldQ30'],
                                          rawQuality= project['RAW_Quality'], PF_Yield= project['PF_Yield'],
                                          PF_YieldQ30= project['PF_YieldQ30'], PF_Quality =project ['PF_Quality'],
                                          barcodeCount= project['barcodeCount'], perfectBarcodeCount= project['perfectBarcodeCount'],
                                          projectName=project)
        raw_stats_xml.save()

def process_xml_stats(stats_projects, run_id):
    for project in stats_projects:
        print('do something')

def process_run_in_samplesent_state (process_list):
     # prepare a dictionary with key as run_name and value the RunID
     for run_item in process_list:
        run_be_processed_id=RunProcess.objects.get(runName__exact=run_item).id
        run_Id_for_searching=RunningParameters.objects.get(runName_id= run_be_processed_id)
        update_state(run_be_processed_id, 'Process Running')
        
def process_run_in_proc_running_state (process_list):
    for run_item in process_list:
        run_be_processed_id=RunProcess.objects.get(runName__exact=run_item).id
        run_Id_for_searching=RunningParameters.objects.get(runName_id= run_be_processed_id)
        # check if runCompletion is avalilable
        file_list = conn.listPath( share_folder_name, run_Id_for_searching)
        for sfh in file_list:
            if sfh.filename=='runCompletion.xml' :
                update_state(run_be_processed_id, 'Bcl2Fq Executed')
    # close samba connection 
    conn.close()


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
   

def process_run_in_bcl2F_q_executed_state (process_list):
    
    for run_item in process_list:
        # change the state to Running Stats
        run_be_processed_id=RunProcess.objects.get(runName__exact=run_item).id
        run_Id_for_searching=RunningParameters.objects.get(runName_id= run_be_processed_id)
        update_state(run_be_processed_id, 'Running Stats')
        # get the directory of samba to fetch the files
        
        local_dir_samba= 'wetlab/tmp/processing'
        demux_file=os.path.join(local_dir_samba,'DemultiplexingStats.xml')
        conversion_file=os.path.join(local_dir_samba,'ConversionStats.xml')
        #copy the demultiplexingStats.xml file to wetlab/tmp/processing
        conn=open_samba_connection()
        with open(demux_file ,'wb') as demux_fp :
            samba_demux_file=os.path.join(run_Id_for_searching,'Stats/DemultiplexingStats.xml')
                #conn.retrieveFile('share', '/path/to/remote_file', fp)
            conn.retrieveFile(share_folder_name, samba_demux_file, demux_fp)
        #copy the demultiplexingStats.xml file to wetlab/tmp/processing
        with open(conversion_file ,'wb') as conv_fp :
            samba_conversion_file=os.path.join(run_Id_for_searching,'Stats/ConversionStats.xml')
                #conn.retrieveFile('share', '/path/to/remote_file', fp)
            conn.retrieveFile(share_folder_name, samba_conversion_file, conv_fp)
        # copy all binary files in interop folder to local  wetlab/tmp/processing/interop  
        interop_local_dir_samba= 'wetlab/tmp/processing/InterOp'
        remote_interop_dir=os.path.join(run_Id_for_searching,'InterOP')
        file_list = conn.listPath( share_folder_name, remote_interop_dir)
        for sh in file_list:
            if sh.isFile:
                remote_file_name=sh.filename
                copy_file=os.path.join(interop_local_dir_samba, remote_file_name)
                with open(copy_file ,'wb') as cp_fp :
                    remote_file=os.path.join(remote_interop_dir,)
                    conn.retrieveFile(share_folder_name, remote_file, cp_fp)
            # close samba connection 
        conn.close() 
        # parsing the files to get the xml Stats
        xml_stats=parsing_statistics_xml(demux_file, conversion_file)
        store_raw_xml_stats(xml_stats,run_be_processed_id)
        process_xml_stats(xml_stats,run_be_processed_id)
 

def find_state_and_save_data(run_name,run_folder):
    run_file='RunInfo.xml'
    run_parameter='RunParameters.xml'

    try:
        rn_found = RunProcess.objects.get(runName__exact=run_name)
    except:
        #os.chdir('wetlab/tmp/logs')
        with open('wetlab/tmp/logs/wetlab.log', 'a') as log_file:
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
        
def found_not_completed_run ():
    working_list={}
    state_list = ['Sample Sent','Process Running','Bcl2Fastq Executed']
    # get the run that are not completed
    for state in state_list:
        if RunProcess.objects.filter(runState__exact = state).exists():
            working_list[state].append(RunProcess.objects.filter(runState__exact = state))
    try:
        conn=open_samba_connection()
        print('Sucessful connection to Flavia \n')
        # clossing the connection after it was checked
        conn.close()
    except:
        return('Error')
    for state in state_list:
        if state in working_list:
            if state == 'Sample Sent':
                process_run_in_samplesent_state(working_list['Sample Sent'])
            elif state == 'Process Running':
                process_run_in_processrunning_state(working_list['Process Running'])
            else:
                process_run_in_bcl2F_q_executed_state(working_list['Bcl2Fastq Executed'])
    return (working_list)
    













#process_run_in_recorded_state()





print ('completed')

    
'''
    array_line=[[] for i in range(5)]
    
    for key, value in xml_statistics[project].items():
        print (key, value )
        array_line[0].append(key)
        count=1
        for val in value:
            array_line[count].append(val)
            count+= 1
    project_stats[project]=array_line
'''
    
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
'''
