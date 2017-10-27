#!/usr/bin/env python3

import sys, os, re
import xml.etree.ElementTree as ET
import time
import shutil
import locale

from  ..models import *
from .interop_statistics import *

from smb.SMBConnection import SMBConnection

def open_samba_connection():
    ## open samba connection
    # There will be some mechanism to capture userID, password, client_machine_name, server_name and server_ip
    # client_machine_name can be an arbitary ASCII string
    # server_name should match the remote machine name, or else the connection will be rejected

    conn=SMBConnection('Luigi', 'Apple123', 'NGS_Data_test', 'LUIGI-PC', use_ntlm_v2=True)
    conn.connect('192.168.1.3', 139)    
    #conn=SMBConnection('bioinfocifs', 'bioinfocifs', 'NGS_Data_test', 'barbarroja', use_ntlm_v2=True)
    #conn.connect('10.15.60.54', 139)
    '''
    conn = SMBConnection(userid, password, client_machine_name, remote_machine_name, use_ntlm_v2 = True)
    conn.connect(server_ip, 139)
    '''
    return conn


def save_run_info(run_info, run_parameter, run_id, logger):
    running_data={}
    image_channel=[]
    #################################################
    ## parsing RunInfo.xml file
    #################################################
    run_data=ET.parse(run_info)
    run_root=run_data.getroot()
    logger.info('Processing the runInfo.xml file')
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
    logger.info('Processing the runParameter.xml file')
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
    
    logger.debug('running_data information', running_data)
    ###########################################
    ## saving data into database
    ###########################################
    logger.info ('Saving to database  the run id ', run_id)
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
    

def fetch_exp_name_from_run_info (local_run_parameter_file):

    ## look for   <ExperimentName>NextSeq_CNM_041</ExperimentName> in RunParameters.xml file
    fh=open(local_run_parameter_file ,'r')
    for line in fh:
        exp_name=re.search('^\s+<ExperimentName>(.*)</ExperimentName>',line)
        if exp_name:
            fh.close()
            return exp_name.group(1)


           
def process_run_in_recorded_state(logger):
    try:
        conn=open_samba_connection()
        logger.info('Sucessfully  SAMBA connection for the process_run_in_recorded_state')
    except:
        return ('Error')
    processed_run_file, runlist = [] , []
    share_folder_name='NGS_Data_test'
    base_directory = 'documents/wetlab/tmp'
    recorded_dir = os.path.join(base_directory, 'recorded')

    local_run_parameter_file = os.path.join(base_directory, 'RunParameters.xml')
    local_run_info_file = os.path.join(base_directory, 'RunInfo.xml')
    process_run_file = os.path.join(base_directory, 'processed_run_file')
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
            if run_dir in processed_run:
                logger.debug('run id %s already processed', run_dir)
                continue
            else:
                #copy the runParameter.xml file to wetlab/tmp/tmp
                logger.info ('Found a new run  %s ,that was not in the processed run file',run_dir)
                with open(local_run_parameter_file ,'wb') as r_par_fp :
                    samba_run_parameters_file=os.path.join(run_dir,'RunParameters.xml')
                    conn.retrieveFile(share_folder_name, samba_run_parameters_file, r_par_fp)
                    logger.debug('Retrieving the remote RunParameter.xml file for %s', run_dir)
                    logger.debug('Local dir for getting RunParameters is %s', base_directory)
                    
                # look for the experience name  for the new run folder. Then find the run_id valued for it
                exp_name=fetch_exp_name_from_run_info(local_run_parameter_file)
                logger.debug('found the experiment name  , %s', exp_name)
                if  RunProcess.objects.filter(runName__icontains = exp_name, runState__exact = 'Recorded').exists():
                    exp_name_id=str(RunProcess.objects.get(runName__exact = exp_name).id)
                    logger.debug('Matching the experimental name id %s with database ', exp_name_id)
                    
                    sample_sheet_tmp_file=os.path.join(recorded_dir,exp_name_id,'samplesheet.csv')
                    if os.path.exists(sample_sheet_tmp_file):
                        # copy Sample heet file to samba directory
                        logger.info('Found run directory %s for the experiment name %s', run_dir, exp_name_id)
                        with open(sample_sheet_tmp_file ,'rb') as  sample_samba_fp:
                            samba_sample_file= os.path.join(run_dir,'samplesheet.csv')
                            logger.debug('Local dir for Shample Sheet %s', sample_sheet_tmp_file)
                            logger.debug('Remote dir to copy Shample Sheet  is %s', run_dir)
                            conn.storeFile(share_folder_name, samba_sample_file, sample_samba_fp)
                            logger.info('Samplesheet.csv file has been copied on the remote node')
                        # retrieve the runInfo.xml file from samba directory
                    else:
                        logger.error('ERROR ---No sample Sheet will be copied to remote dir. Local Directory %s was not found ', sample_sheet_tmp_dir)
                    # get the runIfnfo.xml to collect the  information for this run    
                    with open(local_run_info_file ,'wb') as r_info_fp :
                        samba_run_info_file=os.path.join(run_dir,'RunInfo.xml')
                        logger.debug('Local dir for RunInfo.xml %s', local_run_info_file)
                        logger.debug('Remote file of RunInfo.xml is located in %s', samba_run_info_file)
                        conn.retrieveFile(share_folder_name, samba_run_info_file, r_info_fp)
                    logger.info('copied to local the RunInfo.xml and start the parsing for RunInfo and RunParameter files')
                    save_run_info (local_run_info_file, local_run_parameter_file, exp_name_id, logger)
                        # delete the copy of the run files 
                    os.remove(local_run_info_file)
                    os.remove(local_run_parameter_file)
                    logger.debug('Deleted runInfo and RunParameter files on local server')
                        # delete the file and folder for the sample sheet processed
                    shutil.rmtree(os.path.join(recorded_dir, exp_name_id))
                    logger.debug('Deleted the recorded folder ')
                        # change the run  to SampleSent state
                    update_run_state(exp_name_id, 'Sample Sent', logger)
                    update_project_state(exp_name_id, 'Sample Sent', logger)
                        # add the run_dir inside the processed_run file
                    processed_run.append(run_dir)
                    run_names_processed.append(exp_name)
                    logger.info('SampleSheet.csv file for %s has been sucessful sent to remote server', run_dir)
                else:
                    # error in log file
                    logger.warn ('The run ID %s does not match any run in the RunProcess object.', run_dir)
                    os.remove(local_run_parameter_file) 
                    continue
    conn.close()
    fh =open (process_run_file,'w')
    # update the process_run_file with the new run 
    for process in processed_run:
        fh.write(process)
        fh.write('\n')
    fh.close()
    # check if all run process file are handled
 
    list_dir=os.listdir(recorded_dir)
    if list_dir:
        print ('Warning: There are runs in Recorded state that were not matching on the remote server \n')
        print ('These runs are waiting that the folder exists on the remote server:  \n')
        for item in list_dir:
            print('directory for the run Id : ', item, '\n')
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

def parsing_statistics_xml(demux_file, conversion_file, logger):
    total_p_b_count=[0,0,0,0] 
    stats_result={}

    demux_stat=ET.parse(demux_file)
    root=demux_stat.getroot()
    projects=[]
    logger.info('Starting conversion for demux file')
    for child in root.iter('Project'):
        projects.append(child.attrib['name'])
    total_samples = 0
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
             for  fill in range(4):
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
        if projects[i] != 'default':
            total_samples += len(samples)
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
    for project in stats_projects:
        if project == 'TopUnknownBarcodes':
            continue
        logger.info('processing project %s with rund_id = %s', project, run_id)
        if project == 'all' or project == 'default':
            logger.info('Found project %s setting the project_id to NULL', project)
            project_id= None
            default_all = project
        else:
            p_name_id=Projects.objects.get(projectName__exact = project).id
            project_id= Projects.objects.get(pk=p_name_id)
            default_all = None
           
        raw_stats_xml = RawStatisticsXml (runprocess_id=RunProcess.objects.get(pk=run_id),
                                          project_id = project_id, defaultAll = default_all,
                                          rawYield= stats_projects[project]['RAW_Yield'], rawYieldQ30= stats_projects[project]['RAW_YieldQ30'],
                                          rawQuality= stats_projects[project]['RAW_QualityScore'], PF_Yield= stats_projects[project]['PF_Yield'],
                                          PF_YieldQ30= stats_projects[project]['PF_YieldQ30'], PF_QualityScore =stats_projects[project]['PF_QualityScore'])
        
        logger.info('saving raw stats for %s project', project)
        raw_stats_xml.save()
    logger.info('Raw XML data have been stored for all projects ')
    
    
def process_xml_stats(stats_projects, run_id, logger):
    # get the total number of read per lane
    M_BASE=1.004361/1000000
    logger.debug('starting the process_xml_stats method')
    total_cluster_lane=(stats_projects['all']['PerfectBarcodeCount'])
    logger.info('processing flowcell stats for %s ', run_id)
    for project in stats_projects:
        if project == 'TopUnknownBarcodes':
            continue
        flow_raw_cluster, flow_pf_cluster, flow_yield_mb = 0, 0, 0
        for fl_item in range(4):
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
        
        for i in range (4):
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
            for un_lane in range(4) :
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
                    

def parsing_sample_project_xml(demux_file, conversion_file, logger):
    total_p_b_count=[0,0,0,0] 
    sample_result_dict={}
    #demux_file='example.xml'
    demux_stat=ET.parse(demux_file)
    root=demux_stat.getroot()
    projects=[]
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
            
            for l_index in range(4):
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
        logger.info('Stored sample for the project $s', project)

        #store in database
        logger.info('Processed information sample %s', project)


                                           
def process_run_in_samplesent_state (process_list, logger):
     # prepare a dictionary with key as run_name and value the RunID
     processed_run=[]
     for run_item in process_list:
        logger.info ('running the process sample sent state for %s', run_item)
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
        logger.info('check the Sucessful connection to NGS_Data_test before starting processing runing state method')

    except:
        return('Error')
    
    share_folder_name='NGS_Data_test'
    for run_item in process_list:
        logger.debug ('processing the run %s in process running state' , run_item)
        run_be_processed_id=RunProcess.objects.get(runName__exact=run_item).id
        run_Id_used=str(RunningParameters.objects.get(runName_id= run_be_processed_id))
        logger.debug ('found the run ID  %s' , run_Id_used )
        run_folder=os.path.join('/',run_Id_used,'Data/Intensities/BaseCalls')
        # check if runCompletion is avalilable
        file_list = conn.listPath( share_folder_name, run_folder)
        for sh in file_list:
            if sh.filename =='Reports' :
                logger.info('bcl2fastq has been completed for run %s', run_Id_used)
                processed_run.append(run_Id_used)
                update_run_state(run_be_processed_id, 'Bcl2Fastq Executed', logger)
                update_project_state(run_be_processed_id, 'B2FqExecuted', logger)
                break
            else:
                logger.debug('Report directory not found in file_list %s ', sh.filename)
        processed_run.append(run_Id_used)
            
    # close samba connection 
    conn.close()
    return processed_run



def process_run_in_bcl2F_q_executed_state (process_list, logger):
    processed_run=[]
    logger.debug('Executing process_run_in_bcl2F_q_executed_state method')
    for run_item in process_list:
        # change the state to Running Stats
        logger.info ('Processing the process on bcl2F_q for the run %s', run_item)
        run_processing_id=RunProcess.objects.get(runName__exact=run_item).id
        run_Id_used=str(RunningParameters.objects.get(runName_id= run_processing_id))
        plot_dir='documents/wetlab/images_plot'
        
        update_run_state(run_processing_id, 'Running Stats', logger)
        
        # get the directory of samba to fetch the files
        share_folder_name ='NGS_Data_test'
        local_dir_samba= 'documents/wetlab/tmp/processing'
        demux_file=os.path.join(local_dir_samba,'DemultiplexingStats.xml')
        conversion_file=os.path.join(local_dir_samba,'ConversionStats.xml')
        run_info_file=os.path.join(local_dir_samba, 'RunInfo.xml')
        #copy the demultiplexingStats.xml file to wetlab/tmp/processing
        
        try:
            conn=open_samba_connection()
            logger.info('Successful connection for updating run on bcl2F_q' )
        except:
            return 'Error'
        remote_stats_dir= 'Data/Intensities/BaseCalls/Stats/'
        samba_demux_file=os.path.join('/',run_Id_used,remote_stats_dir, 'DemultiplexingStats.xml')
        logger.debug('path to fetch demultiplexingStats is %s',  samba_demux_file)
        with open(demux_file ,'wb') as demux_fp :
            conn.retrieveFile(share_folder_name, samba_demux_file, demux_fp)
        logger.info('Fetched the DemultiplexingStats.xml')
        #copy the demultiplexingStats.xml file to wetlab/tmp/processing
        samba_conversion_file=os.path.join('/', run_Id_used,remote_stats_dir,'ConversionStats.xml')
        with open(conversion_file ,'wb') as conv_fp :
            conn.retrieveFile(share_folder_name, samba_conversion_file, conv_fp)
        logger.info('Fetched the conversionStats.xml')
        # copy RunInfo.xml  file to process the interop files
        with open(run_info_file ,'wb') as runinfo_fp :
            samba_conversion_file=os.path.join('/', run_Id_used,'RunInfo.xml')
                #conn.retrieveFile('share', '/path/to/remote_file', fp)
            conn.retrieveFile(share_folder_name, samba_conversion_file, runinfo_fp)
        logger.info('Fetched the RunInfo.xml')
        # copy all binary files in interop folder to local  documents/wetlab/tmp/processing/interop  
        interop_local_dir_samba= os.path.join(local_dir_samba, 'InterOp')
        remote_interop_dir=os.path.join('/',run_Id_used,'InterOp')
        file_list = conn.listPath( share_folder_name, remote_interop_dir)
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
                    return ('Error')
            # close samba connection 
        conn.close()
        
        # parsing the files to get the xml Stats
        logger.info('processing the XML files')
        xml_stats=parsing_statistics_xml(demux_file, conversion_file, logger)
        store_raw_xml_stats(xml_stats,run_processing_id, logger)
        process_xml_stats(xml_stats,run_processing_id, logger)
        
        # parsing and processing the project samples
        sample_project_stats = parsing_sample_project_xml (demux_file, conversion_file, logger)
        store_samples_projects (sample_project_stats, run_processing_id, logger)
        
        logger.info('processing interop files')
        # processing information for the interop files

        process_binStats(local_dir_samba, run_processing_id, logger)
        graphic_dir=os.path.join(plot_dir,run_Id_used)
        create_graphics(local_dir_samba, run_processing_id, run_Id_used, logger)
        processed_run.append(run_Id_used)
        logger.info('run id %s is now on Completed state', run_Id_used)
        update_run_state(run_processing_id, 'Completed', logger)
        update_project_state(run_processing_id, 'Completed', logger)
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
        logger.info ('Start processing the run founding runs for state %s', state)
        if state == 'Sample Sent':
            logger.debug ('found sample sent in state %s ')
            processed_run[state]=process_run_in_samplesent_state(working_list['Sample Sent'], logger)
        elif state == 'Process Running':
            logger.debug('Found runs for Process running %s', working_list['Process Running'])
            processed_run[state]=process_run_in_processrunning_state(working_list['Process Running'], logger)
            
        else:
            logger.debug('Found runs for Bcl2Fastq %s', working_list['Bcl2Fastq Executed'])
            processed_run[state]=process_run_in_bcl2F_q_executed_state(working_list['Bcl2Fastq Executed'], logger)
    
    return (processed_run)
    
print( 'executing the parsing_run_info.py')


    
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
