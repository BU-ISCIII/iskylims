import os
from iSkyLIMS_wetlab import wetlab_config
from iSkyLIMS_wetlab.models import *
from .generic_functions import *
from .run_metric_functions import *
import xml.etree.ElementTree as ET

'''
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

'''
'''
def process_run_in_bcl2Fq_executed_state (process_list, logger):
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



def get_bcl2fastq_output_files (conn, run_folder):
    '''
    Description:
        The function will collect files created by bcl2fastq process.
    Input:
        conn # Connection samba object
        run_folder   # folder run to fetch the remote files
    Constant:
        CONVERSION_STATS_FILE
        DEMULTIPLEXION_STATS_FILE
        REPORT_FOLDER
        RUN_TEMP_DIRECTORY
        STATISTICS_FOLDER
    Variables:
        l_conversion # local copy of ConversionStats file
        l_demux      # local copy of DemultiplexingStats file   
        l_metric_folder # local folder to copy the run metrics files
        l_run_info  # local copy of runInfo file
        
        run_folder      # run folder on the remote server
        s_conversion_stats # path for the conversionStats file
        s_conversion # path of the ConversionStats file
        s_demux      # path of the DemultiplexingStats file   
        stats_bcl2fastq_folder # statistics folder on the remote server 
                            for the bcl2fastq
    Return:
        copied_files
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function get_bcl2fastq_output_files')
    # Fetch Demuxtiplexing and conversion files 
    stats_bcl2fastq_folder = os.path.join(wetlab_config.DEMULTIPLEXION_BCL2FASTQ_FOLDER, wetlab_config.STATS_FOLDER)
    # conversion stats file
    l_conversion = os.path.join(wetlab_config.RUN_TEMP_DIRECTORY, wetlab_config.CONVERSION_STATS_FILE)
    s_conversion = os.path.join(wetlab_config.SAMBA_APPLICATION_FOLDER_NAME, run_folder, stats_bcl2fastq_folder, wetlab_config.CONVERSION_STATS_FILE)
    # demultiplexion stats file
    l_demux = os.path.join(wetlab_config.RUN_TEMP_DIRECTORY, wetlab_config.DEMULTIPLEXION_STATS_FILE)
    s_demux = os.path.join(wetlab_config.SAMBA_APPLICATION_FOLDER_NAME, run_folder, stats_bcl2fastq_folder, wetlab_config.DEMULTIPLEXION_STATS_FILE)
    
    try:
        l_conversion = fetch_remote_file (conn, run_folder, s_conversion, l_conversion)
        logger.info('Sucessfully fetch of ConversionStats file')
        
        l_demux = fetch_remote_file (conn, run_folder, s_demux, l_demux)
        logger.info('Sucessfully fetch of Demustiplexing file')

    except:
        string_message = "cannot copy files for getting run metrics"  
        logging_errors(logger,string_message, True , True)
        logger.info('Deleting temporary files')
        os.remove(l_conversion)
        logger.debug ('End function manage_run_in_processed_bcl2fast2_run with error')
        raise
    
    logger.debug ('End function get_bcl2fastq_output_files')
    return  l_demux, l_conversion


def check_run_metrics_processed (run_object_name) :
    '''
    Description:
        The function will check if exists the data stored on StatsRunSummary
        for this run 
        If exists returns True.  
    Input:
        run_object_name   # runProcess object
    Variables:
        run_metric_processed # True or False if there are some rows in 
                            StatsRunSummary for this run 
    Return:
        experiment_name if the run is updated. Empty if not
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function check_run_metrics_processed')
    run_metric_processed = StatsRunSummary.objects.filter(runprocess_id = run_object_name).exists()
    logger.info('Run metrics processed is %s', run_metric_processed)
    logger.debug('End function check_run_metrics_processed')
    return run_metric_processed

def cleanup_demux_tables_if_error (run_object_name):
    '''
    Description:
        The function will check if exists information stored on database
        for the tables : 
        -RawDemuxStats
        
        If they exist, they will be deleted to have a clean and valid 
        information in database
    Input:
        run_object_name   # runProcess object
    Variables:
        
    Return:
        True .
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function cleanup_demux_tables_if_error')
    
    if RawDemuxStats.objects.filter(runprocess_id__exact = run_object_name).exists():
        raw_stats_to_delete = RawDemuxStats.objects.filter(runprocess_id__exact = run_object_name)
        logger.info('Deleting RawDemuxStats rows')
        for raw_stats in raw_stats_to_delete :
            raw_stats.delete()
    projects = Projects.objects.filter(runprocess_id = run_object_name)
    for project in projects:
        if SamplesInProject.objects.filter(project_id = project).exists():
            samples_to_delete = SamplesInProject.objects.filter(project_id = project)
            logger.info('Deleting SampleProjects rows')
            for sample in samples_to_delete:
                sample.delete()

    if StatsFlSummary.objects.filter(runprocess_id = run_object_name).exists():
        fl_summary_to_delete = StatsFlSummary.objects.filter(runprocess_id = run_object_name)
        logger.info('Deleting StatsFlSummary rows')
        for fl_summary in fl_summary_to_delete:
            fl_summary.delete()
        
    if StatsLaneSummary.objects.filter(runprocess_id = run_object_name).exists():
        lane_summary_to_delete = StatsLaneSummary.objects.filter(runprocess_id = run_object_name)
        logger.info('Deleting StatsFlSummary rows')
        for lane_summary in lane_summary_to_delete:
            lane_summary.delete()
    
    logger.debug ('End function cleanup_demux_tables_if_error')
    return True

def cleanup_run_metrics_tables_if_error (run_object_name):
    '''
    Description:
        The function will check if exists information stored on database
        for the tables : 
        -StatsRunSummary
        
        If they exist, they will be deleted to have a clean and valid 
        information in database
    Input:
        run_object_name   # runProcess object
    Variables:
        
    Return:
        True .
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function cleanup_run_metrics_tables_if_error')
    if StatsRunSummary.objects.filter(runprocess_id__exact = run_object_name).exists():
        run_summary_to_delete = RawDemuxStats.objects.filter(runprocess_id__exact = run_object_name)
        logger.info('Deleting StatsRunSummary rows')
        for run_summary in run_summary_to_delete :
            run_summary.delete()
    if StatsRunRead.objects.filter(runprocess_id__exact = run_object_name).exists():
        run_read_to_delete = StatsRunRead.objects.filter(runprocess_id__exact = run_object_name)
        logger.info('Deleting StatsRunRead rows')
        for run_read in run_read_to_delete :
            run_read.delete()
    logger.debug ('End function cleanup_run_metrics_tables_if_error')
    return True

def parsing_demux_and_conversion_files( demux_file, conversion_file, number_of_lanes):
    '''
    Description:
        The function will parse demultiplexing and the conversion files
        that are created during the bcl2fastq process
    Input:
        demux_file   # local copy of the demultiplexting file
        conversion_file # local copy of the conversion file
    Functions:
        get_machine_lanes
    Variables:
        barcodeCount # value of the barcodeCount fetched in the parsing 
        one_mismatch_count # count of number of one mismatch
        p_temp      # pointer in the parse object to position it in each
                    project xml structure
        parsed result  # dictionnary where the parsing information is collected
        perfectBarcodeCount  # value of the perfectBarcodeCount fetched 
                            in the parsing 
        projects
        project_parsed_information # dictionnary to fetch the parsed information
                                    for each project included in the file
        total_p_b_count
        samples         # number of the samples in the run 
        sample_all_index # 
        
        
        run_metric_processed # True or False if there are some rows in 
                            StatsRunSummary for this run 
    Return:
        parsed_result
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function parsing_demux_and_conversion_files')
    
    total_p_b_count = [0,0,0,0]
    parsed_result={}
    projects=[]
    demux_parse = ET.parse(demux_file)
    root=demux_parse.getroot()
    
    logger.info('Start parsing demux file')
    for child in root.iter('Project'):
        projects.append(child.attrib['name'])
    total_samples = 0
    #number_of_lanes = get_machine_lanes(run_id)
    for i in range(len(projects)):
        p_temp=root[0][i]
        barcodeCount ,perfectBarcodeCount, b_count =[], [] ,[]
        p_b_count, one_mismatch_count =[], []
        project_parsed_information = {}
        
        samples=p_temp.findall('Sample')
        sample_all_index=len(samples)-1

        for c in p_temp[sample_all_index].iter('BarcodeCount'):
            barcodeCount.append(c.text)
            
        for c in p_temp[sample_all_index].iter('PerfectBarcodeCount'):
            p_b_count.append(c.text)

        # look for One mismatch barcode
        if p_temp[sample_all_index].find('OneMismatchBarcodeCount') == None:
             for  fill in range(number_of_lanes):
                one_mismatch_count.append('NaN')
        else:
            for c in p_temp[sample_all_index].iter('OneMismatchBarcodeCount'):
                one_mismatch_count.append(c.text)

        #one_mismatch_count.append(one_m_count)
        project_parsed_information['BarcodeCount'] = barcodeCount
        project_parsed_information['PerfectBarcodeCount'] = p_b_count
        project_parsed_information['sampleNumber'] = len(samples) -1
        project_parsed_information['OneMismatchBarcodeCount'] = one_mismatch_count
        parsed_result[projects[i]] = project_parsed_information
        if projects[i] != 'default' and projects[i] != 'all':
            total_samples += len(samples) -1
        logger.info('Completed parsing from demux file for project %s', projects[i])
    # overwrite the value for total samples
    parsed_result['all']['sampleNumber']=total_samples

    conversion_stat=ET.parse(conversion_file)
    root_conv=conversion_stat.getroot()
    projects=[]
    logger.info('Starting parsing for conversion file')
    for child in root_conv.iter('Project'):
        projects.append(child.attrib['name'])
    for i in range(len(projects)):
        p_temp=root_conv[0][i]
        samples=p_temp.findall('Sample')
        sample_all_index=len(samples)-1
        tiles=p_temp[sample_all_index][0][0].findall('Tile')
        tiles_index=len(tiles)-1
        list_raw_yield , list_raw_yield_q30 = [] , []
        list_raw_qualityscore , list_pf_yield = [] , []
        list_pf_yield_q30, list_pf_qualityscore = [], []

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

        parsed_result[projects[i]]['RAW_Yield']=list_raw_yield
        parsed_result[projects[i]]['RAW_YieldQ30']=list_raw_yield_q30
        parsed_result[projects[i]]['RAW_QualityScore']=list_raw_qualityscore
        parsed_result[projects[i]]['PF_Yield']=list_pf_yield
        parsed_result[projects[i]]['PF_YieldQ30']=list_pf_yield_q30
        parsed_result[projects[i]]['PF_QualityScore']=list_pf_qualityscore
        logger.info('Completed parsing for conversion stats for project %s', projects[i])

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
    parsed_result['TopUnknownBarcodes']= unknow_lanes
    logger.debug ('End function parsing_demux_and_conversion_files')

    return parsed_result


def parsing_demux_sample_project(demux_file, conversion_file, number_of_lanes):
    '''
    Description:
        The function will parse demultiplexing and the conversion files
        to get the samples inside the project
    Input:
        demux_file   # local copy of the demultiplexting file
        conversion_file # local copy of the conversion file
        number_of_lanes # numer of lane to be fetched, based on the sequencer type
    Functions:
        
    Variables:
        barcodeCount # value of the barcodeCount fetched in the parsing 
        one_mismatch_count # count of number of one mismatch
        p_temp      # pointer in the parse object to position it in each
                    project xml structure
        parsed result  # dictionnary where the parsing information is collected
        perfectBarcodeCount  # value of the perfectBarcodeCount fetched 
                            in the parsing 
        projects
        project_parsed_information # dictionnary to fetch the parsed information
                                    for each project included in the file
        total_p_b_count
        samples         # number of the samples in the run 
        sample_all_index # 
        
        
        run_metric_processed # True or False if there are some rows in 
                            StatsRunSummary for this run 
    Return:
        parsed_result
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function parsing_demux_and_conversion_files')
    
    total_p_b_count=[0,0,0,0]
    parsed_result = {}
    demux_stat=ET.parse(demux_file)
    root=demux_stat.getroot()
    projects=[]
    #number_of_lanes=get_machine_lanes(run_id)
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

            parsed_result[projects[i]]=sample_dict
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

            parsed_result[projects[i]][sample_name]['RAW_Yield']=raw_yield_value
            parsed_result[projects[i]][sample_name]['RAW_YieldQ30']=raw_yield_q30_value
            parsed_result[projects[i]][sample_name]['RAW_QualityScore']=raw_quality_value
            parsed_result[projects[i]][sample_name]['PF_Yield']=pf_yield_value
            parsed_result[projects[i]][sample_name]['PF_YieldQ30']=pf_yield_q30_value
            parsed_result[projects[i]][sample_name]['PF_QualityScore']=pf_quality_value
        logger.info('completed parsing for xml stats for project %s', projects[i])


    logger.info('Complete XML parsing  for getting Samples')

    return parsed_result


def process_fl_summary_stats (stats_projects, run_object_name):
    '''
    Description:
        The function will process the parsed data, to prepare the 
        Flowcell summary information to store in database 
    Input:
        stats_projects     # dictionnary with parsed data to be processed
        run_object_name   # runProcess object
    Constants:
        M_BASE      # to have values in MB 
    Variables:
        processed_fl_summary_data # dictionary list that contains all processed
                            information.
        project_flowcell    # dictionnary to collect the processed data.
                            After collecting it will be appended to 
                            processed_fl_summary_data.
    Return:
        processed_fl_summary_data.
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function process_fl_summary_stats')
    M_BASE=1.004361/1000000
    total_cluster_lane=(stats_projects['all']['PerfectBarcodeCount'])
    processed_fl_summary_data = []
    number_of_lanes = run_object_name.get_machine_lanes()
    for project in stats_projects:
        project_flowcell = {}
        if project == 'TopUnknownBarcodes':
            continue

        logger.info('Start processing flow Summary for project %s', project)
        flow_raw_cluster, flow_pf_cluster, flow_yield_mb = 0, 0, 0
        for fl_item in range(number_of_lanes):
             # make the calculation for Flowcell
            flow_raw_cluster +=int(stats_projects[project]['BarcodeCount'][fl_item])
            flow_pf_cluster +=int(stats_projects[project]['PerfectBarcodeCount'][fl_item])
            flow_yield_mb +=float(stats_projects[project]['PF_Yield'][fl_item])*M_BASE

        project_flowcell['flowRawCluster'] = '{0:,}'.format(flow_raw_cluster)
        project_flowcell['flowPfCluster'] ='{0:,}'.format(flow_pf_cluster)
        project_flowcell['flowYieldMb'] = '{0:,}'.format(round(flow_yield_mb))
        project_flowcell['sampleNumber'] = stats_projects[project]['sampleNumber']

        if project == 'all' or project == 'default' :
            project_flowcell['project_id'] = None
            project_flowcell['defaultAll'] = project
        else:
            project_flowcell['project_id'] = Projects.objects.get(projectName__exact = project)
            project_flowcell['defaultAll'] = None
        project_flowcell['runprocess_id'] = run_object_name
        #store in database
        logger.info('End processing flow Summary for project %s', project)
        '''
        ns_fl_summary = StatsFlSummary(runprocess_id=RunProcess.objects.get(pk=run_id),
                                project_id=project_id, defaultAll=default_all, flowRawCluster=flow_raw_cluster,
                                flowPfCluster=flow_pf_cluster, flowYieldMb= flow_yield_mb,
                                sampleNumber= sample_number)


        ns_fl_summary.save()
        '''
        processed_fl_summary_data.append(project_flowcell)
    
    logger.debug ('End function process_fl_summary_stats')
    return processed_fl_summary_data

def process_lane_summary_stats (stats_projects, run_object_name):
    '''
    Description:
        The function will process the parsed data, to prepare the 
        Flowcell summary information to store in database 
    Input:
        stats_projects     # dictionnary with parsed data to be processed
        run_object_name   # runProcess object
    Constants:
        M_BASE      # to have values in MB 
    Variables:
        processed_lane_summary_data # dictionary list that contains all processed
                            information.
        project_lane      # dictionnary to collect the processed data.
                            After collecting it will be appended to 
                            processed_fl_summary_data.
    Return:
        processed_lane_summary_data.
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function process_lane_summary_stats')
    M_BASE=1.004361/1000000
    processed_lane_summary_data = []
    number_of_lanes = run_object_name.get_machine_lanes()
    total_cluster_lane=(stats_projects['all']['PerfectBarcodeCount'])
    for project in stats_projects:
        if project == 'TopUnknownBarcodes':
            continue
        logger.info('processing lane stats for %s', project)
        for i in range (number_of_lanes):
            # get the lane information
            project_lane = {}
            project_lane['lane'] = str(i + 1)
            pf_cluster_int=(int(stats_projects[project]['PerfectBarcodeCount'][i]))
            project_lane['pfCluster'] = '{0:,}'.format(pf_cluster_int)
            project_lane['perfectBarcode'] = (format(int(stats_projects[project]['PerfectBarcodeCount'][i])*100/int(stats_projects[project]['BarcodeCount'][i]),'.3f'))
            project_lane['percentLane'] = format(float(int(pf_cluster_int)/int(total_cluster_lane[i]))*100, '.3f')
            project_lane['oneMismatch'] = stats_projects[project]['OneMismatchBarcodeCount'][i]
            project_lane['yieldMb'] = '{0:,}'.format(round(float(stats_projects[project]['PF_Yield'][i])* M_BASE))
            project_lane['biggerQ30'] = format(float(stats_projects[project]['PF_YieldQ30'][i])*100/float( stats_projects[project]['PF_Yield'][i]),'.3f')
            project_lane['meanQuality'] = format(float(stats_projects[project]['PF_QualityScore'][i])/float(stats_projects[project]['PF_Yield'][i]),'.3f')
            
            if project == 'all' or project == 'default':
                project_lane['project_id'] = None
                project_lane['defaultAll'] = project
            else:
                project_lane['project_id'] = Projects.objects.get(projectName__exact = project)
                project_lane['defaultAll'] =  None

            project_lane['runprocess_id'] = run_object_name
            '''
            ns_lane_summary = NextSeqStatsLaneSummary(runprocess_id=RunProcess.objects.get(pk=run_id),
                                                 project_id=project_id, defaultAll = default_all, lane = lane_number,
                                                 pfCluster=pf_cluster, percentLane=percent_lane, perfectBarcode=perfect_barcode,
                                                 oneMismatch= one_mismatch, yieldMb=yield_mb,
                                                 biggerQ30=bigger_q30, meanQuality=mean_quality )

            ns_lane_summary.save()

            '''
            processed_lane_summary_data.append(project_lane)
        logger.info('Processed information for project %s', project)
    logger.debug ('End function process_lane_summary_stats')
    return processed_lane_summary_data


def process_raw_demux_stats(stats_projects, run_object_name):
    '''
    Description:
        The function will process the parsed data, fetched previously and
        it will prepare them to save later in database.
        All the information is saved in a list "processed_raw_data" that
        contains a dictionnary with the values to stored on RawDemuxStats

        If all data are sucessfuly processed it returns processed_raw_data.
        An exeception is raised in case of error while processing
    Input:
        stats_projects     # dictionnary with parsed data to be processed
        run_object_name   # runProcess object
    Variables:
        processed_raw_data # dictionary list that contains all processed
                            information.
        project_list    # list of project for checking if all of them are
                        defined on database
        project_raw_data    # dictionnary to collect the processed data.
                            After collecting it will be appended to 
                            processed_raw_data.
    Return:
        processed_raw_data if information is correctly processed.
        An exception is raise if any project is not defined on database
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function store_raw_xml_stats')
    processed_raw_data = []
    #error_found = False
    project_list = []
    # check first if all project found in the parsing are already defined
    # when sample sheet was processed. 
    for project in stats_projects:
        if project == 'TopUnknownBarcodes'  or project == 'all' or project == 'default':
            continue
        else:
            project_list.append(project)

    if not check_all_projects_exists (project_list) :
        # if there some projects not defined set the run to error state
        string_message = "Found some projects not defined in database"
        logging_errors(logger, string_message, False, True)
        experiment_name = run_object_name.get_run_name()
        handling_errors_in_run(experiment_name, '10')
        raise
    for project in stats_projects:
        project_raw_data = {}
        # ignore on this step the unknowbarcode parsed information
        if project == 'TopUnknownBarcodes':
            continue
        logger.info('processing project %s ', project)
        if project == 'all' or project == 'default':
            logger.info('Found project %s setting the project_id to NULL', project)
            project_raw_data['project_id']= None
            project_raw_data['defaultAll'] = project
        else:
            project_raw_data['project_id']= Projects.objects.get(projectName__exact = project)
            project_raw_data['defaultAll'] = None
        
        project_raw_data['rawYield'] = stats_projects[project]['RAW_Yield']
        project_raw_data['rawYieldQ30'] =  stats_projects[project]['RAW_YieldQ30']
        project_raw_data['rawQuality'] = stats_projects[project]['RAW_QualityScore']
        project_raw_data['PF_Yield'] = stats_projects[project]['PF_Yield']
        project_raw_data['PF_YieldQ30'] = stats_projects[project]['PF_YieldQ30']
        project_raw_data['PF_QualityScore'] = stats_projects[project]['PF_QualityScore']
        # append project_raw_data to list
        processed_raw_data.append(project_raw_data)
        '''
        # delete all information stored on database
        if RawStatisticsXml.objects.filter (runprocess_id__exact = run_id).exists():
            projects_to_delete = RawStatisticsXml.objects.filter (runprocess_id__exact = run_id)
            logger.info('Deleting stored RawStatisticsXml information ')
            for project in projects_to_delete :
                project.delete()
                logger.debug('Deleted the RawStatisticsXml for project %s', project)
        return 'ERROR'
        '''
    logger.info('End function store_raw_xml_stats')
    return processed_raw_data


def process_samples_projects(stats_projects, run_object_name ):
    '''
    Description:
        The function will parse demultiplexing and the conversion files
        that are created during the bcl2fastq process
    Input:
        stats_projects     # dictionnary with parsed data to be processed
        run_object_name   # runProcess object
    Variables:
        processed_sample_data # dictionary list that contains all processed
                            information.
        project_list    # list of project for checking if all of them are
                        defined on database
        project_sample_data    # dictionnary to collect the processed data.
                            After collecting it will be appended to 
                            processed_raw_data.
    Return:
        processed_raw_data if information is correctly processed.
        An exception is raise if any project is not defined on database
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function process_samples_projects')
    # get the total number of read per lane
    M_BASE=1.004361/1000000
    processed_sample_data =[]
    
    for project in stats_projects:
        # find the total number of PerfectBarcodeCount in the procjec to make percent calculations
        total_perfect_barcode_count = 0
        for sample in stats_projects[project]:
            total_perfect_barcode_count += stats_projects[project][sample] ['PerfectBarcodeCount']

        for sample in stats_projects[project]:
            project_sample_data = {}
            project_sample_data['sampleName'] = sample
            project_sample_data['barcodeName'] = stats_projects[project][sample]['barcodeName']
            perfect_barcode = int(stats_projects[project][sample] ['PerfectBarcodeCount'])
            project_sample_data['pfClusters'] =  '{0:,}'.format(perfect_barcode)
            project_sample_data['percentInProject'] = format (float(perfect_barcode) *100 /total_perfect_barcode_count,'.2f')

            project_sample_data['yieldMb'] = '{0:,}'.format(round(float(stats_projects[project][sample] ['PF_Yield'])*M_BASE))
            if stats_projects[project][sample] ['PF_Yield'] > 0:
                bigger_q30=format(float(stats_projects[project][sample]['PF_YieldQ30'])*100/float( stats_projects[project][sample]['PF_Yield']),'.3f')
                mean_quality=format(float(stats_projects[project][sample]['PF_QualityScore'])/float(stats_projects[project][sample]['PF_Yield']),'.3f')
            else:
                bigger_q30 = 0
                mean_quality =0

            project_sample_data['qualityQ30'] = bigger_q30
            project_sample_data['meanQuality'] = mean_quality
            project_sample_data['project_id'] = Projects.objects.get(projectName__exact = project)
            # p_name_id=Projects.objects.get(projectName__exact = project).id
            # project_id= Projects.objects.get(pk=p_name_id)
            '''
            sample_to_store = SamplesInProject (project_id = project_id, sampleName = sample_name,
                            barcodeName = barcode_name, pfClusters = perfect_barcode,
                            percentInProject = percent_in_project , yieldMb = yield_mb,
                            qualityQ30 = bigger_q30, meanQuality = mean_quality)

            sample_to_store.save()
            '''
            processed_sample_data.append(project_sample_data)
        logger.info('Stored sample for the project %s', project)
    logger.debug('End function process_samples_projects')
    return processed_sample_data


def process_unknow_barcode_stats (stats_projects, run_object_name):
    '''
    Description:
        The function will process the parsed data, to prepare the 
        Flowcell summary information to store in database 
    Input:
        stats_projects     # dictionnary with parsed data to be processed
        run_object_name   # runProcess object
    Constants:
        M_BASE      # to have values in MB 
    Variables:
        processed_barcode_data # dictionary list that contains all processed
                            information.
        unknow_barcode    # dictionnary to collect the processed data.
                            After collecting it will be appended to 
                            processed_fl_summary_data.
    Return:
        processed_barcode_data.
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function process_unknow_barcode_stats')
    number_of_lanes = run_object_name.get_machine_lanes()
    processed_barcode_data = []
    for project in stats_projects:
        if project == 'TopUnknownBarcodes':
            for un_lane in range(number_of_lanes) :
                logger.info('Processing lane %s for TopUnknownBarcodes', un_lane)
                count_top=0
                top_number =1
                for barcode_line in stats_projects[project][un_lane]:
                    unknow_barcode = {}
                    unknow_barcode['runprocess_id'] = run_object_name
                    unknow_barcode['lane_number'] = lane_number=str(un_lane + 1)
                    unknow_barcode['top_number'] = str(top_number)
                    unknow_barcode['count'] =  '{0:,}'.format(int(barcode_line['count']))
                    unknow_barcode['sequence'] =  barcode_line['sequence']
                    top_number +=1
                    '''
                    raw_unknow_barcode = RawTopUnknowBarcodes(runprocess_id=RunProcess.objects.get(pk=run_id),
                                                             lane_number = lane_number, top_number=str(top_number),
                                                             count=barcode_count, sequence=barcode_sequence)
                    '''
                    processed_barcode_data.append(unknow_barcode)

    logger.debug ('End function process_unknow_barcode_stats')
    return processed_barcode_data

def manage_run_in_processed_run (conn, run_object_name):
    '''
    Description:
        The function will check if exists the folder that are created
        during the bcl2fastq conversion 
        If exists the run state will move to Processing Bcl2fastq  
    Input:
        conn # Connectio samba object
        run_object_name   # runProcess object
    Functions:
        get_run_metric_files    # Located at utils.run_metics_functions
    Variables:
        experiment_name # Name of the run
        count_file_size # partial size for the subfolder
        run_folder      # run folder on the remote server
        statistics_folder # statistics folder on the remote server 
    Return:
        experiment_name if the run is updated. Empty if not
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function manage_run_in_processed_run')
    experiment_name = run_object_name.get_run_name()
    run_folder = RunningParameters.objects.get(runName_id = run_object_name).get_run_folder()
    statistics_folder = os.path.join(wetlab_config.SAMBA_APPLICATION_FOLDER_NAME, run_folder, wetlab_config.DEMULTIPLEXION_BCL2FASTQ_FOLDER)
    if 'Processed Run' == run_object_name.get_state() :
        if not check_run_metrics_processed (run_object_name) :
            run_state = run_object_name.set_run_state('Processing Metrics')
            run_id = run_object_name.get_run_id()
            try:
                # connect to statistics folder to fetch the run metrics files
                run_metric_files = get_run_metric_files (conn, run_folder)
                
                parsed_run_stats_summary, parsed_run_stats_read = parsing_run_metrics(wetlab_config.RUN_TEMP_DIRECTORY_PROCESSING, run_object_name)
            except:
                logger.debug ('End function manage_run_in_processed_run with error')
                raise

                  
            try:
                for run_stat_summary in parsed_run_stats_summary :
                    saved_run_stat_summary = StatsRunSummary.objects.create_stats_run_summary(run_stat_summary, experiment_name)
                logger.info('run metrics summary data saved to database')
                for run_stat_read in parsed_run_stats_read :
                    saved_run_stat_read = StatsRunRead.objects.create_stats_run_read(run_stat_read, experiment_name)
                logger.info('run metrics read data saved to database')
            except:
                string_message = 'Run metrics data cannot be saved'
                logging_errors(logger, string_message, True, True)
                # delete run metrics tables
                cleanup_run_metrics_tables_if_error (run_object_name)
                delete_run_metric_files (run_metric_files)
                logger.info('Deleted run metrics tables')
                handling_errors_in_run(experiment_name)
                logger.debug ('End function manage_run_in_processed_run with error')
                raise
            
            # create run graphics
            run_graphics = create_graphics (wetlab_config.RUN_TEMP_DIRECTORY_PROCESSING, run_object_name)
            logger.info('run metrics graphics processed and copied to plot image folder')
            # deleting temporary run metrics files
            delete_run_metric_files (run_metric_files)
            # return the state to Processed Run
            run_state = run_object_name.set_run_state('Processed Run')
            
            set_state_in_all_projects(experiment_name, 'Processed Run')
        # Check if Bcl2fastq has started
        try:  
            file_list = conn.listPath( wetlab_config.SAMBA_SHARED_FOLDER_NAME, statistics_folder)
        except: 
            logger.info('No statistics folder were found for %s', experiment_name) 
            logger.debug ('End function manage_run_in_processed_run')
            return experiment_name
        
        run_state = run_object_name.set_run_state('Processing Bcl2fastq')
        set_state_in_all_projects(experiment_name, 'Processing Bcl2fastq')
    
    else:
        string_message = 'Invalid state when calling to ' + experiment_name 
        logging_errors(logger, string_message , False , False)
        logger.debug ('End function manage_run_in_processed_run with error')
        return
    logger.debug ('End function manage_run_in_processed_run')
    return experiment_name


def manage_run_in_processing_bcl2fastq (conn, run_object_name):
    '''
    Description:
        The function will check if report floder exists. Then the bcl2fastq
        conversion is completed
        If exists the run state will move to Processed Bcl2fastq  
    Input:
        conn # Connectio samba object
        run_object_name   # runProcess object
    Constant:
        REPORT_FOLDER
    Functions:
        set_state_in_all_projects # located at utils.generic_functions
    Variables:
        experiment_name # Name of the run
        count_file_size # partial size for the subfolder
        run_folder      # run folder on the remote server
        s_conversion_stats # path for the conversionStats file
        statistics_folder # statistics folder on the remote server 
    Return:
        experiment_name if the run is updated. Empty if not
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function manage_run_in_processing_bcl2fast2_run')
    experiment_name = run_object_name.get_run_name()
    logger.info('Start handling run %s', experiment_name)
    run_folder = RunningParameters.objects.get(runName_id = run_object_name).get_run_folder()
    statistics_folder = os.path.join(wetlab_config.SAMBA_APPLICATION_FOLDER_NAME, run_folder, wetlab_config.DEMULTIPLEXION_BCL2FASTQ_FOLDER)
    
    
    if 'Processing Bcl2fastq' == run_object_name.get_state() :
        # connect to statistics report folder to check if bcl2fastq process is ended
        try:  
            file_list = conn.listPath( wetlab_config.SAMBA_SHARED_FOLDER_NAME, statistics_folder)
        except: 
            string_message = 'Folder statistics have been deleted for ' + experiment_name
            logging_errors(logger, string_message, False, False)
            logger.debug ('End function manage_run_in_processing_bcl2fast2 with error')
            return ''
        for sh in file_list:
            if sh.filename == wetlab_config.REPORT_FOLDER :
                
                logger.info('bcl2fastq has been completed for  %s', experiment_name)
                # Get the time when  the Bcl2Fastq process is ending
                
                s_conversion_stats = os.path.join (statistics_folder, wetlab_config.STATS_FOLDER, wetlab_config.CONVERSION_STATS_FILE)
                conversion_attributes = conn.getAttributes(wetlab_config.SAMBA_SHARED_FOLDER_NAME ,s_conversion_stats)
                bcl2fastq_finish_date = datetime.datetime.fromtimestamp(int(conversion_attributes.create_time)).strftime('%Y-%m-%d %H:%M:%S')
                bcl2fastq_finish_date = run_object_name.set_run_bcl2fastq_finished_date(bcl2fastq_finish_date)
                run_object_name.set_run_state('Processed Bcl2fastq')
                set_state_in_all_projects(experiment_name, 'Processed Bcl2fastq')

                logger.info ('Updated the Bcl2Fastq time in the run %s', experiment_name)
                break
         
    else:
        string_message = 'Invalid state when calling to ' + experiment_name 
        logging_errors(logger, string_message , False , False)
        logger.debug ('End function manage_run_in_processing_bcl2fast2 with error')
        return ''
    
    logger.debug ('End function manage_run_in_processing_bcl2fast2')
    return experiment_name


def manage_run_in_processed_bcl2fastq (conn, run_object_name):
    '''
    Description:
        The function will collect files created by bcl2fastq process.
        Saves on database after performing parsing and processing these
        data.
        When starting the function the run state is changed to Running Stats
        to prevent unconsistance presenting data when database is populated.
        After the sucessful executing of this function the run is set 
        to completed state
    Input:
        conn # Connection samba object
        run_object_name   # runProcess object
    Constant:
        REPORT_FOLDER
    Functions:
        get_bcl2fastq_output_files # located at this file
        parsing_demux_and_conversion_files
        parsing_demux_sample_project
        process_raw_demux_stats
        set_state_in_all_projects # located at utils.generic_functions
    Variables:
        experiment_name # Name of the run
        count_file_size # partial size for the subfolder
        l_metric_folder # local folder to copy the run metrics files
        copied_files    # dictionnary of the temporary folder to be removed
                            in case that any of the remote files could 
                            not be fetched
        parsed_result   # dictionnary with the parsed information from
                        demuxtiplexing and conversion files
        
        processed_samples_stats # list of dictionnary that contains the
                                processed sample information to save in 
                                SamplesInProject database 
        
        processed_raw_stats     list of dictionnary that contains the
                                processed demultiplexing information to 
                                save in RawDemuxStats database 
        run_folder      # run folder on the remote server
        run_metrics_file_name # name of the run metric file. The value
                            is updated for each of the run metric files
                            in the folder
        s_conversion_stats # path for the conversionStats file
        s_metric_folder     # run metric folder on the remote server
        statistics_folder # statistics folder on the remote server 

    Return:
        experiment_name if the run is sucessfully updated. 
        Exception if any error occurs 
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function manage_run_in_processed_bcl2fast2_run')
    experiment_name = run_object_name.get_run_name()
    run_folder = RunningParameters.objects.get(runName_id = run_object_name).get_run_folder()
    number_of_lanes = run_object_name.get_machine_lanes()
    #statistics_folder = os.path.join(run_folder, wetlab_config.CONVERSION_STATS_FOLDER)

    if 'Processed Bcl2fastq' == run_object_name.get_state() :
        run_state = run_object_name.set_run_state('Processing Demultiplexing')
        
        try:
            l_demux , l_conversion= get_bcl2fastq_output_files (conn, run_folder)
        except:
            string_message = 'Unable to fetch stats files for ' + experiment_name
            logging_errors(logger, string_message, False, False)
            handling_errors_in_run(experiment_name, '11')
            logger.debug ('End function manage_run_in_processed_bcl2fast2 with error')
            raise
        # parsing the files to get the xml Stats
        logger.info('Start parsing  demultiplexing files')
        parsed_result = parsing_demux_and_conversion_files(l_demux, l_conversion, number_of_lanes)
        
        # parsing and processing the project samples
        logger.info('Start parsing  samples demultiplexing')
        sample_parsed_result = parsing_demux_sample_project (l_demux, l_conversion, number_of_lanes)
        # clean up the fetched files in the local temporary folder
        os.remove(l_conversion)
        os.remove(l_demux)
        logger.info ('Deleted temporary demultiplexing and conversion files')
        
        processed_raw_stats = process_raw_demux_stats(parsed_result, run_object_name)
        try:
            for raw_stats in processed_raw_stats :
                new_raw_stats = RawDemuxStats.objects.create_stats_run_read(raw_stats, run_object_name)
            logger.info('Saved information to RawDemuxStats')
        except:
            string_message = 'Unable to save raw stats for ' + experiment_name
            logging_errors(logger, string_message, False, False)
            handling_errors_in_run (experiment_name, '11' )
            cleanup_demux_tables_if_error(run_object_name)
            logger.debug('End function manage_run_in_processed_bcl2fast2_run with error')
            raise
        
        processed_samples_stats = process_samples_projects(sample_parsed_result, run_object_name)
        try:
            for sample_stats in processed_samples_stats :
                new_sample_stats = SamplesInProject.objects.create_sample_project(sample_stats)
            logger.info('Saved information to SamplesInProject')
        except:
            string_message = 'Unable to save sample stats saving for ' + experiment_name
            logging_errors(logger, string_message, False, False)
            handling_errors_in_run (experiment_name, '13' )
            cleanup_demux_tables_if_error(run_object_name)
            logger.debug('End function manage_run_in_processed_bcl2fast2_run with error')
            raise
        processed_fl_summary = process_fl_summary_stats(parsed_result, run_object_name)
        try:
            for fl_summary in processed_fl_summary :
                new_fl_summary = StatsFlSummary.objects.create_fl_summary(fl_summary)
            logger.info('Saved information to StatsFlSummary')
        except:
            string_message = 'Unable to save FL Summary  stats for ' + experiment_name
            logging_errors(logger, string_message, False, False)
            handling_errors_in_run (experiment_name, '14' )
            cleanup_demux_tables_if_error(run_object_name)
            logger.debug('End function manage_run_in_processed_bcl2fast2_run with error')
            raise
        
        processed_lane_summary = process_lane_summary_stats(parsed_result, run_object_name)
        try:
            for lane_summary in processed_lane_summary :
                new_lane_summary = StatsLaneSummary.objects.create_lane_summary(lane_summary)
            logger.info('Saved information to StatsLaneSummary')
        except:
            string_message = 'Unable to save Lane Summary stats for ' + experiment_name
            logging_errors(logger, string_message, False, False)
            handling_errors_in_run (experiment_name, '15' )
            cleanup_demux_tables_if_error(run_object_name)
            logger.debug('End function manage_run_in_processed_bcl2fast2_run with error')
            raise
        
        processed_unknow_barcode = process_unknow_barcode_stats(parsed_result, run_object_name)
        try:
            for unknow_barcode in processed_unknow_barcode :
                new_unknow_barcode = RawTopUnknowBarcodes.objects.create_unknow_barcode(unknow_barcode)
            logger.info('Saved information to RawTopUnknowBarcodes')
        except:
            string_message = 'Unable to Unknow Barcode stats for ' + experiment_name
            logging_errors(logger, string_message, False, False)
            handling_errors_in_run (experiment_name, '16' )
            cleanup_demux_tables_if_error(run_object_name)
            logger.debug('End function manage_run_in_processed_bcl2fast2_run with error')
            raise
        
        ## Get the disk space utilization for this run
        import pdb; pdb.set_trace()
        try:
            disk_utilization = get_run_disk_utilization (conn, run_folder)
        except:
            string_message = 'Error when fetching the log file for the run ' + new_run
            logging_errors (logger, string_message, False, False)
            handling_errors_in_run (experiment_name, '17' )
            cleanup_demux_tables_if_error(run_object_name)
            logger.debug('End function manage_run_in_processed_bcl2fast2_run with error')
            raise 
        
        result_store_usage = run_object_name.set_used_space (disk_utilization)
        
        completion_date = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        result_set_completion_date = run_object_name.set_run_completion_date(completion_date)
        # Update the run state to completed
        run_state = run_object_name.set_run_state('Completed')
        projects_state = set_state_in_all_projects(experiment_name, 'Completed')

    else: 
        string_message = 'Invalid state when calling to ' + experiment_name 
        logging_errors(logger, string_message , False , False)
        logger.debug ('End function manage_run_in_processed_bcl2fast2 with error')
        return ''
    logger.debug ('End function manage_run_in_processed_bcl2fast2 ')
    return experiment_name

    
    
