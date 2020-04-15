import os
from iSkyLIMS_wetlab import wetlab_config
from iSkyLIMS_wetlab.models import *
from .generic_functions import *
from .run_metric_functions import *
import xml.etree.ElementTree as ET



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
        logging_errors(tring_message, True , True)
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
    experiment_name = run_object_name.get_run_name()
    logger.debug ('%s : Starting function cleanup_demux_tables_if_error', experiment_name)

    if RawDemuxStats.objects.filter(runprocess_id__exact = run_object_name).exists():
        raw_stats_to_delete = RawDemuxStats.objects.filter(runprocess_id__exact = run_object_name)
        logger.info('%s : Deleting RawDemuxStats rows', experiment_name)
        for raw_stats in raw_stats_to_delete :
            raw_stats.delete()
    projects = Projects.objects.filter(runprocess_id = run_object_name)
    for project in projects:
        if SamplesInProject.objects.filter(project_id = project).exists():
            samples_to_delete = SamplesInProject.objects.filter(project_id = project)
            logger.info('%s : Deleting SampleProjects rows', experiment_name)
            for sample in samples_to_delete:
                sample.delete()

    if StatsFlSummary.objects.filter(runprocess_id = run_object_name).exists():
        fl_summary_to_delete = StatsFlSummary.objects.filter(runprocess_id = run_object_name)
        logger.info('%s : Deleting StatsFlSummary rows', experiment_name)
        for fl_summary in fl_summary_to_delete:
            fl_summary.delete()

    if StatsLaneSummary.objects.filter(runprocess_id = run_object_name).exists():
        lane_summary_to_delete = StatsLaneSummary.objects.filter(runprocess_id = run_object_name)
        logger.info('%s : Deleting StatsFlSummary rows', experiment_name)
        for lane_summary in lane_summary_to_delete:
            lane_summary.delete()

    logger.debug ('%s : End function cleanup_demux_tables_if_error', experiment_name)
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
    experiment_name = run_object_name.get_run_name()
    logger.debug ('%s : Starting function cleanup_run_metrics_tables_if_error', experiment_name)
    if StatsRunSummary.objects.filter(runprocess_id__exact = run_object_name).exists():
        run_summary_to_delete = RawDemuxStats.objects.filter(runprocess_id__exact = run_object_name)
        logger.info('%s : Deleting StatsRunSummary rows', experiment_name)
        for run_summary in run_summary_to_delete :
            run_summary.delete()
    if StatsRunRead.objects.filter(runprocess_id__exact = run_object_name).exists():
        run_read_to_delete = StatsRunRead.objects.filter(runprocess_id__exact = run_object_name)
        logger.info('%s : Deleting StatsRunRead rows', experiment_name)
        for run_read in run_read_to_delete :
            run_read.delete()
    logger.debug ('%s : End function cleanup_run_metrics_tables_if_error', experiment_name)
    return True

def parsing_demux_and_conversion_files( demux_file, conversion_file, number_of_lanes):
    '''
    Description:
        The function will parse demultiplexing and the conversion files
        that are created during the bcl2fastq process
    Input:
        demux_file   # local copy of the demultiplexting file
        conversion_file # local copy of the conversion file
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

        if len(p_b_count) == 0 :
            p_b_count = ['0']*len(barcodeCount)

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
    experiment_name = run_object_name.get_run_name()
    logger.debug ('%s : Starting function process_fl_summary_stats', experiment_name)
    M_BASE=1.004361/1000000
    total_cluster_lane=(stats_projects['all']['PerfectBarcodeCount'])
    processed_fl_summary_data = []
    number_of_lanes = run_object_name.get_sequencing_lanes()
    for project in stats_projects:
        project_flowcell = {}
        if project == 'TopUnknownBarcodes':
            continue

        logger.info('%s : Start processing flow Summary for project %s',experiment_name, project)
        flow_raw_cluster, flow_pf_cluster, flow_yield_mb = 0, 0, 0
        for fl_item in range(number_of_lanes):
             # make the calculation for Flowcell

            flow_raw_cluster +=int(stats_projects[project]['BarcodeCount'][fl_item])
            flow_pf_cluster +=int(stats_projects[project]['PerfectBarcodeCount'][fl_item])
            flow_yield_mb +=float(stats_projects[project]['PF_Yield'][fl_item])*M_BASE
        logger.debug('%s : flow_pf_cluster value is %s',experiment_name, flow_pf_cluster)
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
        logger.info('%s : End processing flow Summary for project %s',experiment_name, project)
        '''
        ns_fl_summary = StatsFlSummary(runprocess_id=RunProcess.objects.get(pk=run_id),
                                project_id=project_id, defaultAll=default_all, flowRawCluster=flow_raw_cluster,
                                flowPfCluster=flow_pf_cluster, flowYieldMb= flow_yield_mb,
                                sampleNumber= sample_number)


        ns_fl_summary.save()
        '''
        processed_fl_summary_data.append(project_flowcell)

    logger.debug ('%s : End function process_fl_summary_stats', experiment_name)
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
    experiment_name = run_object_name.get_run_name()
    logger.debug ('%s : Starting function process_lane_summary_stats',experiment_name)
    M_BASE=1.004361/1000000
    processed_lane_summary_data = []
    number_of_lanes = run_object_name.get_sequencing_lanes()
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
            try:
                project_lane['perfectBarcode'] = (format(int(stats_projects[project]['PerfectBarcodeCount'][i])*100/int(stats_projects[project]['BarcodeCount'][i]),'.3f'))
            except:
                project_lane['perfectBarcode'] = '0'
            try:
                project_lane['percentLane'] = format(float(int(pf_cluster_int)/int(total_cluster_lane[i]))*100, '.3f')
            except:
                project_lane['percentLane'] = '0'
            project_lane['oneMismatch'] = stats_projects[project]['OneMismatchBarcodeCount'][i]
            project_lane['yieldMb'] = '{0:,}'.format(round(float(stats_projects[project]['PF_Yield'][i])* M_BASE))
            try:
                project_lane['biggerQ30'] = format(float(stats_projects[project]['PF_YieldQ30'][i])*100/float( stats_projects[project]['PF_Yield'][i]),'.3f')
            except:
                project_lane['biggerQ30'] = '0'
            try:
                project_lane['meanQuality'] = format(float(stats_projects[project]['PF_QualityScore'][i])/float(stats_projects[project]['PF_Yield'][i]),'.3f')
            except:
                project_lane['meanQuality'] = '0'

            if project == 'all' or project == 'default':
                project_lane['project_id'] = None
                project_lane['defaultAll'] = project
            else:
                project_lane['project_id'] = Projects.objects.get(projectName__exact = project)
                project_lane['defaultAll'] =  None

            project_lane['runprocess_id'] = run_object_name

            processed_lane_summary_data.append(project_lane)
        logger.info('%s : Processed information for project %s',experiment_name, project)

    logger.debug ('%s : End function process_lane_summary_stats', experiment_name)
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
    experiment_name = run_object_name.get_run_name()
    logger.debug ('%s : Starting function store_raw_xml_stats', experiment_name)
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
        string_message = experiment_name +  " : Found some projects not defined in database"
        logging_errors(string_message, False, True)
        handling_errors_in_run(experiment_name, '10')
        raise
    logger.info('%s : Processing demultiplexing project data', experiment_name)
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
    logger.debug('%s : End function store_raw_xml_stats', experiment_name)
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
    experiment_name = run_object_name.get_run_name()
    logger.debug ('%s : Starting function process_samples_projects', experiment_name)
    # get the total number of read per lane
    M_BASE=1.004361/1000000
    processed_sample_data =[]

    for project in stats_projects:
        # find the total number of PerfectBarcodeCount in the procjec to make percent calculations
        logger.info('%s : Starting fetching sample project %s',experiment_name, project)
        total_perfect_barcode_count = 0
        for sample in stats_projects[project]:
            total_perfect_barcode_count += stats_projects[project][sample] ['PerfectBarcodeCount']

        for sample in stats_projects[project]:
            project_sample_data = {}
            project_sample_data['sampleName'] = sample
            project_sample_data['barcodeName'] = stats_projects[project][sample]['barcodeName']
            perfect_barcode = int(stats_projects[project][sample] ['PerfectBarcodeCount'])
            project_sample_data['pfClusters'] =  '{0:,}'.format(perfect_barcode)
            try:
            	project_sample_data['percentInProject'] = format (float(perfect_barcode) *100 /total_perfect_barcode_count,'.2f')
            except:
                project_sample_data['percentInProject'] = '0'
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
        logger.info('%s : Collected  sample data for the project %s',experiment_name, project)
    logger.debug('%s : End function process_samples_projects', experiment_name)
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
    experiment_name = run_object_name.get_run_name()
    logger.debug ('%s : Starting function process_unknow_barcode_stats', experiment_name)
    number_of_lanes = run_object_name.get_sequencing_lanes()
    processed_barcode_data = []
    for project in stats_projects:
        if project == 'TopUnknownBarcodes':
            logger.info ('%s : Collecting TopUnknownBarcodes data ', experiment_name)
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

                    processed_barcode_data.append(unknow_barcode)

    logger.debug ('%s : End function process_unknow_barcode_stats', experiment_name)
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
        check_run_metrics_processed # Located in this file
        create_graphics         # Located at utils.run_metrics_functions
        delete_run_metric_files # Located at utils.run_metrics_functions
        get_run_metric_files    # Located at utils.run_metrics_functions
        handling_errors_in_run  # Located at utils.generic_functions
        parsing_run_metrics     # Located at utils.run_metrics_functions
        set_state_in_all_projects  # Located at utils.generic_functions
    Variables:
        experiment_name # Name of the run
        count_file_size # partial size for the subfolder
        run_folder      # run folder on the remote server
        statistics_folder # statistics folder on the remote server
    Return:
        experiment_name if the run is updated. Empty if not
    '''
    logger = logging.getLogger(__name__)
    experiment_name = run_object_name.get_run_name()
    logger.debug ('%s : Starting function manage_run_in_processed_run', experiment_name)
    experiment_name = run_object_name.get_run_name()
    run_folder = RunningParameters.objects.get(runName_id = run_object_name).get_run_folder()
    statistics_folder = os.path.join(wetlab_config.SAMBA_APPLICATION_FOLDER_NAME, run_folder, wetlab_config.DEMULTIPLEXION_BCL2FASTQ_FOLDER)
    if 'Processed Run' == run_object_name.get_state() :
        if not check_run_metrics_processed (run_object_name) :
            run_state = run_object_name.set_run_state('Processing Metrics')
            run_id = run_object_name.get_run_id()
            try:
                logger.info ('%s : Starts collecting the run metrics data', experiment_name)
                # connect to statistics folder to fetch the run metrics files
                run_metric_files = get_run_metric_files (conn, run_folder, experiment_name)

                parsed_run_stats_summary, parsed_run_stats_read = parsing_run_metrics(wetlab_config.RUN_TEMP_DIRECTORY_PROCESSING, run_object_name)
            except Exception as e:
                string_message = experiment_name + ' : Run metrics data cannot be saved'
                logging_errors(string_message, True, True)
                #print(e)
                logger.debug ('%s : End function manage_run_in_processed_run with error', experiment_name)
                raise

            try:
                for run_stat_summary in parsed_run_stats_summary :
                    saved_run_stat_summary = StatsRunSummary.objects.create_stats_run_summary(run_stat_summary, experiment_name)
                logger.info('%s : run metrics summary data saved to database', experiment_name)
                for run_stat_read in parsed_run_stats_read :
                    saved_run_stat_read = StatsRunRead.objects.create_stats_run_read(run_stat_read, experiment_name)
                logger.info('%s :run metrics read data saved to database',experiment_name)
            except:
                string_message = experiment_name + ' : Run metrics data cannot be saved'
                logging_errors(string_message, True, True)
                # delete run metrics tables
                cleanup_run_metrics_tables_if_error (run_object_name)
                delete_run_metric_files (run_metric_files, experiment_name)
                logger.info('%s : Deleted run metrics tables',experiment_name)
                handling_errors_in_run(experiment_name)
                logger.debug ('%s : End function manage_run_in_processed_run with error',experiment_name)
                raise

            # create run graphics
            run_graphics = create_graphics (wetlab_config.RUN_TEMP_DIRECTORY_PROCESSING, run_object_name)
            logger.info('%s : run metrics graphics processed and copied to plot image folder',experiment_name)
            # deleting temporary run metrics files
            delete_run_metric_files (run_metric_files, experiment_name)
            # return the state to Processed Run
            run_state = run_object_name.set_run_state('Processed Run')

            #set_state_in_all_projects(experiment_name, 'Processed Run')
        # Check if Bcl2fastq has started
        try:
            file_list = conn.listPath( wetlab_config.SAMBA_SHARED_FOLDER_NAME, statistics_folder)
        except:
            logger.info('%s : No statistics folder were found ', experiment_name)
            logger.debug ('%s : End function manage_run_in_processed_run',experiment_name)
            return experiment_name

        run_state = run_object_name.set_run_state('Processing Bcl2fastq')
        # Do not update the projects state. They will keep in Processed Run

    else:
        string_message = experiment_name + ' : Invalid state when calling Processed run function'
        logging_errors(string_message , False , False)
        logger.debug ('%s : End function manage_run_in_processed_run with error',experiment_name)
        return
    logger.debug ('%s : End function manage_run_in_processed_run',experiment_name)
    return experiment_name


def manage_run_in_processing_bcl2fastq (conn, run_object_name):
    '''
    Description:
        The function will check if report folder exists. Then the bcl2fastq
        conversion is completed
        If exists the run state will move to Processed Bcl2fastq
    Input:
        conn # Connectio samba object
        run_object_name   # runProcess object
    Constant:
        REPORT_FOLDER
        SAMBA_APPLICATION_FOLDER_NAME
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
    experiment_name = run_object_name.get_run_name()
    logger.debug ('%s : Starting function manage_run_in_processing_bcl2fast2_run', experiment_name)
    logger.info('%s : Start handling in Processing_Bcl2fastq state', experiment_name)
    run_folder = RunningParameters.objects.get(runName_id = run_object_name).get_run_folder()
    statistics_folder = os.path.join(wetlab_config.SAMBA_APPLICATION_FOLDER_NAME, run_folder, wetlab_config.DEMULTIPLEXION_BCL2FASTQ_FOLDER)

    if 'Processing Bcl2fastq' == run_object_name.get_state() :
        # connect to report folder to check if bcl2fastq process is ended
        try:
            file_list = conn.listPath( wetlab_config.SAMBA_SHARED_FOLDER_NAME, statistics_folder)
        except Exception as e:
            string_message = experiment_name + ' : Folder statistics does not exists '
            logging_errors(string_message, False, False)
            logger.debug ('%s :End function manage_run_in_processing_bcl2fast2 with error',experiment_name)
            return ''
        for sh in file_list:
            if sh.filename == wetlab_config.REPORT_FOLDER :

                logger.info('%s : bcl2fastq has been completed . Extracting data', experiment_name)
                # Get the time when  the Bcl2Fastq process is ending
                s_conversion_stats = os.path.join (statistics_folder, wetlab_config.STATS_FOLDER, wetlab_config.CONVERSION_STATS_FILE)
                conversion_attributes = conn.getAttributes(wetlab_config.SAMBA_SHARED_FOLDER_NAME ,s_conversion_stats)
                bcl2fastq_finish_date = datetime.fromtimestamp(int(conversion_attributes.create_time)).strftime('%Y-%m-%d %H:%M:%S')
                bcl2fastq_finish_date = run_object_name.set_run_bcl2fastq_finished_date(bcl2fastq_finish_date)
                run_object_name.set_run_state('Processed Bcl2fastq')
                # Do not update the projects state. They will keep in Processed Run

                logger.info ('%s : Updated the Bcl2Fastq time in the run ', experiment_name)
                break

    else:
        string_message = experiment_name + ' : Invalid state when calling to '
        logging_errors(string_message , False , False)
        logger.debug ('%s : End function manage_run_in_processing_bcl2fast2 with error',experiment_name)
        return ''

    logger.debug ('%s : End function manage_run_in_processing_bcl2fast2', experiment_name)
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
    experiment_name = run_object_name.get_run_name()
    logger.debug ('%s : Starting function manage_run_in_processed_bcl2fast2_run',experiment_name)

    run_folder = RunningParameters.objects.get(runName_id = run_object_name).get_run_folder()
    number_of_lanes = run_object_name.get_sequencing_lanes()

    if 'Processed Bcl2fastq' == run_object_name.get_state() :
        run_state = run_object_name.set_run_state('Processing Demultiplexing')
        #projects_state = set_state_in_all_projects(experiment_name, 'Processing Demultiplexing')

        try:
            l_demux , l_conversion= get_bcl2fastq_output_files (conn, run_folder)
        except:
            string_message = experiment_name + ' : Unable to fetch stats files '
            logging_errors(string_message, False, False)
            handling_errors_in_run(experiment_name, '11')
            logger.debug ('%s : End function manage_run_in_processed_bcl2fast2 with error',experiment_name)
            raise
        # parsing the files to get the xml Stats
        logger.info('%s : Start parsing  demultiplexing files',experiment_name)
        parsed_result = parsing_demux_and_conversion_files(l_demux, l_conversion, number_of_lanes)

        # parsing and processing the project samples
        logger.info('%s : Start parsing  samples demultiplexing',experiment_name)
        sample_parsed_result = parsing_demux_sample_project (l_demux, l_conversion, number_of_lanes)
        # clean up the fetched files in the local temporary folder
        os.remove(l_conversion)
        os.remove(l_demux)
        logger.info ('%s : Deleted temporary demultiplexing and conversion files', experiment_name)

        processed_raw_stats = process_raw_demux_stats(parsed_result, run_object_name)
        try:
            for raw_stats in processed_raw_stats :
                new_raw_stats = RawDemuxStats.objects.create_stats_run_read(raw_stats, run_object_name)
            logger.info('%s : Saved information to RawDemuxStats', experiment_name)
        except:
            string_message = 'Unable to save raw stats for ' + experiment_name
            logging_errors(string_message, False, False)
            handling_errors_in_run (experiment_name, '11' )
            cleanup_demux_tables_if_error(run_object_name)
            logger.debug('%s : End function manage_run_in_processed_bcl2fast2_run with error', experiment_name)
            raise

        processed_samples_stats = process_samples_projects(sample_parsed_result, run_object_name)
        try:
            for sample_stats in processed_samples_stats :
                new_sample_stats = SamplesInProject.objects.create_sample_project(sample_stats)
            logger.info('%s : Saved information to SamplesInProject', experiment_name)
        except:
            string_message = experiment_name + ' : Unable to save sample stats'
            logging_errors(string_message, False, False)
            handling_errors_in_run (experiment_name, '13' )
            cleanup_demux_tables_if_error(run_object_name)
            logger.debug('%s : End function manage_run_in_processed_bcl2fast2_run with error', experiment_name)
            raise
        processed_fl_summary = process_fl_summary_stats(parsed_result, run_object_name)
        try:
            for fl_summary in processed_fl_summary :
                new_fl_summary = StatsFlSummary.objects.create_fl_summary(fl_summary)
            logger.info('%s : Saved information to StatsFlSummary', experiment_name)
        except:
            string_message = experiment_name + ' : Unable to save FL Summary stats'
            logging_errors(string_message, False, False)
            handling_errors_in_run (experiment_name, '14' )
            cleanup_demux_tables_if_error(run_object_name)
            logger.debug('%s : End function manage_run_in_processed_bcl2fast2_run with error', experiment_name)
            raise

        processed_lane_summary = process_lane_summary_stats(parsed_result, run_object_name)
        try:
            for lane_summary in processed_lane_summary :
                new_lane_summary = StatsLaneSummary.objects.create_lane_summary(lane_summary)
            logger.info('%s : Saved information to StatsLaneSummary', experiment_name)
        except:
            string_message = experiment_name + ' : Unable to save Lane Summary stats'
            logging_errors(string_message, False, False)
            handling_errors_in_run (experiment_name, '15' )
            cleanup_demux_tables_if_error(run_object_name)
            logger.debug('%s :End function manage_run_in_processed_bcl2fast2_run with error', experiment_name)
            raise

        processed_unknow_barcode = process_unknow_barcode_stats(parsed_result, run_object_name)
        try:
            for unknow_barcode in processed_unknow_barcode :
                new_unknow_barcode = RawTopUnknowBarcodes.objects.create_unknow_barcode(unknow_barcode)
            logger.info('%s : Saved information to RawTopUnknowBarcodes',experiment_name)
        except:
            string_message = experiment_name + ' : Unable to Unknow Barcode stats'
            logging_errors(string_message, False, False)
            handling_errors_in_run (experiment_name, '16' )
            cleanup_demux_tables_if_error(run_object_name)
            logger.debug('%s : End function manage_run_in_processed_bcl2fast2_run with error', experiment_name)
            raise

        ## Get the disk space utilization for this run
        try:
            disk_utilization = get_run_disk_utilization (conn, run_folder)
        except:
            string_message = experiment_name + ' : Error when fetching the disk utilizaton'
            logging_errors (string_message, False, False)
            handling_errors_in_run (experiment_name, '17' )
            cleanup_demux_tables_if_error(run_object_name)
            logger.debug('%s : End function manage_run_in_processed_bcl2fast2_run with error', experiment_name)
            raise

        result_store_usage = run_object_name.set_used_space (disk_utilization)
        finish_date = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        result_set_finish_date = run_object_name.set_run_finish_date(finish_date)
        # Update the run state to completed
        run_state = run_object_name.set_run_state('Completed')
        logger.info('%s : is Completed',experiment_name)


    else:
        string_message = experiment_name + ' : Invalid state when calling this function'
        logging_errors(string_message , False , False)
        logger.debug ('%s : End function manage_run_in_processed_bcl2fast2 with error', experiment_name)
        return ''
    logger.debug ('%s : End function manage_run_in_processed_bcl2fast2 ', experiment_name)
    return experiment_name
