import os
import xml.etree.ElementTree as ET
from iSkyLIMS_wetlab.wetlab_config import *
from .handling_crontab_common_functions import *


def check_demultiplexing_folder_exists(conn, run_folder, experiment_name):
    '''
    Description:
        The function check if the folder of demultiplexing files exists
    Input:
        conn                # samba connection instance
        run_folder          # run folder on the remote server
        experiment_name     # experiment name to be checked
    Constants:
        DEMULTIPLEXION_BCL2FASTQ_FOLDER
        REPORT_FOLDER
        STATS_FOLDER
        CONVERSION_STATS_FILE
    Functions:
        get_samba_application_shared_folder     # Located at utils/handling_crontab_common_functions.py
        get_samba_shared_folder                 # Located at utils/handling_crontab_common_functions.py
    Return:
        bcl2fastq_finish_date
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting function check_demultiplexing_folder_exists', experiment_name)
    bcl2fastq_finish_date = ''
    statistics_folder = os.path.join(get_samba_application_shared_folder(), run_folder, DEMULTIPLEXION_BCL2FASTQ_FOLDER)

    try:
        file_list = conn.listPath( get_samba_shared_folder(), statistics_folder)
    except:
        string_message = experiment_name + ' : Unable to fetch folder demultiplexing at  ' + statistics_folder
        logging_errors(string_message, True, True)
        logger.debug ('%s : End function check_demultiplexing_folder_exists with error', experiment_name)
        return {'ERROR':29}

    logger.info('%s : bcl2fastq has been completed . Collecting date when finish the bcl2fastq', experiment_name)

    s_conversion_stats = os.path.join (statistics_folder, STATS_FOLDER, CONVERSION_STATS_FILE)
    try:
        conversion_attributes = conn.getAttributes(get_samba_shared_folder() ,s_conversion_stats)
    except:
        string_message = experiment_name + ' : Unable to fetch ' + CONVERSION_STATS_FILE
        logging_errors(string_message, True, True)
        logger.debug ('%s : End function check_demultiplexing_folder_exists with error', experiment_name)
        return {'ERROR':31}
    bcl2fastq_finish_date = datetime.fromtimestamp(int(conversion_attributes.create_time)).strftime('%Y-%m-%d %H:%M:%S')
    logger.info ('%s : Collected Bcl2Fastq time ', experiment_name)

    logger.debug ('%s : End function waiting_time_expired', experiment_name)
    return bcl2fastq_finish_date

def delete_existing_bcl2fastq_table_processed(run_process_obj, experiment_name):
    '''
    Description:
        The function check if the exists tables created for bcl2fatsq state are already defined.
        If so they are deleted to avoid wrong duplication information
    Input:
        run_process_obj     # run process instance
        experiment_name     # experiment name to be checked
    Return:
        True
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting function delete_existing_bcl2fastq_table_processed', experiment_name)
    if RawDemuxStats.objects.filter(runprocess_id = run_process_obj).exists():
        raw_demux_objs = RawDemuxStats.objects.filter(runprocess_id = run_process_obj)
        for raw_demux_obj in raw_demux_objs:
            raw_demux_obj.delete()
        logger.info('%s : Deleted RawDemuxStats tables ', experiment_name)

    if StatsFlSummary.objects.filter(runprocess_id = run_process_obj).exists():
        stats_fl_summary_objs = StatsFlSummary.objects.filter(runprocess_id = run_process_obj)
        for stats_fl_summary_obj in stats_fl_summary_objs:
            stats_fl_summary_obj.delete()
        logger.info('%s : Deleted StatsFlSummary tables ', experiment_name)

    if StatsLaneSummary.objects.filter(runprocess_id = run_process_obj).exists():
        stats_lane_summary_objs = StatsLaneSummary.objects.filter(runprocess_id = run_process_obj)
        for stats_lane_summary_obj in stats_lane_summary_objs:
            stats_lane_summary_obj.delete()
        logger.info('%s : Deleted StatsLaneSummary tables ', experiment_name)

    logger.debug ('%s : End function delete_existing_bcl2fastq_table_processed', experiment_name)
    return True


def get_demultiplexing_files(conn, run_folder, experiment_name):
    '''
    Description:
        The function check if the folder of demultiplexing files exists
    Input:
        conn                # samba connection instance
        run_folder          # run folder on the remote server
        experiment_name     # experiment name to be checked
    Constants:
        DEMULTIPLEXION_BCL2FASTQ_FOLDER
        REPORT_FOLDER
        STATS_FOLDER
        CONVERSION_STATS_FILE
    Functions:
        get_samba_application_shared_folder     # Located at utils/handling_crontab_common_functions.py
        get_samba_shared_folder                 # Located at utils/handling_crontab_common_functions.py
        fetch_remote_file                       # Located at utils/handling_crontab_common_functions.py
    Return:
        bcl2fastq_finish_date
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting function get_demultiplexing_files', experiment_name)
    statistics_folder = os.path.join(get_samba_application_shared_folder(), run_folder, DEMULTIPLEXION_BCL2FASTQ_FOLDER)
    try:
        file_list = conn.listPath( get_samba_shared_folder(), statistics_folder)
    except:
        string_message = experiment_name + ' : Unable to fetch folder demultiplexing at  ' + statistics_folder
        logging_errors(string_message, True, True)
        logger.debug ('%s : End function get_demultiplexing_files with error', experiment_name)
        return {'ERROR':29}
    # conversion stats file
    l_conversion_stats = os.path.join(RUN_TEMP_DIRECTORY, CONVERSION_STATS_FILE)
    s_conversion_stats = os.path.join (statistics_folder, STATS_FOLDER, CONVERSION_STATS_FILE)
    # demultiplexion stats file
    l_demux_stats = os.path.join(RUN_TEMP_DIRECTORY, DEMULTIPLEXION_STATS_FILE)
    s_demux_stats = os.path.join (statistics_folder, STATS_FOLDER, DEMULTIPLEXION_STATS_FILE)
    demux_files = {}
    try:
        demux_files['conversion_stats'] = fetch_remote_file (conn, run_folder, s_conversion_stats, l_conversion_stats)
        logger.info('%s : Sucessfully fetch of ConversionStats file', experiment_name)
    except:
        string_message = experiment_name + ' : cannot copy or fetch the ' + CONVERSION_STATS_FILE +' file'
        logging_errors(string_message, True , True)
        logger.debug ('%s : End function get_demultiplexing_files with execption', experiment_name)
        return {'ERROR':31}
    try:
        demux_files['demux_stats'] = fetch_remote_file (conn, run_folder, s_demux_stats, l_demux_stats)
        logger.info('%s : Sucessfully fetch of demultiplexing file', experiment_name)
    except:
        string_message = experiment_name + ' : cannot copy or fetch the ' + DEMULTIPLEXION_STATS_FILE +' file'
        logging_errors(string_message, True , True)
        os.remove(l_conversion_stats)
        logger.info('%s : deleted %s file', experiment_name, l_conversion_stats)
        logger.debug ('%s : End function get_demultiplexing_files with execption', experiment_name)
        return {'ERROR':31}
    logger.debug ('%s : End function get_demultiplexing_files ', experiment_name)
    return demux_files



def parsing_demux_and_conversion_files(demux_files, number_of_lanes, experiment_name):
    '''
    Description:
        The function will parse demultiplexing and the conversion files that are created during the bcl2fastq process
    Input:
        demux_files             # local copy of the demultiplexting  and conversion files
        conversion_file # local copy of the conversion file
        experiment_name         # experiment name to be checked
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
    logger.debug ('%s : Starting function parsing_demux_and_conversion_files', experiment_name)

    total_p_b_count = [0,0,0,0]
    parsed_result={}
    projects=[]
    demux_parse = ET.parse(demux_files['demux_stats'])
    root=demux_parse.getroot()

    logger.info('%s : Start parsing demux file', experiment_name)
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

        ## Fill with zeroes if not PerfectBarcodeCount is written in file
        if len(barcodeCount) != len(p_b_count) :
            for i_bar in range(len(barcodeCount)):
                if barcodeCount[i_bar] == '0':
                    p_b_count.insert(i_bar,'0')

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
        logger.info('%s : Completed parsing from demux file for project %s',experiment_name ,projects[i])
    # overwrite the value for total samples
    logger.info('%s : Completed parsed information for all projects in stats files', experiment_name)
    parsed_result['all']['sampleNumber'] = total_samples

    conversion_stat=ET.parse(demux_files['conversion_stats'])
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

    if len(unknow_lanes) != number_of_lanes:
        for index_bar in range(len(barcodeCount)):
            if barcodeCount[index_bar] == '0':
                empty_data ={'count':'0', 'sequence':'Not Applicable'}
                fill_data = [empty_data]*10
                unknow_lanes.insert(index_bar,fill_data)
    parsed_result['TopUnknownBarcodes']= unknow_lanes
    logger.debug ('End function parsing_demux_and_conversion_files')

    return parsed_result
