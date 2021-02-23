import os
import xml.etree.ElementTree as ET
from iSkyLIMS_wetlab.wetlab_config import *
from iSkyLIMS_wetlab.models import Projects, RawDemuxStats, StatsFlSummary, StatsLaneSummary
from .handling_crontab_common_functions import *


def check_demultiplexing_folder_exists(conn, run_folder, experiment_name):
    '''
    Description:
        The function check if the folder of demultiplexing files exists
    Input:
        conn                # samba connection instance
        run_folder          # run folder on the remote server
        experiment_name     # Experiment name
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
        experiment_name     # Experiment name
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
        experiment_name     # Experiment name
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
        demux_files
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
        conversion_file         # local copy of the conversion file
        experiment_name         # Experiment name
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
    logger.info('%s : Starting parsing for conversion file', experiment_name)
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
        logger.info('%s : Completed parsing for conversion stats for project %s', experiment_name, projects[i])

    unknow_lanes  = []
    unknow_barcode_start_index= len(projects)
    counter=0
    logger.info('%s : Collecting the Top Unknow Barcodes', experiment_name)
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
    logger.debug ('%s : End function parsing_demux_and_conversion_files', experiment_name)

    return parsed_result



def parsing_demux_sample_project(demux_files, number_of_lanes, experiment_name):
    '''
    Description:
        The function will parse demultiplexing and the conversion files
        to get the samples inside the project
    Input:
        demux_files         # local copy of the demultiplexting  and conversion files
        number_of_lanes     # numer of lane to be fetched, based on the sequencer type
        experiment_name     # Experiment name
    Return:
        parsed_result
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting function parsing_demux_and_conversion_files', experiment_name)

    total_p_b_count=[0,0,0,0]
    parsed_result = {}
    demux_stat=ET.parse(demux_files['demux_stats'])
    root=demux_stat.getroot()
    projects=[]

    logger.info('%s : Starting parsing DemultiplexingStats.XML for getting Sample information', experiment_name)
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
    logger.info('%s : Complete parsing from demux file for sample and for project %s', experiment_name , projects[i])

    conversion_stat=ET.parse(demux_files['conversion_stats'])
    root_conv=conversion_stat.getroot()
    projects=[]
    logger.info('%s : Starting conversion for conversion file', experiment_name)
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
        logger.info('%s : Completed parsing for xml stats for project %s', experiment_name, projects[i])


    logger.debug('%s : End function parsing_demux_sample_project', experiment_name)

    return parsed_result



def process_and_store_fl_summary_data(parsed_data, run_process_obj, number_of_lanes,  experiment_name ):
    '''
    Description:
        The function manages the parsed data to save them in StatsFlSummary
    Input:
        parsed_data         # dictionnary with parsed data
        run_process_obj     # runProcess object
        experiment_name     # Experiment name
    Return:
        None
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting function process_and_store_fl_summary_data', experiment_name)
    M_BASE=1.004361/1000000
    total_cluster_lane=(parsed_data['all']['PerfectBarcodeCount'])
    processed_fl_summary_data = []
    number_of_lanes = run_process_obj.get_sequencing_lanes()
    for project in parsed_data:
        project_flowcell = {}
        if project == 'TopUnknownBarcodes':
            continue

        logger.info('%s : Start processing flow Summary for project %s',experiment_name, project)
        flow_raw_cluster, flow_pf_cluster, flow_yield_mb = 0, 0, 0

        for fl_item in range(number_of_lanes):
             # make the calculation for Flowcell

            flow_raw_cluster +=int(parsed_data[project]['BarcodeCount'][fl_item])
            flow_pf_cluster +=int(parsed_data[project]['PerfectBarcodeCount'][fl_item])
            flow_yield_mb +=float(parsed_data[project]['PF_Yield'][fl_item])*M_BASE
        logger.debug('%s : flow_pf_cluster value is %s',experiment_name, flow_pf_cluster)
        project_flowcell['flowRawCluster'] = '{0:,}'.format(flow_raw_cluster)
        project_flowcell['flowPfCluster'] ='{0:,}'.format(flow_pf_cluster)
        project_flowcell['flowYieldMb'] = '{0:,}'.format(round(flow_yield_mb))
        project_flowcell['sampleNumber'] = parsed_data[project]['sampleNumber']

        if project == 'all' or project == 'default' :
            project_flowcell['project_id'] = None
            project_flowcell['defaultAll'] = project
        else:
            project_flowcell['project_id'] = Projects.objects.get(projectName__exact = project)
            project_flowcell['defaultAll'] = None
        project_flowcell['runprocess_id'] = run_process_obj
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




def process_and_store_raw_demux_project_data(parsed_data, run_process_obj, experiment_name):
    '''
    Description:
        The function manages the parsed data to save them in RawDemuxStats
    Input:
        parsed_data         # dictionnary with parsed data
        run_process_obj     # runProcess object
        experiment_name     # Experiment name
    Functions:

    Return:
        None
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting function store_raw_xml_stats', experiment_name)
    processed_raw_data = []
    project_list = []
    # check first if all project found in the parsing are already defined.
    if run_process_obj.projects_set.all().exists():
        project_objs = run_process_obj.projects_set.all()
        for project_obj in project_objs:
            project_list.append(project_obj.get_project_name())

    for project in parsed_data.keys():
        if project == 'TopUnknownBarcodes'  or project == 'all' or project == 'default':
            continue
        if not project in project_list:
            # create project and link to the run
            new_project_obj = create_new_empty_project({'user_id':None, 'projectName':project})
            new_project_obj.runProcess.add(run_process_obj)
            string_message = experiment_name +  ' : Created  project name ' + project + 'Because it was not store'
            logging_warnings(string_message, True)

    logger.info('%s : Processing demultiplexing raw project data', experiment_name)
    for project in parsed_data.keys():
        project_raw_data = {}
        # ignore on this step the unknowbarcode parsed information
        if project == 'TopUnknownBarcodes':
            continue
        logger.info('%s : Start processing project %s for raw parsed data', experiment_name, project)
        if project == 'all' or project == 'default':
            project_raw_data['project_id']= None
            project_raw_data['defaultAll'] = project
        else:
            project_raw_data['project_id']= Projects.objects.filter(projectName__exact = project).last()
            project_raw_data['defaultAll'] = None

        project_raw_data['rawYield'] = parsed_data[project]['RAW_Yield']
        project_raw_data['rawYieldQ30'] =  parsed_data[project]['RAW_YieldQ30']
        project_raw_data['rawQuality'] = parsed_data[project]['RAW_QualityScore']
        project_raw_data['PF_Yield'] = parsed_data[project]['PF_Yield']
        project_raw_data['PF_YieldQ30'] = parsed_data[project]['PF_YieldQ30']
        project_raw_data['PF_QualityScore'] = parsed_data[project]['PF_QualityScore']

        import pdb; pdb.set_trace()
        new_raw_obj = RawDemuxStats.objects.create_stats_run_read(project_raw_data, run_process_obj)
        logger.info('%s : End processing project %s for raw parsed data', experiment_name, project)
    logger.info('%s : Completed processing demultiplexing raw project data', experiment_name)
    logger.debug('%s : End function store_raw_xml_stats', experiment_name)
    return


def process_and_store_samples_projects_data(parsed_data, run_process_obj, experiment_name ):
    '''
    Description:
        The function will parse demultiplexing and the conversion files
        that are created during the bcl2fastq process
    Input:
        parsed_data      # dictionnary with parsed data
        run_process_obj     # runProcess object
        experiment_name     # Experiment name
    Return:
        None
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting function process_and_store_samples_projects_data', experiment_name)
    # get the total number of read per lane
    M_BASE=1.004361/1000000
    processed_sample_data =[]

    for project in parsed_data:
        # find the total number of PerfectBarcodeCount in the procjec to make percent calculations
        logger.info('%s : Starting fetching sample project %s',experiment_name, project)
        total_perfect_barcode_count = 0
        for sample in parsed_data[project]:
            total_perfect_barcode_count += parsed_data[project][sample] ['PerfectBarcodeCount']

        for sample in parsed_data[project]:
            project_sample_data = {}
            project_sample_data['sampleName'] = sample
            project_sample_data['barcodeName'] = parsed_data[project][sample]['barcodeName']
            perfect_barcode = int(parsed_data[project][sample] ['PerfectBarcodeCount'])
            project_sample_data['pfClusters'] =  '{0:,}'.format(perfect_barcode)
            try:
            	project_sample_data['percentInProject'] = format (float(perfect_barcode) *100 /total_perfect_barcode_count,'.2f')
            except:
                project_sample_data['percentInProject'] = '0'
            project_sample_data['yieldMb'] = '{0:,}'.format(round(float(parsed_data[project][sample] ['PF_Yield'])*M_BASE))
            if parsed_data[project][sample] ['PF_Yield'] > 0:
                bigger_q30=format(float(parsed_data[project][sample]['PF_YieldQ30'])*100/float( parsed_data[project][sample]['PF_Yield']),'.3f')
                mean_quality=format(float(parsed_data[project][sample]['PF_QualityScore'])/float(parsed_data[project][sample]['PF_Yield']),'.3f')
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
            new_sample_stats = SamplesInProject.objects.create_sample_project(sample_stats)

            processed_sample_data.append(project_sample_data)
        logger.info('%s : Collected  sample data for the project %s',experiment_name, project)
    logger.debug('%s : End function process_samples_projects', experiment_name)
    return processed_sample_data
