#!/usr/bin/env python3
#from  ..models import *

import sys, os, re
import xml.etree.ElementTree as ET
import time
import shutil
import logging
from logging.handlers import RotatingFileHandler
from interop import py_interop_run_metrics, py_interop_run, py_interop_summary, py_interop_plot

def open_log():
    
    LOG_FILENAME = 'testing.log'
    #def create_log ():
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    #create the file handler
    handler = logging.handlers.RotatingFileHandler(LOG_FILENAME, maxBytes=20000, backupCount=5)
    handler.setLevel(logging.DEBUG)
    
    #create a Logging format
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    #add the handlers to the logger
    logger.addHandler(handler)
    
    return logger


def test():
    print ('esto es una prueba de import')
    
def parsing_statistics_xml(demux_file, conversion_file, logger):
    total_p_b_count=[0,0,0,0] 
    stats_result={}
    #demux_file='example.xml'
    demux_stat=ET.parse(demux_file)
    root=demux_stat.getroot()
    projects=[]
    logger.info('Starting conversion for demux file')
    for child in root.iter('Project'):
        projects.append(child.attrib['name'])
    
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
        dict_stats['sampleNumber']=len(samples)
        dict_stats['OneMismatchBarcodeCount']=one_mismatch_count
        stats_result[projects[i]]=dict_stats
        logger.info('Complete parsing from demux file for project %s', projects[i])
    
    
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


def process_xml_stats(stats_projects, run_id, logger):
    
    for project in stats_projects:
        if project == 'TopUnknownBarcodes':
            for un_lane in range(4) :
                logger.info('Processing lane %s for TopUnknownBarcodes')
                count_top=0
                lane_number=str(un_lane + 1)
                for barcode_line in stats_projects[project][un_lane]:
                    #for count, sequence barcode_line.iteritems()
                    barcode_count= barcode_line['count']
                    barcode_sequence= barcode_line['sequence']
    
    
    return stats_result

def process_binStats(run_folder, run_id, logger):
    
    logger.info('starting analyzing the binary statistics ')
    run_metrics = py_interop_run_metrics.run_metrics()
    #run_folder = run_metrics.read(run_folder)
    valid_to_load = py_interop_run.uchar_vector(py_interop_run.MetricCount, 0)
    py_interop_run_metrics.list_summary_metrics_to_load(valid_to_load)
    run_folder = run_metrics.read(run_folder)
    #run_folder = run_metrics.read(run_folder, valid_to_load)

    summary = py_interop_summary.run_summary()
    py_interop_summary.summarize_run_metrics(run_metrics, summary)
    '''
    # get the Run Summary for Read 1 to 4
    for read_level in range(4):
    
        # summary yield total
        read_summary_yield_g=format(summary.at(read_level).summary().yield_g(),'.3f')
        # summary projected total yield
        read_summary_projected_yield_g=format(summary.at(read_level).summary().projected_yield_g(),'.3f')
        # 
        # percent yield
        read_summary_percent_aligned=format(summary.at(read_level).summary().percent_aligned(),'.3f')
        # Error rate
        read_summary_error_rate=format(summary.at(read_level).summary().error_rate(),'.3f')
        # intensity cycle 1
        read_summary_first_cycle_intensity=str(round(summary.at(read_level).summary().first_cycle_intensity()))
        # Q30
        read_percent_gt_q30=format(summary.at(read_level).summary().percent_gt_q30(),'.3f')
        logger.info('Store the Run Summary  for read %s', read_level)
        ns_bin_run_summary = NextSeqStatsBinRunSummary (runprocess_id=RunProcess.objects.get(pk=run_id),
                                                  level=str(read_level), yieldTotal = read_summary_yield_g,
                                                  projectedTotalYield= read_summary_projected_yield_g,
                                                  aligned= read_summary_percent_aligned, errorRate= read_summary_error_rate,
                                                  intensityCycle= read_summary_first_cycle_intensity, biggerQ30= read_percent_gt_q30)
        # ns_bin_run_summary.save()
        
    # get the run summary for Total
    # total summary
    total_s_yield_g=format(summary.total_summary().yield_g(),'.3f')
    # total projected_yield_g
    total_s_projected_yield_g=format(summary.total_summary().projected_yield_g(),'.3f')
    #
    # total percent aligned
    total_s_percent_aligned=format(summary.total_summary().percent_aligned(),'.3f')
    # total error rate
    total_s_error_rate=format(summary.total_summary().error_rate(),'.3f')
    # total intensity cycle
    total_s_first_cycle_intensity=str(round(summary.total_summary().first_cycle_intensity()))
    # total Q 30
    total_s_percent_gt_q30=format(summary.total_summary().percent_gt_q30(),'.3f')
    logger.info('Store the Run Summary  for Total lane')
    ns_bin_run_summary = NextSeqStatsBinRunSummary (runprocess_id=RunProcess.objects.get(pk=run_id),
                                                  level='Total', yieldTotal = total_s_yield_g,
                                                  projectedTotalYield= total_s_projected_yield_g,
                                                  aligned= total_s_percent_aligned, errorRate= total_s_error_rate,
                                                  intensityCycle= total_s_first_cycle_intensity, biggerQ30= total_s_percent_gt_q30)
     # ns_bin_run_summary.save()
     
     # get the run summary for non index
     # non index yield
    nonindex_s_yield_g=format(summary.nonindex_summary().yield_g(),'.3f')
    #  non index projected yield
    nonindex_s_projected_yield_g=format(summary.nonindex_summary().projected_yield_g(),'.3f')
    
    # non index percent aligned
    nonindex_s_percent_aligned=format(summary.nonindex_summary().percent_aligned(),'.3f')
    # non index percent error rate
    nonindex_s_error_rate=format(summary.nonindex_summary().error_rate(),'.3f')
    # non index intensity cycle
    nonindex_s_first_cycle_intensity=str(round(summary.nonindex_summary().first_cycle_intensity()))
    # non index Q 30
    nonindex_s_percent_gt_q30=format(summary.nonindex_summary().percent_gt_q30(),'.3f')
    logger.info('Store the Run Summary  for Non Index lane')
    ns_bin_run_summary = NextSeqStatsBinRunSummary (runprocess_id=RunProcess.objects.get(pk=run_id),
                                                  level='Non Index', yieldTotal = nonindex_s_yield_g,
                                                  projectedTotalYield= nonindex_s_projected_yield_g,
                                                  aligned= nonindex_s_percent_aligned, errorRate= nonindex_s_error_rate,
                                                  intensityCycle= nonindex_s_first_cycle_intensity, biggerQ30= nonindex_s_percent_gt_q30)   

     # ns_bin_run_summary.save()
    '''
    ### information per reads
    
    #lan_summary= py_interop_summary.lane_summary()
    # Tiles
    for read_number in range(4):
        logger.info('Processing bin stats for Read %s', read_number)
        for lane_number in range(4):
            
            read_lane_tiles=str(int(summary.at(0).at(0).tile_count() )*2)
            # Density (k/mm2) divide the value by 1000 to have it K/mm2
            # get the +/- with the steddev
            read_lane_density_mean=str(round(float(summary.at(read_number).at(lane_number).density().mean())/1000))
            read_lane_density_stddev=str(round(float(summary.at(read_number).at(lane_number).density().stddev())/1000))
            # cluster _pf  in % 
            read_lane_percent_pf_mean=format(summary.at(read_number).at(lane_number).percent_pf().mean(),'.3f')
            read_lane_percent_pf_stddev=format(summary.at(read_number).at(lane_number).percent_pf().stddev(),'.3f')  
            read_lane_percent_pf_field=read_lane_percent_pf_mean + '  ' + chr(177) + '  ' +read_lane_percent_pf_stddev
            # phas/ prepas in %
            read_lane_phasing_mean=format(summary.at(read_number).at(lane_number).phasing().mean(),'.3f')
            read_lane_phasing_dev=format(summary.at(read_number).at(lane_number).phasing().stddev(),'.1f')
            read_lane_prephasing_mean=format(summary.at(read_number).at(lane_number).prephasing().mean(),'.3f')
            read_lane_prephasing_stddev=format(summary.at(read_number).at(lane_number).prephasing().stddev(),'.3f')
            read_lane_phas_prephas_field=read_lane_phasing_mean + '  ' + chr(177) + '  ' + read_lane_phasing_dev + '  /  ' + read_lane_prephasing_mean + '  ' + chr(177) + '  ' + read_lane_prephasing_stddev
            # reads (M)
            read_lane_reads=format(float(summary.at(read_number).at(lane_number).reads())/1000000,'.3f')
            #reads PF (M)
            read_lane_reads_pf=format(float(summary.at(read_number).at(lane_number).reads_pf())/1000000,'.3f')
            # percent q30
            read_lane_percent_gt_q30=format(summary.at(read_number).at(lane_number).percent_gt_q30(),'.3f')
            # yield _g
            read_lane_yield_g=format(summary.at(read_number).at(lane_number).yield_g(),'.3f')
            # cycles err Rate
            read_lane_cycles_error_rate=str(summary.at(read_number).at(lane_number).cycle_state().error_cycle_range().first_cycle())
            #percent_aligned
            read_lane_percent_aligned_mean=format(summary.at(read_number).at(lane_number).percent_aligned().mean(),'.3f')
            read_lane_percent_aligned_stddev=format(summary.at(read_number).at(lane_number).percent_aligned().stddev(),'3f')
            read_lane_percent_aligned_field= read_lane_percent_aligned_mean + '  ' + chr(177) + '  ' + read_lane_percent_aligned_stddev
            #error rate
            read_lane_error_rate_mean=format(summary.at(read_number).at(lane_number).error_rate().mean(),'.3f')
            read_lane_error_rate_stddev=format(summary.at(read_number).at(lane_number).error_rate().stddev(),'.3f')
            read_lane_error_rate_field=read_lane_error_rate_mean+ '  ' + chr(177) + '  ' + read_lane_error_rate_stddev
            #error rate_35
            read_lane_error_rate_35_mean=format(summary.at(read_number).at(lane_number).error_rate_35().mean(),'.3f')
            read_lane_error_rate_35_stddev=format(summary.at(read_number).at(lane_number).error_rate_35().stddev(),'.3f')
            read_lane_error_rate_35_field=read_lane_error_rate_35_mean + '  ' + chr(177) + '  ' + read_lane_error_rate_35_stddev
            #error rate 50
            read_lane_error_rate_50_mean=format(summary.at(read_number).at(lane_number).error_rate_50().mean(),'.3f')
            read_lane_error_rate_50_stddev=format(summary.at(read_number).at(lane_number).error_rate_50().stddev(),'.3f')
            read_lane_error_rate_50_field=read_lane_error_rate_50_mean + '  ' + chr(177) + '  ' + read_lane_error_rate_50_stddev
            #error rate 75
            read_lane_error_rate_75_mean=format(summary.at(read_number).at(lane_number).error_rate_75().mean(),'.3f')
            read_lane_error_rate_75_stddev=format(summary.at(read_number).at(lane_number).error_rate_75().stddev(),'.3f')
            read_lane_error_rate_75_field=read_lane_error_rate_75_mean + '  ' + chr(177) + '  ' + read_lane_error_rate_75_stddev
            #error rate 100
            read_lane_error_rate_100_mean=format(summary.at(read_number).at(lane_number).error_rate_100().mean(),'.3f')
            read_lane_error_rate_100_stddev=format(summary.at(read_number).at(lane_number).error_rate_100().stddev(),'.3f')
            read_lane_error_rate_100_field=read_lane_error_rate_100_mean + '  ' + chr(177) + '  ' + read_lane_error_rate_100_stddev
            # intensity cycle 1
            read_lane_intensity_cycle_mean=format(summary.at(read_number).at(lane_number).first_cycle_intensity().mean(),'.3f') # get tiles for read 1 and lane 1
            read_lane_intensity_cycle_stddev=format(summary.at(read_number).at(lane_number).first_cycle_intensity().stddev(),'.3f')
            read_lane_intensity_cycle_field=read_lane_intensity_cycle_mean + '  ' + chr(177) + '  ' + read_lane_intensity_cycle_stddev
    

    logger.info ('Exiting the binary stats ')

local_dir_samba= '../tmp/processing'
logger=open_log()
demux_file='../tmp/processing/DemultiplexingStats.xml'
conversion_file='../tmp/processing/ConversionStats.xml'
run_processing_id=2
#xml_stats=parsing_statistics_xml(demux_file, conversion_file, logger)
#process_xml_stats(xml_stats,run_processing_id, logger)
process_binStats(local_dir_samba, run_processing_id, logger)


print ('completed')
