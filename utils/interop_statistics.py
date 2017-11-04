#!/usr/bin/env python3

from interop import py_interop_run_metrics, py_interop_run, py_interop_summary, py_interop_plot
#import pandas as pd
import os, shutil
from  ..models import *

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
                                                  level=str(read_level+1), yieldTotal = read_summary_yield_g,
                                                  projectedTotalYield= read_summary_projected_yield_g,
                                                  aligned= read_summary_percent_aligned, errorRate= read_summary_error_rate,
                                                  intensityCycle= read_summary_first_cycle_intensity, biggerQ30= read_percent_gt_q30)
        ns_bin_run_summary.save()
        
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
    ns_bin_run_summary.save()
     
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

    ns_bin_run_summary.save()
    
    ### information per reads
    
    #lan_summary= py_interop_summary.lane_summary()
    # Tiles
    for read_number in range(4):
        logger.info('Processing bin stats for Read %s', read_number)
        for lane_number in range(4):
            
            read_lane_tiles=str(int(summary.at(read_number).at(lane_number).tile_count() )*2)
            # Density (k/mm2) divide the value by 1000 to have it K/mm2
            # get the +/- with the steddev
            read_lane_density_mean=str(round(float(summary.at(read_number).at(lane_number).density().mean())/1000))
            read_lane_density_stddev=str(round(float(summary.at(read_number).at(lane_number).density().stddev())/1000))
            read_lane_density_field=read_lane_density_mean + '  ' + chr(177) + '  ' +read_lane_density_stddev
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
    
            ns_bin_read_lane = NextSeqStatsBinRunRead (runprocess_id=RunProcess.objects.get(pk=run_id),
                                                       read= str(read_number+1), lane = str(lane_number+1),
                                                       tiles= read_lane_tiles, density= read_lane_density_field,
                                                       cluster_PF= read_lane_percent_pf_field, phas_prephas= read_lane_phas_prephas_field,
                                                       reads= read_lane_reads, reads_PF= read_lane_reads_pf,
                                                       q30= read_lane_percent_gt_q30, yields= read_lane_yield_g,
                                                       cyclesErrRated= read_lane_cycles_error_rate, aligned= read_lane_percent_aligned_field,
                                                       errorRate= read_lane_error_rate_field, errorRate35= read_lane_error_rate_35_field,
                                                       errorRate50= read_lane_error_rate_50_field, errorRate75= read_lane_error_rate_75_field,
                                                       errorRate100= read_lane_error_rate_100_field, intensityCycle= read_lane_intensity_cycle_field)
            ns_bin_read_lane.save()

    logger.info ('Exiting the binary stats ')
    
def create_graphics(run_folder,run_id, graphic_dir, logger):
    graphic_list=['plot_by_cycle  ', 'plot_by_lane  ', 'plot_flowcell  ', 'plot_qscore_histogram  ',
                  'plot_qscore_heatmap  ', 'plot_sample_qc  ' ]

    # create the graphics
    logger.info('Creating plot graphics for run id %s', run_id)
    for item_graphic in graphic_list:
        plot_command= item_graphic + run_folder + '  | gnuplot'
        os.system(plot_command)
        
    run_graphic_dir=os.path.join('documents/wetlab/images_plot', graphic_dir)
    if not os.path.exists(run_graphic_dir):
        os.mkdir(run_graphic_dir)
        logger.info('created new directory %s', run_graphic_dir)
    #move the graphic files to wetlab directory
    source = os.listdir("./")
    for files in source:
        if files.endswith(".png"):
            shutil.move(files,run_graphic_dir)
            
    #removing the processing_ character in the file names
    source = os.listdir(run_graphic_dir)
    logger.info('Renaming the graphic files')
    for files in source:
        move_file=files.split('_')
        old_file=os.path.join(run_graphic_dir,files)
        new_file=os.path.join(run_graphic_dir,move_file[1])
        shutil.move(old_file,new_file)
    
    folder_run_graphic=graphic_dir.replace('../','')
    # saving the graphic location in database
    ns_graphic_stats= NextSeqGraphicsStats (runprocess_id=RunProcess.objects.get(pk=run_id),
                                            folderRunGraphic= folder_run_graphic, cluserCountGraph = 'ClusterCount-by-lane.png',
                                            flowCellGraph= 'flowcell-Intensity.png', intensityByCycleGraph = 'Intensity-by-cycle.png',
                                            heatMapGraph= 'q-heat-map.png', histogramGraph= 'q-histogram.png',
                                            sampleQcGraph= 'sample-qc.png')
    ns_graphic_stats.save()
    logger.info('Graphic plots saved on database')



'''
columns = ( ('Yield Total (G)', 'yield_g'), ('Projected Yield (G)', 'projected_yield_g'), ('% Aligned', 'percent_aligned'))
rows = [('Non-Indexed Total', summary.nonindex_summary()), ('Total', summary.total_summary())]
d = []
for label, func in columns:
    d.append( (label, pd.Series([getattr(r[1], func)() for r in rows], index=[r[0] for r in rows])))
df = pd.DataFrame.from_items(d)
print(df)
'''



#run_folder = r"/home/bioinfo/Documentos/practicas-CarlosIII/bcl2fastq/InterOp/ejemplo-interop/MiSeqDemo"
#run_folder = r"/home/bioinfo/Documentos/practicas-CarlosIII/bcl2fastq/InterOp/nextSeq"
#process_binStats(run_folder)



'''
run = py_interop_run_metrics.run_metrics()


dataBuffer = numpy.zeros(bufferSize, dtype=numpy.float32)
idBuffer = numpy.zeros(bufferSize, dtype=numpy.uint32)
data = py_interop_plot.flowcell_data()

py_interop_plot.plot_flowcell_map2(run, py_interop_run.Intensity, options, data, dataBuffer, idBuffer)

read_index=0
yield_g = summary.at(read_index).summary().yield_g()
error_rate = summary.at(read_index).summary().error_rate()

first_cycle = summary.at(read_index).summary().first_cycle_intensity()
percent_aligned = summary.at(read_index).summary().percent_aligned()
percent_q30 = summary.at(read_index).summary().percent_gt_q30()

sum_yield = summary.total_summary().yield_g()
sum_non_index_yield = summary.nonindex_summary().yield_g()



print('yield_g = ', yield_g, '  error_rate =  ', error_rate)
value_1 = float(8)
value_2 = float(3)
line_sum = py_interop_summary.lane_summary()
#line_aligned = line_sum.at(read_index).summary().percent_aligned()
metric_stat = py_interop_summary.metric_stat()
#stat_sum = py_interop_summary.stat_summary()
#mean= metric_stat.mean(value_1, value_2)
print('hello')


####codigo de ejemplo
from interop.py_interop_run_metrics import run_metrics as RunMetrics
 
run_dir = r"/home/bioinfo/Documentos/practicas-CarlosIII/bcl2fastq/InterOp/nextSeq"
run_metrics = RunMetrics()
run_metrics.read(run_dir)
 
# Run metrics contains 
# corrected_intensity_metric_set
# error_metric_set
# extraction_metric_set
# image_metric_set
# index_metric_set
# q_by_lane_metric_set
# q_collapsed_metric_set
# q_metric_set
 
q_metric_set = run_metrics.q_metric_set()
# q_metric_set contains metrics but is not a python list 
# it has a size() that give the length of the list
# it contains one metric for each lane/tile/cycle
 
metrics = q_metric_set.metrics_for_cycle(1)
metrics = q_metric_set.metrics_for_lane(1)
metrics = q_metric_set.metrics()
# These function return whole or sub set of metrics which can be subscripted as follow
metrics[0]
# <interop.py_interop_metrics.q_metric; proxy of <Swig Object of type 'std::vector< illumina::interop::model::metrics::q_metric >::value_type *' at 0x2aaab987f1e0> >
# this is a metrics object described here http://illumina.github.io/interop/group__q__metric.html
metrics[0].qscore_hist()
# (0, 109812, 0, 165009, 4582040, 0, 0)
# return number of cluster for each quality score bin
metrics[0].lane()
# 1
# the lane number for this metric
metrics[0].tile()
# 1103
# the tile id for this metric
metrics[0].cycle()
# 1
# the cycle for this metric

############ otro ejemplo de  summarise_metrics_per_tile.py 

from interop.py_interop_run_metrics import run_metrics as RunMetrics
run_dir=r"/home/bioinfo/Documentos/practicas-CarlosIII/bcl2fastq/InterOp/nextSeq"

def calc_avg_q(q_hist):
	q8, q12, q22, q27, q32, q37, q41 = q_hist
	d = q8 * 8 + q12 * 12 + q22 * 22 + q27 * 27 + q32 * 32 + q37 * 37 + q41 * 41
	return d/sum(q_hist)

def sum_q_hist(q_hist1, q_hist2):
  q8, q12, q22, q27, q32, q37, q41 = q_hist1
  tq8, tq12, tq22, tq27, tq32, tq37, tq41 = q_hist2
  return (q8+tq8, q12+tq12, q22+tq22, q27+tq27, q32+tq32, q37+tq37, q41+tq41)

def get_run_metrics_per_tile(run_dir, start_cycle=0, end_cycle=310):
  run_metrics = RunMetrics()
  run_metrics.read(run_dir)
  q_metric_set = run_metrics.q_metric_set()
  metrics = q_metric_set.metrics()

  all_lanes = {}
  for lane in range(1,9):
    all_lanes[lane] = {}
    
  for i in range(metrics.size()):
    tile_set = all_lanes[metrics[i].lane()]
    if metrics[i].cycle() > start_cycle and metrics[i].cycle() < end_cycle:
      q_hist1 = tile_set.get(metrics[i].tile(), (0, 0, 0, 0, 0, 0, 0))
      q_hist2 = metrics[i].qscore_hist()
      tile_set[metrics[i].tile()] = sum_q_hist(q_hist1, q_hist2)
  return all_lanes

run_dir = r"/home/bioinfo/Documentos/practicas-CarlosIII/bcl2fastq/InterOp/nextSeq"
all_lanes = get_run_metrics_per_tile(run_dir)
lane = 4 
for tile in sorted(all_lanes[lane]):
    print(lane, tile, calc_avg_q(all_lanes[lane][tile]))
    
 ###############################################################   

from interop import py_interop_run_metrics
from interop import py_interop_summary
from interop import py_interop_run
import numpy
import logging
import sys
import os
  
  
def main():
    """ Retrieve run folder paths from the command line
      Ensure only metrics required for summary are loaded
      Load the run metrics
      Calculate the summary metrics
      Display error by lane, read
      """
    logging.basicConfig(level=logging.INFO)
 
    run_metrics = py_interop_run_metrics.run_metrics()
    summary = py_interop_summary.run_summary()

    valid_to_load = py_interop_run.uchar_vector(py_interop_run.MetricCount, 0)
    py_interop_run_metrics.list_summary_metrics_to_load(valid_to_load)
 
    for run_folder_path in sys.argv[1:]:
        run_folder = os.path.basename(run_folder_path)
        try:
            run_metrics.read(run_folder_path, valid_to_load)
        except Exception:
            logging.warn("Skipping - cannot read RunInfo.xml: %s - %s"%(run_folder, str(ex)))
            continue
        py_interop_summary.summarize_run_metrics(run_metrics, summary)
        error_rate_read_lane_surface = numpy.zeros((summary.size(), summary.lane_count(), summary.surface_count()))
        for read_index in xrange(summary.size()):
            for lane_index in xrange(summary.lane_count()):
                for surface_index in xrange(summary.surface_count()):
                    error_rate_read_lane_surface[read_index, lane_index, surface_index] = \
                        summary.at(read_index).at(lane_index).at(surface_index).error_rate().mean()
        logging.info("Run Folder: "+run_folder)
        for read_index in xrange(summary.size()):
            read_summary = summary.at(read_index)
            logging.info("Read "+str(read_summary.read().number())+" - Top Surface Mean Error: "+str(error_rate_read_lane_surface[read_index, :, 0].mean()))
 
if __name__ == '__main__':
    main()




    
    options = py_interop_plot.filter_options(run_metrics.run_info().flowcell().naming_method())
    bufferSize = py_interop_plot.calculate_flowcell_buffer_size(run_metrics, options)
       
    print(options)
    print(bufferSize)
'''
