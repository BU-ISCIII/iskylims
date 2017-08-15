#!/usr/bin/env python3

from interop import py_interop_run_metrics, py_interop_run, py_interop_summary, py_interop_plot
#import pandas as pd
import numpy
#run_folder = r"/home/bioinfo/Documentos/practicas-CarlosIII/bcl2fastq/InterOp/ejemplo-interop/MiSeqDemo"
run_folder = r"/home/bioinfo/Documentos/practicas-CarlosIII/bcl2fastq/InterOp/nextSeq"
run_metrics = py_interop_run_metrics.run_metrics()
run_folder = run_metrics.read(run_folder)
valid_to_load = py_interop_run.uchar_vector(py_interop_run.MetricCount, 0)
py_interop_run_metrics.list_summary_metrics_to_load(valid_to_load)

#run_folder = run_metrics.read(run_folder, valid_to_load)

summary = py_interop_summary.run_summary()
py_interop_summary.summarize_run_metrics(run_metrics, summary)
'''
columns = ( ('Yield Total (G)', 'yield_g'), ('Projected Yield (G)', 'projected_yield_g'), ('% Aligned', 'percent_aligned'))
rows = [('Non-Indexed Total', summary.nonindex_summary()), ('Total', summary.total_summary())]
d = []
for label, func in columns:
    d.append( (label, pd.Series([getattr(r[1], func)() for r in rows], index=[r[0] for r in rows])))
df = pd.DataFrame.from_items(d)
print(df)
'''



run = py_interop_run_metrics.run_metrics()
options = py_interop_plot.filter_options(run.run_info().flowcell().naming_method())
bufferSize = py_interop_plot.calculate_flowcell_buffer_size(run, options)
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

