from interop import py_interop_run_metrics, py_interop_run, py_interop_summary, py_interop_plot
import os, shutil
import logging
from  iSkyLIMS_wetlab.models import RunProcess, RunningParameters, StatsRunSummary, StatsRunRead, GraphicsStats

from django.conf import settings
from iSkyLIMS_wetlab import wetlab_config
from .generic_functions import fetch_remote_file, logging_errors

def get_run_metric_files (conn, run_folder):
    '''
    Description:
        The function will collect the run metric files created by 
        sequencer as part of the run process.
    Input:
        conn # Connection samba object
        run_folder   # folder run to fetch the remote files
    Constant:
        RUN_INFO
        RUN_METRIC_FOLDER
        RUN_TEMP_DIRECTORY
        RUN_PARAMETER_NEXTSEQ
        STATISTICS_FOLDER
    Variables:
        l_metric_folder # local folder to copy the run metrics files
        l_run_info  # local copy of runInfo file
        l_run_parameter # local copy of runParamenter file
        copied_files # dictionnary of the temporary files that are copied 
        run_folder      # run folder on the remote server
        run_metrics_file_name # name of the run metric file. The value
                            is updated for each of the run metric files
                            in the folder
        s_interop_folder # local temporary folder to store the metric files
        s_metric_folder # path of the run metrics files at remote server
        s_run_info  # path of the runInfo file
        s_run_parameter # path of the runParamenter file
        statistics_folder # statistics folder on the remote server 
    Return:
        copied_files
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function get_run_metric_files')
    # runInfo needed for run metrics stats
    l_run_info = os.path.join(wetlab_config.RUN_TEMP_DIRECTORY_PROCESSING, wetlab_config.RUN_INFO)
    s_run_info = os.path.join(wetlab_config.SAMBA_APPLICATION_FOLDER_NAME, run_folder,wetlab_config.RUN_INFO)
    # runParameters needed for run metrics stats
    l_run_parameter = os.path.join(wetlab_config.RUN_TEMP_DIRECTORY_PROCESSING, wetlab_config.RUN_PARAMETER_NEXTSEQ)
    s_run_parameter = os.path.join(wetlab_config.SAMBA_APPLICATION_FOLDER_NAME, run_folder,wetlab_config.RUN_PARAMETER_NEXTSEQ)
    l_metric_folder = os.path.join(wetlab_config.RUN_TEMP_DIRECTORY_PROCESSING, wetlab_config.RUN_METRIC_FOLDER)
    s_metric_folder = os.path.join(wetlab_config.SAMBA_APPLICATION_FOLDER_NAME, run_folder, wetlab_config.RUN_METRIC_FOLDER)
    copied_files = {}
    if not os.path.exists(l_metric_folder) :
        try:
            os.makedirs(l_metric_folder)
            logger.info ('%$ was created' , wetlab_config.RUN_METRIC_FOLDER)
        except:
            string_message = "cannot create the folder" + wetlab_config.RUN_METRIC_FOLDER 
            logging_errors(logger,string_message, False , True)
            logger.debug ('End function manage_run_in_processed_bcl2fast2_run with error')
            raise 
    try:
        l_run_info = fetch_remote_file (conn, run_folder, s_run_info, l_run_info)
        logger.info('Sucessfully fetch of RunInfo file')
        
        copied_files[wetlab_config.RUN_INFO] = l_run_info
        
        l_run_parameter = fetch_remote_file (conn, run_folder, s_run_parameter, l_run_parameter)
        logger.info('Sucessfully fetch of RunParameter file')
        copied_files[wetlab_config.RUN_PARAMETER_NEXTSEQ] = l_run_parameter
        
        
        file_list = conn.listPath( wetlab_config.SAMBA_SHARED_FOLDER_NAME, s_metric_folder)
        logger.info('InterOp folder found at  %s', run_folder)

        # copy all binary files in interop folder to local  documents/wetlab/tmp/processing/interop
        copied_files[wetlab_config.RUN_METRIC_FOLDER] = []
        for sh in file_list:
            if sh.isDirectory:
                continue
            else:
                run_metrics_file_name=sh.filename
                s_run_metric_file = os.path.join(s_metric_folder, run_metrics_file_name)
                l_run_metric_file = os.path.join(l_metric_folder, run_metrics_file_name)
                l_run_metric_file = fetch_remote_file (conn, run_folder, s_run_metric_file, l_run_metric_file)
                copied_files[wetlab_config.RUN_METRIC_FOLDER].append(l_run_metric_file)
    except:
            string_message = "cannot copy files for getting run metrics"  
            logging_errors(logger,string_message, False , True)
            logger.info('Deleting temporary files')
            for key in copied_files.keys():
                if key == wetlab_config.RUN_METRIC_FOLDER :
                    for metric_file in copied_files[key]:
                        os.remove(metric_file)
                else:
                    os.remove(copied_files[key])
            logger.debug ('End function manage_run_in_processed_bcl2fast2_run with error')
            raise 
    
    logger.debug ('End function get_run_metric_files')
    return copied_files   
    


def parsing_run_metrics(run_metric_folder, run_object_name):
    '''
    Description:
        The function parse the information from the run metric files
    Input:
        run_metric_folder   # local folder with the run metric files
        run_object_name     RunProcess object for this run
    Import:
        py_interop_run
        py_interop_run_metrics
    Variables:
        bin_run_stats_summary_list # list of dictionnary with the summary
                                    information
        run_stats_read_list  # list of dictionnary with the read information
    Return:
        bin_run_stats_summary_list, run_stats_read_list
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function parsing_run_metrics')
    # get the number of lanes for the sequencer
    number_of_lanes = run_object_name.get_machine_lanes()
    # get run folder
    run_folder = RunningParameters.objects.get(runName_id = run_object_name).get_run_folder()
    # get number of reads for the run
    num_of_reads = RunningParameters.objects.get(runName_id = run_object_name).get_number_of_reads()
    logger.info('Fetched run information  needed for running metrics')
    
    
    run_metrics = py_interop_run_metrics.run_metrics()
    valid_to_load = py_interop_run.uchar_vector(py_interop_run.MetricCount, 0)
    py_interop_run_metrics.list_summary_metrics_to_load(valid_to_load)
    run_metric_folder = run_metrics.read(run_metric_folder)

    summary = py_interop_summary.run_summary()
 
    py_interop_summary.summarize_run_metrics(run_metrics, summary)

    bin_run_stats_summary_list = []
    # get the Run Summary for each Read
    for read_level in range(num_of_reads):
        run_summary_stats_level = {}
        # summary yield total
        run_summary_stats_level['yieldTotal'] = format(summary.at(read_level).summary().yield_g(),'.3f')
        # summary projected total yield
        run_summary_stats_level['projectedTotalYield'] = format(summary.at(read_level).summary().projected_yield_g(),'.3f')
        #
        # percent yield
        run_summary_stats_level['aligned'] = format(summary.at(read_level).summary().percent_aligned(),'.3f')
        # Error rate
        run_summary_stats_level['errorRate'] = format(summary.at(read_level).summary().error_rate(),'.3f')
        # intensity cycle 1
        run_summary_stats_level['intensityCycle'] = str(round(summary.at(read_level).summary().first_cycle_intensity()))
        # Q30
        run_summary_stats_level['biggerQ30'] = format(summary.at(read_level).summary().percent_gt_q30(),'.3f')
        
        run_summary_stats_level['level'] = str(read_level+1)

        bin_run_stats_summary_list.append(run_summary_stats_level)
    logger.info('Parsed run Metrics on summary level ')
    
    # get the run summary for Total
    run_summary_stats_level = {}
    # total summary
    run_summary_stats_level['yieldTotal'] = format(summary.total_summary().yield_g(),'.3f')
    # total projected_yield_g
    run_summary_stats_level['projectedTotalYield'] = format(summary.total_summary().projected_yield_g(),'.3f')
    # total percent aligned
    run_summary_stats_level['aligned'] = format(summary.total_summary().percent_aligned(),'.3f')
    # total error rate
    run_summary_stats_level['errorRate'] = format(summary.total_summary().error_rate(),'.3f')
    # total intensity cycle
    run_summary_stats_level['intensityCycle'] = str(round(summary.total_summary().first_cycle_intensity()))
    # total Q 30
    run_summary_stats_level['biggerQ30'] = format(summary.total_summary().percent_gt_q30(),'.3f')
    
    run_summary_stats_level['level'] = 'Total'
    logger.info('Parsed run Metrics on Total lane')

    bin_run_stats_summary_list.append(run_summary_stats_level)

     # get the run summary for non index
    run_summary_stats_level = {}
     # non index yield
    run_summary_stats_level['yieldTotal'] = format(summary.nonindex_summary().yield_g(),'.3f')
    #  non index projected yield
    run_summary_stats_level['projectedTotalYield'] =format(summary.nonindex_summary().projected_yield_g(),'.3f')

    # non index percent aligned
    run_summary_stats_level['aligned'] = format(summary.nonindex_summary().percent_aligned(),'.3f')
    # non index percent error rate
    run_summary_stats_level['errorRate'] = format(summary.nonindex_summary().error_rate(),'.3f')
    # non index intensity cycle
    run_summary_stats_level['intensityCycle'] = str(round(summary.nonindex_summary().first_cycle_intensity()))
    # non index Q 30
    run_summary_stats_level['biggerQ30'] = format(summary.nonindex_summary().percent_gt_q30(),'.3f')
    
    run_summary_stats_level['level'] = 'Non Index'
    logger.info('Parsed run metric for Non Index lane')

    bin_run_stats_summary_list.append(run_summary_stats_level)
    
    ### information per reads
    run_stats_read_list = []
    #lan_summary= py_interop_summary.lane_summary()
    # Tiles
    for read_number in range(num_of_reads):
        logger.info('Processing run metrics stats on Read %s', read_number)
        for lane_number in range(number_of_lanes):
            run_read_stats_level ={}
            run_read_stats_level['tiles'] = str(int(summary.at(read_number).at(lane_number).tile_count() )*2)
            # Density (k/mm2) divide the value by 1000 to have it K/mm2
            # get the +/- with the steddev
            read_lane_density_mean=str(round(float(summary.at(read_number).at(lane_number).density().mean())/1000))
            read_lane_density_stddev=str(round(float(summary.at(read_number).at(lane_number).density().stddev())/1000))
            run_read_stats_level['density'] = read_lane_density_mean + '  ' + chr(177) + '  ' +read_lane_density_stddev
            # cluster _pf  in %
            read_lane_percent_pf_mean=format(summary.at(read_number).at(lane_number).percent_pf().mean(),'.3f')
            read_lane_percent_pf_stddev=format(summary.at(read_number).at(lane_number).percent_pf().stddev(),'.3f')
            run_read_stats_level['cluster_PF']  = read_lane_percent_pf_mean + '  ' + chr(177) + '  ' +read_lane_percent_pf_stddev
            # phas/ prepas in %
            read_lane_phasing_mean=format(summary.at(read_number).at(lane_number).phasing().mean(),'.3f')
            read_lane_phasing_dev=format(summary.at(read_number).at(lane_number).phasing().stddev(),'.1f')
            read_lane_prephasing_mean=format(summary.at(read_number).at(lane_number).prephasing().mean(),'.3f')
            read_lane_prephasing_stddev=format(summary.at(read_number).at(lane_number).prephasing().stddev(),'.3f')
            run_read_stats_level['phas_prephas']  = read_lane_phasing_mean + '  ' + chr(177) + '  ' + read_lane_phasing_dev + '  /  ' + read_lane_prephasing_mean + '  ' + chr(177) + '  ' + read_lane_prephasing_stddev
            # reads (M)
            run_read_stats_level['reads'] = format(float(summary.at(read_number).at(lane_number).reads())/1000000,'.3f')
            #reads PF (M)
            run_read_stats_level['reads_PF'] = format(float(summary.at(read_number).at(lane_number).reads_pf())/1000000,'.3f')
            # percent q30
            run_read_stats_level['q30'] = format(summary.at(read_number).at(lane_number).percent_gt_q30(),'.3f')
            # yield _g
            run_read_stats_level['yields'] = format(summary.at(read_number).at(lane_number).yield_g(),'.3f')
            # cycles err Rate
            run_read_stats_level['cyclesErrRated'] = str(summary.at(read_number).at(lane_number).cycle_state().error_cycle_range().first_cycle())
            #percent_aligned
            read_lane_percent_aligned_mean=format(summary.at(read_number).at(lane_number).percent_aligned().mean(),'.3f')
            read_lane_percent_aligned_stddev=format(summary.at(read_number).at(lane_number).percent_aligned().stddev(),'3f')
            run_read_stats_level['aligned'] = read_lane_percent_aligned_mean + '  ' + chr(177) + '  ' + read_lane_percent_aligned_stddev
            #error rate
            read_lane_error_rate_mean=format(summary.at(read_number).at(lane_number).error_rate().mean(),'.3f')
            read_lane_error_rate_stddev=format(summary.at(read_number).at(lane_number).error_rate().stddev(),'.3f')
            run_read_stats_level['errorRate']  = read_lane_error_rate_mean+ '  ' + chr(177) + '  ' + read_lane_error_rate_stddev
            #error rate_35
            read_lane_error_rate_35_mean=format(summary.at(read_number).at(lane_number).error_rate_35().mean(),'.3f')
            read_lane_error_rate_35_stddev=format(summary.at(read_number).at(lane_number).error_rate_35().stddev(),'.3f')
            run_read_stats_level['errorRate35']  = read_lane_error_rate_35_mean + '  ' + chr(177) + '  ' + read_lane_error_rate_35_stddev
            #error rate 50
            read_lane_error_rate_50_mean=format(summary.at(read_number).at(lane_number).error_rate_50().mean(),'.3f')
            read_lane_error_rate_50_stddev=format(summary.at(read_number).at(lane_number).error_rate_50().stddev(),'.3f')
            run_read_stats_level['errorRate50']  = read_lane_error_rate_50_mean + '  ' + chr(177) + '  ' + read_lane_error_rate_50_stddev
            #error rate 75
            read_lane_error_rate_75_mean=format(summary.at(read_number).at(lane_number).error_rate_75().mean(),'.3f')
            read_lane_error_rate_75_stddev=format(summary.at(read_number).at(lane_number).error_rate_75().stddev(),'.3f')
            run_read_stats_level['errorRate75']  = read_lane_error_rate_75_mean + '  ' + chr(177) + '  ' + read_lane_error_rate_75_stddev
            #error rate 100
            read_lane_error_rate_100_mean=format(summary.at(read_number).at(lane_number).error_rate_100().mean(),'.3f')
            read_lane_error_rate_100_stddev=format(summary.at(read_number).at(lane_number).error_rate_100().stddev(),'.3f')
            run_read_stats_level['errorRate100']  = read_lane_error_rate_100_mean + '  ' + chr(177) + '  ' + read_lane_error_rate_100_stddev
            # intensity cycle 1
            read_lane_intensity_cycle_mean=format(summary.at(read_number).at(lane_number).first_cycle_intensity().mean(),'.3f') # get tiles for read 1 and lane 1
            read_lane_intensity_cycle_stddev=format(summary.at(read_number).at(lane_number).first_cycle_intensity().stddev(),'.3f')
            run_read_stats_level['intensityCycle'] = read_lane_intensity_cycle_mean + '  ' + chr(177) + '  ' + read_lane_intensity_cycle_stddev

            run_read_stats_level['read'] = str(read_number+1)
            run_read_stats_level['lane'] = str(lane_number+1)
            # append run_read_stats_level information to run_stats_read_list 
            run_stats_read_list.append(run_read_stats_level)
    logger.info ('End function parsing_run_metrics')
    return bin_run_stats_summary_list, run_stats_read_list

def create_graphics(run_metric_folder,run_object_name):
    '''
    Description:
        The function create an entry on database with the run graphics 
        by using the run metrics files
    Input:
        run_metric_folder   # local folder with the run metric files
        run_object_name     RunProcess object for this run
    Import:
        py_interop_run
        py_interop_run_metrics
    Variables:
        
    Return:
        True
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function create_graphics')
    experiment_name = run_object_name.get_run_name()
    run_folder = RunningParameters.objects.get(runName_id = run_object_name).get_run_folder()
    
    graphic_list=['plot_by_cycle  ', 'plot_by_lane  ', 'plot_flowcell  ', 
                    'plot_qscore_histogram  ', 'plot_qscore_heatmap  ', 'plot_sample_qc  ' ]
    # create the graphics
    logger.info('Creating plot graphics for run id ')
    
    for graphic in graphic_list:
        graphic_command = os.path.join(wetlab_config.INTEROP_PATH, graphic )
        plot_command= graphic_command + wetlab_config.RUN_TEMP_DIRECTORY_PROCESSING + '  | gnuplot'
        logger.debug('command used to create graphic is : %s', plot_command)
        os.system(plot_command)
    run_graphic_dir=os.path.join(settings.MEDIA_ROOT,wetlab_config.RUN_IMAGES_DIRECTORY, run_folder)
    if os.path.exists(run_graphic_dir):
        shutil.rmtree(run_graphic_dir)
    os.mkdir(run_graphic_dir)
    logger.info('created new directory %s', run_graphic_dir)

    #move the graphic files to wetlab directory
    source = os.listdir()
    for files in source:
        if files.endswith(wetlab_config.PLOT_EXTENSION):
            logger.debug('moving file %s', files)
            shutil.move(files,run_graphic_dir)

    #removing the processing_ character in the file names
    graphic_files = os.listdir(run_graphic_dir)
    logger.info('Renaming the graphic files')

    for graphic_file in graphic_files :
        old_file_name = os.path.join(run_graphic_dir, graphic_file)
        split_file_name = graphic_file.split('_')
        if not split_file_name[1].endswith(wetlab_config.PLOT_EXTENSION):
            split_file_name[1] = split_file_name[1] + wetlab_config.PLOT_EXTENSION
        new_file_name = os.path.join(run_graphic_dir, split_file_name[1])
        os.rename(old_file_name , new_file_name)
        logger.debug('Renamed file from %s to %s', old_file_name, new_file_name)


    # saving the graphic location in database
    ns_graphic_stats= GraphicsStats (runprocess_id = RunProcess.objects.get(runName__exact = experiment_name),
                                            folderRunGraphic= run_folder, cluserCountGraph = 'ClusterCount-by-lane.png',
                                            flowCellGraph= 'flowcell-Intensity.png', intensityByCycleGraph = 'Intensity-by-cycle.png',
                                            heatMapGraph= 'q-heat-map.png', histogramGraph= 'q-histogram.png',
                                            sampleQcGraph= 'sample-qc.png')
    ns_graphic_stats.save()

    logger.info('Store Graphic plots in database')
    logger.debug ('End function create_graphics')
    return True


