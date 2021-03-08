from interop import py_interop_run_metrics, py_interop_run, py_interop_summary, py_interop_plot
import os, shutil
import logging
from  iSkyLIMS_wetlab.models import RunProcess, RunningParameters, StatsRunSummary, StatsRunRead, GraphicsStats

from django.conf import settings
from iSkyLIMS_wetlab.wetlab_config import *

from .handling_crontab_common_functions import *


def create_run_metric_graphics(run_metric_folder,run_process_obj, run_folder,experiment_name):
    '''
    Description:
        The function create an entry on database with the run graphics
        by using the run metrics files
    Input:
        run_metric_folder   # local folder with the run metric files
        run_process_obj     # RunProcess object for this run
        run_folder          # run folder to store the figures
        experiment_name     # Experiment name
    Constants:
        RUN_METRIC_GRAPHIC_COMMANDS
        INTEROP_PATH
        RUN_IMAGES_DIRECTORY
        MEDIA_ROOT
    Return:
        graphic_stats_obj
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting create_run_metric_graphics', experiment_name)
    present_working_dir = os.getcwd()
    run_graphic_dir=os.path.join(settings.MEDIA_ROOT,RUN_IMAGES_DIRECTORY, run_folder)
    try:
        if os.path.exists(run_graphic_dir):
            shutil.rmtree(run_graphic_dir)
        os.mkdir(run_graphic_dir)
        logger.info('%s : created new directory %s',experiment_name, run_graphic_dir)
    except:
        string_message = experiment_name + ' : Unable to create folder to store graphics on ' + run_graphic_dir
        logging_errors(string_message, True, True)
        logger.debug ('%s : End create_run_metric_graphics with exception', experiment_name)
        return {'ERROR':28}
    os.chdir(run_graphic_dir)
    logger.info('%s : Changed working direcory to copy run metric graphics', experiment_name)
    # create the graphics
    logger.info('%s : Creating plot graphics for run id ',experiment_name)

    full_path_run_processing_tmp = os.path.join(present_working_dir, RUN_TEMP_DIRECTORY_PROCESSING)
    for graphic in RUN_METRIC_GRAPHIC_COMMANDS:
        graphic_command = os.path.join(INTEROP_PATH, graphic )
        plot_command= graphic_command + full_path_run_processing_tmp + '  | gnuplot'
        logger.debug('%s : command used to create graphic is : %s',experiment_name, plot_command)
        os.system(plot_command)

    os.chdir(present_working_dir)
    logger.info('%s : Returning back the working directory', experiment_name)

    #removing the processing_ character in the file names
    graphic_files = os.listdir(run_graphic_dir)
    logger.info('%s : Renaming the graphic files',experiment_name)

    for graphic_file in graphic_files :
        old_file_name = os.path.join(run_graphic_dir, graphic_file)
        split_file_name = graphic_file.split('_')
        if not split_file_name[1].endswith(PLOT_EXTENSION):
            split_file_name[1] = split_file_name[1] + PLOT_EXTENSION
        new_file_name = os.path.join(run_graphic_dir, split_file_name[1])
        os.rename(old_file_name , new_file_name)
        logger.debug('%s : Renamed file from %s to %s', experiment_name, old_file_name, new_file_name)

    # saving the graphic location in database
    graphic_stats_obj= GraphicsStats.objects.create_graphic_run_metrics(run_process_obj,run_folder )

    logger.info('%s : Store Graphic plots in database',experiment_name)
    logger.debug ('%s : End function create_graphics',experiment_name)
    return {'graph_stats_obj':graphic_stats_obj}

def delete_existing_run_metrics_table_processed (run_process_obj, experiment_name) :
    '''
    Description:
        The function will check if exists data stored on StatsRunSummary  and/or
        StatsRunRead for the run
    Input:
        run_process_obj     # runProcess object
        experiment_name     # experiment name
    Return:
        None
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting function check_run_metrics_processed', experiment_name)
    if StatsRunSummary.objects.filter(runprocess_id = run_process_obj).exists():
        run_summary_objs = StatsRunSummary.objects.filter(runprocess_id = run_process_obj)
        for run_summary_obj in run_summary_objs:
            run_summary_obj.delete()
        logger.info('%s : Deleted rows on StatsRunSummary table', experiment_name )

    if StatsRunRead.objects.filter(runprocess_id = run_process_obj).exists():
        run_read_objs = StatsRunRead.objects.filter(runprocess_id = run_process_obj)
        for run_read_obj in run_read_objs:
            run_read_obj.delete()
        logger.info('%s : Deleted rows on StatsRunSummary table', experiment_name )

    if GraphicsStats.objects.filter(runprocess_id = run_process_obj).exists():
        graph_stats_objs = GraphicsStats.objects.filter(runprocess_id = run_process_obj)
        for graph_stats_obj in graph_stats_objs:
            graph_stats_obj.delete()
        logger.info('%s : Deleted rows on GraphicsStats table', experiment_name )
    logger.debug('%s : End function check_run_metrics_processed', experiment_name)
    return None


def delete_run_metric_files (experiment_name):
    '''
    Description:
        The function delete the files used for collecting the run metrics
    Input:
        experiment_name     # experiment name
    Return:
        None
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting function delete_run_metric_files', experiment_name)
    local_metric_folder = os.path.join(RUN_TEMP_DIRECTORY_PROCESSING, RUN_METRIC_FOLDER)
    l_run_parameter = os.path.join(RUN_TEMP_DIRECTORY_PROCESSING, RUN_PARAMETER_FILE)
    l_run_info = os.path.join(RUN_TEMP_DIRECTORY_PROCESSING, RUN_INFO)
    local_files = [l_run_parameter, l_run_info]
    for local_file in local_files:
        if os.path.exists(local_file):
            try:
                os.remove(local_file)
            except:
                string_message = experiment_name + ' : Unable to delete ' + local_file
                logging_errors(string_message, True, True)
                continue
    logger.info('%s : Deleted temporary files',experiment_name)
    if os.path.exists(local_metric_folder):
        try:
            shutil.rmtree(local_metric_folder)
            logger.info('%s : Deleted Folder %s',experiment_name, local_metric_folder)
        except:
            string_message = experiment_name + ' : Unable to delete  folder ' + local_metric_folder
            logging_errors(string_message, True, True)

    logger.debug ('%s : End function delete_run_metric_files',experiment_name)
    return


def get_run_metric_files (conn, run_folder, experiment_name):
    '''
    Description:
        The function will collect the run metric files created by sequencer as part of the run process.
    Input:
        conn # Connection samba object
        run_folder   # folder run to fetch the remote files
        experiment_name # experiment name
    Constant:
        RUN_INFO
        RUN_METRIC_FOLDER
        RUN_TEMP_DIRECTORY
        RUN_PARAMETER_FILE
        STATISTICS_FOLDER
    Functions:
        get_samba_application_shared_folder     # Located at utils/handling_crontab_common_functions.py
        fetch_remote_file                       # Located at utils.handling_crontab_common_functions.py
    Return:
        copied_files
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting function get_run_metric_files', experiment_name)

    # runInfo needed for run metrics stats
    l_run_info = os.path.join(RUN_TEMP_DIRECTORY_PROCESSING, RUN_INFO)
    s_run_info = os.path.join(get_samba_application_shared_folder(), run_folder, RUN_INFO)
    # runParameters needed for run metrics stats
    l_run_parameter = os.path.join(RUN_TEMP_DIRECTORY_PROCESSING, RUN_PARAMETER_FILE)
    s_run_parameter = os.path.join(get_samba_application_shared_folder(), run_folder,RUN_PARAMETER_FILE)
    l_metric_folder = os.path.join(RUN_TEMP_DIRECTORY_PROCESSING, RUN_METRIC_FOLDER)
    s_metric_folder = os.path.join(get_samba_application_shared_folder(), run_folder, RUN_METRIC_FOLDER)
    copied_files = {}

    if not os.path.exists(l_metric_folder) :
        try:
            os.makedirs(l_metric_folder)
            logger.info ('%s : Created folder %s' , experiment_name, l_metric_folder)
        except:
            string_message = experiment_name + " : cannot create folder on " + l_metric_folder
            logging_errors(string_message, True , True)
            logger.debug ('%s : End function get_run_metric_files with error',experiment_name)
            return {'ERROR':26}

    try:
        l_run_info = fetch_remote_file (conn, run_folder, s_run_info, l_run_info)
        logger.info('%s : Sucessfully fetch of RunInfo file',experiment_name)
    except:
        string_message = experiment_name + ' : Unable to fetch ' + s_run_info
        logging_errors(string_message, True , True)
        logger.debug ('%s : End function get_run_metric_files with error',experiment_name)
        return {'ERROR':20}
    copied_files[RUN_INFO] = l_run_info

    try:
        l_run_parameter = fetch_remote_file (conn, run_folder, s_run_parameter, l_run_parameter)
        logger.info('%s : Sucessfully fetch of RunParameter file',experiment_name)
    except:
        string_message = experiment_name + ' : Unable to fetch ' + s_run_parameter
        logging_errors(string_message, True , True)
        logger.debug ('%s : End function get_run_metric_files with error',experiment_name)
        return {'ERROR':21}
    copied_files[RUN_PARAMETER_FILE] = l_run_parameter

    try:
        file_list = conn.listPath( get_samba_shared_folder(), s_metric_folder)
        logger.info('%s : InterOp folder found at  %s', experiment_name, s_metric_folder)
    except:
        string_message = experiment_name + ' : Unable to fetch ' + s_run_parameter
        logging_errors(string_message, True , True)
        shutil.rmtree(l_metric_folder)
        logger.debug ('%s : End function get_run_metric_files with error',experiment_name)
        return {'ERROR':27}
    # copy all binary files in interop folder to local  documents/wetlab/tmp/processing/interop
    copied_files[RUN_METRIC_FOLDER] = []
    try:
        for sh in file_list:
            if sh.isDirectory:
                continue
            else:
                run_metrics_file_name=sh.filename
                s_run_metric_file = os.path.join(s_metric_folder, run_metrics_file_name)
                l_run_metric_file = os.path.join(l_metric_folder, run_metrics_file_name)
                l_run_metric_file = fetch_remote_file (conn, run_folder, s_run_metric_file, l_run_metric_file)
                # copied_files[RUN_METRIC_FOLDER].append(l_run_metric_file)
    except:
        string_message = experiment_name + ' : Unable to fetch ' + s_run_metric_file
        logging_errors(string_message, True , True)
        shutil.rmtree(l_metric_folder)
        logger.debug ('%s : End function get_run_metric_files with error',experiment_name)
        return {'ERROR':27}

    logger.debug ('%s : End function get_run_metric_files', experiment_name)
    return copied_files


def parsing_run_metrics_files(local_run_metric_folder, run_process_obj, experiment_name):
    '''
    Description:
        The function parse the information from the run metric files
    Input:
        local_run_metric_folder   # local folder with the run metric files
        run_process_obj           # RunProcess object for this run
        experiment_name           # experiment name
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
    logger.debug ('%s : Starting function parsing_run_metrics',experiment_name)
    # get the number of lanes for the sequencer
    number_of_lanes = run_process_obj.get_sequencing_lanes()
    # get number of reads for the run
    num_of_reads = RunningParameters.objects.get(runName_id = run_process_obj).get_number_of_reads()
    logger.info('%s : Fetched run information  needed for running metrics',experiment_name)

    run_metrics = py_interop_run_metrics.run_metrics()
    valid_to_load = py_interop_run.uchar_vector(py_interop_run.MetricCount, 0)
    py_interop_run_metrics.list_summary_metrics_to_load(valid_to_load)
    run_metric_folder = run_metrics.read(local_run_metric_folder)

    summary = py_interop_summary.run_summary()

    py_interop_summary.summarize_run_metrics(run_metrics, summary)

    bin_run_stats_summary_list = []
    logger.info('%s : Starts collecting data for run metric ', experiment_name)
    # get the Run Summary for each Read
    for read_level in range(num_of_reads):
        run_summary_stats_level = {}
        # summary yield total
        run_summary_stats_level['yieldTotal'] = format(summary.at(read_level).summary().yield_g(),'.3f')
        # summary projected total yield
        run_summary_stats_level['projectedTotalYield'] = format(summary.at(read_level).summary().projected_yield_g(),'.3f')

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
    logger.info('%s : Parsed run Metrics on summary level ', experiment_name)

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

    logger.info('%s : Parsed run Metrics on Total lane',experiment_name)

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
    logger.info('%s : Parsed run metric for Non Index lane', experiment_name)

    bin_run_stats_summary_list.append(run_summary_stats_level)

    ### information per reads
    run_stats_read_list = []
    #lan_summary= py_interop_summary.lane_summary()
    # Tiles
    for read_number in range(num_of_reads):
        logger.info('%s : Processing run metrics stats on Read %s',experiment_name, read_number)
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

    logger.debug ('%s : End function parsing_run_metrics',experiment_name)
    return bin_run_stats_summary_list, run_stats_read_list
