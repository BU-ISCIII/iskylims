
import statistics
import os
from datetime import datetime
from ..fusioncharts.fusioncharts import FusionCharts

from django.conf import settings
from django.contrib.auth.models import Group
from iSkyLIMS_core.models import SequencerInLab
from iSkyLIMS_wetlab.models import *
from iSkyLIMS_wetlab.wetlab_config import *
from .generic_functions import normalized_data, get_run_in_same_year_to_compare
from iSkyLIMS_wetlab.utils.sample_functions import get_sample_in_project_obj_from_id
from .stats_graphics import *


def get_boxplot_comparation_runs (run_object):
    '''
    Description:
        The function collect information run to compare with the run
        executed on the same year and with the same chemistry
    Input:
        run_object      # contains the runProcess object
    functions:
        normalized_data # located at utils.generic_functions
        bloxplot_graphic # located at utils.
        get_run_in_same_year_to_compare # located at utils.generic_functions
    Variables:
        categories          # category list of the data to display
        chem_high_mid       # chemistry value of the run to compare
                            runs with the same value
        data_source         # data in json format
        run_year            # contains the year of the selected run
        same_run_in_year    # contains all the runs in the year of the
                            selected run excluding those which do not
                            have the same chemistry value
        same_runs_in_year_list # contains the id list for all runs that match
        series              # tupla list to assign comparation color in the graphic
        start_date          # it is the first of january of the run_year
        end_date            # it is the 31st of december of the run_year
        ################################################################
        q_30_run_value        # q30 values of the run
        q_30_run_value_float  # q30 _value in float format
        q_30_all_value_float  # q30 values for all runs in the comparations
        q_30_run_normalized   # q30 normalized value for run
        q_30_all_normalized   # q30 normalized value for all runs
        q30_run_str     # q_30 normalized value in string format
        q_30_all_str    # q_30_all_normalized in string format
        ################################################################
        #### same structure for mean_value, yield_mb and cluster_pf
        ################################################################
        cluster_pf_run_value_float
        cluster_pf_all_value_float

        mean_value_list
        mean_run_value
        mean_run_value_float
        mean_all_value_float

        yield_mb_run_value_float
        yield_mb_all_value_float
        ################################################################
    return:
        FusionCharts object with the graphic data
    '''
    # fetch Q>30 , mean_q and yield mb for all projects per lane to create the boxplot
    if not StatsLaneSummary.objects.filter(runprocess_id__exact =run_object ).exclude(defaultAll__isnull = False).exists():
        # return empty information . No information for stats stored on database
        return
    run_lane_summary = StatsLaneSummary.objects.filter(runprocess_id__exact =run_object ).exclude(defaultAll__isnull = False)
    q_30_run_value , mean_run_value = [] , []
    q_30_run_value_float , mean_run_value_float , yield_mb_run_value_float, cluster_pf_run_value_float = [] , [], [] , []
    q_30_all_value_float , mean_all_value_float , yield_mb_all_value_float , cluster_pf_all_value_float = [] , [] , [], []

    for item in run_lane_summary:
        q_30_value, mean_value , yield_mb_value , cluster_pf_value, = item.get_stats_info()

        q_30_run_value_float.append(float(q_30_value))
        mean_run_value_float.append(float(mean_value))
        yield_mb_run_value_float.append(float(yield_mb_value.replace(',','')))
        cluster_pf_run_value_float.append(float(cluster_pf_value.replace(',','')))

    # get the chemistry type for the run, that will be used to compare runs with the same chemistry value
    '''
    chem_high_mid = RunningParameters.objects.get(runName_id__exact = run_object).Chemistry
    run_different_chemistry = RunningParameters.objects.all(). exclude(Chemistry__exact = chem_high_mid)
    run_year = run_object.run_date.timetuple().tm_year

    start_date = str(run_year) + '-1-1'
    end_date = str(run_year) +'-12-31'
    same_run_in_year = RunProcess.objects.filter(run_date__range=(start_date, end_date)).exclude(runName__in = run_different_chemistry)
    '''
    same_run_in_year = get_run_in_same_year_to_compare(run_object)
    same_runs_in_year_list = []
    for run in same_run_in_year :
        same_runs_in_year_list.append(run.get_run_id())

    all_lane_summary = StatsLaneSummary.objects.filter(runprocess_id__in = same_runs_in_year_list).exclude(defaultAll__isnull = False).exclude(runprocess_id__exact =run_object)
    if len(all_lane_summary) == 0 :
        # It is the first run in the year. Then include it until more than one run was stored
        all_lane_summary = StatsLaneSummary.objects.filter(runprocess_id__in = same_runs_in_year_list).exclude(defaultAll__isnull = False)
    for item in all_lane_summary:
        q_30_value, mean_value , yield_mb_value , cluster_pf_value = item.get_stats_info()

        q_30_all_value_float.append(float(q_30_value))
        mean_all_value_float.append(float(mean_value))
        yield_mb_all_value_float.append(float(yield_mb_value.replace(',','')))
        cluster_pf_all_value_float.append(float(cluster_pf_value.replace(',','')))

    # normalized data to display in the graphic
    q_30_run_normalized , q_30_all_normalized = normalized_data (q_30_run_value_float, q_30_all_value_float)
    mean_run_normalized , mean_all_normalized = normalized_data (mean_run_value_float, mean_all_value_float)
    yield_mb_run_normalized , yield_mb_all_normalized = normalized_data (yield_mb_run_value_float, yield_mb_all_value_float)
    cluster_pf_run_normalized , cluster_pf_all_normalized = normalized_data (cluster_pf_run_value_float, cluster_pf_all_value_float)
    q30_run_str = ','.join(q_30_run_normalized)
    mean_run_str = ','.join(mean_run_normalized)
    yield_mb_run_str = ','.join(yield_mb_run_normalized)
    cluster_pf_run_str = ','.join(cluster_pf_run_normalized)

    q_30_all_str = ','.join(q_30_all_normalized)
    mean_all_str = ','.join(mean_all_normalized)
    yield_mb_all_str = ','.join(yield_mb_all_normalized)
    cluster_pf_all_str = ','.join(cluster_pf_all_normalized)

    # prepare the graphic
    run_year = run_object.get_run_year()
    heading =  run_object.get_run_name() +' versus runs executed on '
    sub_caption = str( 'year ' + str(run_year))
    theme = 'fint'
    x_axis_name = 'Quatilty measures (normalized data)'
    y_axis_name = 'Normalized values '
    series = [[run_object.runName,'#0075c2', '#1aaf5d'],['All runs','#f45b00','#f2c500']]
    data = [[q30_run_str,mean_run_str, yield_mb_run_str, cluster_pf_run_str],[q_30_all_str, mean_all_str, yield_mb_all_str, cluster_pf_all_str]]

    categories = ['Q > 30', 'Mean Quality Score', 'Yield MB', 'Cluster PF']
    data_source = bloxplot_graphic(heading, sub_caption, x_axis_name, y_axis_name, theme, categories, series, data)

    return FusionCharts("boxandwhisker2d", "box1" , "800", "400", "box_chart1", "json", data_source).render()


def graphics_state (state):
    '''
    Description:
        The function will match the state to return the color and the percentage
        value to display.
    Input:
        state      # contains the state value of the run
    Variables:
        g_color     # contains the color list
        g_value     # contains the percentage value
        state_list     # list with the possible states to display the run
                    information
    return:
        the index value of g_value and g_color [index]
    '''
    state_list = ['Error',  'Recorded',  'Sample Sent', 'Processing Run',
                'Processed Run', 'Processing Bcl2fastq', 'Processed Bcl2fastq',
                'Completed', 'Cancelled' ]
    g_value = [ 10, 15, 30, 45, 60, 75, 90, 100, 10]
    g_color = [ 'red', 'violet',  'pink', 'brown', 'orange', 'yellow',
                'yellow', 'green',  'red']
    index = state_list.index(state)

    return g_value[index], g_color [index]


def get_run_graphics (run_object) :
    '''
    Description:
        The function will get the run graphics.
    Input:
        run_object      # contains the run object
    Constants:
        RUN_IMAGES_DIRECTORY
        MEDIA_URL
    Variables:
        folder_graphic  # contains the folder where are the run graphics
        run_graphics_object  # contains GraphicsStats object
        run_graphics     # contain the tupla list with the graphic name
                            and the grafic path
    return:
        run_graphics
    '''
    # prepare the graphics for the run
    run_graphics = []
    if not GraphicsStats.objects.filter(runprocess_id__exact = run_object).exists():
        return run_graphics
    run_graphics_object = GraphicsStats.objects.get(runprocess_id__exact = run_object)
    folder_graphic = os.path.join( settings.MEDIA_URL, wetlab_config.RUN_IMAGES_DIRECTORY,
                            run_graphics_object.get_folder_graphic() )
    graphics = run_graphics_object.get_graphics().split(';')

    graphic_text= ['Data By Lane','Flow Cell Chart','Data By Cycle','QScore Heatmap','QScore Distribution','Indexing QC']


    for index_graph in range (len(graphics)):
        run_graphics.append([graphic_text[index_graph], os.path.join(folder_graphic, graphics[index_graph])])
    return  run_graphics


def get_run_read_data(run_object, num_of_reads, number_of_lanes) :
    '''
    Description:
        The function will get the run metrics per read information.
    Input:
        run_object      # contains the run object
        num_of_reads    # number of read used in the run
        number_of_lanes # number of lanes used in the run
    Variables:
        read_summary_values # contains the information per read

        index_run_summary   # contains the number of index used in the run
        rp_list     # list with the items to display in the running parameter
                    information
        rp_data     # get the list of the running parameter data
        run_parameters # running parameter object to fetch the number of reads

    return:
        info_reads_dict
    '''

    info_reads = {}
    info_reads_dict ={}

    for read_number in range (1, num_of_reads +1) :
        read_summary_values=[]
        for lane_number in range(1, number_of_lanes+1):
            read_lane_id= StatsRunRead.objects.filter(runprocess_id__exact =run_object, read__exact = read_number, lane__exact = lane_number)
            lane_values=read_lane_id[0].get_bin_run_read().split(';')
            read_summary_values.append(lane_values)

        read_number_index2 = str('Read '+str(read_number))
        info_reads_dict[read_number_index2] = read_summary_values

    #info_dict['reads']= info_reads_dict
    return info_reads_dict


def get_run_summary_data(run_object, num_of_reads) :
    '''
    Description:
        The function will get the summary run metrics information.
    Input:
        run_object      # contains the run object
        num_of_reads    # number of read used in the run
    Variables:
        index_run_summary   # contains the number of index used in the run
        rp_list     # list with the items to display in the running parameter
                    information
        rp_data     # get the list of the running parameter data
        run_parameters # running parameter object to fetch the number of reads

    return:
        run_summary_values and per_read_run_summary     # contain a tupla list with text and its value
    '''

    run_parameters = RunningParameters.objects.get(runName_id__exact = run_object)
    num_of_reads = run_parameters.get_number_of_reads ()
    index_run_summary = [i +1 for i in range(num_of_reads)]
    index_run_summary.append('Non Index')
    index_run_summary.append('Total')
    line_description = ['Read ' +str( i+1) for i in range (num_of_reads)]
    line_description.append('Non Index')
    line_description.append('Totals')

    per_read_run_summary = []

    for index in range (len(index_run_summary)):
        run_summary_id = StatsRunSummary.objects.filter(runprocess_id__exact =run_object , level__exact = index_run_summary[index])
        run_summary_values = run_summary_id[0].get_bin_run_summary().split(';')
        run_summary_values.insert(0, line_description[index])
        if index_run_summary[index] == 'Total':
            total_run_summary = run_summary_values
        else:
            per_read_run_summary.append(run_summary_values)

    return total_run_summary, per_read_run_summary


def get_running_parameters (run_object) :
    '''
    Description:
        The function will get the information from a specific run requested
        on the input parameter.
    Input:
        run_object      # contains the run object
    Variables:
        rp_list     # list with the items to display in the running parameter
                    information
        rp_data     # get the list of the running parameter data
        running_parameter_object # running parameter object to fetch the
                                information
    return:
        rp_info     # contain a tupla list with text and its value
    '''
    rp_info = []
    # finding the running parameters index for the run
    running_parameter_object = RunningParameters.objects.get(runName_id = run_object)

    # Adding the Run Parameters information
    rp_list=['Run ID','Experiment Name ','RTA version ','System Suite Version','Library ID ',
            'Chemistry','Run Start Date', 'Analysis Work Flow Type','Run Management Type',
            'Planned Read1 Cycles', 'Planned Read2 Cycles','Planned Index1 Read Cycles',
            'Planned Index2 Read Cycles','Application Version','Num Tiles per Swatch',
            'Image Channel', 'Flowcel','Image Dimensions', 'Flowcell Layout']
    #rp_data=running_parameter_object.get_run_parameters_info().split(';')
    rp_data = running_parameter_object.get_run_parameters_info()
    for i in range (len(rp_list)):
        rp_info.append([rp_list[i], rp_data[i]])

    return rp_info


def match_unkownbarcodes_with_index (unknow_dict) :
    '''
    Description:
        The function will match the unknow barcodes found in the run aginst
        the one already stored in the library kits. Unmatch or
    Input:
        unknow_dict      # dictionary with the unknow barcodes
    Variables:
        libraries_using_base   # contains the indexLibrary objects for
                            index that matches
        g_value     # contains the percentage value
        state_list     # list with the possible states to display the run
                    information
    return:
        index_match_list
    '''
    # prepare the data to match unknow barcodes against the index base sequence
    index_match_list = []
    for key , value in unknow_dict.items():
        found_unknow_index = []
        found_unknow_index.append(key)
        index_temp = ''
        library_info = []
        #
        if '+' in key:
            split_base = key.split('+')

            if CollectionIndexValues.objects.filter(i_7_seq__exact = split_base[0]).exists():
                libraries_using_base = CollectionIndexValues.objects.filter(i_7_seq__exact = split_base[0])
                index_temp = split_base[0]
                for library in libraries_using_base :
                    library_info.append([library.index_7,library.get_collection_index_name()])


            if CollectionIndexValues.objects.filter(i_5_seq__exact = split_base[1]).exists():
                if len(index_temp) == 1:
                    index_temp += (str (' + ' + split_base[1]))
                else:
                    index_temp = split_base[1]
                libraries_using_base = CollectionIndexValues.objects.filter(i_5_seq__exact = split_base[1])
                for library in libraries_using_base :
                    library_info.append([library.index_5,library.get_collection_index_name()])

        else:
            if CollectionIndexValues.objects.filter(i_7_seq__exact = key).exists():
                found_unknow_index.append(key)
                libraries_using_base = CollectionIndexValues.objects.filter(i_7_seq__exact = key)
                for library in libraries_using_base :
                    library_info.append([library.index_7,library.get_collection_index_name()])

        if len (index_temp) == 0 :
            index_temp= 'Index not match '
        if len (library_info) == 0 :
            library_info = ['Index bases not found in library']

        found_unknow_index.append(index_temp)
        found_unknow_index.append(library_info)
        index_match_list.append(found_unknow_index)

    return index_match_list



def get_information_for_incompleted_run():
    '''
    Description:
        The function will get the information from the runs that are incompleted.
        It creates a 4 groups. Recorded, Error, Canceled and rest of the states.
        Creates a pie graphic
    Functions:
        graphic_3D_pie      # located at iSkyLIMS_wetlab/utils/stats_graphics
    Return:
        run_information with information collected for each run
    '''
    run_information = {}
    today = datetime.date.today()
    if RunProcess.objects.filter(state__runStateName = 'Recorded').exists():
        run_information['recorded'] = []
        run_objs = RunProcess.objects.filter(state__runStateName = 'Recorded').order_by('runName')
        for run_obj in run_objs :
            data = []
            data.append(run_obj.get_run_id())
            data.append(run_obj.get_run_name())
            recorded_date = run_obj.get_recorded_date_no_format().date()
            data.append(recorded_date.strftime("%d %B  %Y"))
            data.append(str((today - recorded_date).days))
            run_information['recorded'].append(data)


    if RunProcess.objects.filter(state__runStateName = 'Error').exists():
        run_information['error'] = []
        run_objs = RunProcess.objects.filter(state__runStateName = 'Error').order_by('runName')
        for run_obj in run_objs :
            data = []
            data.append(run_obj.get_run_id())
            data.append(run_obj.get_run_name())
            run_date = run_obj.get_run_date_no_format()
            if run_date == None:
                # if no value stored on run date use the recorded date instead
                run_date = run_obj.get_recorded_date_no_format().date()
            data.append(run_date.strftime("%d %B  %Y"))
            data.append(run_obj.get_state_before_error())
            data.append(run_obj.get_error_text())
            data.append(str((today - run_date).days))
            run_information['error'].append(data)

    if RunProcess.objects.filter(state__runStateName = 'Cancelled').exists():
        run_information['cancelled'] = []
        run_objs = RunProcess.objects.filter(state__runStateName = 'Cancelled').order_by('runName')
        for run_obj in run_objs :
            data = []
            data.append(run_obj.get_run_id())
            data.append(run_obj.get_run_name())
            #import pdb; pdb.set_trace()
            run_date = run_obj.get_run_date_no_format()
            if run_date == None:
                # if no value stored on run date use the recorded date instead
                run_date = run_obj.get_recorded_date_no_format().date()
            data.append(run_date.strftime("%d %B %Y"))

            data.append(str((today - run_date).days))
            run_information['cancelled'].append(data)

    exclude_state = ['Recorded', 'Error', 'Cancelled', 'Completed','Pre-Recorded']

    if RunProcess.objects.all().exclude(state__runStateName__in = exclude_state).exists():
        run_information['other'] = []
        run_objs = RunProcess.objects.all().exclude(state__runStateName__in = exclude_state).order_by('state__runStateName')
        for run_obj in run_objs :
            data = []
            data.append(run_obj.get_run_id())
            data.append(run_obj.get_run_name())
            run_date = run_obj.get_run_date_no_format()
            if run_date == None:
                # if no value stored on run date use the recorded date instead
                run_date = run_obj.get_recorded_date_no_format().date()
            data.append(run_date.strftime("%d %B %Y"))
            data.append(run_obj.get_state())
            data.append(str((today - run_date).days))
            run_information['other'].append(data)

    if len(run_information) > 0:
        runs_in_state = {}
        for key in run_information.keys():
            runs_in_state[key] = len(run_information[key])
    #graphic_3D_pie (heading, sub_title, axis_x_description, axis_y_description, theme, source_data)
    data_source = graphic_3D_pie('Incomplete Runs', '', '', '', 'fint',runs_in_state)

    run_information['incompleted_graphic'] = FusionCharts("pie3d", "ex1" , "550", "400", "chart-1", "json", data_source).render()
    return run_information


def get_information_run(run_object):
    '''
    Description:
        The function will get the information from a specific run requested
        on the input parameter.
    Input:
        run_object      # contains the run object
        run_id              # contains the run id of the run_object
    Functions:
        graphics_state      # located at this file
        get_running_parameters # located at this file
        normalized_data     # imported from wetlab_misc_utilities
        check_run_in_same_year # imported from generic_functions
    Constants:
        RUN_IMAGES_DIRECTORY
        MEDIA_URL
    Variables:
        all_lane_summary # contains all LaneSummary objects from the year
                           of run excluding the default lane
        chem_high_mid   # has the chemistry of the run
        d_list      # Contains the list of parameters to be fetched from database.
                        2 different set of values, based on the run state
        data_source # reused variable to have the json data format for
                        displaying the graphics
        info_dict   # dictionary where collect all the run information that
                        will be returned
        p_list      # list of projects objects in the run
        p_info      # Tupla containing the project names and their index in Database

        rp_data     # run parameters information

        rp_list     # list of run parameters to be fetched from database.

        run_info  # Tupla list contating the string from d_list and
                        the value get from run_info_data
        run_data # list having the run information get from run object
                        function "get_info_process()"
        run_different_chemistry # contains the object runs that are different
                                that the selected run
        run_lane_summary # Lane summary object for the run
        run_state   # state of the run
    Return:
        info_dict with all information collected in the function
    '''
    info_dict={}
    run_info=[]
    p_info=[]
    p_library_kits = []
    ## collect the state to get the valid information of run that matches the run name
    run_state=run_object.get_state()
    info_dict['run_state'] = run_state
    # if run is processing data to insert in table show a message that
    # going back again after some minutes .
    no_valid_information = ['Processing Demultiplexing', 'Processing test', 'None', 'Pre-Recorded','Processing Metrics']
    if run_state in no_valid_information :
        info_dict['no_stable_data'] = [[run_object.get_run_name(), run_state]]
        return info_dict

    # allow to change the run name in case that run state was recorded or Sample Sent
    if run_state == 'Recorded' or run_state == 'Sample Sent':
        info_dict['change_run_name'] = [[run_object.get_run_name(), run_object.get_run_id()]]

    d_list=['Run name','State of the Run is','Run was requested by','Run was recorded on date', 'Run date', 'Run Completion Date']
    if run_state == 'Completed':
        d_list += ['Bcl2Fastq finish Date', 'Run Finish Date', 'Disk space used for Images(in MB)',
                'Disk space used for Fasta Files(in MB)', 'Disk space used for other Files(in MB)',]
    run_data=run_object.get_info_process().split(';')

    for i in range (len (d_list)):
        run_info.append([d_list[i],run_data[i]])

    info_dict['run_name'] = run_object.get_run_name()
    info_dict['data']=run_info

    info_dict['Sample_Sheet'] = [['Sample Sheet File', run_object.get_sample_file()]]
    info_dict['graphic_value'], info_dict['graphic_color'] = graphics_state(run_state)

    if run_state == 'Error' :
        # get the state before the error to present run information
        run_state = run_object.get_state_before_error()
        info_dict['error_run'] = [run_object.get_run_name(), run_state, run_object.get_error_text(), run_object.get_run_id()]

    p_list= Projects.objects.filter(runProcess=run_object)
    if p_list !='':
        #p_info = []
        #p_library_kit_list = []
        for p in range (len(p_list)):
            p_info.append([p_list[p].projectName,p_list[p].id])
            # get information about the library kits used for this run
            lib_kit = p_list[p].get_index_library_name()
            if not lib_kit in p_library_kits :
                p_library_kits.append(lib_kit)
        info_dict['projects']=p_info
        info_dict['library_kit'] = p_library_kits
        info_dict['run_id'] = run_object.get_run_id()
    ## get information up on state level
    if RunningParameters.objects.filter(runName_id = run_object).exists():
        # Adding the Run Parameters information
        info_dict['running_parameters'] = get_running_parameters(run_object)

    # get the run metric  statistics if they are already processed
    if StatsRunSummary.objects.filter(runprocess_id__exact =run_object).exists():
        run_parameters = RunningParameters.objects.get(runName_id__exact = run_object)
        num_of_reads = run_parameters.get_number_of_reads ()
        number_of_lanes= int(run_parameters.get_number_of_lanes())

        # prepare data for Run Binary summary stats
        info_dict ['runSummaryHeading'] = ['Level','Yield','Projected Yield','Aligned (%)','Error Rate (%)','Intensity Cycle 1','Quality >=30 (%)']
        info_dict ['runSummaryTotal'] , info_dict ['runSummary'] = get_run_summary_data(run_object, num_of_reads)

        # prepare the data for Reads Binary summary stats
        info_dict ['laneSummaryHeading']= ['Lane','Tiles','Density (K/mm2)','Cluster PF (%)','Phas/Prephas (%)',
                    'Reads (M)','Reads PF (M)','%>= Q30','Tield (G)','Cycles Err Rate',
                    'Aligned (%)','Error Rate (%)','Error Rate 35 cycle (%)',
                    'Error Rate 50 cycle (%)','Error Rate 75 cycle (%)',
                    'Error Rate 100 cycle (%)','Intensity Cycle 1']
        info_dict ['reads'] = get_run_read_data(run_object, num_of_reads, number_of_lanes)

        info_dict['runGraphic'] = get_run_graphics (run_object)

    ## get the stats information if run is completed
    if run_state == 'Completed':
        # check if there are more runs in the same year to make the box comparation
        if len(get_run_in_same_year_to_compare(run_object)) > 1 :
            # prepare the data for run comparations
            info_dict ['boxplot'] = get_boxplot_comparation_runs (run_object)
        else:
            info_dict['not_boxplot'] = "True"
        percent_projects = {}

        # get the demultiplexion information for projects included in the run

        for project_demultiplexion in p_list :
            if not StatsLaneSummary.objects.filter(runprocess_id__exact = run_object, project_id = project_demultiplexion.id ).exists():
                # continue with the stats for the next proyect in the run
                continue
            lanes_for_percent_graphic = StatsLaneSummary.objects.filter(runprocess_id__exact = run_object, project_id = project_demultiplexion.id )
            percent_lane = []
            for lane in lanes_for_percent_graphic :
                percent_lane.append(float(lane.percentLane))
            percent_projects[project_demultiplexion.projectName] =format(statistics.mean(percent_lane),'2f')
            #series.append(project_demultiplexion.projectName)

        # get the demultiplexion information for the default

        percent_default_lane = []

        default_lanes_for_percent_graphic = StatsLaneSummary.objects.filter(runprocess_id__exact = run_object, defaultAll__exact = 'default')
        for default_lane in default_lanes_for_percent_graphic :
            percent_default_lane.append(float(default_lane.percentLane))

        percent_projects['Unable to identify the project'] = format(statistics.mean(percent_default_lane),'2f')
        heading = 'Percentage of each project in the Run'
        sub_caption = ''
        theme = 'fint'
        #x_axis_name = 'Lanes'
        x_axis_name = 'Projects names'
        y_axis_name = 'Percentage '
        #categories = ['Lane 1', 'Lane 2', 'Lane 3','Lane 4']
        data_source = column_graphic_simple (heading, sub_caption, x_axis_name, y_axis_name, theme, percent_projects)
        #data_source = column_graphic_with_categories(heading, sub_caption, x_axis_name, y_axis_name, theme, categories, series, data)
        #
        info_dict ['run_project_comparation'] = FusionCharts("column3d", "column1" , "600", "400", "column_chart1", "json", data_source).render()

        fl_data_display=[]

        fl_summary_id = StatsFlSummary.objects.filter(runprocess_id__exact =run_object , project_id__isnull=True, defaultAll='all')
        fl_list = ['Cluster (Raw)', 'Cluster (PF)', 'Yield (MBases)', 'Number of Samples']
        fl_data_display.append(fl_list)
        fl_values = fl_summary_id[0].get_fl_summary()
        fl_data_display.append(fl_values)
        info_dict['fl_summary']=fl_data_display

        # prepare the data for Lane Summary
        lane_data_display = []
        lane_summary_id = StatsLaneSummary.objects.filter(runprocess_id__exact =run_object , project_id__isnull=True, defaultAll='all')
        lane_list = ['Lane', 'PF Clusters', '% of the lane','% Perfect barcode',
                    '% One mismatch barcode','Yield (Mbases)','% >= Q30 bases',
                    'Mean Quality Score']
        lane_data_display.append(lane_list)
        for lane_sum in lane_summary_id:
            lane_values = lane_sum.get_lane_summary()
            lane_data_display.append(lane_values)
        info_dict['lane_summary'] = lane_data_display

        # prepare the data for default Flowcell summary
        default_fl_data_display=[]

        default_fl_summary_id = StatsFlSummary.objects.filter(runprocess_id__exact =run_object , project_id__isnull=True, defaultAll='default')
        default_fl_data_display.append(fl_list)
        default_fl_values = default_fl_summary_id[0].get_fl_summary()
        default_fl_data_display.append(default_fl_values)
        info_dict['default_fl_summary']=default_fl_data_display

        # prepare the data for default Lane Summary
        default_lane_data_display = []
        default_lane_summary_id = StatsLaneSummary.objects.filter(runprocess_id__exact =run_object , project_id__isnull=True, defaultAll='default')
        default_lane_data_display.append(lane_list)
        for default_lane_sum in default_lane_summary_id:
            default_lane_values = default_lane_sum.get_lane_summary()
            default_lane_data_display.append(default_lane_values)
        info_dict['default_lane_summary'] = default_lane_data_display

        # prepare the data for top unknown barcode
        unknow_dict = {}
        for lane_un in range (number_of_lanes):
            lane_unknow_barcode = []
            lane_number=str(lane_un +1)

            unknow_bar_id = RawTopUnknowBarcodes.objects.filter(runprocess_id__exact =run_object , lane_number__exact = lane_number)
            #
            for item_id in unknow_bar_id:
                #
                unknow_values = item_id.get_unknow_barcodes().split(';')
                lane_unknow_barcode.append(unknow_values)
                unknow_bar_value = int(unknow_values[0].replace(',',''))
                if unknow_values[1] in unknow_dict:
                    unknow_dict [unknow_values[1]] += unknow_bar_value
                else:
                    unknow_dict [unknow_values[1]] = unknow_bar_value

            lane_number=str('unknow_bar_'+ str(lane_un))
            info_dict[lane_number] = lane_unknow_barcode
            #
            # keep the top 10 unknow bar by deleting the lowest values
            unknow_dict_len = len (unknow_dict)

        # create chart with the top unknown barcode in the run
        data_source = json_unknow_barcode_graphic('Unknow Sequence', unknow_dict)
        #data_source = json_unknow_barcode_graphic('Unknow Sequence', list(unknow_dict.keys()),list(unknow_dict.values()))
        unknow_pie3d = FusionCharts("pie3d", "ex1" , "600", "400", "chart-1", "json", data_source)

        info_dict ['unknow_pie3d'] = unknow_pie3d.render()

        info_dict['match_unknows']= match_unkownbarcodes_with_index(unknow_dict)


    return info_dict


def get_information_project (project_id, request):
    '''
    Description:
        The function will get the information from a specific project requested
        on the input parameter.
    Input:
        project_id      # contains the project id
        request         # contains gjango request to be used to identify
                        the user group the run id of the run_name_found
    Constants:
        WETLAB_MANAGER  #
        HEADING_FOR_PROJECT_DATES
        HEADING_SINGLE_PROJECT_FL_SUMMARY
        HEADING_SINGLE_PROJECT_STATS_LANE
        HEADING_SINGLE_PROJECT_SAMPLES
    Return:
        project_info_dict with all information collected in the function
    '''
    project_info_dict = {}
    run_names = []
    #project_info_dict['project_id'] = project_id.id
    project_info_dict['p_name'] = project_id.get_project_name()
    project_info_dict ['user_id'] = project_id.get_user_name()
    run_objs = project_id.runProcess.all()
    for run_obj in run_objs:
        run_names.append([run_obj.get_run_id(),run_obj.get_run_name()])
    #run_name = project_id.runprocess_id.runName
    groups = Group.objects.get(name = WETLAB_MANAGER)

    if groups in request.user.groups.all():
        project_info_dict['manager_display_runs'] = True

    project_info_dict['run_name'] = run_names
    if len(run_objs) == 1 :
        run_obj = run_objs[0]
        project_info_dict['collection_index'] = project_id.get_index_library_name()
        # project_info_dict['base_space_file'] = project_id.get_base_space_file()
        project_info_dict['dates'] = list(zip(HEADING_FOR_PROJECT_DATES, project_id.get_project_dates() ))

        # prepare the data for Flowcell Summary
        if StatsFlSummary.objects.filter(project_id__exact = project_id, runprocess_id = run_obj).exists():
            fl_summary_obj = StatsFlSummary.objects.filter(project_id__exact = project_id, runprocess_id = run_obj).last()
            project_info_dict['fl_summary_heading'] = wetlab_config.HEADING_SINGLE_PROJECT_FL_SUMMARY
            #fl_data_display.append(fl_list)
            project_info_dict['fl_summary_data'] = fl_summary_obj.get_fl_summary()

        # prepare the data for Lane Summary
        lane_data_display = []
        if StatsLaneSummary.objects.filter(project_id__exact = project_id , runprocess_id = run_obj).exists():
            lane_summary_obj = StatsLaneSummary.objects.filter(project_id__exact = project_id , runprocess_id = run_obj)
            project_info_dict['lane_summary_heading'] = wetlab_config.HEADING_SINGLE_PROJECT_STATS_LANE
            project_info_dict['lane_summary_data'] = []
            for lane_sum in lane_summary_obj:
                project_info_dict['lane_summary_data'].append(lane_sum.get_lane_summary())

        # prepare the data for sample information
        if SamplesInProject.objects.filter(project_id__exact = project_id, runProcess_id = run_obj).exists():
            sample_objs = SamplesInProject.objects.filter(project_id__exact = project_id, runProcess_id = run_obj)
            project_info_dict['sample_heading'] = wetlab_config.HEADING_SINGLE_PROJECT_SAMPLES
            project_info_dict['sample_data'] = []
            for sample_obj in sample_objs :
                project_info_dict['sample_data'].append([sample_obj.get_sample_id(),  sample_obj.get_sample_information()])

    return project_info_dict



def get_info_sample_in_run (sample_run_obj):
    '''
    Description:
        The function will get the information from a specific sample
        requested on the input parameter.
    Input:
        sample_run_obj           # contains the sample obj of the requested sample
    Variables:
        data_source # reused variable to have the json data format for
                        displaying the graphics
        sample_info_dict   # dictionary where collect all the run information that
                        will be returned
        quality_sample  # contains the quality value of the sample
        quality_sample_angular  # contains the graphic script and data for
                        # displaying the quality graphic
        percentage_chart         # contains the graphic script and data for
                        # showing the sample percentage in the project
    Return:
        sample_info_dict with all information collected in the function
    '''

    sample_info_dict ={}
    #collect the general information from the Sample
    sample_info_dict['sample_name'] = sample_run_obj.get_sample_name()
    sample_info_dict['project_name'] = sample_run_obj.get_project_name()
    project_id= sample_run_obj.project_id.id
    sample_info_dict['project_id'] = project_id
    #
    sample_info_dict['run_id'] = sample_run_obj.get_run_id()
    sample_info_dict['run_name'] = sample_run_obj.get_run_name()
    user_name_id = sample_run_obj.project_id.user_id
    sample_info_dict['investigator_name'] = user_name_id.username
    # collect the Sample information
    sample_info_dict['heading_samples_info'] = ['Sample', 'Index used', 'PF Cluster', '% of Project','Yield (Mbases)','>= Q30 bases','Mean Quality Score']
    sample_info_dict['data_samples_info'] = sample_run_obj.get_sample_information()
    # Quality graphic
    quality_sample = sample_run_obj.get_quality_sample()
    heading_chart_quality = 'Quality for the Sample ' + sample_run_obj.sampleName
    data_source = graphic_for_quality_angular(heading_chart_quality, quality_sample)
    quality_sample_angular = FusionCharts("angulargauge", "ex1" , "350", "200", "chart-1", "json", data_source)
    sample_info_dict['quality_chart1'] = quality_sample_angular.render()
    percentage_in_project ={}
    samples_in_project = SamplesInProject.objects.filter(project_id__exact = project_id)
    for sample in samples_in_project:
        percentage_in_project[sample.sampleName] = sample.percentInProject
    heading_samples_in_project = 'Samples belonging to the same project'
    sub_caption = ''
    x_axis_name = 'Samples in project ' + sample_run_obj.get_project_name()
    y_axis_name = '% of Project '
    theme = 'fint'
    data_source = column_graphic_one_column_highligthed (heading_samples_in_project, sub_caption, x_axis_name, y_axis_name, theme, percentage_in_project, sample_run_obj.sampleName)
    #
    percentage_chart = FusionCharts("column3d", 'samplesProject' , "750", "300", 'samples-chart-2', "json", data_source)
    sample_info_dict['percentage_chart'] = percentage_chart.render()

    return sample_info_dict

def get_last_runs_by_sequencer():
    '''
    Description:
        The function will get the last run executed for each sequencer
    Return:
        runs_using_sequencer
    '''
    last_runs = []
    for sequencer in get_sequencer_installed_names():
        if RunProcess.objects.filter(usedSequencer__sequencerName__exact = sequencer).exists():
            last_run_obj = RunProcess.objects.filter(usedSequencer__sequencerName__exact = sequencer).last()
            last_runs.append([sequencer, last_run_obj.get_run_platform(), last_run_obj.get_run_id(),last_run_obj.get_run_name() ])
    return last_runs

def get_sequencers_run_from_time_interval(sequencer,start_date, end_date):
    '''
    Description:
        The function will get the runs executed in the sequencer, during a period of time
        divided in completed runs and not completed runs
    Return:
        runs_using_sequencer
    '''
    runs_using_sequencer = {}
    if SequencerInLab.objects.filter(sequencerName__exact = sequencer).exists():
        sequencer_obj = SequencerInLab.objects.get(sequencerName__exact = sequencer)
        if RunProcess.objects.filter(usedSequencer = sequencer_obj, run_date__range=(start_date, end_date)).exists():
            runs_using_sequencer['completed_run_objs'] = []
            runs_using_sequencer['not_completed_run_objs'] = []
            run_objs_found = RunProcess.objects.filter(usedSequencer = sequencer_obj, run_date__range=(start_date, end_date))
            query_completed = run_objs_found.filter(state__runStateName__exact = 'Completed')
            for item in query_completed :
                runs_using_sequencer['completed_run_objs'].append(item)
            query_not_completed = run_objs_found.exclude(state__runStateName__exact = 'Completed')
            for item in query_not_completed :
                runs_using_sequencer['not_completed_run_objs'].append(item)

    return runs_using_sequencer



def get_sequencer_installed_names():
    '''
    Description:
        The function will get list names of the sequencers
    Return:
        sample_info_dict with all information collected in the function
    '''
    sequencer_list = []
    if SequencerInLab.objects.all().exists():
        sequencers = SequencerInLab.objects.all().order_by('platformID')
        for sequencer in sequencers:
            sequencer_list.append(sequencer.get_sequencer_name())
    return sequencer_list

def get_number_of_runs_per_sequencer():
    '''
    Description:
        The function will get a dictionary with the sequencer name and the number
        of runs executed in each
    Functions:
        get_sequencer_installed_names   # located at this file
    Return:
        number_of_runs
    '''
    number_of_runs = {}
    sequencer_list = get_sequencer_installed_names()
    for sequencer in sequencer_list:
        if RunProcess.objects.filter(usedSequencer__sequencerName__exact = sequencer).exists():
            number_of_runs[sequencer]= RunProcess.objects.filter(usedSequencer__sequencerName__exact = sequencer).count()
    return number_of_runs


def get_sequencer_names_from_platform(platform):
    '''
    Description:
        The function will get list of the sequencers using the platform value in the input data
    Input:
        platform           # contains the platform name
    Return:
        sample_info_dict with all information collected in the function
    '''
    sequencer_list = []
    if SequencerInLab.objects.filter(platformID__platformName__exact = platform).exists():
        sequencers = SequencerInLab.objects.filter(platformID__platformName__exact = platform)
        for sequencer in sequencers:
            sequencer_list.append(sequencer.get_sequencer_name())
    return sequencer_list

def  get_stats_sequencer_data_from_selected_runs (runs_using_sequencer, sequencer,start_date, end_date):
    '''
    Description:
        The function will get the sequencer statistic information from the runs selected
        in runs_using_sequencer
    Input:
        runs_using_sequencer  # contains the dictionary with completed and not completed runs
        sequencer       # sequencer name
        start_date
        end_date
    Constants:
        HEADING_FOR_STATISTICS_RUNS_BASIC_DATA

    Return:
        sequencer_data        # with all information collected in the function
    '''
    sequencer_data = {}
    sequencer_data['sequencer_name'] = sequencer
    sequencer_data['run_completed'] = []
    sequencer_data['not_run_completed'] = []
    sequencer_data['run_name_heading'] = HEADING_FOR_STATISTICS_RUNS_BASIC_DATA
    # get the run completed and not completed for the sequencer
    for run_in_seq in runs_using_sequencer['completed_run_objs']:
        sequencer_data['run_completed'].append([run_in_seq.get_run_name(),run_in_seq.get_run_date(),run_in_seq.get_run_id()])
    for run_in_seq in runs_using_sequencer['not_completed_run_objs'] :
        sequencer_data['not_run_completed'].append([run_in_seq.get_run_name(),run_in_seq.get_run_date(),run_in_seq.get_run_id()])

    # get the data for comparation sequencer vs all sequencers
    sequencer_runs = get_number_of_runs_per_sequencer()
    heading = 'Runs executed in each sequencer in (%)'
    data_source = pie_graphic_standard (heading, '', 'ocean', sequencer_runs)
    run_pie_graph = FusionCharts("pie3d", "run_pie_graph" , "500", "400", "run_pie_chart", "json", data_source).render()
    sequencer_data ['run_pie_graph'] = run_pie_graph

    # get data to display the time between run stats and finish
    # mean of all data is done
    time_diference =[]
    for run_in_seq in runs_using_sequencer['completed_run_objs']:

        time_diference.append((run_in_seq.get_run_finish_date_no_format().date()- run_in_seq.get_run_date_no_format()).days)
    sequencer_data['mean_time_duration'] = format(statistics.mean(time_diference), '.2f')

    # get the data for run executed in the sequencer per months
    run_time_dict ={}
    for run_in_seq in runs_using_sequencer['completed_run_objs']:
        run_year_month = run_in_seq.get_run_date_no_format().strftime("%Y_%m")
        if run_year_month in run_time_dict:
            run_time_dict[run_year_month] +=1
        else:
            run_time_dict[run_year_month] = 1
    run_time_tupla =[]
    for key in sorted(run_time_dict.keys()):
        run_time_tupla.append([key,run_time_dict[key]])
    heading = 'Runs executed per months in the sequencer'
    data_source = column_graphic_tupla (heading, '', '', 'Number of Runs','ocean', run_time_tupla, None)
    sequencer_data ['sequencer_runs_per_month_graph'] = FusionCharts("column3d", "run_per_month_graph" , "500", "400", "chart_seq_month", "json", data_source).render()

    # get the data for run executed in other sequencers per months
    if RunProcess.objects.filter(state__runStateName__exact = 'Completed').exclude(usedSequencer__sequencerName__exact = sequencer).exists():
        runs_in_other_sequencers = RunProcess.objects.filter(state__runStateName__exact = 'Completed').exclude(usedSequencer__sequencerName__exact = sequencer)
        run_time_dict ={}
        for run_in_other_seq in runs_in_other_sequencers:
            run_year_month = run_in_other_seq.get_run_date_no_format().strftime("%Y_%m")
            if run_year_month in run_time_dict:
                run_time_dict[run_year_month] +=1
            else:
                run_time_dict[run_year_month] = 1
        run_time_tupla =[]
        for key in sorted(run_time_dict.keys()):
            run_time_tupla.append([key,run_time_dict[key]])
        heading = 'Runs executed per months in the rest of the sequencers'
        data_source = column_graphic_tupla (heading, '', '', 'Number of Runs','fint', run_time_tupla, None)
        sequencer_data ['other_sequencers_runs_per_month_graph'] = FusionCharts("column3d", "other_run_per_month_graph" , "500", "400", "chart_other_seq_month", "json", data_source).render()

    return sequencer_data
