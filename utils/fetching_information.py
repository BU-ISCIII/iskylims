
import statistics
import os
from ..fusioncharts.fusioncharts import FusionCharts

from django.conf import settings
from django.contrib.auth.models import Group

from iSkyLIMS_wetlab.models import *
from iSkyLIMS_wetlab.wetlab_config import RUN_IMAGES_DIRECTORY, WETLAB_MANAGER
from .wetlab_misc_utilities import normalized_data
from .stats_graphics import *

def get_information_run(run_name_found,run_id):
    '''
    Description:
        The function will get the information from a specific run requested
        on the input parameter. 
    Input:
        run_name_found      # contains the run object 
        run_id              # contains the run id of the run_name_found
    Functions:
        get_machine_lanes   # imported from parsing_run_info
        normalized_data     # imported from wetlab_misc_utilities
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
        rp_info     # Tupla containing the list of run parameters and their value
        rp_list     # list of run parameters to be fetched from database.

        run_info  # Tupla list contating the string from d_list and
                        the value get from run_info_data
        run_data # list having the run information get from run object
                        function "get_info_process()"
        run_different_chemistry # contains the object runs that are different
                                that the selected run
        run_lane_summary # Lane summary object for the run
        run_state   # state of the run
        run_year    # contains the year of the selected run
        same_run_in_year    # contains all the runs in the year of the 
                            selected run excluding those which do not
                            have the same chemistry value
        same_runs_in_year_list  # contains the id list of the same_runs_in_year
        start_date  # it is the first of january of the run_year
        end_date    # it is the 31st of december of the run_year
        
        q_30_value      # string value containing the q_30 value from lane
        q_30_run_value_float  # q30 _value in float format
        q_30_run_normalized # q_30 normalized value for run
        q_30_all_normalized # q_30 normalized value for all runs
        q30_run_str     # q_30 normalized value in string format
        q_30_all_str    # q_30_all_normalized in string format
        ################################################################
        #### same structure for mean_value, yield_mb and cluster_pf
        ################################################################
        q_30_value_list # list containing the q30 values in float format
        q_30_run_value
        q_30_run_value_float
        q_30_all_value_float
        
    Return:
        info_dict with all information collected in the function
    '''
    info_dict={}
    rp_info=[]
    p_info=[]
    p_library_kits = []
    ## collect the state to get the valid information of run that matches the run name
    run_state=run_name_found.get_state()

    # allow to change the run name in case that run state was recorded or Sample Sent
    if run_state == 'Recorded' or run_state == 'Sample Sent':
        info_dict['change_run_name'] = [[run_name_found.runName, run_id]]
    if (run_state != 'Completed'):
        d_list=['Run name','State of the Run is','Run was requested by','Run was recorded on date', 'Run date', 'Run Finish Date','RunID']
    else:
        number_of_lanes=run_name_found.get_machine_lanes()
        d_list=['Run name','State of the Run is','Run was requested by',
                'Disk space used for Images(in MB)','Disk space used for Fasta Files(in MB)',
                'Disk space used for other Files(in MB)','Run recorded date','Run date', 'Run Finish Date',
                'Bcl2Fastq finish Date','Run Completion Date']
    run_data=run_name_found.get_info_process().split(';')
    info_dict['Sample_Sheet'] = [['Sample Sheet File', run_name_found.get_sample_file()]]
    run_info=[]
    for i in range (len (d_list)):
        run_info.append([d_list[i],run_data[i]])
    info_dict['data']=run_info
    info_dict['run_state'] = run_state
    if (run_state.startswith('ERROR')):
        info_dict['graphic_value']=10
        info_dict['graphic_color']= 'red'
    if (run_state == 'Recorded'):
        info_dict['graphic_value']=25
        info_dict['graphic_color']= 'violet'
    if (run_state == 'Sample Sent'):
        info_dict['graphic_value']=40
        info_dict['graphic_color']= 'pink'
    if (run_state == 'Process Running'):
        info_dict['graphic_value']=50
        info_dict['graphic_color']= 'brown'
    if (run_state == 'Bcl2Fastq Executed'):
        info_dict['graphic_value']=60
        info_dict['graphic_color']= 'orange'
    if (run_state == 'Running Stats'):
        info_dict['graphic_value']=75
        info_dict['graphic_color']= 'yellow'
    if (run_state == 'Completed'):
        info_dict['graphic_value']=100
        info_dict['graphic_color']= 'green'

    p_list= Projects.objects.filter(runprocess_id=run_id)
    if p_list !='':
        #p_info = []
        #p_library_kit_list = []
        for p in range (len(p_list)):
            p_info.append([p_list[p].projectName,p_list[p].id])
            # get information about the library kits used for this run
            lib_kit = p_list[p].libraryKit
            if not lib_kit in p_library_kits :
                p_library_kits.append(lib_kit)
        info_dict['projects']=p_info
        info_dict['library_kit'] = p_library_kits
        info_dict['run_id'] = run_id

    ## get the stats information if run is completed
    if run_state == 'Completed':
        # finding the running parameters index for the run
        run_parameter_object = RunningParameters.objects.get(pk=run_id)

        # Adding the Run Parameters information
        rp_list=['Run ID','Experiment Name ','RTA version ','System Suite Version','Library ID ','Chemistry','Run Start Date', 'Analysis Work Flow Type','Run Management Type','Planned Read1 Cycles',
                'Planned Read2 Cycles','Planned Index1 Read Cycles','Planned Index2 Read Cycles','Application Version','Num Tiles per Swatch','Image Channel',
                'Flowcel','Image Dimensions', 'Flowcell Layout']
        rp_data=run_parameter_object.get_run_parameters_info().split(';')
        
        for i in range (len(rp_list)):
            if i == 'Image Channel':
                img_data_list=rp_data[i].split(',')
                rp_info.append([rp_list[i],[img_data_list]])
            else:
                rp_info.append([rp_list[i], rp_data[i]])
        info_dict['parameters']=rp_info
        
        # prepare the data for q-means
        # fetch Q>30 , mean_q and yield mb for all projects per lane to create the boxplot
        run_lane_summary = NextSeqStatsLaneSummary.objects.filter(runprocess_id__exact =run_id ).exclude(defaultAll__isnull = False)
        q_30_value_list, mean_value_list = [] , []
        q_30_run_value , mean_run_value = [] , []
        q_30_run_value_float , mean_run_value_float , yield_mb_run_value_float, cluster_pf_run_value_float = [] , [], [] , []
        q_30_all_value_float , mean_all_value_float , yield_mb_all_value_float , cluster_pf_all_value_float = [] , [] , [], []

        for item in run_lane_summary:
            q_30_value, mean_value , yield_mb_value , cluster_pf_value, = item.get_stats_info().split(';')

            q_30_run_value_float.append(float(q_30_value))
            mean_run_value_float.append(float(mean_value))
            yield_mb_run_value_float.append(float(yield_mb_value.replace(',','')))
            cluster_pf_run_value_float.append(float(cluster_pf_value.replace(',','')))

        # get the chemistry type for the run, that will be used to compare runs with the same chemistry value
        chem_high_mid = RunningParameters.objects.get(runName_id__exact = run_id).Chemistry
        run_different_chemistry = RunningParameters.objects.all(). exclude(Chemistry__exact = chem_high_mid)
        run_year = run_name_found.run_date.timetuple().tm_year

        start_date = str(run_year) + '-1-1'
        end_date = str(run_year) +'-12-31'
        same_run_in_year = RunProcess.objects.filter(run_date__range=(start_date, end_date)).exclude(runName__in = run_different_chemistry)

        same_runs_in_year_list = []
        for run in same_run_in_year :
            same_runs_in_year_list.append(run.id)


        all_lane_summary = NextSeqStatsLaneSummary.objects.filter(runprocess_id__in = same_runs_in_year_list).exclude(defaultAll__isnull = False).exclude(runprocess_id__exact =run_id)
        for item in all_lane_summary:
            q_30_value, mean_value , yield_mb_value , cluster_pf_value = item.get_stats_info().split(';')

            q_30_all_value_float.append(float(q_30_value))
            mean_all_value_float.append(float(mean_value))
            yield_mb_all_value_float.append(float(yield_mb_value.replace(',','')))
            cluster_pf_all_value_float.append(float(cluster_pf_value.replace(',','')))

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



        heading =  run_name_found.runName +' versus runs executed on '
        sub_caption = str( 'year ' + str(run_year))
        theme = 'fint'
        x_axis_name = 'Quatilty measures (normalized data)'
        y_axis_name = 'Normalized values '
        series = [[run_name_found.runName,'#0075c2', '#1aaf5d'],['All runs','#f45b00','#f2c500']]
        data = [[q30_run_str,mean_run_str, yield_mb_run_str, cluster_pf_run_str],[q_30_all_str, mean_all_str, yield_mb_all_str, cluster_pf_all_str]]


        categories = ['Q > 30', 'Mean Quality Score', 'Yield MB', 'Cluster PF']
        data_source = bloxplot_graphic(heading, sub_caption, x_axis_name, y_axis_name, theme, categories, series, data)
        info_dict ['boxplot'] = FusionCharts("boxandwhisker2d", "box1" , "800", "400", "box_chart1", "json", data_source).render()

        percent_projects = {}
        # get the demultiplexion information for projects included in the run
        percent_lane = []
        for project_demultiplexion in p_list :
            lanes_for_percent_graphic = NextSeqStatsLaneSummary.objects.filter(runprocess_id__exact = run_id, project_id = project_demultiplexion.id )
            for lane in lanes_for_percent_graphic :
                percent_lane.append(float(lane.percentLane))
            percent_projects[project_demultiplexion.projectName] =format(statistics.mean(percent_lane),'2f')
            #series.append(project_demultiplexion.projectName)

        # get the demultiplexion information for the default

        percent_default_lane = []

        default_lanes_for_percent_graphic = NextSeqStatsLaneSummary.objects.filter(runprocess_id__exact = run_id, defaultAll__exact = 'default')
        for default_lane in default_lanes_for_percent_graphic :
            percent_default_lane.append(float(default_lane.percentLane))

        #series.append('Unable to identify the project')
        #data.append(statistics.mean(percent_default_lane))
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

        fl_summary_id = NextSeqStatsFlSummary.objects.filter(runprocess_id__exact =run_id , project_id__isnull=True, defaultAll='all')
        fl_list = ['Cluster (Raw)', 'Cluster (PF)', 'Yield (MBases)', 'Number of Samples']
        fl_data_display.append(fl_list)
        fl_values = fl_summary_id[0].get_fl_summary().split(';')
        fl_data_display.append(fl_values)
        info_dict['fl_summary']=fl_data_display

        # prepare the data for Lane Summary
        lane_data_display = []
        lane_summary_id = NextSeqStatsLaneSummary.objects.filter(runprocess_id__exact =run_id , project_id__isnull=True, defaultAll='all')
        lane_list = ['Lane', 'PF Clusters', '% of the lane','% Perfect barcode',
                    '% One mismatch barcode','Yield (Mbases)','% >= Q30 bases',
                    'Mean Quality Score']
        lane_data_display.append(lane_list)
        for lane_sum in lane_summary_id:
            lane_values = lane_sum.get_lane_summary().split(';')
            lane_data_display.append(lane_values)
        info_dict['lane_summary'] = lane_data_display

        # prepare the data for default Flowcell summary
        default_fl_data_display=[]

        default_fl_summary_id = NextSeqStatsFlSummary.objects.filter(runprocess_id__exact =run_id , project_id__isnull=True, defaultAll='default')
        default_fl_data_display.append(fl_list)
        default_fl_values = default_fl_summary_id[0].get_fl_summary().split(';')
        default_fl_data_display.append(default_fl_values)
        info_dict['default_fl_summary']=default_fl_data_display

        # prepare the data for default Lane Summary
        default_lane_data_display = []
        default_lane_summary_id = NextSeqStatsLaneSummary.objects.filter(runprocess_id__exact =run_id , project_id__isnull=True, defaultAll='default')
        default_lane_data_display.append(lane_list)
        for default_lane_sum in default_lane_summary_id:
            default_lane_values = default_lane_sum.get_lane_summary().split(';')
            default_lane_data_display.append(default_lane_values)
        info_dict['default_lane_summary'] = default_lane_data_display

        # prepare the data for top unknown barcode
        unknow_dict = {}
        for lane_un in range (number_of_lanes):
            lane_unknow_barcode = []
            lane_number=str(lane_un +1)

            unknow_bar_id = RawTopUnknowBarcodes.objects.filter(runprocess_id__exact =run_id , lane_number__exact = lane_number)
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

                if IndexLibraryValues.objects.filter(indexBase__exact = split_base[0]).exists():
                    libraries_using_base = IndexLibraryValues.objects.filter(indexBase__exact = split_base[0])
                    index_temp = split_base[0]
                    for library in libraries_using_base :
                        library_info.append([library.indexName,library.indexLibraryKit_id.indexLibraryName])


                if IndexLibraryValues.objects.filter(indexBase__exact = split_base[1]).exists():
                    if len(index_temp) == 1:
                        index_temp += (str (' + ' + split_base[1]))
                    else:
                        index_temp = split_base[1]
                    libraries_using_base = IndexLibraryValues.objects.filter(indexBase__exact = split_base[1])
                    for library in libraries_using_base :
                        library_info.append([library.indexName,library.indexLibraryKit_id.indexLibraryName])

            else:
                if IndexLibraryValues.objects.filter(indexBase__exact = key).exists():
                    found_unknow_index.append(key)
                    libraries_using_base = IndexLibraryValues.objects.filter(indexBase__exact = key)
                    for library in libraries_using_base :
                        library_info.append([library.indexName,library.indexLibraryKit_id.indexLibraryName])

            if len (index_temp) == 0 :
                index_temp= 'Index not match '
            if len (library_info) == 0 :
                library_info = ['Index bases not found in library']
                #libraries_using_base = IndexLibraryValues.objects.filter(indexBase__exact = split_base[1])
                #for library in libraries_using_base :
                #    found_unknow_index.append(library.indexBase)
            #index_item = 5
            found_unknow_index.append(index_temp)
            found_unknow_index.append(library_info)
            index_match_list.append(found_unknow_index)
        #

        info_dict['match_unknows']= index_match_list

        # prepare data for Run Binary summary stats

        run_parameters = RunningParameters.objects.get(runName_id__exact = run_id)
        num_of_reads = run_parameters.get_number_of_reads ()
        index_run_summary = [i +1 for i in range(num_of_reads)]
        index_run_summary.append('Non Index')
        index_run_summary.append('Total')
        #index_run_summary = ['1','2','3','4', 'Non Index', 'Total']
        info_dict ['runSummaryHeading']= ['Level','Yield','Projected Yield','Aligned (%)','Error Rate (%)','Intensity Cycle 1','Quality >=30 (%)']
        line_description = ['Read ' +str( i+1) for i in range (num_of_reads)]
        line_description.append('Non Index')
        line_description.append('Totals')

        #line_description=['Read 1','Read 2','Read 3','Read 4','Non Index','Totals']
        line_run_summary = []
        for index in range (len(index_run_summary)):
            #
            run_summary_id = NextSeqStatsBinRunSummary.objects.filter(runprocess_id__exact =run_id , level__exact = index_run_summary[index])
            run_summary_values = run_summary_id[0].get_bin_run_summary().split(';')
            run_summary_values.insert(0, line_description[index])
            if index_run_summary[index] == 'Total':
                info_dict ['runSummaryTotal'] = run_summary_values
            else:
                line_run_summary.append(run_summary_values)
        #
        info_dict ['runSummary'] = line_run_summary

        # prepare the data for Reads Binary summary stats
        info_dict ['laneSummaryHeading']= ['Lane','Tiles','Density (K/mm2)','Cluster PF (%)','Phas/Prephas (%)',
                    'Reads (M)','Reads PF (M)','%>= Q30','Tield (G)','Cycles Err Rate',
                    'Aligned (%)','Error Rate (%)','Error Rate 35 cycle (%)',
                    'Error Rate 50 cycle (%)','Error Rate 75 cycle (%)',
                    'Error Rate 100 cycle (%)','Intensity Cycle 1']
        info_reads_dict ={}
        for read_number in range (1, num_of_reads +1) :
            read_summary_values=[]
            for lane_number in range(1, number_of_lanes+1):
                read_lane_id= NextSeqStatsBinRunRead.objects.filter(runprocess_id__exact =run_id, read__exact = read_number, lane__exact = lane_number)
                lane_values=read_lane_id[0].get_bin_run_read().split(';')
                read_summary_values.append(lane_values)
            #read_number_index = str('laneSummary'+str(read_number))
            read_number_index2 = str('Read '+str(read_number))
            info_reads_dict[read_number_index2] = read_summary_values
            #info_dict[read_number_index] = read_summary_values
        info_dict['reads']= info_reads_dict

        # prepare the graphics for the run
        folder_for_plot='/documents/wetlab/images_plot/'

        run_graphics_id = NextSeqGraphicsStats.objects.filter(runprocess_id__exact =run_id)
        folder_graphic = os.path.join( settings.MEDIA_URL, RUN_IMAGES_DIRECTORY, run_graphics_id[0].get_folder_graphic() )
        graphics = run_graphics_id[0].get_graphics().split(';')
        graphic_text= ['Data By Lane','Flow Cell Chart','Data By Cycle','QScore Heatmap','QScore Distribution','Indexing QC']
        for index_graph in range (len(graphics)):
            tmp_value = graphics[index_graph]
            graphics[index_graph] = [graphic_text[index_graph], os.path.join(folder_graphic, tmp_value)]

        info_dict['runGraphic'] = graphics
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
    Functions:
        get_machine_lanes   # imported from parsing_run_info
        normalized_data     # imported from wetlab_misc_utilities
    Constants:
        WETLAB_MANAGER  # 
    Variables:
        groups          # get the group objects to check if requested user
                        belongs to wetlab manager group
        
        fl_data_display # contains the list of the flowcell values
        fl_values       # Tupla containing the flowcell summary heading 
                        and their values
        lane_values    # contains the list of the lanes values
        lane_data_display   # Tupla containing the lanes summary heading 
                        and their values
        p_data          # Tupla containing the project names and their
                        index in Database
        p_state         # project state
        project_info_dict   # dictionary where collect all the project
                        information that will be returned
        project_values  # contains the information retuned by get_project_info
        run_name        # contain the run name for the requested project
    Return:
        project_info_dict with all information collected in the function
    '''
    project_info_dict = {}
    p_data = []
    project_info_dict['project_id'] = project_id.id
    project_info_text = ['Project Name','Library Kit','File to upload to BaseSpace','Project Recorder date', 'Project date','Run name']
    project_values = project_id.get_project_info().split(';')
    run_name = project_id.runprocess_id.runName
    groups = Group.objects.get(name = WETLAB_MANAGER)

    if groups not in request.user.groups.all():
        project_info_dict['run_id'] = ''
    else:
        project_info_dict['run_id'] = project_id.runprocess_id.id
    project_values.append(run_name )
    for item in range(len(project_info_text)):
        p_data.append([project_info_text[item], project_values[item]])
    project_info_dict['p_data'] = p_data

    project_info_dict ['user_name'] = project_id.get_user_name()
    p_state = project_id.get_state()
    project_info_dict['state'] = p_state
    if p_state.startswith('ERROR'):
        project_info_dict['graphic_value']=10
        project_info_dict['graphic_color']='red'
    if p_state == 'Recorded':
        project_info_dict['graphic_value']=25
        project_info_dict['graphic_color']='violet'
    if p_state == 'Sample Sent':
        project_info_dict['graphic_value']= 50
        project_info_dict['graphic_color']='brown'
    if p_state == 'B2FqExecuted':
        project_info_dict['graphic_value']= 75
        project_info_dict['graphic_color']='yellow'
    if p_state == 'Completed':
        project_info_dict['graphic_value']= 100
        project_info_dict['graphic_color']='green'
        fl_data_display=[]

        # prepare the data for Flowcell Summary
        fl_summary_id = NextSeqStatsFlSummary.objects.get(project_id__exact = project_id)
        fl_list = ['Cluster (Raw)', 'Cluster (PF)', 'Yield (MBases)', 'Number of Samples']
        fl_data_display.append(fl_list)
        fl_values = fl_summary_id.get_fl_summary().split(';')
        fl_data_display.append(fl_values)
        project_info_dict['fl_summary']=fl_data_display

        # prepare the data for Lane Summary
        lane_data_display = []
        lane_summary_id = NextSeqStatsLaneSummary.objects.filter(project_id__exact = project_id)
        lane_list = ['Lane', 'PF Clusters', '% of the lane','% Perfect barcode',
                    '% One mismatch barcode','Yield (Mbases)','% >= Q30 bases',
                    'Mean Quality Score']
        lane_data_display.append(lane_list)
        for lane_sum in lane_summary_id:
            lane_values = lane_sum.get_lane_summary().split(';')
            lane_data_display.append(lane_values)
        project_info_dict['lane_summary'] = lane_data_display

        # prepare the data for sample information
        sample_found_list = SamplesInProject.objects.filter(project_id__exact = project_id)
        sample_heading_list = ['Sample','Barcode','PF Clusters','Percent of Project', 'Yield (Mbases)','% >= Q30 bases', 'Mean Quality Score']
        project_info_dict['sample_heading'] = sample_heading_list
        sample_list ={}
        for sample_item in sample_found_list :
            sample_line = sample_item.get_sample_information().split(';')
            sample_list[sample_item.id] = [sample_line]

        project_info_dict['sample_table'] = sample_list
    return project_info_dict



def get_info_sample (sample_id):
    '''
    Description:
        The function will get the information from a specific sample 
        requested on the input parameter. 
    Input:
        sample_id           # contains the sample id of the requested sample
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
    sample_info_dict['sample_name'] = sample_id.sampleName
    sample_info_dict['project_name'] = sample_id.get_project_name()
    project_id= sample_id.project_id.id
    sample_info_dict['project_id'] = project_id
    #
    sample_info_dict['run_id'] = sample_id.project_id.runprocess_id.id
    sample_info_dict['run_name'] = sample_id.project_id.runprocess_id.runName
    user_name_id = sample_id.project_id.user_id
    sample_info_dict['investigator_name'] = user_name_id.username
    # collect the Sample information
    sample_info_dict['heading_samples_info'] = ['Sample', 'Barcode', 'PF Cluster', '% of Project','Yield (Mbases)','>= Q30 bases','Mean Quality Score']
    sample_info_dict['data_samples_info'] = sample_id.get_sample_information().split(';')
    # Quality graphic
    quality_sample = sample_id.get_quality_sample()
    heading_chart_quality = 'Quality for the Sample ' + sample_id.sampleName
    data_source = graphic_for_quality_angular(heading_chart_quality, quality_sample)
    quality_sample_angular = FusionCharts("angulargauge", "ex1" , "350", "200", "chart-1", "json", data_source)
    sample_info_dict['quality_chart1'] = quality_sample_angular.render()
    percentage_in_project ={}
    samples_in_project = SamplesInProject.objects.filter(project_id__exact = project_id)
    for sample in samples_in_project:
        percentage_in_project[sample.sampleName] = sample.percentInProject
    heading_samples_in_project = 'Samples belonging to the same project'
    sub_caption = ''
    x_axis_name = 'Samples in project ' + sample_id.get_project_name()
    y_axis_name = '% of Project '
    theme = 'fint'
    data_source = column_graphic_one_column_highligthed (heading_samples_in_project, sub_caption, x_axis_name, y_axis_name, theme, percentage_in_project, sample_id.sampleName)
    #
    percentage_chart = FusionCharts("column3d", 'samplesProject' , "750", "300", 'samples-chart-2', "json", data_source)
    sample_info_dict['percentage_chart'] = percentage_chart.render()
    return sample_info_dict
