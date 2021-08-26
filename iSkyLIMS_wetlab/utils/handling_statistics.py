import statistics
from iSkyLIMS_wetlab.fusioncharts.fusioncharts import FusionCharts
from iSkyLIMS_wetlab.models import *
from iSkyLIMS_wetlab.wetlab_config import *
from iSkyLIMS_wetlab.utils.generic_functions import check_valid_date_format, get_configuration_value
from iSkyLIMS_wetlab.utils.fetching_information import get_sequencer_installed_names
from iSkyLIMS_wetlab.utils.stats_graphics import column_graphic_simple, column_graphic_tupla

def get_min_mean_and_max_values(values_data, reference_data, number_to_split):
    '''
    Description:
        The function get the lower and the higher values that are inside the "values "
        variable. but only the number of items defined in the "number_to_split"
        are considered. The middle values are join into and a mean value is returned
        In case that number of items is smaller than 2 times the number_to_split
        they are sorted and no mean value is calculated .
    Input:
        values_data          # list of values to extract the information
        reference_data  # list of data to get the releated value (run_name)
        number_to_split  # number of items to get
    Return:
        reference_query_values a tupla list with te reference name and value
    '''
    reference_query_values = []
    value_sort = values_data.copy()
    value_sort.sort()
    index_used_in_reference = []
    if len(values_data) <= 2*number_to_split:
        list_of_values = value_sort
    else:
        list_of_values = value_sort[0:number_to_split] + value_sort[-number_to_split:]
    for val in list_of_values:
        tmp_run_list_index = [index for index, value in enumerate(values_data) if value == val]

        for tmp_index in tmp_run_list_index:
            if tmp_index in index_used_in_reference:
                continue
            else:
                break
        reference_query_values.append([reference_data[tmp_index], val])
        index_used_in_reference.append(tmp_index)
    if len(values_data)  > 2*number_to_split:
        # Add the median value
        reference_query_values.insert(number_to_split + 1,['Median values',round(statistics.median(value_sort[number_to_split:-number_to_split]),2)])
    return reference_query_values



def get_researcher_statistics(researcher_name, start_date, end_date):
    '''
    Description:
        The function get the statitics information for the researcher
    Input:
        researcher_name     # name of the researcher
        start_date          # start date to collect statistics
        end_date            # end date for collecting statistics
    Return:
        ERROR if no match user conditions, or researcher_statistics
    '''
    researcher_statistics = {}
    if not User.objects.filter(username__icontains = researcher_name).exists():
        researcher_statistics['ERROR'] = wetlab_config.ERROR_NO_MATCHES_FOR_INPUT_CONDITIONS
        return researcher_statistics

    user_objs = User.objects.filter(username__icontains = researcher_name)

    if len(user_objs) > 1:
        researcher_statistics['ERROR'] = wetlab_config.ERROR_MANY_USER_MATCHES_FOR_INPUT_CONDITIONS
        return researcher_statistics
    # over write the user name with the full user id
    researcher_statistics['researcher_name'] = user_objs[0].username

    if start_date != '' and not check_valid_date_format(start_date):
        errror_message = wetlab_config.ERROR_INVALID_FORMAT_FOR_DATES
        return researcher_statistics
    if end_date != '' and not check_valid_date_format(start_date):
        errror_message = wetlab_config.ERROR_INVALID_FORMAT_FOR_DATES
        return researcher_statistics

    if SamplesInProject.objects.filter(user_id = user_objs[0]).exists():
        sample_objs = SamplesInProject.objects.filter(user_id = user_objs[0]).order_by('runProcess_id')
        # check if start and end date are present in the form
        if start_date != '' and end_date !='':
            sample_objs= sample_objs.filter(runProcess_id__state__runStateName__range=(start_date, end_date))
        elif start_date != '':
            sample_objs = sample_objs.filter(runProcess_id__state__runStateName__gte = start_date)
        elif end_date != '':
            sample_objs = sample_objs.filter(runProcess_id__state__runStateName__lte = end_date)
        if len(sample_objs) == 0:
            researcher_statistics['ERROR'] = wetlab_config.ERROR_NO_MATCHES_FOR_INPUT_CONDITIONS
            return researcher_statistics


        # Get data from researcher projects
        researcher_statistics ['sample_researcher_heading'] = wetlab_config.RESEARCHER_SAMPLE_HEADING_STATISTICS

        installed_sequencers = get_sequencer_installed_names()
        researcher_sample_data = {}
        for installed_sequencer in installed_sequencers:
            researcher_sample_data[installed_sequencer] = []
        for sample_obj in sample_objs:
            sample_data = sample_obj.get_sample_information_with_project_run()
            researcher_sample_data[sample_obj.get_sequencer_used()].append(sample_data)
        researcher_statistics['researcher_sample_data'] = researcher_sample_data
        # return researcher_statistics

        # collect number of saples for run and for projects
        nun_sample_in_run = {}
        num_sample_in_project = {}
        for sequencer, sample_in_sequencer in researcher_statistics['researcher_sample_data'].items():
            if len(sample_in_sequencer) == 0:
                continue
            for sample_values in sample_in_sequencer:
                run_name = sample_values[2]
                project_name = sample_values[1]
                if run_name not in nun_sample_in_run:
                    nun_sample_in_run[run_name] = 0
                nun_sample_in_run[run_name] += 1
                if project_name not in num_sample_in_project :
                    num_sample_in_project[project_name] = 0
                num_sample_in_project[project_name] += 1

        # Create run grapic
        heading = 'Graphics for number of samples in runs'
        data_source = column_graphic_simple (heading, '', 'Runs',  'Number of Samples', 'ocean', nun_sample_in_run)
        researcher_statistics['run_graphic'] = FusionCharts("column3d", 'run_graph' , "500", "350",'run_chart' , "json", data_source).render()
        # Create projert grapic
        heading = 'Graphics for number of samples in projects'
        data_source = column_graphic_simple (heading, '', 'Projects',  'Number of Samples', 'ocean', num_sample_in_project)
        researcher_statistics['project_graphic'] = FusionCharts("column3d", 'project_graph' , "500", "350",'project_chart' , "json", data_source).render()

        # collect Q> 30  and mean data for each sequencer used
        runs_index_sample = []
        q30_sample_value = []
        mean_sample_value =[]
        for sequencer, sample_in_sequencer in researcher_statistics['researcher_sample_data'].items():
            if len(sample_in_sequencer) == 0:
                continue
            for sample_values in sample_in_sequencer:
                runs_index_sample.append(sample_values[2])
                q30_sample_value.append(float(sample_values[6]))
                mean_sample_value.append(float(sample_values[7]))

        q30_data_in_run = get_min_mean_and_max_values(q30_sample_value,runs_index_sample, wetlab_config.NUMBER_OF_VALUES_TO_FETCH_FROM_RESEARCHER)
        mean_data_in_run = get_min_mean_and_max_values(mean_sample_value,runs_index_sample, wetlab_config.NUMBER_OF_VALUES_TO_FETCH_FROM_RESEARCHER)
        # create the graphic for q30 quality
        heading = 'Graphics for Q > 30'
        data_source = column_graphic_tupla (heading, '', 'Runs',  'Q > 30 value', 'ocean', q30_data_in_run, 'Median values')
        researcher_statistics['q30_graphic'] = FusionCharts("column3d", 'q30_graph' , "550", "350",'q30_chart' , "json", data_source).render()
        # create the graphic for mean quality
        heading = 'Graphics for Mean Quality'
        data_source = column_graphic_tupla (heading, '', 'Runs',  'Mean value', 'ocean', mean_data_in_run, 'Median values')
        researcher_statistics['mean_graphic'] = FusionCharts("column3d", 'mean_graph' , "550", "350",'mean_chart' , "json", data_source).render()
        return researcher_statistics
    else:
        # check the setting for having user name in the description column in sample sheet
        if get_configuration_value('DESCRIPTION_IN_SAMPLE_SHEET_MUST_HAVE_USERNAME') :
            researcher_statistics['ERROR'] = wetlab_config.ERROR_NOT_SAMPLES_FOR_USER_FOUND_BECAUSE_OF_CONFIGURATION_SETTINGS
        else:
            researcher_statistics['ERROR'] = wetlab_config.ERROR_USER_DOES_NOT_HAVE_ANY_SAMPLE
        researcher_statistics['ERROR'].append(researcher_name)
        return researcher_statistics
        '''


        # Create the table with projects executed by the researcher
        researcher_seq_graphs, researcher_graphs = [], []
        for sequencer, projects_name_list in projects_name_dict.items() :
            sequencer_proj = {}
            proj_data =[]
            for project_name in projects_name_list :
                proj_data.append([project_name, p_researcher_date[sequencer][project_name], p_researcher_lib_kit[sequencer][project_name],
                        p_researcher_num_sample[sequencer][project_name], '{0:,}'.format(int(p_researcher_cluster_pf_dict[sequencer][project_name])),
                        '{0:,}'.format(int(p_researcher_yield_mb_dict[sequencer][project_name])), p_researcher_q30_dict[sequencer][project_name],
                        p_researcher_mean_dict[sequencer][project_name], p_researcher_sequencer[sequencer][project_name]])
            sequencer_proj[sequencer] = proj_data
            projects_data.append(sequencer_proj)
            # create the graphic for q30 quality
            theme = 'ocean'
            heading = 'Graphics for Q > 30 for investigator ' + r_name
            sub_caption = 'Sequencer ' + sequencer
            x_axis_name = 'Projects'
            y_axis_name = 'Q 30 (in %)'

            data_source = column_graphic_simple (heading, sub_caption, x_axis_name, y_axis_name, theme, p_researcher_q30_dict[sequencer])
            seq_chart = sequencer + 'q30_chart'
            seq_graph = sequencer + 'q30_graph'
            q30_researcher_seq_graph = FusionCharts("column3d", seq_graph , "500", "350",seq_chart , "json", data_source).render()

            researcher_seq_graphs.append([seq_chart, q30_researcher_seq_graph])
            # create the graphic for mean quality
            theme = 'carbon'
            heading = 'Graphics for Mean quality for investigator ' + r_name
            sub_caption = 'Sequencer ' + sequencer
            x_axis_name = 'Projects'
            y_axis_name = 'Mean Quality'
            data_source = column_graphic_simple (heading, sub_caption, x_axis_name, y_axis_name, theme, p_researcher_mean_dict[sequencer])
            seq_chart = sequencer + 'mean_q_chart'
            seq_graph = sequencer + 'mean_q_graph'
            mean_q_researcher_seq_graph = FusionCharts("column3d", seq_graph , "500", "350", seq_chart, "json", data_source).render()
            researcher_seq_graphs.append([seq_chart, mean_q_researcher_seq_graph])

            # create the graphic for yield Mb
            theme = 'zune'
            heading = 'Graphics for Yield Mb for investigator ' + r_name
            sub_caption = 'Sequencer ' + sequencer
            x_axis_name = 'Projects'
            y_axis_name = 'Yield Mb'
            data_source = column_graphic_simple (heading, sub_caption, x_axis_name, y_axis_name, theme, p_researcher_yield_mb_dict[sequencer])
            seq_chart = sequencer + 'yield_mb_chart'
            seq_graph = sequencer + 'yield_mb_graph'
            yield_mb_researcher_graph = FusionCharts("column3d", seq_graph , "500", "350", seq_chart, "json", data_source).render()
            researcher_seq_graphs.append([seq_chart, yield_mb_researcher_graph])
            # create the graphic for cluster Pf
            theme = 'ocean'
            heading = 'Graphics for Cluster Pf for investigator ' + r_name
            sub_caption = 'Sequencer ' + sequencer
            x_axis_name = 'Projects'
            y_axis_name = 'Cluster Pf'
            data_source = column_graphic_simple (heading, sub_caption, x_axis_name, y_axis_name, theme, p_researcher_cluster_pf_dict[sequencer])
            seq_chart = sequencer + 'cluster_pf_chart'
            seq_graph = sequencer + 'cluster_pf_graph'
            cluster_pf_researcher_graph = FusionCharts("column3d", seq_graph , "500", "350", seq_chart, "json", data_source).render()
            researcher_seq_graphs.append([seq_chart, cluster_pf_researcher_graph])


            researcher_graphs.append(researcher_seq_graphs)

        researcher_statistics ['researcher_graph'] = researcher_graphs
        researcher_statistics ['researcher_name'] = r_name
        researcher_statistics['projects_data'] = projects_data


        #collecting data for comparation graphics

        # Calculating the mean for all projects performed by researcher
        comp_q30_dict, comp_mean_q_dict = {} , {}
        comp_yield_mb_dict, comp_cluster_pf_dict = {} , {}

        q30_val = p_researcher_q30_dict.values()
        mean_q_val = p_researcher_mean_dict.values()
        yield_mb_val = p_researcher_yield_mb_dict.values()
        cluster_pf = p_researcher_cluster_pf_dict.values()

        for sequencer in projects_name_dict.keys() :
            #sequencer_proj = {}
            #proj_data =[]
            comp_q30_dict[sequencer] , comp_mean_q_dict [sequencer]= {} , {}
            comp_yield_mb_dict[sequencer], comp_cluster_pf_dict [sequencer] = {}, {}
            comp_q30_dict [sequencer][r_name] = format(statistics.mean( [float(x) for x in list(p_researcher_q30_dict[sequencer].values())]),'.2f')
            comp_mean_q_dict[sequencer] [r_name] = format(statistics.mean( [float(x) for x in list(p_researcher_mean_dict[sequencer].values())]),'.2f')

            comp_yield_mb_dict[sequencer] [r_name] = sum(list(p_researcher_yield_mb_dict[sequencer].values()))
            comp_cluster_pf_dict[sequencer] [r_name] = sum(list(p_researcher_cluster_pf_dict[sequencer].values()))

        total_q_30_list, total_mean_q_list = [] , []
        total_yield_mb_list, total_cluster_pf_list = [] , []
        total_lanes_summary = {}

        for sequencer in projects_name_dict.keys() :
            runs_sequencer = RunProcess.objects.filter(usedSequencer__sequencerName__exact = sequencer)
            run_sequencer_id_list = []
            for run in runs_sequencer :
                run_sequencer_id_list.append(run.pk)

            if StatsLaneSummary.objects.filter(runprocess_id__in  = run_sequencer_id_list).exclude(defaultAll__isnull = False).exclude(project_id__in = projects_id_list[sequencer]).exists():
                total_lanes_summary[sequencer] = StatsLaneSummary.objects.filter(runprocess_id__in  = run_sequencer_id_list).exclude(defaultAll__isnull = False).exclude(project_id__in = projects_id_list[sequencer])
            else:
                total_lanes_summary[sequencer] = ''

        if len(total_lanes_summary) > 0:
            comp_graphs, comp_seq_graphs = [] , []
            for sequencer in projects_name_dict.keys() :
                for lane_summary in total_lanes_summary[sequencer] :
                    q_30_value, mean_q_value , yield_mb_value , cluster_pf_value = lane_summary.get_stats_info()
                    total_q_30_list.append(float(q_30_value))
                    total_mean_q_list.append(float(mean_q_value))
                    total_yield_mb_list.append(int(yield_mb_value.replace(',','')))
                    total_cluster_pf_list.append(int(cluster_pf_value.replace(',','')))
                comp_q30_dict[sequencer]['Other investigators']= format(statistics.mean(total_q_30_list), '.2f')
                comp_mean_q_dict[sequencer]['Other investigators'] = format(statistics.mean(total_mean_q_list), '.2f')
                comp_yield_mb_dict[sequencer]['Other investigators'] = sum(total_yield_mb_list)
                comp_cluster_pf_dict[sequencer]['Other investigators'] = sum(total_cluster_pf_list)
                # create the graphic for q30 quality

                theme = ''
                heading = 'Comparation graphics for Q > 30 for investigator ' + r_name
                sub_caption = ''
                x_axis_name = r_name + ' versus other investigators'
                y_axis_name = 'Q 30 (in %)'

                data_source = column_graphic_simple (heading, sub_caption, x_axis_name, y_axis_name, theme, comp_q30_dict[sequencer])
                seq_chart = sequencer + 'comparation_q30_chart'
                seq_graph = sequencer + 'comparation_q30_graph'
                comp_q30_seq_graph = FusionCharts("column3d", seq_graph , "500", "350",seq_chart , "json", data_source).render()
                comp_seq_graphs.append([seq_chart, comp_q30_seq_graph])

                theme = ''
                heading = 'Comparation graphics for Mean Quality for investigator ' + r_name
                sub_caption = ''
                x_axis_name = r_name + ' versus other investigators'
                y_axis_name = 'Mean Quality'
                data_source = column_graphic_simple (heading, sub_caption, x_axis_name, y_axis_name, theme, comp_mean_q_dict[sequencer])
                seq_chart = sequencer + 'comparation_mean_q_chart'
                seq_graph = sequencer + 'comparation_mean_q_graph'
                comp_mean_q_seq_graph = FusionCharts("column3d", seq_graph , "500", "350",seq_chart , "json", data_source).render()
                comp_seq_graphs.append([seq_chart, comp_mean_q_seq_graph])

                theme = ''
                heading = 'Comparation graphics for Yield (Mb) for investigator ' + r_name
                sub_caption = ''
                x_axis_name = r_name + ' versus other investigators'
                y_axis_name = '(Mb)'
                data_source = column_graphic_simple (heading, sub_caption, x_axis_name, y_axis_name, theme, comp_yield_mb_dict[sequencer])
                seq_chart = sequencer + 'comparation_yield_mb_chart'
                seq_graph = sequencer + 'comparation_yield_mb_graph'
                comp_yield_mb_seq_graph = FusionCharts("column3d", seq_graph , "500", "350",seq_chart , "json", data_source).render()
                comp_seq_graphs.append([seq_chart, comp_yield_mb_seq_graph])

                theme = ''
                heading = 'Comparation graphics for Cluster PF for investigator ' + r_name
                sub_caption = ''
                x_axis_name = r_name + ' versus other investigators'
                y_axis_name = 'Cluster pf'
                data_source = column_graphic_simple (heading, sub_caption, x_axis_name, y_axis_name, theme, comp_cluster_pf_dict[sequencer])
                seq_chart = sequencer + 'comparation_cluster_pf_chart'
                seq_graph = sequencer + 'comparation_cluster_pf_graph'
                comp_cluster_pf_seq_graph = FusionCharts("column3d", seq_graph , "500", "350",seq_chart , "json", data_source).render()
                comp_seq_graphs.append([seq_chart, comp_cluster_pf_seq_graph])
                comp_graphs.append(comp_seq_graphs)

            researcher_statistics ['comp_graphs'] = comp_graphs

        # Sequencer graphic utilization
        sequencer_used = {}
        for sequencer in projects_name_dict.keys() :
            sequencer_used[sequencer] = len( projects_name_dict[sequencer])

        theme = 'ocean'
        heading = 'Sequencer utilization for investigator ' + r_name
        sub_caption = ''
        data_source = pie_graphic_standard (heading, sub_caption, theme, sequencer_used)
        sequencer_pie_graph = FusionCharts("pie3d", "sequencer_pie_graph" , "500", "400", "sequencer_pie_chart", "json", data_source).render()
        researcher_statistics ['sequencer_pie_graph'] = sequencer_pie_graph
        '''
