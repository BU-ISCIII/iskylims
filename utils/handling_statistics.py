from iSkyLIMS_wetlab.models import *
from iSkyLIMS_wetlab.wetlab_config import *

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

    if SamplesInProject.objects.filter(user_id = name_objs[0], runprocess_id__state__runStateName__exact = 'Completed').exists():
        sample_objs = SamplesInProject.objects.filter(user_id = name_objs[0], runprocess_id__state__runStateName__exact = 'Completed').order_by('runProcess_id')
        # check if start and end date are present in the form
        if start_date != '' and end_date !='':
            sample_objs= sample_objs.filter(runprocess_id__state__runStateName__range=(start_date, end_date))
        elif start_date != '':
            sample_objs = sample_objs.filter(runprocess_id__state__runStateName__gte = start_date)
        elif end_date != '':
            sample_objs = sample_objs.filter(runprocess_id__state__runStateName__lte = end_date)
        if len(sample_objs) == 0:
            researcher_statistics['ERROR'] = wetlab_config.ERROR_NO_MATCHES_FOR_INPUT_CONDITIONS
            return researcher_statistics

        # Get data from researcher projects

        researcher_statistics ['projects_heading'] = ['Project name', 'Date', 'Libraty Kit','Samples', 'Cluster PF', 'Yield Mb', '% Q> 30', 'Mean','Sequencer ID']
        projects_data =[]
        p_researcher_date , p_researcher_num_sample = {} , {}
        p_researcher_lib_kit, p_researcher_sequencer = {}, {}
        p_researcher_q30_dict, p_researcher_mean_dict = {} , {}
        p_researcher_yield_mb_dict, p_researcher_cluster_pf_dict ={} , {}
        projects_name_dict , projects_id_list = {} , {}


        for project_researcher in r_project_by_researcher:
            q_30_list , mean_q_list = [] , []
            yield_mb_list,  cluster_pf_list = [], []
            p_name = project_researcher.get_project_name()

            sequencer_in_project = project_researcher.runprocess_id.get_run_used_sequencer()
            if not sequencer_in_project in projects_name_dict :
                p_researcher_num_sample[sequencer_in_project] ={}
                p_researcher_sequencer[sequencer_in_project] ={}
                p_researcher_date[sequencer_in_project] ={}
                p_researcher_lib_kit[sequencer_in_project] ={}
                p_researcher_q30_dict[sequencer_in_project] ={}
                p_researcher_mean_dict[sequencer_in_project] ={}
                p_researcher_yield_mb_dict[sequencer_in_project] ={}
                p_researcher_cluster_pf_dict[sequencer_in_project] ={}
                projects_name_dict[sequencer_in_project] = []
                projects_id_list[sequencer_in_project] =  []
            projects_name_dict[sequencer_in_project].append(p_name)
            r_project_id = project_researcher.id
            projects_id_list[sequencer_in_project].append(r_project_id)
            p_researcher_num_sample[sequencer_in_project][p_name] = StatsFlSummary.objects.get(project_id__exact = r_project_id).sampleNumber
            p_researcher_date [sequencer_in_project][p_name] = project_researcher.get_date()
            p_researcher_lib_kit[sequencer_in_project][p_name]= project_researcher.get_index_library_name()
            #
            p_researcher_sequencer[sequencer_in_project][p_name] = str(project_researcher.runprocess_id.sequencerModel)
            lanes_in_project = StatsLaneSummary.objects.filter( project_id__exact = r_project_id)
            for lane in lanes_in_project :
                q_30_value, mean_q_value , yield_mb_value , cluster_pf_value = lane.get_stats_info()
                q_30_list.append(float(q_30_value))
                mean_q_list.append(float(mean_q_value))
                yield_mb_list.append(float(yield_mb_value.replace(',','')))
                cluster_pf_list.append(float(cluster_pf_value.replace(',','')))
            p_researcher_q30_dict[sequencer_in_project] [p_name]= format(statistics.mean(q_30_list), '.2f')
            p_researcher_mean_dict[sequencer_in_project][p_name] = format(statistics.mean(mean_q_list), '.2f')
            p_researcher_yield_mb_dict[sequencer_in_project][p_name] = round(sum(yield_mb_list))
            p_researcher_cluster_pf_dict[sequencer_in_project][p_name] = round(sum(cluster_pf_list))

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
