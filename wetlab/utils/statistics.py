import statistics
from django.contrib.auth.models import User
from django.db.models import Avg, F
import core.fusioncharts.fusioncharts
import core.models
import wetlab.models
import wetlab.utils.common
# from wetlab.utils.fetch_info import *
import core.utils.graphics # import (column_graphic_simple,
                                   #               column_graphic_tupla)
import wetlab.config


def get_min_mean_and_max_values(values_data, reference_data, number_to_split):
    """
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
    """
    reference_query_values = []
    value_sort = values_data.copy()
    value_sort.sort()
    index_used_in_reference = []
    if len(values_data) <= 2 * number_to_split:
        list_of_values = value_sort
    else:
        list_of_values = value_sort[0:number_to_split] + value_sort[-number_to_split:]
    for val in list_of_values:
        tmp_run_list_index = [
            index for index, value in enumerate(values_data) if value == val
        ]

        for tmp_index in tmp_run_list_index:
            if tmp_index in index_used_in_reference:
                continue
            else:
                break
        reference_query_values.append([reference_data[tmp_index], val])
        index_used_in_reference.append(tmp_index)
    if len(values_data) > 2 * number_to_split:
        # Add the median value
        reference_query_values.insert(
            number_to_split + 1,
            [
                "Median values",
                round(
                    statistics.median(value_sort[number_to_split:-number_to_split]), 2
                ),
            ],
        )
    return reference_query_values


def get_researcher_statistics(researcher_name, start_date, end_date):
    """_summary_

    Parameters
    ----------
    researcher_name : str
        _description_
    start_date : str
        _description_
    end_date : str, 
        _description_, by default None

    Returns
    -------
    _type_
        _description_
    """
    researcher_statistics = {}
    if not User.objects.filter(username__icontains=researcher_name).exists():
        researcher_statistics[
            "ERROR"
        ] = wetlab.config.ERROR_NO_MATCHES_FOR_INPUT_CONDITIONS
        return researcher_statistics

    user_objs = User.objects.filter(username__icontains=researcher_name)

    if len(user_objs) > 1:
        researcher_statistics[
            "ERROR"
        ] = wetlab.config.ERROR_MANY_USER_MATCHES_FOR_INPUT_CONDITIONS
        return researcher_statistics
    researcher_name = user_objs[0].username
    # validate date format
    if start_date != "" and not wetlab.utils.common.check_valid_date_format(start_date):
        researcher_statistics[
            "ERROR"
        ] = wetlab.config.ERROR_INVALID_FORMAT_FOR_DATES
        return researcher_statistics
    if end_date != "" and not wetlab.utils.common.check_valid_date_format(start_date):
        researcher_statistics[
            "ERROR"
        ] = wetlab.config.ERROR_INVALID_FORMAT_FOR_DATES
        return researcher_statistics

    # check if start and end date are present in the form
    if start_date != "" and end_date != "":
        sample_objs = wetlab.models.SamplesInProject.objects.filter(
            run_process_id__run_date__range=(start_date, end_date)
        )
    elif start_date != "":
        sample_objs = wetlab.models.SamplesInProject.objects.filter(run_process_id__run_date__gte=start_date)
    elif end_date != "":
        sample_objs = wetlab.models.SamplesInProject.objects.filter(run_process_id__run_date__lte=end_date)
    else:
        sample_objs = wetlab.models.SamplesInProject.objects.all()
            
    other_user_sample_objs = sample_objs.exclude(user_id=user_objs[0])
    user_sample_objs = sample_objs.filter(user_id=user_objs[0])
    if len(user_sample_objs) == 0:
        researcher_statistics[
            "ERROR"
        ] = wetlab.config.ERROR_NO_MATCHES_FOR_INPUT_CONDITIONS
        return researcher_statistics
    # sample table
    researcher_statistics["samples"] = user_sample_objs.values_list("sample_name", "project_id__project_name",  "run_process_id__run_name", "run_process_id__used_sequencer__sequencer_name")
    researcher_statistics["table_heading"]  = wetlab.config.RESEARCHER_SAMPLE_HEADING_STATISTICS
    
    # pie graph percentage researcher vs others 
    per_data_user = {}
    per_data_user[researcher_name] = user_sample_objs.count()
    per_data_user["all researchers"] = other_user_sample_objs.count()
    # heading, sub_title, axis_x_description, axis_y_description, theme, source_data
    g_data = core.utils.graphics.preparation_3D_pie("Percentage of samples", "Research vs all", "ocean", per_data_user)
    
    researcher_statistics["research_vs_other_graphic"] = core.fusioncharts.fusioncharts.FusionCharts(
            "pie3d", "research_vs_other_graph", "600", "300", "research_vs_other_chart", "json", g_data
        ).render()
    
    # pie graph for sequencers used
    # #############################
    seq_objs = core.models.SequencerInLab.objects.all()
    sample_per_sequencer = {}
    for seq_obj in seq_objs:
        sample_per_sequencer[seq_obj.get_sequencer_name()] = user_sample_objs.filter(run_process_id__used_sequencer=seq_obj).count()
    g_data = core.utils.graphics.preparation_3D_pie("Sequencer usage", "", "ocean", sample_per_sequencer)
    researcher_statistics["research_usage_sequencer_graphic"] = core.fusioncharts.fusioncharts.FusionCharts(
        "pie3d", "research_usage_sequencer_graph", "600", "300", "research_usage_sequencer_chart", "json", g_data
    ).render()
    
    # chart graph for runs
    # ####################
    researcher_runs = {}
    runs = list(user_sample_objs.values_list("run_process_id__run_name", flat=True).distinct())
    for run in runs:
        researcher_runs[run] = user_sample_objs.filter(run_process_id__run_name__exact=run).count()

    g_data = core.utils.graphics.preparation_bar_column("Number of samples per run", "", "Run name", "Number of samples", "ocean", researcher_runs)
    researcher_statistics["research_run_graphic"] = core.fusioncharts.fusioncharts.FusionCharts(
        "column3d", "research_run_graph", "550", "350", "research_run_chart", "json", g_data
    ).render()
        
    # chart graph for projects
    # ######################## 
    researcher_projects = {}
    projects = list(user_sample_objs.values_list("project_id__project_name", flat=True).distinct())
    for project in projects:
        researcher_projects[project] = user_sample_objs.filter(project_id__project_name__exact=project).count()
    g_data = core.utils.graphics.preparation_bar_column("Number of samples per project", "", "Project name", "Number of samples", "ocean", researcher_projects)
    researcher_statistics["research_project_graphic"] = core.fusioncharts.fusioncharts.FusionCharts(
        "column3d", "research_project_graph", "550", "350", "research_project_chart", "json", g_data
    ).render()
        
    # chart graph for Q > 30 on runs
    # ##############################
    researcher_q_30 = list(user_sample_objs.values(run_name=F("run_process_id__run_name")).annotate(q_30_value=Avg("quality_q30")).order_by("run_process_id__run_name"))
    
    g_data = core.utils.graphics.preparation_bar_column("Percentage of samples with Q > 30", "", "Run name", "Percentage of Q>30", "ocean", researcher_q_30, "run_name", "q_30_value")
    researcher_statistics["research_q_30_graphic"] = core.fusioncharts.fusioncharts.FusionCharts(
        "column3d", "research_q_30_graph", "550", "350", "research_q_30_chart", "json", g_data
    ).render()
    
    # chart graph for mean
    # ####################
    researcher_mean = list(user_sample_objs.values(run_name=F("run_process_id__run_name")).annotate(mean_value=Avg("mean_quality")).order_by("run_process_id__run_name"))
    g_data = core.utils.graphics.preparation_bar_column("Qualiy mean of samples per run", "", "Run name", "Quality mean", "ocean", researcher_mean, "run_name", "mean_quality")
    researcher_statistics["research_mean_graphic"] = core.fusioncharts.fusioncharts.FusionCharts(
        "column3d", "research_mean_graph", "550", "350", "research_mean_chart", "json", g_data
    ).render()
    # import pdb; pdb.set_trace()
    researcher_statistics["researcher_name"] = researcher_name
    """
        # collect number of saples for run and for projects
        nun_sample_in_run = {}
        num_sample_in_project = {}
        for sequencer, sample_in_sequencer in researcher_statistics[
            "researcher_sample_data"
        ].items():
            if len(sample_in_sequencer) == 0:
                continue
            for sample_values in sample_in_sequencer:
                run_name = sample_values[2]
                project_name = sample_values[1]
                if run_name not in nun_sample_in_run:
                    nun_sample_in_run[run_name] = 0
                nun_sample_in_run[run_name] += 1
                if project_name not in num_sample_in_project:
                    num_sample_in_project[project_name] = 0
                num_sample_in_project[project_name] += 1
        
        # Create run grapic
        heading = "Graphics for number of samples in runs"
        data_source = column_graphic_simple(
            heading, "", "Runs", "Number of Samples", "ocean", nun_sample_in_run
        )
        researcher_statistics["run_graphic"] = FusionCharts(
            "column3d", "run_graph", "500", "350", "run_chart", "json", data_source
        ).render()
        # Create projert grapic
        heading = "Graphics for number of samples in projects"
        data_source = column_graphic_simple(
            heading, "", "Projects", "Number of Samples", "ocean", num_sample_in_project
        )
        researcher_statistics["project_graphic"] = FusionCharts(
            "column3d",
            "project_graph",
            "500",
            "350",
            "project_chart",
            "json",
            data_source,
        ).render()

        # collect Q> 30  and mean data for each sequencer used
        runs_index_sample = []
        q30_sample_value = []
        mean_sample_value = []
        for sequencer, sample_in_sequencer in researcher_statistics[
            "researcher_sample_data"
        ].items():
            if len(sample_in_sequencer) == 0:
                continue
            for sample_values in sample_in_sequencer:
                runs_index_sample.append(sample_values[2])
                q30_sample_value.append(float(sample_values[6]))
                mean_sample_value.append(float(sample_values[7]))

        q30_data_in_run = get_min_mean_and_max_values(
            q30_sample_value,
            runs_index_sample,
            config.NUMBER_OF_VALUES_TO_FETCH_FROM_RESEARCHER,
        )
        mean_data_in_run = get_min_mean_and_max_values(
            mean_sample_value,
            runs_index_sample,
            config.NUMBER_OF_VALUES_TO_FETCH_FROM_RESEARCHER,
        )
        # create the graphic for q30 quality
        heading = "Graphics for Q > 30"
        data_source = column_graphic_tupla(
            heading,
            "",
            "Runs",
            "Q > 30 value",
            "ocean",
            q30_data_in_run,
            "Median values",
        )
        researcher_statistics["q30_graphic"] = FusionCharts(
            "column3d", "q30_graph", "550", "350", "q30_chart", "json", data_source
        ).render()
        # create the graphic for mean quality
        heading = "Graphics for Mean Quality"
        data_source = column_graphic_tupla(
            heading,
            "",
            "Runs",
            "Mean value",
            "ocean",
            mean_data_in_run,
            "Median values",
        )
        researcher_statistics["mean_graphic"] = FusionCharts(
            "column3d", "mean_graph", "550", "350", "mean_chart", "json", data_source
        ).render()
        return researcher_statistics
    else:
        # check the setting for having user name in the description column in sample sheet
        if get_configuration_value("DESCRIPTION_IN_SAMPLE_SHEET_MUST_HAVE_USERNAME"):
            researcher_statistics[
                "ERROR"
            ] = (
                config.ERROR_NOT_SAMPLES_FOR_USER_FOUND_BECAUSE_OF_CONFIGURATION_SETTINGS
            )
        else:
            researcher_statistics[
                "ERROR"
            ] = config.ERROR_USER_DOES_NOT_HAVE_ANY_SAMPLE
        researcher_statistics["ERROR"].append(researcher_name)
    """
    return researcher_statistics
