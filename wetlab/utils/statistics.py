# Generic imports
from django.contrib.auth.models import User
from django.db.models import Avg, F, Count, Func, Value, CharField
from django.db.models.functions import ExtractWeek, ExtractYear

# Local imports
import core.fusioncharts.fusioncharts
import core.models
import wetlab.models
import wetlab.utils.common
import core.utils.graphics
import core.utils.common
import wetlab.config


def get_per_time_statistics(start_date, end_date):
    """_summary_

    Parameters
    ----------
    start_date : str
        Date from starting the statistics
    end_date : str
        Date from the statistics ends

    Returns
    -------
    dict

    """
    per_time_statistics = {}
    # validate date format
    if start_date != "" and not wetlab.utils.common.check_valid_date_format(start_date):
        per_time_statistics["ERROR"] = wetlab.config.ERROR_INVALID_FORMAT_FOR_DATES
        return per_time_statistics
    if end_date != "" and not wetlab.utils.common.check_valid_date_format(start_date):
        per_time_statistics["ERROR"] = wetlab.config.ERROR_INVALID_FORMAT_FOR_DATES
        return per_time_statistics
    run_objs = wetlab.models.RunProcess.objects.filter(
        run_date__range=(start_date, end_date)
    )
    if len(run_objs) == 0:
        per_time_statistics[
            "ERROR"
        ] = wetlab.config.ERROR_NOT_RUNS_FOUND_IN_SELECTED_PERIOD
    project_objs = wetlab.models.Projects.objects.filter(run_process__in=run_objs)
    sample_objs = wetlab.models.SamplesInProject.objects.filter(
        run_process_id__in=run_objs
    )

    # Graphic chart for run and states
    run_states = list(
        run_objs.values(sum_state=F("state__run_state_name")).annotate(
            value=Count("run_name")
        )
    )
    g_data = core.utils.graphics.preparation_graphic_data(
        "Run states", "", "", "", "ocean", run_states, "sum_state", "value"
    )
    per_time_statistics[
        "time_state_graphic"
    ] = core.fusioncharts.fusioncharts.FusionCharts(
        "pie3d", "time_state_graph", "600", "300", "time_state_chart", "json", g_data
    ).render()
    # Graphic chart for number of runs per weeks
    run_per_date = (
        run_objs.annotate(year=ExtractYear("run_date"))
        .annotate(week=ExtractWeek("run_date"))
        .values("year", "week")
        .annotate(value=Count("run_name"))
    )
    format_run_per_date = core.utils.common.week_month_number_to_date(
        run_per_date, "week", "value", "%Y-%m-%d"
    )
    g_data = core.utils.graphics.preparation_graphic_data(
        "Run per week",
        "",
        "Date",
        "Number of runs per week",
        "ocean",
        format_run_per_date,
    )
    per_time_statistics[
        "time_run_weeks_graphic"
    ] = core.fusioncharts.fusioncharts.FusionCharts(
        "column3d",
        "time_run_weeks_graph",
        "600",
        "350",
        "time_run_weeks_chart",
        "json",
        g_data,
    ).render()
    # Graphic chart for projects
    project_per_date = (
        project_objs.annotate(year=ExtractYear("run_process__run_date"))
        .annotate(week=ExtractWeek("run_process__run_date"))
        .values("year", "week")
        .annotate(value=Count("project_name"))
    )
    format_project_per_date = core.utils.common.week_month_number_to_date(
        project_per_date, "week", "value", "%Y-%m-%d"
    )
    g_data = core.utils.graphics.preparation_graphic_data(
        "Projects per week",
        "",
        "Date",
        "Number of runs per week",
        "zune",
        format_project_per_date,
    )
    per_time_statistics[
        "time_project_weeks_graphic"
    ] = core.fusioncharts.fusioncharts.FusionCharts(
        "column3d",
        "time_project_weeks_graph",
        "600",
        "350",
        "time_project_weeks_chart",
        "json",
        g_data,
    ).render()

    # Graphic chart for samples per researcher
    sample_per_researcher = list(
        sample_objs.values(Researcher=F("user_id__username")).annotate(
            value=Count("sample_name")
        )
    )
    g_data = core.utils.graphics.preparation_graphic_data(
        "Samples per researcher",
        "",
        "",
        "",
        "flint",
        sample_per_researcher,
        "Researcher",
        "value",
    )
    per_time_statistics[
        "time_researcher_graphic"
    ] = core.fusioncharts.fusioncharts.FusionCharts(
        "pie3d",
        "time_researcher_graph",
        "600",
        "300",
        "time_researcher_chart",
        "json",
        g_data,
    ).render()

    # Graphic chart for unknown barcodes
    barcode_objs = wetlab.models.RawTopUnknowBarcodes.objects.filter(
        runprocess_id__in=run_objs
    )

    # chart graph for Q > 30 based on runs
    # ##############################
    researcher_q_30 = list(
        sample_objs.values(run_name=F("run_process_id__run_name"))
        .annotate(q_30_value=Avg("quality_q30"))
        .order_by("run_process_id__run_name")
    )

    g_data = core.utils.graphics.preparation_graphic_data(
        "Quality Q > 30 for each run",
        "",
        "Run name",
        "Percentage of Q>30",
        "ocean",
        researcher_q_30,
        "run_name",
        "q_30_value",
    )
    per_time_statistics[
        "time_q_30_graphic"
    ] = core.fusioncharts.fusioncharts.FusionCharts(
        "column3d", "time_q_30_graph", "550", "350", "time_q_30_chart", "json", g_data
    ).render()

    # chart graph for mean based on runs
    # ####################
    researcher_mean = list(
        sample_objs.values(run_name=F("run_process_id__run_name"))
        .annotate(mean_value=Avg("mean_quality"))
        .order_by("run_process_id__run_name")
    )
    g_data = core.utils.graphics.preparation_graphic_data(
        "Qualiy mean for each run",
        "",
        "Run name",
        "Quality mean",
        "ocean",
        researcher_mean,
        "run_name",
        "mean_value",
    )
    per_time_statistics[
        "time_mean_graphic"
    ] = core.fusioncharts.fusioncharts.FusionCharts(
        "column3d", "time_mean_graph", "550", "350", "time_mean_chart", "json", g_data
    ).render()

    # chart graph for Q > 30 based on researcher
    # ##############################
    researcher_q_30 = list(
        sample_objs.values(run_name=F("user_id__username"))
        .annotate(q_30_value=Avg("quality_q30"))
        .order_by("user_id__username")
    )

    g_data = core.utils.graphics.preparation_graphic_data(
        "Quality Q > 30 for reseacher",
        "",
        "Run name",
        "Percentage of Q>30",
        "zune",
        researcher_q_30,
        "run_name",
        "q_30_value",
    )
    per_time_statistics[
        "time_researcher_q_30_graphic"
    ] = core.fusioncharts.fusioncharts.FusionCharts(
        "column3d",
        "time_researcher_q_30_graph",
        "550",
        "350",
        "time_researcher_q_30_chart",
        "json",
        g_data,
    ).render()

    # chart graph for mean based on researcher
    # ####################
    researcher_mean = list(
        sample_objs.values(run_name=F("user_id__username"))
        .annotate(mean_value=Avg("mean_quality"))
        .order_by("user_id__username")
    )
    g_data = core.utils.graphics.preparation_graphic_data(
        "Qualiy mean for researcher",
        "",
        "Run name",
        "Quality mean",
        "zune",
        researcher_mean,
        "run_name",
        "mean_value",
    )
    per_time_statistics[
        "time_researcher_mean_graphic"
    ] = core.fusioncharts.fusioncharts.FusionCharts(
        "column3d",
        "time_researcher_mean_graph",
        "550",
        "350",
        "time_researcher_mean_chart",
        "json",
        g_data,
    ).render()

    # Table information for run data
    per_time_statistics["run_data"] = list(
        run_objs.values_list(
            "pk", "run_name", "state__run_state_name", "used_sequencer__sequencer_name"
        ).annotate(
            formated_date=Func(
                F("run_date"),
                Value("%Y-%m-%d"),
                function="DATE_FORMAT",
                output_field=CharField(),
            )
        )
    )
    per_time_statistics[
        "run_table_heading"
    ] = wetlab.config.HEADING_STATISTICS_FOR_TIME_RUN

    # Table information for sample data
    per_time_statistics["sample_data"] = list(
        sample_objs.values_list(
            "pk",
            "sample_name",
            "user_id__username",
            "project_id__project_name",
            "run_process_id__run_name",
            "barcode_name",
        )
    )
    per_time_statistics[
        "sample_table_heading"
    ] = wetlab.config.HEADING_STATISTICS_FOR_TIME_SAMPLE
    per_time_statistics["start_date"] = start_date
    per_time_statistics["end_date"] = end_date
    per_time_statistics["num_runs"] = len(run_objs)
    per_time_statistics["num_projects"] = len(project_objs)
    # import pdb; pdb.set_trace()
    return per_time_statistics


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
        researcher_statistics["ERROR"] = wetlab.config.ERROR_INVALID_FORMAT_FOR_DATES
        return researcher_statistics
    if end_date != "" and not wetlab.utils.common.check_valid_date_format(start_date):
        researcher_statistics["ERROR"] = wetlab.config.ERROR_INVALID_FORMAT_FOR_DATES
        return researcher_statistics

    # check if start and end date are present in the form
    if start_date != "" and end_date != "":
        sample_objs = wetlab.models.SamplesInProject.objects.filter(
            run_process_id__run_date__range=(start_date, end_date)
        )
    elif start_date != "":
        sample_objs = wetlab.models.SamplesInProject.objects.filter(
            run_process_id__run_date__gte=start_date
        )
    elif end_date != "":
        sample_objs = wetlab.models.SamplesInProject.objects.filter(
            run_process_id__run_date__lte=end_date
        )
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
    researcher_statistics["samples"] = user_sample_objs.values_list(
        "sample_name",
        "project_id__project_name",
        "run_process_id__run_name",
        "run_process_id__used_sequencer__sequencer_name",
    )
    researcher_statistics[
        "table_heading"
    ] = wetlab.config.HEADING_STATISTICS_FOR_RESEARCHER_SAMPLE

    # pie graph percentage researcher vs others
    per_data_user = {}
    per_data_user[researcher_name] = user_sample_objs.count()
    per_data_user["all researchers"] = other_user_sample_objs.count()
    # heading, sub_title, axis_x_description, axis_y_description, theme, source_data
    g_data = core.utils.graphics.preparation_3D_pie(
        "Percentage of samples", "Research vs all", "ocean", per_data_user
    )

    researcher_statistics[
        "research_vs_other_graphic"
    ] = core.fusioncharts.fusioncharts.FusionCharts(
        "pie3d",
        "research_vs_other_graph",
        "600",
        "300",
        "research_vs_other_chart",
        "json",
        g_data,
    ).render()

    # pie graph for sequencers used
    # #############################
    seq_objs = core.models.SequencerInLab.objects.all()
    sample_per_sequencer = {}
    for seq_obj in seq_objs:
        sample_per_sequencer[seq_obj.get_sequencer_name()] = user_sample_objs.filter(
            run_process_id__used_sequencer=seq_obj
        ).count()
    g_data = core.utils.graphics.preparation_3D_pie(
        "Sequencer usage", "", "ocean", sample_per_sequencer
    )
    researcher_statistics[
        "research_usage_sequencer_graphic"
    ] = core.fusioncharts.fusioncharts.FusionCharts(
        "pie3d",
        "research_usage_sequencer_graph",
        "600",
        "300",
        "research_usage_sequencer_chart",
        "json",
        g_data,
    ).render()

    # chart graph for runs
    # ####################
    researcher_runs = {}
    runs = list(
        user_sample_objs.values_list("run_process_id__run_name", flat=True).distinct()
    )
    for run in runs:
        researcher_runs[run] = user_sample_objs.filter(
            run_process_id__run_name__exact=run
        ).count()

    g_data = core.utils.graphics.preparation_graphic_data(
        "Number of samples per run",
        "",
        "Run name",
        "Number of samples",
        "ocean",
        researcher_runs,
    )
    researcher_statistics[
        "research_run_graphic"
    ] = core.fusioncharts.fusioncharts.FusionCharts(
        "column3d",
        "research_run_graph",
        "550",
        "350",
        "research_run_chart",
        "json",
        g_data,
    ).render()

    # chart graph for projects
    # ########################
    researcher_projects = {}
    projects = list(
        user_sample_objs.values_list("project_id__project_name", flat=True).distinct()
    )
    for project in projects:
        researcher_projects[project] = user_sample_objs.filter(
            project_id__project_name__exact=project
        ).count()
    g_data = core.utils.graphics.preparation_graphic_data(
        "Number of samples per project",
        "",
        "Project name",
        "Number of samples",
        "ocean",
        researcher_projects,
    )
    researcher_statistics[
        "research_project_graphic"
    ] = core.fusioncharts.fusioncharts.FusionCharts(
        "column3d",
        "research_project_graph",
        "550",
        "350",
        "research_project_chart",
        "json",
        g_data,
    ).render()

    # chart graph for Q > 30 on runs
    # ##############################
    researcher_q_30 = list(
        user_sample_objs.values(run_name=F("run_process_id__run_name"))
        .annotate(q_30_value=Avg("quality_q30"))
        .order_by("run_process_id__run_name")
    )

    g_data = core.utils.graphics.preparation_graphic_data(
        "Percentage of samples with Q > 30",
        "",
        "Run name",
        "Percentage of Q>30",
        "ocean",
        researcher_q_30,
        "run_name",
        "q_30_value",
    )
    researcher_statistics[
        "research_q_30_graphic"
    ] = core.fusioncharts.fusioncharts.FusionCharts(
        "column3d",
        "research_q_30_graph",
        "550",
        "350",
        "research_q_30_chart",
        "json",
        g_data,
    ).render()

    # chart graph for mean
    # ####################
    researcher_mean = list(
        user_sample_objs.values(run_name=F("run_process_id__run_name"))
        .annotate(mean_value=Avg("mean_quality"))
        .order_by("run_process_id__run_name")
    )
    g_data = core.utils.graphics.preparation_graphic_data(
        "Qualiy mean of samples per run",
        "",
        "Run name",
        "Quality mean",
        "ocean",
        researcher_mean,
        "run_name",
        "mean_value",
    )
    researcher_statistics[
        "research_mean_graphic"
    ] = core.fusioncharts.fusioncharts.FusionCharts(
        "column3d",
        "research_mean_graph",
        "550",
        "350",
        "research_mean_chart",
        "json",
        g_data,
    ).render()
    researcher_statistics["researcher_name"] = researcher_name

    return researcher_statistics


def get_pending_graphic_data(
    pend_data, heading, theme, ex_value, width, heigth, chart_value
):

    g_data = core.utils.graphics.preparation_graphic_data(
        heading, "", "", "", theme, pend_data
    )
    return core.fusioncharts.fusioncharts.FusionCharts(
        "pie3d", ex_value, width, heigth, chart_value, "json", g_data
    ).render()


def get_sequencer_statistics(sequencer_name, start_date, end_date):
    """_summary_

    Parameters
    ----------
    sequencer_name : str
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
    sequencer_statistics = {}
    if not core.models.SequencerInLab.objects.filter(
        sequencer_name__iexact=sequencer_name
    ).exists():
        sequencer_statistics[
            "ERROR"
        ] = wetlab.config.ERROR_NO_MATCHES_FOR_INPUT_CONDITIONS
        return sequencer_statistics
    # validate date format
    if start_date != "" and not wetlab.utils.common.check_valid_date_format(start_date):
        sequencer_statistics["ERROR"] = wetlab.config.ERROR_INVALID_FORMAT_FOR_DATES
        return sequencer_statistics
    if end_date != "" and not wetlab.utils.common.check_valid_date_format(start_date):
        sequencer_statistics["ERROR"] = wetlab.config.ERROR_INVALID_FORMAT_FOR_DATES
        return sequencer_statistics

    if start_date != "" and end_date != "":
        run_objs = wetlab.models.RunProcess.objects.filter(
            run_date__range=(start_date, end_date)
        )
    elif start_date != "":
        run_objs = wetlab.models.RunProcess.objects.filter(run_date__gte=start_date)
    elif end_date != "":
        run_objs = wetlab.models.RunProcess.objects.filter(run_date__lte=end_date)
    else:
        run_objs = wetlab.models.RunProcess.objects.all()
    all_run_objs = run_objs
    other_run_objs = run_objs.exclude(
        used_sequencer__sequencer_name__iexact=sequencer_name
    )
    run_objs = run_objs.filter(used_sequencer__sequencer_name__iexact=sequencer_name)

    if len(run_objs) == 0:
        sequencer_statistics[
            "ERROR"
        ] = wetlab.config.ERROR_NO_MATCHES_FOR_INPUT_CONDITIONS
        return sequencer_statistics

    # Chart graphic for run states
    run_states = list(
        run_objs.values(sum_state=F("state__run_state_name")).annotate(
            value=Count("run_name")
        )
    )
    g_data = core.utils.graphics.preparation_graphic_data(
        "Run states handled by sequencer",
        "",
        "",
        "",
        "ocean",
        run_states,
        "sum_state",
        "value",
    )
    sequencer_statistics[
        "sequencer_state_graphic"
    ] = core.fusioncharts.fusioncharts.FusionCharts(
        "pie3d",
        "sequencer_state_graph",
        "600",
        "300",
        "sequencer_state_chart",
        "json",
        g_data,
    ).render()

    # Chart graphic for comparation of sequencers
    sequencers_usage = list(
        all_run_objs.values(
            sequencer_name=F("used_sequencer__sequencer_name")
        ).annotate(value=Count("run_name"))
    )
    g_data = core.utils.graphics.preparation_graphic_data(
        "Sequencer utilization for the same period of time",
        "",
        "",
        "",
        "ocean",
        sequencers_usage,
        "sequencer_name",
        "value",
    )
    sequencer_statistics[
        "sequencer_usage_graphic"
    ] = core.fusioncharts.fusioncharts.FusionCharts(
        "pie3d",
        "sequencer_usage_graph",
        "600",
        "300",
        "sequencer_usage_chart",
        "json",
        g_data,
    ).render()

    # run table
    sequencer_statistics["run_table_data"] = list(
        run_objs.values_list(
            "pk", "run_name", "state__run_state_name", "used_sequencer__sequencer_name"
        ).annotate(
            formated_date=Func(
                F("run_date"),
                Value("%Y-%m-%d"),
                function="DATE_FORMAT",
                output_field=CharField(),
            )
        )
    )
    sequencer_statistics[
        "run_table_heading"
    ] = wetlab.config.HEADING_STATISTICS_FOR_SEQUENCER_RUNS

    sequencer_statistics["sequencer_name"] = sequencer_name
    return sequencer_statistics
