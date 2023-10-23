# Generic imports
from django.db.models import F, Count
from django.db.models.functions import ExtractMonth, ExtractYear

# Local imports
import core.fusioncharts.fusioncharts
import core.models
import core.utils.common
import core.utils.graphics
import wetlab.config
import wetlab.models


def get_annual_report(rep_year):
    """Collect data from the requested year and make comparation again previous years

    Parameters
    ----------
    rep_year : string
        Year number to collect statistics data

    Return
    ------
    report_data : dictionnary
        Contains the graphics and data to display statistics
    """
    report_data = {}
    start_date = rep_year + "-01-01"
    end_date = rep_year + "-12-31"
    run_objs = wetlab.models.RunProcess.objects.filter(
        run_date__range=(start_date, end_date)
    )
    if len(run_objs) == 0:
        report_data["ERROR"] = wetlab.config.ERROR_NOT_RUNS_FOUND_IN_SELECTED_PERIOD

    sample_objs = wetlab.models.SamplesInProject.objects.filter(
        run_process_id__in=run_objs
    )
    report_data["year_num_runs"] = len(run_objs)
    report_data["year_num_samples"] = len(sample_objs)
    # Graphic chart for run and states
    rep_year_run_states = list(
        run_objs.values(sum_state=F("state__run_state_name")).annotate(
            value=Count("run_name")
        )
    )
    g_data = core.utils.graphics.preparation_graphic_data(
        "Run states in " + rep_year,
        "",
        "",
        "",
        "ocean",
        rep_year_run_states,
        "sum_state",
        "value",
    )
    report_data["rep_state_graphic"] = core.fusioncharts.fusioncharts.FusionCharts(
        "pie3d", "rep_state_graph", "600", "300", "rep_state_chart", "json", g_data
    ).render()

    # Graphic for runs per months
    run_per_month = (
        run_objs.annotate(year=ExtractYear("run_date"))
        .annotate(month=ExtractMonth("run_date"))
        .values("year", "month")
        .annotate(value=Count("run_name"))
    )
    format_run_per_month = core.utils.common.week_month_number_to_date(
        run_per_month, "month", "value", "%Y-%m-%d"
    )
    g_data = core.utils.graphics.preparation_graphic_data(
        "Run per months",
        "",
        "Date",
        "Number of runs per months",
        "ocean",
        format_run_per_month,
    )
    report_data["rep_run_month_graphic"] = core.fusioncharts.fusioncharts.FusionCharts(
        "column3d",
        "rep_run_month_graph",
        "600",
        "350",
        "rep_run_month_chart",
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
    report_data["rep_researcher_graphic"] = core.fusioncharts.fusioncharts.FusionCharts(
        "pie3d",
        "rep_researcher_graph",
        "600",
        "300",
        "rep_researcher_chart",
        "json",
        g_data,
    ).render()

    # Graphic chart for samples
    sample_per_month = (
        sample_objs.annotate(year=ExtractYear("run_process_id__run_date"))
        .annotate(month=ExtractMonth("run_process_id__run_date"))
        .values("year", "month")
        .annotate(value=Count("sample_name"))
    )
    format_sample_per_month = core.utils.common.week_month_number_to_date(
        sample_per_month, "month", "value", "%Y-%m-%d"
    )
    g_data = core.utils.graphics.preparation_graphic_data(
        "Samples per month",
        "",
        "Date",
        "Number of runs per month",
        "zune",
        format_sample_per_month,
    )
    report_data[
        "rep_sample_month_graphic"
    ] = core.fusioncharts.fusioncharts.FusionCharts(
        "column3d",
        "rep_sample_month_graph",
        "600",
        "350",
        "rep_sample_month_chart",
        "json",
        g_data,
    ).render()

    # Collecting data for previous years
    # ##################################
    all_run_objs = wetlab.models.RunProcess.objects.filter(run_date__lte=end_date)
    all_sample_objs = wetlab.models.SamplesInProject.objects.filter(
        run_process_id__in=all_run_objs
    )

    # Graphic for runs comparision
    comp_run = (
        all_run_objs.annotate(year=ExtractYear("run_date"))
        .values("year")
        .annotate(values=Count("run_name"))
        .order_by("year")
    )

    g_data = core.utils.graphics.preparation_graphic_data(
        "Runs per year",
        "",
        "Year",
        "Number of runs per year",
        "zune",
        comp_run,
        "year",
        "values",
    )
    report_data["rep_comp_run_graphic"] = core.fusioncharts.fusioncharts.FusionCharts(
        "column3d",
        "rep_comp_run_graph",
        "600",
        "350",
        "rep_comp_run_chart",
        "json",
        g_data,
    ).render()

    # Graphic for sample comparision
    comp_sample = (
        all_sample_objs.annotate(year=ExtractYear("run_process_id__run_date"))
        .values("year")
        .annotate(values=Count("sample_name"))
        .order_by("year")
    )

    g_data = core.utils.graphics.preparation_graphic_data(
        "Samples per year",
        "",
        "Year",
        "Number of samples per year",
        "zune",
        comp_sample,
        "year",
        "values",
    )
    report_data[
        "rep_comp_sample_graphic"
    ] = core.fusioncharts.fusioncharts.FusionCharts(
        "column3d",
        "rep_comp_sample_graph",
        "600",
        "350",
        "rep_comp_sample_chart",
        "json",
        g_data,
    ).render()

    # Graphic chart for samples per researcher
    comp_sample_per_researcher = list(
        all_sample_objs.values(Researcher=F("user_id__username")).annotate(
            value=Count("sample_name")
        )
    )
    g_data = core.utils.graphics.preparation_graphic_data(
        "Samples per researcher",
        "",
        "",
        "",
        "flint",
        comp_sample_per_researcher,
        "Researcher",
        "value",
    )
    report_data[
        "rep_comp_researcher_graphic"
    ] = core.fusioncharts.fusioncharts.FusionCharts(
        "pie3d",
        "rep_comp_researcher_graph",
        "600",
        "350",
        "rep_comp_researcher_chart",
        "json",
        g_data,
    ).render()

    report_data["year"] = rep_year

    return report_data
