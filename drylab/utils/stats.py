# Generic import
import django.contrib.auth.models
from django.db.models import F, Count, Func, Value, CharField
from django.db.models.functions import ExtractWeek, ExtractYear

# Local imports
import core.fusioncharts
import core.utils.graphics
import drylab.config
import drylab.models


def create_statistics_by_user(user_id, start_date=None, end_date=None):
    stats_info = {}
    stats_info["service_by_user"] = []
    stats_info["user_name"] = (
        django.contrib.auth.models.User.objects.filter(pk__exact=user_id)
        .last()
        .username
    )

    # check if start and end date are present in the form
    if start_date != "" and end_date != "":
        service_objs = drylab.models.Service.objects.filter(
            service_created_date__range=(start_date, end_date)
        )
    elif start_date != "":
        service_objs = drylab.models.Service.objects.filter(
            service_created_date__gte=start_date
        )
    elif end_date != "":
        service_objs = drylab.models.Service.objects.filter(
            service_created_date__lte=end_date
        )
    else:
        service_objs = drylab.models.Service.objects.all()

    other_service_objs = service_objs.exclude(service_user_id__pk=user_id)
    service_objs = service_objs.filter(service_user_id__pk__exact=user_id)

    if len(service_objs) == 0:
        stats_info["ERROR"] = (
            drylab.config.ERROR_NO_MATCHES_FOUND_FOR_YOUR_SERVICE_SEARCH
        )
        return stats_info

    # create table of services
    stats_info["service_by_user"] = list(
        service_objs.values_list(
            "pk", "service_request_number", "service_state__state_display"
        )
        .annotate(
            create_date=Func(
                F("service_created_date"),
                Value("%Y-%m-%d"),
                function="DATE_FORMAT",
                output_field=CharField(),
            )
        )
        .annotate(
            approved_date=Func(
                F("service_approved_date"),
                Value("%Y-%m-%d"),
                function="DATE_FORMAT",
                output_field=CharField(),
            )
        )
        .annotate(
            delivery_date=Func(
                F("service_delivered_date"),
                Value("%Y-%m-%d"),
                function="DATE_FORMAT",
                output_field=CharField(),
            )
        )
    )

    # create graphic for user service states
    serv_per_state = list(
        service_objs.values(sum_services=F("service_state__state_display")).annotate(
            value=Count("service_request_number")
        )
    )
    g_data = core.utils.graphics.preparation_graphic_data(
        "Service states", "", "", "", "ocean", serv_per_state, "sum_services", "value"
    )
    stats_info["user_services_graphic"] = core.fusioncharts.fusioncharts.FusionCharts(
        "pie3d",
        "user_services_graph",
        "580",
        "300",
        "user_services_chart",
        "json",
        g_data,
    ).render()

    # create graphic for services requested for all users
    service_user = {}
    service_user[stats_info["user_name"]] = service_objs.count()
    service_user["all researchers"] = other_service_objs.count()
    # heading, sub_title, axis_x_description, axis_y_description, theme, source_data
    g_data = core.utils.graphics.preparation_3D_pie(
        "Percentage of services", "Research vs all", "ocean", service_user
    )

    stats_info["research_vs_other_graphic"] = (
        core.fusioncharts.fusioncharts.FusionCharts(
            "pie3d",
            "research_vs_other_graph",
            "580",
            "300",
            "research_vs_other_chart",
            "json",
            g_data,
        ).render()
    )

    # getting statistics of the created services per week
    service_per_week = (
        service_objs.annotate(year=ExtractYear("service_created_date"))
        .annotate(week=ExtractWeek("service_created_date"))
        .values("year", "week")
        .annotate(value=Count("service_request_number"))
    )
    format_service_per_week = core.utils.common.week_month_number_to_date(
        service_per_week, "week", "value", "%Y-%m-%d"
    )
    g_data = core.utils.graphics.preparation_graphic_data(
        "Services per week",
        "",
        "Date",
        "Number of services per week",
        "ocean",
        format_service_per_week,
    )
    stats_info["research_service_weeks_graphic"] = (
        core.fusioncharts.fusioncharts.FusionCharts(
            "column3d",
            "research_service_weeks_graph",
            "580",
            "400",
            "research_service_weeks_chart",
            "json",
            g_data,
        ).render()
    )

    # create graphic for requested service on level 2
    avail_serv_level_2 = drylab.models.AvailableService.objects.filter(level=2)
    serv_level_2_objs = service_objs.filter(
        service_available_service__in=avail_serv_level_2
    )
    services_level_2 = serv_level_2_objs.values(
        sum_avail_serv=F("service_available_service__avail_service_description")
    ).annotate(value=Count("service_request_number"))
    g_data = core.utils.graphics.preparation_graphic_data(
        "Available services requested",
        "",
        "",
        "",
        "fint",
        services_level_2,
        "sum_avail_serv",
        "value",
    )

    stats_info["research_avail_services_graphic"] = (
        core.fusioncharts.fusioncharts.FusionCharts(
            "column3d",
            "research_avail_services_graph",
            "580",
            "400",
            "research_avail_services_chart",
            "json",
            g_data,
        ).render()
    )

    return stats_info
