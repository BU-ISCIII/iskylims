import clinic.models
import core.fusioncharts.fusioncharts
import core.utils.stats_graphics


def check_empty_fields(row_data):
    for data in row_data:
        if data == "":
            return True
    return False


def get_configuration_from_database(configuration_name):
    """
    Description:
        The function fetch from database the configuration setting value
    Input:
        configuration_name      # configuration settings name
    """
    configuration_value = ""
    if clinic.models.ConfigSetting.objects.filter(
        configuration_name__exact=configuration_name
    ).exists():
        configuration_settings_obj = clinic.models.ConfigSetting.objects.filter(
            configuration_name__exact=configuration_name
        ).last()
        configuration_value = configuration_settings_obj.get_configuration_value()
    return configuration_value


def pending_clinic_samples_for_grafic(pending):
    number_of_pending = {}
    number_of_pending["DEFINED"] = pending["defined"]["length"]
    number_of_pending["PATIENT UPDATE"] = pending["patient_update"]["length"]
    number_of_pending["SEQUENCING"] = pending["sequencing"]["length"]
    number_of_pending["PENDING RESULTS"] = (
        pending["pending_results"]["length"] + pending["pending_protocol"]["length"]
    )

    data_source = core.utils.stats_graphics.graphic_3D_pie(
        "Number of Pending Clinic Samples", "", "", "", "fint", number_of_pending
    )
    graphic_pending_samples = core.fusioncharts.fusioncharts.FusionCharts(
        "pie3d", "ex1", "430", "450", "chart-1", "json", data_source
    )
    return graphic_pending_samples
