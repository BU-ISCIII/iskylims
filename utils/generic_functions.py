from iSkyLIMS_core.fusioncharts.fusioncharts import FusionCharts
from iSkyLIMS_core.utils.stats_graphics import *

def check_empty_fields (row_data):
    for data in row_data:
        if data == '':
            return True
    return False
def pending_clinic_samples_for_grafic(pending):
    number_of_pending = {}
    number_of_pending ['DEFINED'] = pending['defined']['length']
    number_of_pending ['PATIENT UPDATE'] = pending['patient_update']['length']
    number_of_pending ['SEQUENCING'] = pending['sequencing']['length']
    number_of_pending ['PENDING RESULTS'] = pending['pending_results']['length'] + pending['pending_protocol']['length']

    data_source = graphic_3D_pie('Number of Pending Clinic Samples', '', '', '', 'fint',number_of_pending)
    graphic_pending_samples = FusionCharts("pie3d", "ex1" , "430", "450", "chart-1", "json", data_source)
    return graphic_pending_samples
