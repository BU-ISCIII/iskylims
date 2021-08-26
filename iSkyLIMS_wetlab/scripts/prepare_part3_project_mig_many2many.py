from iSkyLIMS_wetlab.models import StatsLaneSummary, StatsFlSummary, RawDemuxStats
'''
    This file is part of the migration preparation for the version 2.0.0.
    It collect the information for the StatsLaneSummary, StatsFlSummary, RawDemuxStats
    tables for having the funcitonality that Project name can be repeated in different runs
'''

def run ():

    stats_lanes_to_update = StatsLaneSummary.objects.all()
    s_to_update = len (stats_lanes_to_update)
    if s_to_update == 0 :
        print('There is no StatsLaneSummary defined in iSkyLIMS\n')
    else:
        with open ('part3_StatsLaneSummary_migration_data.csv', 'w') as fh:
            for stat_lane in stats_lanes_to_update :
                project_id = stat_lane.project_id
                if project_id == None:
                    continue
                fh.write(str(stat_lane.pk) + ',' + str(project_id.pk) + '\n')

    stats_flsummary_to_update = StatsFlSummary.objects.all()
    fl_to_update = len (stats_flsummary_to_update)
    if fl_to_update == 0 :
        print('There is no StatsFlSummary defined in iSkyLIMS\n')
    else:
        with open ('part3_StatsFlSummary_migration_data.csv', 'w') as fh:
            for stats_flsummary in stats_flsummary_to_update :
                project_id = stats_flsummary.project_id
                if project_id == None:
                    continue
                fh.write(str(stats_flsummary.pk) + ',' + str(project_id.pk) + '\n')

    stats_demux_to_update = RawDemuxStats.objects.all()
    demux_to_update = len (stats_demux_to_update)
    if demux_to_update == 0 :
        print('There is no RawDemuxStats defined in iSkyLIMS\n')
    else:
        with open ('part3_RawDemuxStats_migration_data.csv', 'w') as fh:
            for demux in stats_demux_to_update :
                project_id = demux.project_id
                if project_id == None:
                    continue
                fh.write(str(demux.pk) + ',' + str(project_id.pk) + '\n')

    print('Successful fetched information for part 3\n')
    return
