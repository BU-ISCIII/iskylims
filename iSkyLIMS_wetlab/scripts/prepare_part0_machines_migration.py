from iSkyLIMS_wetlab.models import RunProcess
#from iSkyLIMS_core.models import SequencerInLab
#from iSkyLIMS_drylab.models import Machines
#import sys


def run ():
    '''
    The script migration is needed to copy the sequencer information from the table in
    drylab to core.

    '''

    run_objs = RunProcess.objects.all()
    number_runs = str(len(run_objs))
    print('Starting migration ', number_runs , ' runs were found' )
    with open ('part_0_machines_migration.csv', 'w') as fh:
        for run_obj in run_objs :
            if run_obj.sequencerModel != None:
                fh.write(str(str(run_obj.pk) + ',' + run_obj.sequencerModel.machineName + '\n' ))
    print('Successful fetched information for part 0\n')
    return
