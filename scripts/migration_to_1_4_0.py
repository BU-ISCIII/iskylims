from iSkyLIMS_wetlab.models import RunProcess
from iSkyLIMS_core.models import SequencerInLab
from iSkyLIMS_drylab.models import Machines
import sys

def check_migration_is_allow():
    sequencer_objs = get_list_defined_sequencers()
    if not sequencer_objs :
        print('No Sequencers were found. Define them before starts the migration')
        sys.exit (1)
    machines_objs = get_list_machines()
    if not machines_objs :
        print('No Machines were found. There is not need to run migration')
        sys.exit (1)
    sequencer_list = []
    for sequencer in sequencer_objs:
        sequencer_list.append(sequencer.get_sequencer_name())

    for machine in machines_objs :
        if not machine.get_machine_name() in sequencer_list:
            print('Machine ',  machine.get_machine_name() ,' not defined in sequencer in lab database')
            sys.exit (1)
    if not RunProcess.objects.all().exists():
        print('No Runs were found. There is not need to run migration')
        sys.exit (1)
    return

def get_list_defined_sequencers ():
    if SequencerInLab.objects.all().exists():
        sequencers = SequencerInLab.objects.all().order_by('sequencerName')
    else:
        return False
    return sequencers

def get_list_machines ():
    if Machines.objects.all().exists():
        machines = Machines.objects.all().order_by('machineName')
    else:
        return False
    return machines



def run ():
    check_migration_is_allow()
    run_objs = RunProcess.objects.all()
    number_runs = len(run_objs)
    print('Starting migration ', number_runs , ' runs were found' )
    counter = 0
    for run_obj in run_objs :
        machine_name = run_obj.get_run_sequencerModel()
        if machine_name != 'None':
            sequencer_obj = SequencerInLab.objects.get(sequencerName__exact = machine_name)
            run_obj.set_used_sequencer(sequencer_obj)
        counter +=1
        print( 'Updated run ', counter,  'of ' , number_runs , end='\r')
    print('\n Successfully completed')
    return
