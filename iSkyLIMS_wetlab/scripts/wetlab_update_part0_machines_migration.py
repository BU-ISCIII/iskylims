from iSkyLIMS_wetlab.models import RunProcess
from iSkyLIMS_core.models import SequencerInLab

'''
The script migration restore the information collected in the prepare_part0_machines_migration .
Fetching the run_id and the sequencer name.
The csv file contains "run_id, sequencer_name"
'''

def get_list_defined_sequencers ():
    if SequencerInLab.objects.all().exists():
        sequencers = {}
        sequencer_objs = SequencerInLab.objects.all().order_by('sequencerName')
        for sequencer_obj in sequencer_objs:
            sequencers[sequencer_obj.get_sequencer_name()] = sequencer_obj
    else:
        return False
    return sequencers

def run ():
    sequencers = get_list_defined_sequencers()
    if not sequencers :
        print('Unable to restore machine migration\n', 'Define first the machine on SequencerInLab table\n')
        exit(1)
    invalids = 0
    with open ('wetlab_part_0_machines_migration.csv', 'r') as fh:
        for line in fh.readlines():
            line = line.rstrip()
            split_line = line.split(',')
            if len(split_line )!= 2:
                print('invalid line :' , line , '\n')
                continue
            r_id = split_line[0]
            machine_name = split_line[1]
            if RunProcess.objects.filter(pk__exact = r_id).exists():
                run_obj = RunProcess.objects.get(pk__exact = r_id)
                if machine_name in sequencers:
                    run_obj.set_used_sequencer(sequencers[machine_name])
                else:
                    print('Run id ' , r_id ,'Sequencer ', machine_name, ' is not defined on SequencerInLab table\n')
                    invalids +=1
                    continue
            else:
                print('Run id ' , r_id , 'does not longer exists in database\n')
                invalids +=1
                continue
    if invalids == 0:
        print ('Sucessfully for part0 machine migration\n')
    else:
        print('Unable to migrate all machines')
