from iSkyLIMS_wetlab.models import RunProcess


def run ():
    
    runs_to_update = RunProcess.objects.all()
    n_to_update = len (runs_to_update)
    sucessful_counter = 0  
    unsucessful_runs = []
    for run in runs_to_update :
        state_of_run = run.runState
        if run.set_run_state(state_of_run):
            sucessful_counter += 1
        else:
            unsucessful_runs.append(run.get_run_name())
    if sucessful_counter != n_to_update:
        print('The following ', len(unsucessful_runs), ' runs were not updated:\n')
        for run in unsucessful_runs :
            print(run)
    else:   
        print( 'All runs are now updated ')
