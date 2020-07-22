from iSkyLIMS_drylab.models import PreparationPipelineJobs



def check_pending_jobs():
    '''
    Description:
        The function checks if there are jobs in Queued state to start handling
    Return:
        True if there are jobs
    '''
    if PreparationPipelineJobs.objects.filter(jobState__jobStateName__exact = 'Queued').exists():
        return True
    return False

def check_preconditions(preparation_pipeline_obj):
    '''
    Description:
        The function check if there are requirements to fulfil before starting
    Return:
        True if preconditions are fulfilled
    '''
    if preparation_pipeline_obj.get_used_run_folder() == 'True':

    return True

def handling_pending_jobs():
    '''
    Description:
        The function collect the available services and the fields to present information
        in the form
    Return:
        children_services
    '''
    pending_jobs = PreparationPipelineJobs.objects.filter(jobState__jobStateName__exact = 'Queued').order_by('generated_at')
    for pending_job in pending_jobs:
        if not check_preconditions(pending_job):
            continue
