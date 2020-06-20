from iSkyLIMS_drylab.models import PreparationPipelineJobs



def check_pending_jobs():
    if PreparationPipelineJobs.objects.filter(jobState__jobStateName__exact = 'Queued').exists():
        return True
    return False

def handling_pending_jobs():
    pending_jobs = PreparationPipelineJobs.objects.filter(jobState__jobStateName__exact = 'Queued').order_by('generated_at')
    for pending_job in pending_jobs:
        pass
        check_preconditions
