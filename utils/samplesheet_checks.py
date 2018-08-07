## Module containing checks performed on a SampleSheet file for Illumina sequencers

def check_run_name_already_used(run_name, illumina_sequencer):
    ## Function checks whether run_name is already used in the database.
    ## error page is showed if run_name is already  defined
    ##

    #import pdb; pdb.set_trace()
    check_run_name_result='ERROR when executing: check_run_name_result()'
    if (RunProcess.objects.filter(runName = run_name)).exists():
        if RunProcess.objects.filter(runName = run_name, runState__exact ='Pre-Recorded'):
            ## Delete (for NextSeq, the sample sheet file and) the row in database
            delete_run = RunProcess.objects.filter(runName = run_name, runState__exact ='Pre-Recorded')
            delete_run[0].delete()
            if  (illumina_sequencer == "NextSeq"):
                sample_sheet_file = str(delete_run[0].sampleSheet)
                #import pdb; pdb.set_trace()
                full_path_sample_sheet_file = os.path.join(settings.MEDIA_ROOT, sample_sheet_file)
                os.remove(full_path_sample_sheet_file)
            check_run_name_result="OK"
        else:#runState != 'Pre-Recorded'
            check_run_name_result='Run Name is already used. ',
                                  'Run Name must be unique in database.', ' ',
                                  'ADVICE:',
                                  'Change the value in the Sample Sheet  file '

    else: # run_name is new
        check_run_name_result="OK"
    return check_run_name_result




## Check that there are indeed projects within the run_name
### TODO project_list=get_projects_in_run(stored_file)
#
#        if len (project_list) == 0 :
#        .....




def check_existence_run_users(project_list):

    check_existence_run_users_result='ERROR when executing: check_existence_run_usersxx()'

    ...

    return check_existence_run_users_result






def check_existence_run_projects(project_list)
    ## Check that the projects are already defined in the database

    check_existence_run_projects_result='ERROR when executing: check_existence_run_projects()'

    ...

    return check_existence_run_projects_result
