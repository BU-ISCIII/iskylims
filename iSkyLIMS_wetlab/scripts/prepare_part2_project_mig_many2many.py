from iSkyLIMS_wetlab.models import SamplesInProject
'''
    This file is part of the migration preparation for the version 2.0.0.
    It collect the information for the sampleInProject table for having
    the functionality that Project name can be repeated in different runs
    Two csv files are generated:
    - part2_Samples_linked_to_Runs_migration_data.csv having smaple_id, run_id
    - part2_Samples_linked_to_Users_migration_data.csv having sample_id, user_id
'''

def run ():

    samples_to_update = SamplesInProject.objects.all()
    s_to_update = len (samples_to_update)
    if s_to_update == 0 :
        print('There is no SamplesInProjects defined in iSkyLIMS\n')
        print('No upload SamplesInProject script is required to execute')
        return
    with open ('part2_SamplesID_RunsID_migration.csv', 'w') as fh:
        for sample in samples_to_update :
            project_obj = sample.project_id
            run_obj =  project_obj.runprocess_id
            fh.write(str(sample.pk) + ',' + str(run_obj.pk) + '\n')

    with open ('part2_SamplesID_UsersID_migration.csv', 'w') as fh:
        for sample in samples_to_update :
            project_obj = sample.project_id
            user_id =  project_obj.user_id.pk
            fh.write(str(sample.pk) + ',' + str(user_id) + '\n')
    print( 'All sample information is collected\n')


    return
