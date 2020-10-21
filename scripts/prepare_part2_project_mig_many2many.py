from iSkyLIMS_wetlab.models import SamplesInProject
'''
    This file is part of the migration preparation for the version 2.0.0.
    It collect the information for the sampleInProject table for having
    the funcitonality that Project name can be repeated in different runs
'''

def run ():

    samples_to_update = SamplesInProject.objects.all()
    s_to_update = len (samples_to_update)
    if s_to_update == 0 :
        print('There is no SamplesInProjects defined in iSkyLIMS\n')
        print('No upload SamplesInProject script is required to execute')
        return
    with open ('part1_sampleInProject_migration_data.csv', 'w') as fh:
        for sample in samples_to_update :
            project_obj = sample.project_id
            run_objs =  project_obj.runProcess.all()
            fh.write(str(sample.pk) + ',' + str(run_objs[0].pk) + '\n')

    with open ('part2_sampleInProject_migration_data.csv', 'w') as fh:
        for sample in samples_to_update :
            project_obj = sample.project_id
            user_id =  project_obj.user_id.pk
            fh.write(str(sample.pk) + ',' + str(user_id) + '\n')
    print( 'All sample information is collected\n')


    return
