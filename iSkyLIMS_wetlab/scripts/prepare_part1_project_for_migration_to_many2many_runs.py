from iSkyLIMS_wetlab.models import Projects

'''
    This file is part of the migration preparation for the version 2.0.0.
    It collect the information for the projects table for having the
    functionality that Project name can be repeated in different runs
'''
def run ():

    projects_to_update = Projects.objects.all()
    p_to_update = len (projects_to_update)
    if p_to_update == 0 :
        print('There is no projects defined in iSkyLIMS\n')
        print('No upload project script is required to execute')
        return
    with open ('project_migration_data.csv', 'w') as fh:
        for project in projects_to_update :
            fh.write(str(project.pk) + ',' + str(project.runprocess_id.pk) + '\n')

    print( 'All projects information is collected\n')
    return
