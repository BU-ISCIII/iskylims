from iSkyLIMS_wetlab.models import SamplesInProject, RunProcess, Projects
from django.contrib.auth.models import User

def run ():
    invalids = 0

    with open ('part1_sampleInProject_migration_data.csv', 'r') as fh:
        for line in fh.readlines():
            split_line = line.split(',')
            if len(split_line )!= 2:
                print('invalid line :' , line , '\n')
                continue
            p_id = split_line[0]
            r_id = split_line[1]
            if Projects.objects.filter(pk__exact = p_id).exists():
                p_obj = Projects.objects.get(pk__exact = p_id)
            else:
                print('Project id ' , p_id , 'does not longer exists in database\n')
                invalids += 1
                continue
            if RunProcess.objects.filter(pk__exact = r_id).exists():
                r_obj = RunProcess.objects.get(pk__exact = r_id)
            else:
                print('Run id ' , r_id, 'does not longer exists in database\n')
                invalids += 1
                continue
            p_obj.projects.add(r_obj)


    if invalids == 0:
        print( 'All projects have been updated to get many to many relation with run\n')
    else:
        print('Unable to migrate all projects to have many to many relation ')
