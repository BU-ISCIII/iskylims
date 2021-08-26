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
            s_id = split_line[0]
            r_id = split_line[1]
            if SamplesInProject.objects.filter(pk__exact = s_id).exists():
                s_obj = SamplesInProject.objects.get(pk__exact = s_id)
            else:
                print('sample in project id ' , s_id , 'does not longer exists in database\n')
                invalids += 1
                continue
            if RunProcess.objects.filter(pk__exact = r_id).exists():
                r_obj = RunProcess.objects.get(pk__exact = r_id)
            else:
                print('Run id ' , r_id, 'does not longer exists in database\n')
                invalids += 1
                continue
            s_obj.runProcess_id = (r_obj)
            s_obj.save()

    with open ('part1.1_sampleInProject_migration_data.csv', 'r') as fh:
        for line in fh.readlines():
            split_line = line.split(',')
            if len(split_line )!= 2:
                print('invalid line :' , line , '\n')
                continue
            s_id = split_line[0]
            p_id = split_line[1]
            if SamplesInProject.objects.filter(pk__exact = s_id).exists():
                s_obj = SamplesInProject.objects.get(pk__exact = s_id)
            else:
                print('sample in project id ' , s_id , 'does not longer exists in database\n')
                invalids += 1
                continue
            if Projects.objects.filter(pk__exact = r_id).exists():
                p_obj = Projects.objects.get(pk__exact = p_id)
            else:
                print('Project id ' , p_id, 'does not longer exists in database\n')
                invalids += 1
                continue
            s_obj.project_id = (p_obj)
            s_obj.save()

    if invalids == 0:
        print( 'All samples have been updated to get relation with run\n')
    else:
        print('Unable to migrate all Samples to have runProcess relation ')

    invalids = 0

    with open ('part2_sampleInProject_migration_data.csv', 'r') as fh:
        for line in fh.readlines():
            split_line = line.split(',')
            if len(split_line )!= 2:
                print('invalid line :' , line , '\n')
                continue
            s_id = split_line[0]
            u_id = split_line[1]
            if SamplesInProject.objects.filter(pk__exact = s_id).exists():
                s_obj = SamplesInProject.objects.get(pk__exact = s_id)
            else:
                print('sample in project id ' , s_id , 'does not longer exists in database\n')
                invalids += 1
                continue
            if User.objects.filter(pk__exact = u_id).exists():
                u_obj = User.objects.get(pk__exact = u_id)
            else:
                print('User id ' , r_id, 'does not longer exists in database\n')
                invalids += 1
                continue
            s_obj.user_id = (u_obj)
            s_obj.save()

    if invalids == 0:
        print( 'All samples have been updated to get relation with User\n')
    else:
        print('Unable to migrate all Samples to have User relation ')
    return
