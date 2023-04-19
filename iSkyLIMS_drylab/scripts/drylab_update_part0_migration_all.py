from iSkyLIMS_drylab.models import Service, Resolution, UploadServiceFile, RequestedSamplesInServices
from iSkyLIMS_drylab.drylab_config import *
from iSkyLIMS_core.models import SequencingPlatform
from iSkyLIMS_wetlab.models import Projects, SamplesInProject, RunProcess, RunningParameters

from django.conf import settings
import os


def platforms_defined ():
    if SequencingPlatform.objects.all().exists():
        platform_dict = {}
        plat_objs = SequencingPlatform.objects.all()
        for plat_obj in plat_objs:
            platform_dict[plat_obj.platformName] = plat_obj
        return platform_dict
    else:
        return None

def get_samples_info_for_migration(project_id):
    '''
    This function get the sample name and sample id for the samples in the project

    '''
    if Projects.objects.filter(pk__exact = project_id).exists():
        if SamplesInProject.objects.filter(project_id__exact = project_id).exists():
            sample_objs = SamplesInProject.objects.filter(project_id__exact = project_id)
            return [[sample_obj.sample_name, sample_obj.pk] for sample_obj in sample_objs]
    return None

def get_run_information(project_id):
    '''
    This functions use the  SamplesInProject to get the run information.
    Not checking for project it assume that exists
    '''
    sample_objs = SamplesInProject.objects.filter(project_id__exact = project_id)
    # import pdb; pdb.set_trace()
    if len(sample_objs) >0 :
        path = RunningParameters.objects.get(runName_id = sample_objs[0].runProcess_id).RunID
        return sample_objs[0].runProcess_id.runName, sample_objs[0].runProcess_id.pk, path
    else:
        return None , None, None

def run():

    ### restore platform
    invalid = 0
    platforms = platforms_defined ()
    if platforms != None:
        if SequencingPlatform.objects.all().exists():
            with open('drylab_part1_platform_to_migrate.csv', 'r') as fh:
                for line in fh.readlines():
                    line = line.rstrip()
                    split_line = line.split(',')
                    if len(split_line )!= 2:
                        print('invalid line :' , line , '\n')
                        invalid += 1
                        continue
                    s_id = split_line[0]
                    platform_name = split_line[1]
                    if 'Next-Seq' in platform_name:
                        p_obj = platforms['NextSeq 550 Series']
                    else:
                        p_obj = platforms['MiSeq']
                    if Service.objects.filter(pk__exact = s_id).exists():
                        s_obj = Service.objects.get(pk__exact = s_id)
                        s_obj.update_sequencing_platform(p_obj)
            if invalid == 0:
                print('Successfully file creation for platform migration\n')
            else:
                print('migration perform with errors\n')
        else:
            print('No sequencingPlatform migration is needed\n')
    else:
        print('SequencingPlatform has not data. Unable to execute the migration\n')


    ## restore service state to apply to resolution state
    map_state = {'recorded':'Recorded', 'queued':'Recorded', 'rejected':'Rejected', 'in_progress':'In Progress',  'delivered':'Delivery', 'archived':'Delivery'}
    with open('drylab_part2_resolution_state_to_migrate.csv', 'r') as fh:
        for line in fh.readlines():
            line = line.rstrip()
            split_line = line.split(',')
            if len(split_line )!= 2:
                print('invalid line :' , line , '\n')
                continue
            if split_line[1] == 'recorded':
                continue
            s_id = split_line[0]
            if Resolution.objects.filter(resolution_serviceID__pk__exact = s_id).exists():
                resolution_objs = Resolution.objects.filter(resolution_serviceID__pk__exact = s_id)
                for resolution_obj in resolution_objs:
                    resolution_obj.update_resolution_state(map_state[split_line[1]])
        print('Successfully restore resolution state migration')

    ## restore files for multiple upload files
    invalids = 0
    folder_files = os.path.join(settings.MEDIA_ROOT, USER_REQUESTED_SERVICE_FILE_DIRECTORY)
    if not os.path.isdir(folder_files):
        os.makedirs(folder_files)
    with open('drylab_part3_service_files_to_migrate.csv', 'r') as fh:
        for line in fh.readlines():
            line = line.rstrip()
            split_line = line.split(',')
            if len(split_line )!= 2:
                print('invalid line :' , line , '\n')
                continue
            s_id = split_line[0]
            if Service.objects.filter(pk__exact = s_id).exists():
                service_obj  = Service.objects.get(pk__exact = s_id)
                data = {'file':split_line[1], 'file_name': os.path.basename(split_line[1]) }
                #import pdb; pdb.set_trace()
                upload_file_obj = UploadServiceFile.objects.create_upload_file(data)
                upload_file_obj.update_service_id(service_obj)
            else:
                print('Service id', s_id ,' Does not longer exists\n')
                invalids += 1
        if invalids == 0:
            print ('Sucessfully for upload service files migration\n')
        else:
            print('Unable to migrate all service files')

    ### collect the wetlab proyect used in the services

    invalids = 0
    with open('drylab_part4_project_used_in_services_to_migrate.csv', 'r') as fh:
        for line in fh.readlines():
            line = line.rstrip()
            split_line = line.split(',')
            if len(split_line )!= 2:
                print('invalid line :' , line , '\n')
                continue
            serv_id, proj_id = split_line
            if Service.objects.filter(pk__exact = serv_id).exists():
                serv_obj = Service.objects.get(pk__exact = serv_id)
            else:
                print('Service id', serv_id ,' Does not longer exists\n')
                invalids += 1
                continue
            if not Projects.objects.filter(pk__exact = proj_id).exists():
                print('Project id', proj_id ,' Does not longer exists\n')
                invalids += 1
                continue
            data = {}
            project_obj = Projects.objects.get(pk__exact = proj_id)
            data['project_name'] = project_obj.projectName
            data['project_id'] = proj_id
            data['samples_in_service'] = serv_obj
            data['only_recorded'] = False
            data['run_name'], data['run_id'] , data['sample_path'] = get_run_information(proj_id)
            samples_data = get_samples_info_for_migration(proj_id)
            if not samples_data:
                print('Project id', proj_id ,' Does not have any sample\n')
                invalids += 1
                continue
            for sample in samples_data:
                data['sample_name'] , data['sample_id'] = sample
                RequestedSamplesInServices.objects.create_request_sample(data)
    if invalids == 0:
        print ('Sucessfully for samples in service migration\n')
    else:
        print('Unable to migrate all service files')

    return
