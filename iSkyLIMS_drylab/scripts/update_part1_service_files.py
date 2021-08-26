from iSkyLIMS_drylab.models import Service, Resolution, UploadServiceFile



def run():
    map_state = {'recorded':'Recorded','queued':'Recorded', 'in_progress':'In Progress',  'delivered':'Delivery', 'archived':'Delivery'}
    with open('resolution_state_to_migrate.csv', 'r') as fh:
        for line in fh.readlines():
            line = line.rstrip()
            split_line = line.split(',')
            if len(split_line )!= 2:
                print('invalid line :' , line , '\n')
                continue
            if split_line[1] == 'recorded':
                continue
            s_id = split_line[0]
            if Resolution.objects.filter(resolutionServiceID__pk__exact = s_id).exists():
                resolution_objs = Resolution.objects.filter(resolutionServiceID__pk__exact = s_id)
                for resolution_obj in resolution_objs:
                    resolution_obj.update_resolution_state(map_state[split_line[1]])
        print('Successfully restore resolution state migration')
    invalids = 0
    with open('service_files_to_migrate.csv', 'r') as fh:
        for line in fh.readlines():
            line = line.rstrip()
            split_line = line.split(',')
            if len(split_line )!= 2:
                print('invalid line :' , line , '\n')
                continue
            s_id = split_line[0]
            if Service.objects.filter(pk__exact = s_id).exists():
                service_obj  = Service.objects.get(pk__exact = s_id)
                upload_file_obj = UploadServiceFile.objects.create(uploadService = service_obj, uploadFile = split_line[1])
            else:
                print('Service id', s_id ,' Does not longer exists\n')
                invalids += 1
        if invalids == 0:
            print ('Sucessfully for part0 upload service files migration\n')
        else:
            print('Unable to migrate all service files')
