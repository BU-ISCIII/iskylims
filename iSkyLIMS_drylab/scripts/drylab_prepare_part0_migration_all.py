from iSkyLIMS_drylab.models import Service



def run():
    ## migrate the servicefile that it was in the service request for allowing
    ## to add multiple files.

    if Service.objects.all().exists():
        services = Service.objects.all()
        ## dump platform
        with open('drylab_part1_platform_to_migrate.csv', 'w') as fh:
            for service in services:
                if service.servicePlatform != None :  # for doing on master branch
                #if service.serviceSequencingPlatform != None : # for doing on develop
                    fh.write(str(str(service.pk) + ',' + service.servicePlatform.platformName + '\n'))
        print('Successfully file creation for platform migration')
        # dump services states for apliying them to resolution
        with open('drylab_part2_resolution_state_to_migrate.csv', 'w') as fh:
            for item in services:
                fh.write(str(str(item.pk) + ',' + item.serviceStatus + '\n'))
        print('Successfully file creation for resolution state migration')
        # dump projects used in services
        with open('drylab_part4_project_used_in_services_to_migrate.csv', 'w') as fh:
            for item in services:
                project_objs = item.serviceProjectNames.all()
                if len(project_objs) > 0:
                    for project_obj in project_objs:
                        fh.write(str(str(item.pk) + ',' + str(project_obj.pk) + '\n'))
        print('Successfully file creation for project in services migration')

    ## migrate the servicefile that it was in the service request for allowing
    ## to add multiple files.

    if Service.objects.all().exclude(serviceFile = None).exists():
        services = Service.objects.all().exclude(serviceFile = None)
        with open('drylab_part3_service_files_to_migrate.csv', 'w') as fh:
            for service in services:
                if str(service.serviceFile) != '':
                    fh.write(str(str(service.pk) + ',' + str(service.serviceFile) + '\n'))
        print('Successfully file creation for file migration')
    else:
        print('There is no service file. Do not need to perform a migration')
    return

    ##  migrate the projects used in the services
    ##
