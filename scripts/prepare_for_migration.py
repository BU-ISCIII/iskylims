from iSkyLIMS_drylab.models import Service



def run():
    ## migrate the servicefile that it was in the service request for allowing
    ## to add multiple files.

    if Service.objects.all().exists():
        services = Service.objects.all()
        with open('platform_to_migrate.csv', 'w') as fh:
            for service in services:
                if service.servicePlatform != None :
                    fh.write(str(str(service.pk) + ',' + service.servicePlatform.platformName + '\n'))
        print('Successfully file creation for platform migration')

        with open('resolution_state_to_migrate.csv', 'w') as fh:
            for item in services:
                fh.write(str(str(item.pk) + ',' + item.serviceStatus + '\n'))
        print('Successfully file creation for resolution state migration')
        
    if Service.objects.all().exclude(serviceFile = None).exists():
        services = Service.objects.all().exclude(serviceFile = None)
        with open('service_files_to_migrate.csv', 'w') as fh:
            for service in services:
                if str(service.serviceFile) != '':
                    fh.write(str(str(service.pk) + ',' + str(service.serviceFile) + '\n'))
        print('Successfully file creation for file migration')
    else:
        print('There is no service file. Do not need to perform a migration')
    return
