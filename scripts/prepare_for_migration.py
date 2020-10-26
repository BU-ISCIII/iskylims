from iSkyLIMS_drylab.models import Service



def run():
    ## migrate the servicefile that it was in the service request for allowing
    ## to add multiple files.
    if Service.objects.all().exclude(serviceFile = None).exists():

        services = Service.objects.all().exclude(serviceFile = None)
        with open('service_files_to_migrate.csv', 'w') as fh:
            for service in services:
                fh.write(service.pk, + ',' + service.serviceFile)
        print('Successfully file cretion for migration')
    else:
        print('There is no service file. Do not need to perform a migration')
    return
