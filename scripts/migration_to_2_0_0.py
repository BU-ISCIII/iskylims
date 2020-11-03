from iSkyLIMS_drylab.models import Service
from iSkyLIMS_core.models import SequencingPlatform

def platforms_defined ():
    if SequencingPlatform.objects.all().exists():
        platform_dict = {}
        plat_objs = SequencingPlatform.objects.all()
        for plat_obj in plat_objs:
            platform_dict[plat_obj.platformName] = plat_obj
        return platform_dict
    else:
        return None

def run():
    invalid = 0
    platforms = platforms_defined ()
    if platforms != None:
        if SequencingPlatform.objects.all().exists():
            with open('platform_to_migrate.csv', 'r') as fh:
                for line in fh.readlines():
                    line = line.rstrip()
                    split_line = line.split(',')
                    if len(split_line )!= 2:
                        print('invalid line :' , line , '\n')
                        invalid += 1
                        continue
                    s_id = split_line[0]
                    p_name = split_line[1]
                    if 'Next-Seq' in p_name:
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

    return
