from django.conf import settings

def get_installed_apps () :
    app_prefix = 'iSkyLIMS'
    core = 'core'
    apps_list = []
    #import pdb; pdb.set_trace()
    apps = list(apps for apps in settings.INSTALLED_APPS if (app_prefix in apps and not core in apps))
    for app in apps :
        apps_list.append([app,settings.APPS_NAMES[app]])
    return apps_list
     
