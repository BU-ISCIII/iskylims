from django.conf import settings

def get_installed_apps () :
    app_prefix = 'iSkyLIMS'
    core = 'core'
    return list(apps for apps in settings.INSTALLED_APPS if (app_prefix in apps and not core in apps))
     
