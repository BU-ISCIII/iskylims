from django.conf import settings
from django.contrib.auth.models import User

def get_installed_apps () :
    app_prefix = 'iSkyLIMS'
    core = 'core'
    apps_list = []
    #import pdb; pdb.set_trace()
    apps = list(apps for apps in settings.INSTALLED_APPS if (app_prefix in apps and not core in apps))
    for app in apps :
        apps_list.append([app,settings.APPS_NAMES[app]])
    return apps_list

def get_friend_list(user_name):
    friend_list = []
    user_groups = user_name.groups.values_list('name',flat=True)
    if len (user_groups) > 0 :
        for user in user_groups :

            if User.objects.filter(username__exact = user).exists():
                # friend_list.append(User.objects.get(username__exact = user).id)
                friend_list.append(User.objects.get(username__exact = user))

    friend_list.append(user_name)
    return friend_list
