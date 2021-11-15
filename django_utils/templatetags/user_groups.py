from django.contrib.auth.models import Group
from django import template



register = template.Library()


@register.filter(name='has_group')
def has_group(user, group_name) :
    if Group.objects.filter(name=group_name).exists():
        groups = Group.objects.get(name=group_name) 
        return True if groups in user.groups.all() else False
    else:
        return False

