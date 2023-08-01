
from django import template

register = template.Library()

@register.simple_tag 
def get_verbose_name(object): 
    return object._meta.verbose_name


@register.simple_tag 
def get_name(object): 
    return object._meta.name