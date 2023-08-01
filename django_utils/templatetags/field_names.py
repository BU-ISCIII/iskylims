
from django import template

register = template.Library()

@register.filter 
def get_verbose_name(object, field_name): 
    return object._meta.get_field(field_name).verbose_name


@register.filter
def get_name(object, field_name): 
    return object._meta.get_field(field_name).name