from django import template

register = template.Library()


@register.filter
def replace_spaces_with_underscores(value):
    return value.replace(" ", "_")
