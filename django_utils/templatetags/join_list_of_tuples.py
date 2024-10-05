from django import template


register = template.Library()

@register.filter
def join_tuples(value):
    text = ""
    for item in value:
        text += item[0] + "," + str(item[1]) + ";;"
    return text