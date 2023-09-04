import re
from django import template

register = template.Library()


@register.filter
def get(mapping, key):
    return mapping.get(key, "")


@register.filter(name="clean_string")
def clean_string(value):
    # Use regular expression to remove spaces, "-" and "/"
    cleaned_value = re.sub(r"[-/.,;:\s]", "", value)
    return cleaned_value
