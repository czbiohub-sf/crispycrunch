# TODO (gdingle): consider switch to Jinja templating
# TODO (gdingle): consider eval tag: https://djangosnippets.org/snippets/1820/

from django import template

register = template.Library()


@register.simple_tag
def return_item(l, row, col=''):
    """Makes it possible to return an item in a dict by variables. Django templates suck."""
    try:
        return l[str(row) + str(col)]
    except Exception:
        return None
