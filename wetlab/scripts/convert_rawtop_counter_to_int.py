import wetlab.models


def run():
    """The script implemted the issue #158 Unable to convert barcode count to
    integer, to convert the counter that are in a string format, separated
    by ","  to int.
    """
    top_bar_objs = wetlab.models.RawTopUnknowBarcodes.objects.all()
    for top_bar_obj in top_bar_objs:
        try:
            top_bar_obj.count = int(top_bar_obj.count.replace(",", ""))
            top_bar_obj.save()
        except AttributeError:
            continue
    return
