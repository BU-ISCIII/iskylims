import core.models


def get_laboratory_instance(lab_name):
    """return the laboratory instance, or None if not match"""
    if core.models.LabRequest.objects.filter(lab_name__iexact=lab_name).exists():
        return core.models.LabRequest.objects.filter(lab_name__iexact=lab_name).last()
    elif core.models.LabRequest.objects.filter(
        lab_name_coding__iexact=lab_name
    ).exists():
        return core.models.LabRequest.objects.filter(
            lab_name_coding__iexact=lab_name
        ).last()
    return None
