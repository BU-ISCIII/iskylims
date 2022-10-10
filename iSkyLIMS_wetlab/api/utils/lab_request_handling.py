from iSkyLIMS_core.models import LabRequest


def get_laboratory_instance(lab_name):
    """return the laboratory instance, or None if not match"""
    if LabRequest.objects.filter(labName__iexact=lab_name).exists():
        return LabRequest.objects.filter(labName__iexact=lab_name).last()
    return None
