from core.models import LabRequest


def get_laboratory_instance(lab_name):
    """return the laboratory instance, or None if not match"""
    if LabRequest.objects.filter(lab_name__iexact=lab_name).exists():
        return LabRequest.objects.filter(lab_name__iexact=lab_name).last()
    return None
