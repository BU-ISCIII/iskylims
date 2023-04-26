from core.models import *


def get_defined_platforms_and_ids(sequencing_technology):
    """
    Functions:
        The function collect a lis wit all platform defined
    Input;
        sequencing_technology   # Technology for filtering
    Return:
        platform_list with a tuple list
    """
    platform_list = []
    if SequencingPlatform.objects.filter(
        sequencing_technology__exact=sequencing_technology
    ).exists():
        platform_objs = SequencingPlatform.objects.filter(
            sequencing_technology__exact=sequencing_technology
        ).order_by("platform_name")
        for platform_obj in platform_objs:
            platform_list.append(
                [platform_obj.get_platform_id(), platform_obj.get_platform_name()]
            )
    return platform_list
