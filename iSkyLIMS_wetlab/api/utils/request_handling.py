from iSkyLIMS_core.models import SampleProjects, SampleProjectsFields


def get_sample_project_obj(project_name):
    """Check if sampleProyect is defined in database"""
    if SampleProjects.objects.filter(sampleProjectName__exact=project_name).exists():
        return SampleProjects.objects.filter(sampleProjectName__exact=project_name).last()
    return False


def get_project_fields_objs(p_obj):
    """Fetch the fields defined for project """
    fields = []
    if SampleProjectsFields.objects.filter(sampleProjects_id=p_obj).exists():
        p_field_objs = SampleProjectsFields.objects.filter(sampleProjects_id=p_obj)
        for p_field_obj in p_field_objs:
            fields.append(p_field_obj)
    return fields


def split_sample_data(data):
    """Split the json data in 2 dictionary, for having data to create the
        sample an data to create the project related data
    """
    project_obj = get_sample_project_obj(data['project'])
    if not project_obj:
        return False
    split_data = {"s_data": {}, "p_data": {}}
    fields_in_sample = [
        "sampleState",
        "patientCore",
        "labRequest",
        "sampleType",
        "sampleUser",
        "sampleCodeID",
        "uniqueSampleID",
        "species",
        "sampleLocation",
        "sampleEntryDate",
        "uniqueSampleID",
        "sampleCodeID",
    ]
    p_fields = get_project_fields_objs(project_obj)
