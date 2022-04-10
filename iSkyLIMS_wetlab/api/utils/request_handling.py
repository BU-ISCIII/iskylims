from iSkyLIMS_core.models import SampleProjects, SampleProjectsFields, Samples
from iSkyLIMS_core.utils.handling_samples import increase_unique_value


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


def include_additional_sample_data(s_data):
    """Add required coding data for sample"""
    if not Samples.objects.exclude(uniqueSampleID__isnull=True).exists():
        s_data['uniqueSampleID'] = 'AAA-0001'
    else:
        last_unique_value = Samples.objects.exclude(uniqueSampleID__isnull=True).last().uniqueSampleID
        s_data['uniqueSampleID'] = increase_unique_value(last_unique_value)


def split_sample_data(data):
    """Split the json data in 2 dictionary, for having data to create the
        sample an data to create the project related data
    """
    project_obj = get_sample_project_obj(data['project'])
    if not project_obj:
        return False
    split_data = {"s_data": {}, "p_data": {}}
    fields_in_sample = [
        "patientCore",
        "sampleName",
        "labRequest",
        "sampleType",
        "species",
        "collectionSampleDate",
        "sampleEntryDate",
        "sampleLocation",
        "onlyRecorded"
    ]
    p_fields = get_project_fields_objs(project_obj)
    "sampleCodeID",
    "uniqueSampleID",
