from datetime import datetime

from iSkyLIMS_core.models import (
    SampleProjects,
    SampleProjectsFields,
    Samples,
    PatientCore,
    LabRequest,
    SampleType,
    Species,
    StatesForSample,
    OntologyMap,
)
from iSkyLIMS_core.utils.handling_samples import increase_unique_value
from iSkyLIMS_core.core_config import HEADING_FOR_RECORD_SAMPLES

from iSkyLIMS_wetlab.models import SamplesInProject
# HEADING_FOR_RECORD_SAMPLES = ['Patient Code ID', 'Sample Name', 'Lab requested', 'Type of Sample', 'Species', 'Project/Service', 'Date sample reception', 'Collection Sample Date', 'Sample Storage Location', 'Only recorded']


def get_sample_fields(apps_name):
    sample_fields = {
        "Patient Code ID": {"field_name": "patientCore"},
        "Sample Name": {"field_name": "sampleName"},
        "Lab Requested": {"field_name": "labRequest"},
        "Type of Sample": {"field_name": "sampleType"},
        "Species": {"field_name": "species"},
        "Project/Service": {"field_name": "sampleProject"},
        "Date sample reception": {"field_name": "sampleEntryDate"},
        "Collection Sample Date": {"field_name": "collectionSampleDate"},
        "Sample Storage": {"field_name": "sampleLocation"},
        "Only recorded": {"field_name": "onlyRecorded"},
    }
    if SampleType.objects.filter(apps_name__exact=apps_name).exists():
        s_type_objs = SampleType.objects.filter(apps_name__exact=apps_name).order_by(
            "sampleType"
        )
        sample_fields["Type of Sample"]["options"] = []
        for s_type_obj in s_type_objs:
            sample_fields["Type of Sample"]["options"].append(s_type_obj.get_name())
    if Species.objects.filter(apps_name__exact=apps_name).exists():
        sample_fields["Species"]["options"] = []
        species_objs = Species.objects.filter(apps_name__exact=apps_name).order_by(
            "speciesName"
        )
        for species_obj in species_objs:
            sample_fields["Species"]["options"].append(species_obj.get_name())
    if LabRequest.objects.all().exists():
        sample_fields["Lab Requested"]["options"] = []
        lab_request_objs = LabRequest.objects.all().order_by("labName")
        for lab_request_obj in lab_request_objs:
            sample_fields["Lab Requested"]["options"].append(lab_request_obj.get_name())
    if SampleProjects.objects.filter(apps_name__exact=apps_name).exists():
        sample_fields["Project/Service"]["options"] = []
        s_proj_objs = SampleProjects.objects.filter(
            apps_name__exact=apps_name
        ).order_by("sampleProjectName")
        for s_proj_obj in s_proj_objs:
            sample_fields["Project/Service"]["options"].append(
                s_proj_obj.get_sample_project_name()
            )
    for key in sample_fields.keys():
        # import pdb; pdb.set_trace()
        if OntologyMap.objects.filter(
            label__iexact=sample_fields[key]["field_name"]
        ).exists():
            sample_fields[key]["ontology"] = (
                OntologyMap.objects.filter(
                    label__iexact=sample_fields[key]["field_name"]
                )
                .last()
                .get_ontology()
            )

    return sample_fields


def get_sample_project_obj(project_name):
    """Check if sampleProyect is defined in database"""
    if SampleProjects.objects.filter(sampleProjectName__exact=project_name).exists():
        return SampleProjects.objects.filter(
            sampleProjectName__exact=project_name
        ).last()
    return False


def get_project_fields_id_and_name(p_obj):
    """Fetch the fields defined for project"""
    fields = []
    if SampleProjectsFields.objects.filter(sampleProjects_id=p_obj).exists():
        p_field_objs = SampleProjectsFields.objects.filter(sampleProjects_id=p_obj)
        for p_field_obj in p_field_objs:
            fields.append([p_field_obj.get_field_id(), p_field_obj.get_field_name()])
    return fields


def get_patient_obj(patient):
    """Get the patient instance"""
    if PatientCore.objects.filter(patientCode__iexact=patient).exists():
        return PatientCore.objects.filter(patientCode__iexact=patient).last()
    return False


def split_sample_data(data):
    """Split the json data in 2 dictionaries, for having data to create the
    sample and data to create the project related info
    """
    split_data = {"s_data": {}, "p_data": []}
    sample_fields = [
        "patientCore",
        "sampleName",
        "labRequest",
        "sampleType",
        "species",
        "sampleProject",
        "sampleEntryDate",
        "collectionSampleDate",
        "sampleLocation",
        "onlyRecorded",
    ]

    for sample_field in sample_fields:
        try:
            if "date" in sample_field.lower():
                split_data["s_data"][sample_field] = datetime.strptime(
                    data[sample_field], "%Y-%m-%d"
                )
            else:
                split_data["s_data"][sample_field] = data[sample_field]
        except KeyError as e:
            return str(str(e) + " is not defined in your query")
    if split_data["s_data"]["onlyRecorded"] == "Yes":
        split_data["s_data"]["onlyRecorded"] = True
    else:
        split_data["s_data"]["onlyRecorded"] = False
    # fetch project fields
    project_obj = get_sample_project_obj(data["sampleProject"])
    if not project_obj:
        return "Project is not defined"
    project_fields = get_project_fields_id_and_name(project_obj)
    # check fields that are linked to other table
    if len(project_fields) > 0:
        for p_field in project_fields:
            try:
                split_data["p_data"].append(
                    {
                        "sampleProjecttField_id": p_field[0],

                        "sampleProjectFieldValue": data[p_field[1]],
                    }
                )
            except KeyError:
                # if not entry for the field set it to empty
                split_data["p_data"].append(
                    {
                        "sampleProjecttField_id": p_field[0],

                        "sampleProjectFieldValue": "",
                    }
                )
    return split_data


def include_instances_in_sample(data):

    """Fecth the patient instance"""
    if data["patientCore"] == "" or data["patientCore"].lower() == "null":
        data["patientCore"] = None
    else:
        try:
            data["patientCore"] = get_patient_obj(data["patientCore"]).get_patient_id()
        except AttributeError:
            return str(
                "patientCore " + data["patientCore"] + " is not defined in database"
            )
    if LabRequest.objects.filter(labName__exact=data["labRequest"]).exists():
        data["labRequest"] = (
            LabRequest.objects.filter(labName__exact=data["labRequest"]).last().get_id()
        )
    else:
        return str("labRequest " + data["labRequest"] + " is not defined in database")
    if SampleType.objects.filter(sampleType__exact=data["sampleType"]).exists():
        data["sampleType"] = (
            SampleType.objects.filter(sampleType__exact=data["sampleType"])
            .last()
            .get_sample_type_id()
        )
    else:
        return str("sampleType " + data["sampleType"] + " is not defined in database")
    if Species.objects.filter(speciesName__exact=data["species"]).exists():
        data["species"] = (
            Species.objects.filter(speciesName__exact=data["species"]).last().get_id()
        )
    else:
        return str("species " + data["species"] + " is not defined in database")
    if SampleProjects.objects.filter(
        sampleProjectName__exact=data["sampleProject"]
    ).exists():
        data["sampleProject"] = (
            SampleProjects.objects.filter(
                sampleProjectName__exact=data["sampleProject"]
            )
            .last()
            .get_id()
        )
    else:
        return str("sampleProject " + data["sampleProject"] + " is not defined")
    if data["onlyRecorded"]:
        data["sampleState"] = (
            StatesForSample.objects.filter(sampleStateName="Completed").last().get_id()
        )
        data["completedDate"] = datetime.now()
    else:
        data["sampleState"] = (
            StatesForSample.objects.filter(sampleStateName="Defined").last().get_id()
        )
        data["completedDate"] = None
    return data


def include_codding(user_name, sample):
    """Include Unique_id and Code_"""
    c_data = {}
    if not Samples.objects.exclude(uniqueSampleID__isnull=True).exists():
        c_data["uniqueSampleID"] = "AAA-0001"
    else:
        last_unique_value = (
            Samples.objects.exclude(uniqueSampleID__isnull=True).last().uniqueSampleID
        )
        c_data["uniqueSampleID"] = increase_unique_value(last_unique_value)
    c_data["sampleCodeID"] = str(user_name + "_" + sample)
    return c_data
