from datetime import datetime

from core.models import (
    City,
    LabRequest,
    OntologyMap,
    PatientCore,
    SampleProjects,
    SampleProjectsFields,
    SampleProjectsFieldsValue,
    Samples,
    SampleType,
    Species,
    StateInCountry,
    StatesForSample,
)
from core.utils.samples import increase_unique_value
from wetlab.api.serializers import CreateSampleTypeSerializer
from wetlab.config import (
    ERROR_API_NO_SAMPLE_DEFINED,
    ERROR_API_NO_SAMPLE_PROJECT_DEFINED,
    ERROR_API_NO_SAMPLE_PROJECT_FIELD_DEFINED,
    ERROR_API_SAMPLE_STATE_VALUE_IS_NOT_DEFINED,
)


def create_state(state, apps_name):
    """Create state instance"""
    data = {"state": state, "apps_name": apps_name}
    return StateInCountry.objects.create_new_state(data)


def create_city(data, apps_name):
    """Create a City instance"""
    data["state"] = create_state(data["geo_loc_state"], apps_name).get_state_id()
    data["cityName"] = data["geo_loc_city"]
    data["latitude"] = data["geo_loc_latitude"]
    data["longitude"] = data["geo_loc_longitude"]
    data["apps_name"] = apps_name
    return City.objects.create_new_city(data)


def create_new_laboratory(lab_data):
    """Create new laboratory instance with the data collected in the request"""
    if City.objects.filter(city_name__exact=lab_data["geo_loc_city"]).exists():
        city_id = (
            City.objects.filter(city_name__exact=lab_data["geo_loc_city"])
            .last()
            .get_city_id()
        )
    else:
        if StateInCountry.objects.filter(
            state_name__exact=lab_data["geo_loc_state"]
        ).exists():
            lab_data["state"] = (
                StateInCountry.objects.filter(
                    state_name__exact=lab_data["geo_loc_state"]
                )
                .last()
                .get_state_id()
            )
        else:
            lab_data["state"] = create_state(
                lab_data["geo_loc_state"], lab_data["apps_name"]
            ).get_state_id()
        city_id = create_city(lab_data, lab_data["apps_name"]).get_city_id()
    lab_data["city"] = city_id
    return LabRequest.objects.create_lab_request(lab_data)


def create_new_sample_type(sample_type, apps_name):
    """Create new Sample type with optional fields Patient, sample storage
    location
    """
    data = {}
    data["optional_fields"] = "0,8"
    data["apps_name"] = apps_name
    data["sampleType"] = sample_type
    sample_type_serializers = CreateSampleTypeSerializer(data=data)
    if sample_type_serializers.is_valid():
        sample_type_obj = sample_type_serializers.save()
        return sample_type_obj
    return None


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
            "sample_type"
        )
        sample_fields["Type of Sample"]["options"] = []
        for s_type_obj in s_type_objs:
            sample_fields["Type of Sample"]["options"].append(s_type_obj.get_name())
    if Species.objects.filter(apps_name__exact=apps_name).exists():
        sample_fields["Species"]["options"] = []
        species_objs = Species.objects.filter(apps_name__exact=apps_name).order_by(
            "species_name"
        )
        for species_obj in species_objs:
            sample_fields["Species"]["options"].append(species_obj.get_name())
    if LabRequest.objects.all().exists():
        sample_fields["Lab Requested"]["options"] = []
        lab_request_objs = LabRequest.objects.all().order_by("lab_name")
        for lab_request_obj in lab_request_objs:
            sample_fields["Lab Requested"]["options"].append(lab_request_obj.get_name())
    if SampleProjects.objects.filter(apps_name__exact=apps_name).exists():
        sample_fields["Project/Service"]["options"] = []
        s_proj_objs = SampleProjects.objects.filter(
            apps_name__exact=apps_name
        ).order_by("sample_project_name")
        for s_proj_obj in s_proj_objs:
            sample_fields["Project/Service"]["options"].append(
                s_proj_obj.get_sample_project_name()
            )
    for key in sample_fields.keys():
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
    if SampleProjects.objects.filter(sample_project_name__exact=project_name).exists():
        return SampleProjects.objects.filter(
            sample_project_name__exact=project_name
        ).last()
    return False


def include_instances_in_sample(data, lab_data, apps_name):
    """Collect the instances before creating the sample instance
    If laboratory will be created if it is not defined
    """
    if LabRequest.objects.filter(lab_name__exact=data["labRequest"]).exists():
        data["labRequest"] = (
            LabRequest.objects.filter(lab_name__exact=data["labRequest"])
            .last()
            .get_id()
        )
    else:
        lab_data["apps_name"] = apps_name
        data["labRequest"] = create_new_laboratory(lab_data).get_id()
    if SampleType.objects.filter(sample_type__exact=data["sampleType"]).exists():
        data["sampleType"] = (
            SampleType.objects.filter(sample_type__exact=data["sampleType"])
            .last()
            .get_sample_type_id()
        )
    else:
        # create Sample type instance
        sample_type_obj = create_new_sample_type(data["sampleType"], apps_name)
        data["sampleType"] = sample_type_obj.get_sample_type_id()
        # return str("sampleType " + data["sampleType"] + " is not defined in database")
    if Species.objects.filter(species_name__exact=data["species"]).exists():
        data["species"] = (
            Species.objects.filter(species_name__exact=data["species"]).last().get_id()
        )
    else:
        return str("species " + data["species"] + " is not defined in database")
    if data["patientCore"] == "" or data["patientCore"].lower() == "null":
        data["patientCore"] = None
    else:
        try:
            data["patientCore"] = get_patient_obj(data["patientCore"]).get_patient_id()
        except AttributeError:
            return str(
                "patientCore " + data["patientCore"] + " is not defined in database"
            )
    if SampleProjects.objects.filter(
        sample_project_name__exact=data["sampleProject"]
    ).exists():
        data["sampleProject"] = (
            SampleProjects.objects.filter(
                sample_project_name__exact=data["sampleProject"]
            )
            .last()
            .get_id()
        )
    else:
        return str("sampleProject " + data["sampleProject"] + " is not defined")
    if data["onlyRecorded"]:
        data["sampleState"] = (
            StatesForSample.objects.filter(sample_state_name="Completed")
            .last()
            .get_id()
        )
        data["completedDate"] = datetime.now()
    else:
        data["sampleState"] = (
            StatesForSample.objects.filter(sample_state_name="Defined").last().get_id()
        )
        data["completedDate"] = None
    return data


def include_coding(user_name, sample):
    """Include Unique_id and Code_"""
    c_data = {}
    if not Samples.objects.exclude(unique_sample_id__isnull=True).exists():
        c_data["uniqueSampleID"] = "AAA-0001"
    else:
        last_unique_value = (
            Samples.objects.exclude(unique_sample_id__isnull=True)
            .last()
            .unique_sample_id
        )
        c_data["uniqueSampleID"] = increase_unique_value(last_unique_value)
    c_data["sampleCodeID"] = str(user_name + "_" + sample)
    return c_data


def get_project_fields_id_and_name(p_obj):
    """Fetch the fields defined for project"""
    fields = []
    if SampleProjectsFields.objects.filter(sample_projects_id=p_obj).exists():
        p_field_objs = SampleProjectsFields.objects.filter(sample_projects_id=p_obj)
        for p_field_obj in p_field_objs:
            fields.append([p_field_obj.get_field_id(), p_field_obj.get_field_name()])
    return fields


def get_patient_obj(patient):
    """Get the patient instance"""
    if PatientCore.objects.filter(patient_code__iexact=patient).exists():
        return PatientCore.objects.filter(patient_code__iexact=patient).last()
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
                try:
                    # Check if value is in date format with - separation
                    split_data["s_data"][sample_field] = datetime.strptime(
                        data[sample_field], "%Y-%m-%d"
                    )
                except ValueError:
                    try:
                        # Check if no separation is in date format
                        split_data["s_data"][sample_field] = datetime.strptime(
                            data[sample_field], "%Y%m%d"
                        )
                    except ValueError:
                        # Value is not a date. Set to None to allow that serialzer
                        # store it in database.
                        split_data["s_data"][sample_field] = None
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
    # fetched data to define new lab
    lab_data = {}
    lab_data["labName"] = data["labRequest"]
    lab_split = data["labRequest"].strip().split(" ")
    lab_code = ""
    for word in lab_split:
        lab_code += word[0]
    lab_data["labNameCoding"] = lab_code
    lab_data["labUnit"] = ""
    lab_data["labContactName"] = ""
    lab_data["labPhone"] = ""
    lab_data["labEmail"] = data["collecting_institution_email"]
    lab_data["address"] = data["collecting_institution_email"]
    lab_data["geo_loc_city"] = data["geo_loc_city"]
    lab_data["geo_loc_state"] = data["geo_loc_state"]
    lab_data["geo_loc_latitude"] = data["geo_loc_latitude"]
    lab_data["geo_loc_longitude"] = data["geo_loc_longitude"]

    split_data["lab_data"] = lab_data

    return split_data


def summarize_samples(data):
    summarize = {}
    sample_objs = Samples.objects.all()
    # Filter the samples to get the summary information
    if len(sample_objs) == 0:
        return {"ERROR": ERROR_API_NO_SAMPLE_DEFINED}
    if "sampleList" in data:
        sample_objs = sample_objs.filter(sample_name__in=data["sampleList"])
    if "sampleState" in data:
        if not StatesForSample.object.filter(
            sample_state_name__iexact=data["state"]
        ).exists():
            return {"ERROR": ERROR_API_SAMPLE_STATE_VALUE_IS_NOT_DEFINED}
        sample_objs = sample_objs.filter(
            sampleState__sample_state_name__iexact=data["state"]
        )
    if "startDate" in data and "endDate" in data:
        sample_objs = sample_objs.filter(
            generated_at__range=(data["startDate"], data["end_date"])
        )
    elif "startDate" in data:
        sample_objs = sample_objs.filter(generated_at__gte=data["startDate"])
    elif "endDate" in data:
        sample_objs.filter(generated_at__lte=data["endDate"])

    if "region" in data:
        sample_objs = sample_objs.filter(
            lab_request__lab_city__belongs_to_state__state_name__iexact=data["region"]
        )
    if "laboratory" in data:
        sample_objs = sample_objs.filter(
            lab_request__lab_name__iexact=data["laboratory"]
        )
    if "sample_project_name" in data:
        if not SampleProjects.objects.filter(
            sample_project_name__iexact=data["sample_project_name"]
        ).exists():
            return {"ERROR": ERROR_API_NO_SAMPLE_PROJECT_DEFINED}
        s_project_obj = SampleProjects.objects.filter(
            sample_project_name__iexact=data["sample_project_name"]
        ).last()
        s_project_field_objs = SampleProjectsFields.objects.filter(
            sample_projects_id=s_project_obj
        )
        # filter the samples for this  sample project
        sample_objs = sample_objs.filter(sample_project=s_project_obj)
        if "sample_project_field" in data:
            if not SampleProjectsFields.objects.filter(
                sample_projects_id=s_project_obj,
                sample_project_field_description__iexact=data["sample_project_field"],
            ).exists():
                return {"ERROR": ERROR_API_NO_SAMPLE_PROJECT_FIELD_DEFINED}

            s_project_field_objs = SampleProjectsFields.objects.filter(
                sample_projects_id=s_project_obj,
                sample_project_field_description__iexact=data["sample_project_field"],
            ).order_by("sample_project_field_classification_id")
    #
    # get the sumarize infomation for the selected samples
    #

    sample_list = list(sample_objs.values_list("sampleName", flat=True))
    summarize["region"] = {}
    summarize["laboratory"] = {}
    if "region" in data:
        regions = [data["region"]]
    else:
        regions = (
            StateInCountry.objects.all().values_list("stateName", flat=True).distinct()
        )
    for region in regions:
        summarize["region"][region] = sample_objs.filter(
            lab_request__lab_city__belongs_to_state__state_name__iexact=region
        ).count()

    if "laboratory" in data:
        laboratories = [data["laboratory"]]
        summarize["samples"] = sample_list
        summarize.pop("region", None)
    else:
        laboratories = (
            LabRequest.objects.all().values_list("labName", flat=True).distinct()
        )
    for laboratory in laboratories:
        summarize["laboratory"][laboratory] = sample_objs.filter(
            lab_request__lab_name__iexact=laboratory
        ).count()
    # Show only the parameters when the sample project name is given
    if "sample_project_name" in data:
        summarize["parameters"] = {}

        for s_project_field_obj in s_project_field_objs:
            p_name = s_project_field_obj.get_field_name()
            summarize["parameters"][p_name] = {}
            # summarize["parameters"][] =
            # check if sample Proyect fields has options
            # if true then get the used values and get their numbers
            if s_project_field_obj.get_field_type() == "Options List":
                if SampleProjectsFieldsValue.objects.filter(
                    sample_project_field_id=s_project_field_obj
                ).exists():
                    # get the unique values to iter over them

                    f_values = (
                        SampleProjectsFieldsValue.objects.filter(
                            sample_project_field_id=s_project_field_obj
                        )
                        .values_list("sample_project_field_value", flat=True)
                        .distinct()
                    )
                    for f_value in f_values:
                        summarize["parameters"][p_name][
                            f_value
                        ] = SampleProjectsFieldsValue.objects.filter(
                            sample_project_field_id=s_project_field_obj,
                            sample_project_field_value__exact=f_value,
                            sample_id__sample_name__in=sample_list,
                        ).count()
                else:
                    summarize["parameters"][p_name]["value"] = 0
            else:
                summarize["parameters"][p_name][
                    "value"
                ] = SampleProjectsFieldsValue.objects.filter(
                    sample_project_field_id=s_project_field_obj,
                    sample_id__sample_name__in=sample_list,
                ).count()
            summarize["parameters"][p_name][
                "classification"
            ] = s_project_field_obj.get_classification_name()

    summarize["samples_number"] = len(sample_list)
    return summarize


def collect_statistics_information(data):
    """Collect statistics for the fields utilization for the requested project"""

    if "sample_project_name" in data:
        if not SampleProjects.objects.filter(
            sample_project_name__iexact=data["sample_project_name"]
        ).exists():
            return ""
        s_project_obj = SampleProjects.objects.filter(
            sample_project_name__iexact=data["sample_project_name"]
        ).last()
        if "project_field" in data:
            query_params = data["project_field"].split(",")

            if len(query_params) > 2:
                return {"ERROR": ""}
            stats_data = {}
            par1_values = (
                SampleProjectsFieldsValue.objects.filter(
                    sample_project_field_id__sample_project_id=s_project_obj,
                    sample_project_field_id__sample_project_field_name__iexact=query_params[
                        0
                    ],
                )
                .values_list("sampleProjectFieldValue", flat=True)
                .distinct()
            )

            if len(query_params) == 2:
                for par1_val in par1_values:
                    stats_data[par1_val] = {}

                    samples = SampleProjectsFieldsValue.objects.filter(
                        sample_project_field_id__sample_project_id=s_project_obj,
                        sample_project_field_id__sample_project_field_name__iexact=query_params[
                            0
                        ],
                        sample_project_field_value__exact=par1_val,
                    ).values_list("sample_id", flat=True)
                    par2_values = (
                        SampleProjectsFieldsValue.objects.filter(
                            sample_id__in=samples,
                            sample_project_field_id__sample_project_field_name__iexact=query_params[
                                1
                            ],
                        )
                        .values_list("sample_project_field_value", flat=True)
                        .distinct()
                    )
                    for par2_val in par2_values:
                        value = SampleProjectsFieldsValue.objects.filter(
                            sample_id__in=samples,
                            sample_project_field_id__sample_project_field_name=query_params[
                                1
                            ],
                            sample_project_field_value__exact=par2_val,
                        ).count()
                        if value > 0:
                            stats_data[par1_val][par2_val] = value
            else:
                for par1_val in par1_values:
                    stats_data[par1_val] = SampleProjectsFieldsValue.objects.filter(
                        sample_project_field_id__sample_projects_id=s_project_obj,
                        sample_project_field_id__sample_project_field_name__iexact=query_params[
                            0
                        ],
                        sample_project_field_value=par1_val,
                    ).count()

            return stats_data
        else:  # Collect info stats for all fields
            # Collect the fields utilization for sample projects
            stats_data = {
                "always_none": [],
                "fields_norm": {},
                "never_used": [],
                "fields_value": {},
            }

            num_samples = Samples.objects.filter(sample_project=s_project_obj).count()
            s_project_field_objs = SampleProjectsFields.objects.filter(
                sample_project_id=s_project_obj
            )
            for s_project_field_obj in s_project_field_objs:
                f_name = s_project_field_obj.get_field_name()
                if not SampleProjectsFieldsValue.objects.filter(
                    sample_project_field_id=s_project_field_obj
                ).exists():
                    stats_data["never_used"].append(f_name)
                    stats_data["fields_value"][f_name] = 0
                    continue
                count_not_none = (
                    SampleProjectsFieldsValue.objects.filter(
                        sample_project_field_id=s_project_field_obj
                    )
                    .exclude(sample_project_field_value__in=["None", ""])
                    .count()
                )
                stats_data["fields_value"][f_name] = count_not_none
                if count_not_none == 0:
                    stats_data["always_none"].append(f_name)
                    continue
                try:
                    stats_data["fields_norm"][f_name] = count_not_none / num_samples
                except ZeroDivisionError:
                    stats_data["fields_norm"][f_name] = 0
            stats_data["num_samples"] = num_samples
            return stats_data

    return {"ERROR": ""}
