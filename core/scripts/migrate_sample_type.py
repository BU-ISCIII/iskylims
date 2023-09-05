import core.models


"""
    The script is applicable for the upgrade from 2.3.0 to 3.0.0.
    Because the new version changes the value that is stored now is the field
    name and not the number and instead of optional values now are the 
    mandatory fields-
    The impacted tables is:
        - SampleType
    Script reads the table which containg the optional index value and replace
    them for the mandatory name fields.
"""


def run():
    available_fields = [
        "patient_core",
        "sample_name",
        "lab_request",
        "sample_type",
        "species",
        "sample_project",
        "sample_entry_date",
        "collection_sample_date",
        "sample_location",
        "only_recorded" 
    ]
    sample_type_objs = core.models.SampleType.objects.all()
    for sample_type_obj in sample_type_objs:
        s_type_opt_field = sample_type_obj.get_mandatory_values()
        mandatory_fields = []
        try:
            _ = int(s_type_opt_field[0])
        except ValueError:
            # Fields have already in the new format
            continue
        for field in available_fields:
            if str((available_fields.index(field) + 1)) not in s_type_opt_field:
                mandatory_fields.append(field)
        m_fields = ",".join(mandatory_fields)
        sample_type_obj.update_mandatory_fields(m_fields)

    print("sample type migration was completed")

    return