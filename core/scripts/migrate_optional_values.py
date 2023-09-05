import core.models


"""
    The script is applicable for the upgrade from 2.3.0 to 3.0.0.
    Because the new version changes the way that optional values are stored.
    Instead of having a field "sample_project_option_list" now values are
    saved in SamplesProjectsTableOptions table. During the upgrade the existing
    field is still in the SampleProjectsFields, however this field will be
    deleted in the next release.
    The impacted tables is:
        - SampleProjectsFields
        - SamplesProjectsTableOptions
    Script reads the SampleProjectsFields table which containg
    sample_project_option_list and moves these values to SamplesProjectsTableOptions
"""


def run():
    s_project_field_objs = core.models.SampleProjectsFields.objects.all()
    for s_project_field_obj in s_project_field_objs:
        if s_project_field_obj.sample_project_field_type == "Options List":
            if (
                s_project_field_obj.sample_project_option_list == ""
                or s_project_field_obj.sample_project_option_list is None
            ):
                continue
            try:
                opt_values = s_project_field_obj.sample_project_option_list.split(",")
            except:
                import pdb

                pdb.set_trace()
            if len(opt_values) == 0:
                # field already migrated or contains no values
                continue
            # check if SamplesProjectsTableOptions exists for this parameter
            if core.models.SamplesProjectsTableOptions.objects.filter(
                sample_project_field=s_project_field_obj
            ).exists():
                # delete existing values
                s_p_opt_objs = core.models.SamplesProjectsTableOptions.objects.filter(
                    sample_project_field=s_project_field_obj
                )
                for s_p_opt_obj in s_p_opt_objs:
                    s_p_opt_obj.delete()
            for opt_val in opt_values:
                data = {"s_proj_obj": s_project_field_obj, "opt_value": opt_val}
                _ = core.models.SamplesProjectsTableOptions.objects.create_new_s_proj_table_opt(
                    data
                )
            # delete values in sample_project_option_list
            s_project_field_obj.sample_project_option_list = ""
            s_project_field_obj.save()
    print("Sample Projects optional fields migration was completed")

    return
