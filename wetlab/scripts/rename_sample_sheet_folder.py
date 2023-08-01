import wetlab.models


"""
    The script is applicable for the upgrade from 2.3.0 to 3.0.0.
    Because the folder where the sample sheets were recorded have been
    modified to sample_sheet. All the recorded run must update the field
    of sample sheet.
    Script reads RunProcess table and check if sample_sheet contains the
    SampleSheet string and change it to sample_sheet.
"""

def run():
    # update SampleSheets
    run_objs = wetlab.models.RunProcess.objects.filter(
        sample_sheet__contains="SampleSheet"
    )
    for run_obj in run_objs:
        old_name = run_obj.get_sample_file()
        new_name = old_name.replace("/SampleSheets/", "/sample_sheet/")
        run_obj.set_run_sample_sheet(new_name)
    # update SampleSheets4LibPrep
    lib_user_objs = wetlab.models.LibUserSampleSheet.objects.filter(
        sample_sheet__contains="SampleSheets4LibPrep"
    )
    for lib_user_obj in lib_user_objs:
        old_name = lib_user_obj.get_lib_user_sample_sheet()
        new_name = old_name.replace(
            "/SampleSheets4LibPrep/", "/sample_sheets_lib_prep/"
        )
        lib_user_obj.set_lib_user_sample_sheet(new_name)
    return
