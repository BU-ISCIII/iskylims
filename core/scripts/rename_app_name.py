import core.models


"""
    The script is applicable for the upgrade from 2.3.0 to 3.0.0.
    Because the application in iSkylims have been renamed, this required that
    some tables where was indicated the application name must be
    renamed too.
    The impacted tables are:
        - LabRequest
        - MoleculeType
        - ProtocolType
        - SampleType
        - Species
        - PatientProjects
        - SampleProjects
        - MoleculeUsedFor
    Script reads the table which containg the old apps_name and remove the iskylims
    string.
"""


def run():
    tables_to_update = [
        "LabRequest",
        "MoleculeType",
        "ProtocolType",
        "SampleType",
        "Species",
        "PatientProjects",
        "SampleProjects",
        "MoleculeUsedFor",
    ]
    str_to_remove = "iSkyLIMS_"
    for table in tables_to_update:
        table_name = "core.models." + table
        if eval(table_name).objects.filter(apps_name__icontains=str_to_remove).exists():
            eval(table_name).objects.filter(apps_name__icontains=str_to_remove).update(
                apps_name="wetlab"
            )
    # perform this dummy check to use core.models to avoid litin error
    if core.models.MoleculeUsedFor.objects.filter(
        apps_name__icontains=str_to_remove
    ).exists():
        print("ERROR   IN  MIGRATIONS")
    return
