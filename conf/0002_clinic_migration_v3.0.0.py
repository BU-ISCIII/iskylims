# Generated by Django 3.2 on 2023-04-29 10:52

from django.db import migrations


class Migration(migrations.Migration):
    dependencies = [
        ("clinic", "0001_initial"),
    ]

    operations = [
        migrations.RenameField(
            model_name="clinicsamplerequest",
            old_name="clinicSampleState",
            new_name="clinic_sample_state",
        ),
        migrations.RenameField(
            model_name="clinicsamplerequest",
            old_name="confirmationCode",
            new_name="confirmation_code",
        ),
        migrations.RenameField(
            model_name="clinicsamplerequest",
            old_name="orderInEntry",
            new_name="entry_order",
        ),
        migrations.RenameField(
            model_name="clinicsamplerequest",
            old_name="patientCore",
            new_name="patient_core",
        ),
        migrations.RenameField(
            model_name="clinicsamplerequest",
            old_name="sampleCore",
            new_name="sample_core",
        ),
        migrations.RenameField(
            model_name="clinicsamplerequest",
            old_name="sampleRequestUser",
            new_name="sample_request_user",
        ),
        migrations.RenameField(
            model_name="clinicsamplerequest",
            old_name="serviceDate",
            new_name="service_date",
        ),
        migrations.RenameField(
            model_name="clinicsamplerequest",
            old_name="serviceUnit_id",
            new_name="service_unit_id",
        ),
        migrations.RenameField(
            model_name="clinicsamplestate",
            old_name="clinicState",
            new_name="clinic_state",
        ),
        migrations.RenameField(
            model_name="configsetting",
            old_name="configurationName",
            new_name="configuration_name",
        ),
        migrations.RenameField(
            model_name="configsetting",
            old_name="configurationValue",
            new_name="configuration_value",
        ),
        migrations.RenameField(
            model_name="doctor",
            old_name="doctorName",
            new_name="doctor_name",
        ),
        migrations.RenameField(
            model_name="doctor",
            old_name="serviceUnit_id",
            new_name="service_unit_id",
        ),
        migrations.RenameField(
            model_name="family",
            old_name="familyComments",
            new_name="family_comments",
        ),
        migrations.RenameField(
            model_name="family",
            old_name="familyID",
            new_name="family_id",
        ),
        migrations.RenameField(
            model_name="family",
            old_name="FamilyRelatives_id",
            new_name="family_relatives_id",
        ),
        migrations.RenameField(
            model_name="family",
            old_name="patienCore_id",
            new_name="patien_core_id",
        ),
        migrations.RenameField(
            model_name="family",
            old_name="relative1",
            new_name="relative_1",
        ),
        migrations.RenameField(
            model_name="family",
            old_name="relative2",
            new_name="relative_2",
        ),
        migrations.RenameField(
            model_name="familyrelatives",
            old_name="relationShip",
            new_name="relationship",
        ),
        migrations.RenameField(
            model_name="patientdata",
            old_name="notificationPreference",
            new_name="notification_preference",
        ),
        migrations.RenameField(
            model_name="patientdata",
            old_name="patienCore",
            new_name="patien_core",
        ),
        migrations.RenameField(
            model_name="patienthistory",
            old_name="entryDate",
            new_name="entry_date",
        ),
        migrations.RenameField(
            model_name="patienthistory",
            old_name="patientCore",
            new_name="patient_core",
        ),
        migrations.RenameField(
            model_name="serviceunits",
            old_name="serviceUnitName",
            new_name="service_unit_name",
        ),
    ]
