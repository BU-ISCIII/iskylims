from django.contrib import admin

from core.models import (
    City,
    CommercialKits,
    Contact,
    LabRequest,
    MoleculeParameterValue,
    MoleculePreparation,
    MoleculeType,
    MoleculeUsedFor,
    OntologyMap,
    PatientCore,
    PatientProjects,
    PatientSex,
    ProtocolParameters,
    Protocols,
    ProtocolType,
    SampleProjectFieldClassification,
    SampleProjects,
    SampleProjectsFields,
    SampleProjectsFieldsValue,
    Samples,
    SamplesProjectsTableOptions,
    SampleType,
    SequencerInLab,
    SequencingConfiguration,
    SequencingPlatform,
    Species,
    StateInCountry,
    StatesForMolecule,
    StatesForSample,
    UserLotCommercialKits,
)


class ContactAdmin(admin.ModelAdmin):
    list_display = ["contact_name", "contact_mail"]


class StateInCountryAdmin(admin.ModelAdmin):
    list_display = ["state_name", "apps_name"]


class CityAdmin(admin.ModelAdmin):
    list_display = ["city_name", "apps_name"]


class CommercialKitsAdmin(admin.ModelAdmin):
    list_display = ("name", "provider", "platform_kits", "cat_number")


class LabRequestAdmin(admin.ModelAdmin):
    list_display = [
        "lab_name",
        "lab_name_coding",
        "lab_contact_name",
        "lab_phone",
        "lab_email",
    ]


class MoleculeTypeAdmin(admin.ModelAdmin):
    list_display = ("molecule_type",)


class MoleculeParameterValueAdmin(admin.ModelAdmin):
    list_display = ("molecule_parameter_id", "molecule_id", "parameter_value")


class MoleculePreparationAdmin(admin.ModelAdmin):
    list_display = (
        "molecule_code_id",
        "state",
        "sample",
        "molecule_type",
        "extraction_type",
        "protocol_used",
        "molecule_extraction_date",
        "molecule_used_for",
        "reused_number",
    )
    list_filter = ("generated_at",)
    search_fields = ("sample__startswith",)


class MoleculeUsedForAdmin(admin.ModelAdmin):
    list_display = ["used_for", "apps_name", "massive_use"]


class OntologyMapAdmin(admin.ModelAdmin):
    list_display = ["label", "ontology"]


class PatientCoreAdmin(admin.ModelAdmin):
    list_display = ("patient_code", "patient_name", "patient_surname", "patient_sex")
    search_fields = ("patient_code__icontains",)


class PatientSexAdmin(admin.ModelAdmin):
    list_display = ("sex",)


class PatientProjectsAdmin(admin.ModelAdmin):
    list_display = ("project_name", "project_description")


class ProtocolsAdmin(admin.ModelAdmin):
    list_display = ("name", "type", "description")


class ProtocolTypeAdmin(admin.ModelAdmin):
    list_display = ("protocol_type", "molecule", "apps_name")


class ProtocolParametersAdmin(admin.ModelAdmin):
    list_display = (
        "protocol_id",
        "parameter_name",
        "parameter_order",
        "parameter_used",
        "parameter_min_value",
        "parameter_max_value",
        "parameter_description",
    )
    list_filter = ("protocol_id",)
    search_fields = ("parameter_name__startswith",)


class SamplesAdmin(admin.ModelAdmin):
    list_display = (
        "sample_code_id",
        "sample_name",
        "sample_state",
        "lab_request",
        "sample_type",
        "sample_user",
        "species",
        "sample_project",
        "sample_entry_date",
        "unique_sample_id",
        "reused_number",
        "sequencing_date",
    )
    list_filter = ("generated_at",)
    search_fields = ("sample_name__icontains",)


class SampleTypeAdmin(admin.ModelAdmin):
    list_display = ("sample_type", "apps_name", "optional_fields")


class SampleProjectsAdmin(admin.ModelAdmin):
    list_display = [
        "sample_project_name",
        "sample_project_manager",
        "sample_project_contact",
        "sample_project_description",
        "apps_name",
    ]


class SampleProjectsFieldsAdmin(admin.ModelAdmin):
    list_display = [
        "sample_projects_id",
        "sample_project_field_name",
        "sample_project_field_type",
        "sample_project_option_list",
        "sample_project_field_order",
        "sample_project_field_used",
    ]


class SampleProjectFieldClassificationAdmin(admin.ModelAdmin):
    list_display = ["classification_name", "classification_display"]


class SamplesProjectsTableOptionsAdmin(admin.ModelAdmin):
    list_display = ["sample_project_field", "option_value"]


class SampleProjectsFieldsValueAdmin(admin.ModelAdmin):
    list_display = [
        "sample_id",
        "sample_project_field_id",
        "sample_project_field_value",
    ]
    list_filter = ["sample_id"]
    search_fields = (
        "sample_project_field_id__sample_project_field_name",
        "sample_project_field_value",
    )


class SequencingPlatformAdmin(admin.ModelAdmin):
    list_display = ("platform_name", "company_name", "sequencing_technology")


class SequencingConfigurationAdmin(admin.ModelAdmin):
    list_display = ["platform_id", "configuration_name"]


class SequencerInLabAdmin(admin.ModelAdmin):
    list_display = (
        "sequencer_name",
        "platform_id",
        "sequencer_description",
        "sequencer_location",
        "sequencer_serial_number",
        "sequencer_state",
        "sequencer_operation_start",
        "sequencer_operation_end",
        "sequencer_number_lanes",
    )


class SpeciesAdmin(admin.ModelAdmin):
    list_display = (
        "species_name",
        "ref_genome_name",
        "ref_genome_size",
        "ref_genome_id",
    )


class StatesForMoleculeAdmin(admin.ModelAdmin):
    list_display = ("molecule_state_name",)


class StatesForSampleAdmin(admin.ModelAdmin):
    list_display = ("sample_state_name",)


class UserLotCommercialKitsAdmin(admin.ModelAdmin):
    list_display = (
        "user",
        "based_commercial",
        "run_out",
        "uses_number",
        "chip_lot",
        "latest_used_date",
        "expiration_date",
    )


admin.site.register(StateInCountry, StateInCountryAdmin)
admin.site.register(City, CityAdmin)
admin.site.register(CommercialKits, CommercialKitsAdmin)
admin.site.register(LabRequest, LabRequestAdmin)
admin.site.register(MoleculeType, MoleculeTypeAdmin)
admin.site.register(MoleculeUsedFor, MoleculeUsedForAdmin)
admin.site.register(MoleculePreparation, MoleculePreparationAdmin)
admin.site.register(OntologyMap, OntologyMapAdmin)
admin.site.register(PatientCore, PatientCoreAdmin)
admin.site.register(PatientSex, PatientSexAdmin)
admin.site.register(PatientProjects, PatientProjectsAdmin)
admin.site.register(ProtocolType, ProtocolTypeAdmin)
admin.site.register(Protocols, ProtocolsAdmin)
admin.site.register(ProtocolParameters, ProtocolParametersAdmin)
admin.site.register(MoleculeParameterValue, MoleculeParameterValueAdmin)
admin.site.register(SampleType, SampleTypeAdmin)
admin.site.register(Samples, SamplesAdmin)
admin.site.register(SampleProjects, SampleProjectsAdmin)
admin.site.register(
    SampleProjectFieldClassification, SampleProjectFieldClassificationAdmin
)
admin.site.register(SampleProjectsFields, SampleProjectsFieldsAdmin)
admin.site.register(SamplesProjectsTableOptions, SamplesProjectsTableOptionsAdmin)
admin.site.register(SampleProjectsFieldsValue, SampleProjectsFieldsValueAdmin)
admin.site.register(SequencingPlatform, SequencingPlatformAdmin)
admin.site.register(SequencingConfiguration, SequencingConfigurationAdmin)
admin.site.register(SequencerInLab, SequencerInLabAdmin)
admin.site.register(Species, SpeciesAdmin)
admin.site.register(StatesForMolecule, StatesForMoleculeAdmin)
admin.site.register(StatesForSample, StatesForSampleAdmin)
admin.site.register(UserLotCommercialKits, UserLotCommercialKitsAdmin)
admin.site.register(Contact, ContactAdmin)
