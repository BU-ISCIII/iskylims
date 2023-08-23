# Generic imports
import json
import re

from Bio.Seq import Seq

# Local imports
import core.models
import core.utils.commercial_kits
import core.utils.protocols
import core.utils.samples
import wetlab.config
import wetlab.models
import wetlab.utils.common
import wetlab.utils.samplesheet
import wetlab.utils.sequencers
import wetlab.utils.stats_graphs


def check_empty_fields(data):
    """
    Description:
        The function check if row_data contains empty values.
    Input:
        data:       # data to be checked
    Return:
        False is all fields contain data. True if any of them are empty
    """
    for row_data in data:
        for field in row_data:
            if field == "":
                return True
    return False


def analyze_and_store_prot_lib_param_values(form_data):
    """
    Description:
        The function get the user input  for the library preparation parameters and store them
        in database.
    Input:
        form_data   # User form
    Constant:
        HEADING_FIX_FOR_ADDING_LIB_PARAMETERS
        ERROR_EMPTY_VALUES
    Return:
        ERROR message if some of the data are missing, or a list of recorded library prepartion obj
    """
    stored_params = {}
    json_data = json.loads(form_data["protocol_data"])
    if check_empty_fields(json_data):
        stored_params = {}
        stored_params["ERROR"] = wetlab.config.ERROR_EMPTY_VALUES
        return stored_params
    prot_obj = core.utils.protocols.get_protocol_obj_from_id(form_data["protocol_id"])
    parameters = core.utils.protocols.get_protocol_parameters(prot_obj)
    if wetlab.models.AdditionaKitsLibPrepare.objects.filter(
        protocol_id=prot_obj
    ).exists():
        additional_kits = True
    else:
        additional_kits = False
    headings = (
        wetlab.config.HEADING_FIX_FOR_ADDING_LIB_PROT_PARAMETERS
        + ["lib_id"]
        + parameters
    )
    stored_params["data"] = []
    for row_index in range(len(json_data)):
        library_prep_obj = get_lib_prep_obj_from_id(json_data[row_index][2])
        for p_index in range(
            len(wetlab.config.HEADING_FIX_FOR_ADDING_LIB_PROT_PARAMETERS) + 1,
            len(json_data[row_index]),
        ):
            lib_parameter_value = {}
            lib_parameter_value[
                "parameter_id"
            ] = core.models.ProtocolParameters.objects.filter(
                protocol_id=prot_obj,
                parameter_name__exact=headings[p_index],
                parameter_used=True,
            ).last()
            lib_parameter_value["library_id"] = library_prep_obj
            lib_parameter_value["parameterValue"] = json_data[row_index][p_index]
            wetlab.models.LibParameterValue.objects.create_library_parameter_value(
                lib_parameter_value
            )

            if additional_kits:
                library_prep_obj.set_state("Updated parameters")
            else:
                library_prep_obj.set_state("Updated additional kits")
        stored_params["data"].append(
            [library_prep_obj.get_sample_name(), library_prep_obj.get_lib_prep_code()]
        )
    return stored_params


def create_library_preparation_instance(samples_data, user):
    """
    Description:
        The function create a new library preparation instance in database with
        user, moleculeID, sampleID, protocolID
        protocols.
    Input:
        samples_data :  Contains  moleculeID, sampleID, protocolID for each sample
    Constants:
        SAMPLE_NAMES_IN_SAMPLE_SHEET_CONTAIN_PROTOCOL_PREFIX
        PROTOCOL_SEPARATION_IN_SAMPLE_SHEET
    Functions:
        get_library_code_and_unique_id  # located at this file
    Return:
        lib_prep_for_params
    """
    lib_prep_for_params = {}
    prot_in_samples = wetlab.utils.common.get_configuration_value(
        "SAMPLE_NAMES_IN_SAMPLE_SHEET_CONTAIN_PROTOCOL_PREFIX"
    )
    if prot_in_samples == "TRUE":
        protocol_separation = wetlab.utils.common.get_configuration_value(
            "PROTOCOL_SEPARATION_IN_SAMPLE_SHEET"
        )

    for key, values in samples_data.items():
        s_name = core.utils.samples.get_sample_obj_from_id(key).get_sample_name()
        lib_prep_data = {}
        lib_prep_data["sample_id"] = key
        lib_prep_data["molecule_id"] = values["mol_id"]
        lib_prep_data["protocol_obj"] = core.models.Protocols.objects.filter(
            type__protocol_type__exact="Library Preparation",
            name__exact=values["prot_name"],
        ).last()
        lib_prep_data["prefixProtocol"] = values["prot_name"]
        if prot_in_samples == "TRUE":
            lib_prep_data["samplename_in_samplesheet"] = str(
                values["prot_name"] + protocol_separation + s_name
            )
        else:
            lib_prep_data["samplename_in_samplesheet"] = s_name
        lib_prep_data["registerUser"] = user
        lib_prep_data["user_sampleID"] = core.utils.samples.get_sample_obj_from_id(
            key
        ).get_sample_code()
        (
            lib_prep_data["lib_prep_code_id"],
            lib_prep_data["uniqueID"],
        ) = get_library_code_and_unique_id(
            lib_prep_data["sample_id"], lib_prep_data["molecule_id"]
        )
        lib_prep_obj = wetlab.models.LibPrepare.objects.create_lib_preparation(
            lib_prep_data
        )
        if values["prot_name"] not in lib_prep_for_params:
            lib_prep_for_params[values["prot_name"]] = []
        # create the structure data to get later the list of parameters
        lib_prep_for_params[values["prot_name"]].append(key)
        lib_prep_for_params[values["prot_name"]].append(lib_prep_obj.get_id())

    return lib_prep_for_params


def extract_protocol_library_preparation_form(form_data):
    """The function get the user input to assign sample to library preparation
        protocols. Get the lines that contains library protocol and return
        dictionnary with sample and protocol name
    Input:
        form_data :     User input form
    Return:
        extraction_data
    """
    json_data = json.loads(form_data["protocol_data"])
    extraction_data = {}

    for row_index in range(len(json_data)):
        if json_data[row_index][4] == "":
            continue
        extraction_data[json_data[row_index][1]] = {}
        extraction_data[json_data[row_index][1]]["mol_id"] = json_data[row_index][3]
        extraction_data[json_data[row_index][1]]["prot_name"] = json_data[row_index][4]
    return extraction_data


def get_protocol_parameters_for_library_preparation(protocol_libs):
    """
    Description:
        The function get the protocols parameters for the list of library preparation.
        In case that library preparations do not have the same protocol, only the
        ones that matches with the protocol of the first library preparation are
        considered

    Return:
        lib_prep_same_prot_parameters
    """

    protocol_lib_prep_data = []
    for prot_name, lib_prep in protocol_libs.items():
        prot_obj = core.utils.protocols.get_protocol_obj_from_name(prot_name)
        prot_lib_data = {"protocol_name": prot_obj.get_name()}
        prot_lib_data["protocol_id"] = prot_obj.get_protocol_id()
        prot_lib_data["lib_prep_data"] = list(
            wetlab.models.LibPrepare.objects.filter(pk__in=lib_prep).values_list(
                "sample_id__sample_name", "lib_prep_code_id", "pk"
            )
        )
        prot_lib_data[
            "param_heading_type"
        ] = core.utils.protocols.get_protocol_parameters_and_type(prot_obj)
        prot_lib_data[
            "fix_heading"
        ] = wetlab.config.HEADING_FIX_FOR_ADDING_LIB_PROT_PARAMETERS
        protocol_lib_prep_data.append(prot_lib_data)

    return protocol_lib_prep_data


def get_samples_for_library_preparation(user=None, friend_list=False):
    """_summary_

    Parameters
    ----------
    user : User instance, optional
        user instance to filter the samples only to the user, by default None
    friend_list : bool, optional
        Allow to extend the search for friend list, by default False

    Returns
    -------
    _type_
        _description_
    """
    samples_in_lib_prep = {}
    samples_in_lib_prep["avail_samples"] = {}
    s_in_lib_prep = core.utils.samples.get_sample_objs_in_state(
        "Library preparation", user, friend_list
    )

    s_lib_prep_defined = wetlab.models.LibPrepare.objects.filter(
        sample_id__in=s_in_lib_prep, lib_prep_state__lib_prep_state__exact="Defined"
    )
    s_lib_prep_upd_param = wetlab.models.LibPrepare.objects.filter(
        sample_id__in=s_in_lib_prep,
        lib_prep_state__lib_prep_state__exact="Updated parameters",
    )
    s_lib_prep_upd_kit = wetlab.models.LibPrepare.objects.filter(
        sample_id__in=s_in_lib_prep,
        lib_prep_state__lib_prep_state__exact="Updated additional kits",
    )

    if s_lib_prep_defined:
        samples_in_lib_prep["lib_prep_defined"] = list(
            s_lib_prep_defined.values_list(
                "sample_id__sample_name", "lib_prep_code_id", "protocol_id__name", "pk"
            )
        )

    if s_lib_prep_upd_param:
        samples_in_lib_prep["lib_prep_updated_param"] = ""

    if s_lib_prep_upd_kit:
        samples_in_lib_prep["sample_to_be_s_sheet"] = list(
            s_lib_prep_upd_kit.values_list(
                "sample_id__sample_name",
                "sample_id__pk",
                "lib_prep_code_id",
                "protocol_id__name",
            )
        )
        if wetlab.utils.sequencers.configuration_sequencer_exists():
            samples_in_lib_prep.update(
                wetlab.utils.sequencers.get_configuration_sequencers_data()
            )
            samples_in_lib_prep["display_sample_sheet"] = True
            samples_in_lib_prep["avail_samples"][
                "lib_prep_protocols"
            ] = get_protocols_for_library_preparation()

    # find out the samples that have not be defined in LibPrepare
    if len(s_in_lib_prep) > (
        len(s_lib_prep_defined) + len(s_lib_prep_upd_param) + len(s_lib_prep_upd_kit)
    ):
        s_def = list(s_lib_prep_defined.values_list("sample_id__pk", flat=True))
        s_param = list(s_lib_prep_upd_param.values_list("sample_id__pk", flat=True))
        s_kit = list(s_lib_prep_upd_kit.values_list("sample_id__pk", flat=True))
        s_in_lib_prep_ids = list(s_in_lib_prep.values_list("pk", flat=True))
        def_lib = s_def + s_param + s_kit
        s_not_def_ids = list(set(s_in_lib_prep_ids).symmetric_difference(set(def_lib)))
        # s_not_def = wetlab.models.LibPrepare.objects.filter(pk__in=s_not_def_ids)

        molecules_obj = core.models.MoleculePreparation.objects.filter(
            sample__pk__in=s_not_def_ids
        )
        samples_in_lib_prep["avail_samples"]["data"] = list(
            molecules_obj.values_list(
                "sample__sample_name", "sample_id__pk", "molecule_code_id", "pk"
            )
        )
        samples_in_lib_prep["avail_samples"][
            "heading"
        ] = wetlab.config.HEADING_FOR_SAMPLES_TO_DEFINE_PROTOCOL
        samples_in_lib_prep["avail_samples"][
            "lib_prep_protocols"
        ] = get_protocols_for_library_preparation()

    return samples_in_lib_prep


def extract_userids_from_sample_sheet_data(file_read):
    """
    Description:
        The function extract userID and checks if userIDs are defined in database
    Input:
        file_read    # String having the file from IEM
    Functions:
        get_userid_in_user_iem_file # located at utils/sample_sheet_utils.py
        delete_stored_file     # located at utils/sample_sheet_utils.py
    Constant:
        ERROR_UNABLE_TO_DELETE_USER_FILE
        ERROR_SAMPLE_SHEET_DOES_NOT_HAVE_DESCRIPTION_FIELD
        ERROR_SAMPLE_SHEET_WHEN_FETCHING_USERID_NAMES
    Return
        data with the file name and the content of the sample sheet.
        data['Error'] if file is invalid
    """

    user_id_list = wetlab.utils.common.get_userid_list()

    user_ids = wetlab.utils.samplesheet.validate_userid_in_user_iem_file(
        file_read, user_id_list
    )
    if user_ids == [""]:
        user_ids["ERROR"] = wetlab.config.ERROR_SAMPLE_SHEET_WHEN_FETCHING_USERID_NAMES
    return user_ids


def validate_sample_sheet_data(input_data):
    """
    Description:
        The function checks if library preparation samples are defined are they are in Updated additional kits state ,
        and no duplication index exists.
    Input:
        input_data  # contain a dictionary with sample names, heading to match the values and
                    a list with all  information samples
    Functions:
        valid_samples_for_lib_preparation    # located at this file
        find_duplicate_index                 # located at this file
        check_collection_index_exists        # located at wetlab/utils/collection_index_functions.py
    Constant:
        ERROR_SAMPLES_INVALID_STATE_FOR_LIBRARY_PREPARATION
        ERROR_SAMPLE_SHEET_CONTAINS_NOT_DEFINED_SAMPLES
    Return
        sample data objects if all checks are valid or ERROR if file is invalid
    """
    # check that samples are defined and in the right state

    not_defined_samples = []
    invalid_state_samples = []

    for sample in input_data["samples"]:
        if not wetlab.models.LibPrepare.objects.filter(
            samplename_in_samplesheet__exact=sample
        ).exists():
            not_defined_samples.append(sample)
            continue
        if not wetlab.models.LibPrepare.objects.filter(
            samplename_in_samplesheet__exact=sample,
            lib_prep_state__lib_prep_state__exact="Updated additional kits",
        ).exists():
            invalid_state_samples.append(sample)
    if len(not_defined_samples) > 0:
        if (
            wetlab.utils.common.get_configuration_value(
                "SAMPLE_NAMES_IN_SAMPLE_SHEET_CONTAIN_PROTOCOL_PREFIX"
            )
            == "TRUE"
        ):
            not_defined = (
                wetlab.config.ERROR_SAMPLE_SHEET_CONTAINS_NOT_DEFINED_SAMPLES_WITH_PROTOCOL
            )
        else:
            not_defined = wetlab.config.ERROR_SAMPLE_SHEET_CONTAINS_NOT_DEFINED_SAMPLES
        return {"ERROR": not_defined}
    if len(invalid_state_samples) > 0:
        return {
            "ERROR": wetlab.config.ERROR_SAMPLES_INVALID_STATE_FOR_LIBRARY_PREPARATION
        }
    # check if sample sheet has duplicate index
    duplicate_index = find_duplicate_index(
        input_data["sample_data"], input_data["heading"]
    )
    if "ERROR" in duplicate_index:
        return duplicate_index
    return "Validated"


def find_duplicate_index(sample_row_data, heading):
    """
    Description:
        The function get the sample rows from sample sheet to check if the samples have duplicated
        index.
    Input:
        sample_row_data  # contains the row sample list
        heading     # contains the list of heading to get the index colums
    Constant:
        ERROR_SAMPLES_INVALID_DUPLICATED_INDEXES
    Return
        False if not duplication found. error message if exists duplicted index
    """
    index_values = {}
    duplicated_index_sample = []
    i5_index = False
    i7_index = heading.index("I7_Index_ID")
    if "I5_Index_ID" in heading:
        i5_index = heading.index("I5_Index_ID")
    sample_name_index = heading.index("Sample_Name")
    for sample_row in sample_row_data:
        if i5_index:
            indexes_in_sample = str(sample_row[i7_index] + "_" + sample_row[i5_index])
        else:
            indexes_in_sample = sample_row[i7_index]
        if indexes_in_sample not in index_values:
            index_values[indexes_in_sample] = []
        else:
            duplicated_index_sample.append(sample_row[sample_name_index])

        if not sample_row[sample_name_index] in index_values:
            index_values[sample_row[sample_name_index]] = []
        index_values[sample_row[sample_name_index]].append([indexes_in_sample])
    if len(duplicated_index_sample) > 0:
        return {"ERROR": wetlab.config.ERROR_SAMPLES_INVALID_DUPLICATED_INDEXES}
    return "False"


def get_protocols_for_library_preparation():
    protocol_list = []
    if core.models.Protocols.objects.filter(
        type__protocol_type__exact="Library Preparation"
    ).exists():
        protocols = core.models.Protocols.objects.filter(
            type__protocol_type__exact="Library Preparation"
        )
        for protocol in protocols:
            protocol_list.append(protocol.get_name())
    return protocol_list


def get_all_library_information(sample_id):
    """
    Description:
        The function get the library preparation information for sample.
        It return a dictionary with heading and lib_prep_data which is a list
        having index values, parameter heading, parameter values, library preparation id
        and library preparation codeID
    Input:
        sample_id  # sample id
    Constants:
        HEADING_FOR_LIBRARY_PREPARATION_DEFINITION
        HEADING_FOR_DISPLAY_POOL_INFORMATION_IN_SAMPLE_INFO
    Return
        library_information
    """
    library_information = {}
    if wetlab.models.LibPrepare.objects.filter(sample_id__pk__exact=sample_id).exists():
        library_information[
            "library_definition_heading"
        ] = wetlab.config.HEADING_FOR_LIBRARY_PREPARATION_DEFINITION
        library_information["library_definition"] = []
        library_information["pool_information"] = []
        library_preparation_items = wetlab.models.LibPrepare.objects.filter(
            sample_id__pk__exact=sample_id
        ).exclude(lib_prep_state__lib_prep_state__exact="Created for Reuse")
        library_information["lib_prep_param_value"] = []
        library_information["lib_prep_data"] = []
        for library_item in library_preparation_items:
            lib_prep_data = []
            lib_prep_data.append(library_item.get_info_for_display())
            protocol_used_obj = library_item.get_protocol_obj()
            if core.models.ProtocolParameters.objects.filter(
                protocol_id=protocol_used_obj
            ).exists():
                parameter_names = core.models.ProtocolParameters.objects.filter(
                    protocol_id=protocol_used_obj
                ).order_by("parameter_order")
                lib_prep_param_heading = ["Lib Preparation CodeID"]
                lib_prep_param_value = [library_item.get_lib_prep_code()]
                for p_name in parameter_names:
                    lib_prep_param_heading.append(p_name.get_parameter_name())
                    if wetlab.models.LibParameterValue.objects.filter(
                        library_id=library_item
                    ).exists():
                        try:
                            lib_prep_param_value.append(
                                wetlab.models.LibParameterValue.objects.filter(
                                    library_id=library_item, parameter_id=p_name
                                )
                                .last()
                                .get_parameter_information()
                            )
                        except wetlab.models.LibParameterValue.DoesNotExist:
                            lib_prep_param_value.append("")
                lib_prep_data.append(lib_prep_param_heading)
                lib_prep_data.append(lib_prep_param_value)
            else:
                lib_prep_data.append("")
                lib_prep_data.append("")
            lib_prep_data.append(library_item.get_id())
            lib_prep_data.append(library_item.get_lib_prep_code())
            lib_prep_data.append(library_item.get_molecule_code_id())
            lib_prep_data.append(library_item.get_sample_id())
            library_information["lib_prep_data"].append(lib_prep_data)

            if library_item.pools.all().exists():
                pools = library_item.pools.all()
                lib_prep_code_id = library_item.get_lib_prep_code()
                for pool in pools:
                    pool_name = pool.get_pool_name()
                    pool_code = pool.get_pool_code_id()
                    run_name = pool.get_run_name()
                    library_information["pool_information"].append(
                        [
                            lib_prep_code_id,
                            pool_name,
                            pool_code,
                            run_name,
                            library_item.get_id(),
                        ]
                    )

        if library_information["pool_information"]:
            library_information[
                "pool_heading"
            ] = wetlab.config.HEADING_FOR_DISPLAY_POOL_INFORMATION_IN_SAMPLE_INFO

    return library_information


def get_iem_version_for_library_prep_ids(lib_prep_id_list):
    """
    Description:
        The function get the list of library preparation ids and return the IEM version used
        when user create the sample sheet .
    Input:
        lib_prep_id_list        # ID List of the library preparation
    Return:
        iem_version
    """
    versions = []
    for lib_prep_id in lib_prep_id_list:
        lib_prep_obj = get_lib_prep_obj_from_id(lib_prep_id)
        if lib_prep_obj:
            version = lib_prep_obj.get_iem_version()
            if version != "None" and version not in versions:
                versions.append(version)
    return versions


def get_lib_prep_to_add_parameters():
    """
    Description:
        The function will return a list with samples which are needs to add library preparation parameters
    Variables:
        library_prep_information # Dictionary with the heading and the molecule information
    Return:
        lib_prep_parameters.
    """
    lib_prep_parameters = {}
    lib_prep_parameters["length"] = 0
    if wetlab.models.LibPrepare.objects.filter(
        lib_prep_state__lib_prep_state__exact="Defined"
    ).exists():
        samples = wetlab.models.LibPrepare.objects.filter(
            lib_prep_state__lib_prep_state__exact="Defined"
        )
        sample_info = []
        for sample in samples:
            lib_prep_info = []
            lib_prep_info.append(sample.get_lib_prep_code())
            lib_prep_info.append(sample.get_sample_name())
            lib_prep_info.append(sample.get_protocol_used())
            lib_prep_info.append(sample.get_id())
            sample_info.append(lib_prep_info)
        lib_prep_parameters["lib_prep_info"] = sample_info
        lib_prep_parameters[
            "lib_prep_heading"
        ] = wetlab.config.HEADING_FOR_ADD_LIBRARY_PREPARATION_PARAMETERS
        lib_prep_parameters["length"] = len(sample_info)
    return lib_prep_parameters


def get_protocol_from_library_id(library_prep_id):
    """
    Description:
        The function will return a list with samples which are needs to add library preparation parameters
    Input:
        library_prep_id # id to get the protocol name
    Return:
        protocol name or empty if library id does not exists.
    """
    if wetlab.models.LibPrepare.objects.filter(pk__exact=library_prep_id).exists():
        return wetlab.models.LibPrepare.objects.get(
            pk__exact=library_prep_id
        ).get_protocol_used()
    return ""


def get_samples_in_lib_prep_state():
    """
    Description:
        The function will return a list with samples which are in add_library_preparation state.
        Include the ones that are requested to reprocess
    Return:
        samples_in_lib_prep.
    """

    samples_in_lib_prep = {}
    lib_prep_data = []
    states_excluded = ["Completed", "Reused pool"]
    if core.models.Samples.objects.filter(
        sample_state__sample_state_name__exact="Library preparation"
    ).exists():
        samples_obj = (
            core.models.Samples.objects.filter(
                sample_state__sample_state_name__exact="Library preparation"
            )
            .order_by("sample_user")
            .order_by("sample_entry_date")
        )

        for sample in samples_obj:
            if (
                not wetlab.models.LibPrepare.objects.filter(sample_id=sample).exists()
            ) or (
                wetlab.models.LibPrepare.objects.filter(sample_id=sample)
                .exclude(lib_prep_state__lib_prep_state__in=states_excluded)
                .exists()
            ):
                sample_information = sample.get_info_in_defined_state()
                sample_information.append(sample.get_register_user())
                molecule_obj = core.models.MoleculePreparation.objects.filter(
                    sample=sample, state__molecule_state_name="Completed"
                ).last()
                if molecule_obj:
                    molecule_data = molecule_obj.get_molecule_information()
                else:
                    molecule_data = [""] * 4
                lib_prep_data.append(sample_information + molecule_data)

        samples_in_lib_prep["library_information"] = lib_prep_data
        samples_in_lib_prep[
            "lib_prep_heading"
        ] = wetlab.config.HEADING_FOR_LIBRARY_PREPARATION_STATE
        samples_in_lib_prep["length"] = len(lib_prep_data)

        return samples_in_lib_prep
    else:
        samples_in_lib_prep["length"] = 0
        return samples_in_lib_prep


def find_index_sequence_collection_values_kit(sequence):
    """
    Description:
        The function will try to find the sequence by looking on the I7 index and if not matched
        check on the I5 first on the forward and then on the reverse sequence.
    Input:
        sequence        : Sequence for searching
    Return:
        index_found and the sequence
    """

    if wetlab.models.CollectionIndexValues.objects.filter(
        i_7_seq__icontains=sequence
    ).exists():
        return ["I7", sequence]
    if wetlab.models.CollectionIndexValues.objects.filter(
        i_5_seq__icontains=sequence
    ).exists():
        return ["I5", sequence]
    rev_sequence = str(Seq(sequence).reverse_complement())
    if wetlab.models.CollectionIndexValues.objects.filter(
        i_5_seq__icontains=rev_sequence
    ).exists():
        return ["I5", rev_sequence]
    return "None", sequence


def store_confirmation_library_preparation_index(form_data):
    """
    Description:
        The function will fetch the indexes defined in the confirmed sample sheet
        and updated the library preparation sample with index information.
        It updated also the userid field.
    Input:
        form_data               # data included in the form
        user_sample_sheet_obj   # user sample sheet object
    Constant:
        ERROR_USER_SAMPLE_SHEET_NO_LONGER_EXISTS
        ERROR_LIBRARY_PREPARATION_NOT_EXISTS
        MAP_USER_SAMPLE_SHEET_TO_DATABASE_ALL_PLATFORMS
    Return:
        store_result .
    """
    store_result = {}
    unable_store_lib_prep = []
    json_data = json.loads(form_data["index_data"])
    heading = form_data["heading_excel"].split(",")
    store_result = {}
    if not wetlab.models.LibUserSampleSheet.objects.filter(
        pk__exact=form_data["libPrepUserSampleSheetId"]
    ).exists():
        store_result["ERROR"] = wetlab.config.ERROR_USER_SAMPLE_SHEET_NO_LONGER_EXISTS
        return store_result
    user_sample_sheet_obj = wetlab.models.LibUserSampleSheet.objects.get(
        pk__exact=form_data["libPrepUserSampleSheetId"]
    )
    sample_name_index = heading.index("Sample_Name")

    for row_index in range(len(json_data)):
        lib_prep_data = {}
        sample_name = json_data[row_index][sample_name_index]
        if wetlab.models.LibPrepare.objects.filter(
            samplename_in_samplesheet__exact=sample_name,
            lib_prep_state__lib_prep_state__exact="Updated additional kits",
        ).exists():
            lib_prep_obj = wetlab.models.LibPrepare.objects.filter(
                samplename_in_samplesheet__exact=sample_name,
                lib_prep_state__lib_prep_state__exact="Updated additional kits",
            ).last()

            for item in wetlab.config.MAP_USER_SAMPLE_SHEET_TO_DATABASE_ALL_PLATFORMS:
                if item[0] in heading:
                    try:
                        lib_prep_data[item[1]] = json_data[row_index][
                            heading.index(item[0])
                        ]
                    except Exception:
                        lib_prep_data[item[1]] = None
            # if Single reads then set the index 5 to empty
            if "I5_Index_ID" not in heading:
                lib_prep_data["i5IndexID"] = ""
                lib_prep_data["i5Index"] = ""
            lib_prep_data["user_sample_sheet"] = user_sample_sheet_obj
            lib_prep_obj.update_library_preparation_with_indexes(lib_prep_data)
            # Update library preparation and sample state
            lib_prep_obj.set_state("Completed")
            lib_prep_obj.get_sample_obj().set_state("Pool Preparation")
        else:
            # ERROR
            unable_store_lib_prep.append(sample_name)
            continue

    if len(unable_store_lib_prep) > 0:
        store_result["ERROR"] = wetlab.config.ERROR_LIBRARY_PREPARATION_NOT_EXISTS
        store_result["ERROR"].append(unable_store_lib_prep)
    else:
        user_sample_sheet_obj.update_confirm_used(True)
        store_result["Successful"] = True
    return store_result


def store_library_preparation_sample_sheet(
    sample_sheet_data, user, platform, configuration
):
    """
    Description:
        The function will get the extracted data from sample sheet.
        Then store the libraryPreparation data for each sample and update the sample state to "library Preparation"
    Input:
        sample_sheet_data   # extracted data from sample sheet in dictionary format
        user        # user object
        platform    # platform used in the sample sheet
        configuration # configuration used in the sample sheet
    Return:
        new_user_s_sheet_obj .
    """
    sample_sheet_data["user"] = user
    sample_sheet_data["platform"] = platform
    sample_sheet_data["configuration"] = configuration
    new_user_s_sheet_obj = (
        wetlab.models.LibUserSampleSheet.objects.create_lib_prep_user_sample_sheet(
            sample_sheet_data
        )
    )

    return new_user_s_sheet_obj


def get_library_code_and_unique_id(sample_id, molecule_id):
    """
    Description:
        The function will find out the latest library preparation uniqueID", increment the value
        and will return the updated value to use
    Input:
        sample_id        # id of the sample
        molecule_id      # id of the molecule
    Return:
        uniqueID .
    """
    sample_obj = core.utils.samples.get_sample_obj_from_id(sample_id)
    molecule_obj = core.utils.samples.get_molecule_obj_from_id(molecule_id)
    if wetlab.models.LibPrepare.objects.filter(
        sample_id=sample_obj, molecule_id=molecule_obj
    ).exists():
        lib_prep_obj = wetlab.models.LibPrepare.objects.filter(
            sample_id=sample_obj, molecule_id=molecule_obj
        ).last()
        last_lib_prep_code_id = lib_prep_obj.get_lib_prep_code()
        split_code = re.search("(.*_)(\d+)$", last_lib_prep_code_id)
        index_val = int(split_code.group(2))
        new_index = str(index_val + 1).zfill(2)
        lib_prep_code_id = split_code.group(1) + new_index
        s_uniqueID = sample_obj.get_unique_sample_id()
        # count the number that library preparation was used on the same sample
        lib_prep_times = str(
            wetlab.models.LibPrepare.objects.filter(sample_id=sample_obj).count()
        )
        uniqueID = s_uniqueID + "-" + lib_prep_times
    else:
        lib_prep_code_id = molecule_obj.get_molecule_code_id() + "_LIB_01"
        uniqueID = sample_obj.get_unique_sample_id() + "-1"
    return lib_prep_code_id, uniqueID


def format_sample_sheet_to_display_in_form(sample_sheet_data):
    """
    Description:
        The function gets the information to display the library preparation to assing
        the protocol parameter values sample names, extracted data from sample sheet, index, .
        Then store the libraryPreparation data for each sample and update the sample state to
        "library Preparation"
    Input:
        lib_prep_ids   # Library preparation ids
        user        # user object
        user_sample_sheet_obj # user_sample_sheet object for assigning to each library preparation
    Constamt:
        HEADING_MAIN_DATA_SAMPLE_SHEET
        HEADING_SUMMARY_DATA_SAMPLE_SHEET
    Return:
        display_data
    """
    display_data = {}
    main_data_heading = wetlab.config.HEADING_MAIN_DATA_SAMPLE_SHEET.copy()
    extract_values = [
        "application",
        "instrument type",
        "assay",
        "index_adapters",
        "reads",
        "adapter1",
        "adapter2",
    ]
    if "" == sample_sheet_data["adapter2"]:
        extract_values.pop()
        main_data_heading.pop()
        display_data["adapter2"] = False
    else:
        display_data["adapter2"] = True

    display_data["sample_data"] = sample_sheet_data["sample_data"]
    display_data["heading"] = sample_sheet_data["heading"]
    main_values = []
    for item in extract_values:
        main_values.append(sample_sheet_data[item])
    summary_values = []
    summary_values.append(len(sample_sheet_data["samples"]))
    summary_values.append(sample_sheet_data["proyects"])
    summary_values.append(sample_sheet_data["userid_names"])
    display_data["main_data"] = list(zip(main_data_heading, main_values))
    display_data["summary_data"] = list(
        zip(wetlab.config.HEADING_SUMMARY_DATA_SAMPLE_SHEET, summary_values)
    )
    display_data["heading_excel"] = ",".join(sample_sheet_data["heading"])
    # if len(sample_sheet_data['userid_names']) == 0:
    #    display_data['no_user_defined'] = True

    return display_data


def get_lib_prep_obj_from_id(library_preparation_id):
    """
    Description:
        The function gets the library preparation id and it returns the object instance
    Input:

    Return:
        library_preparation_obj or None if not match
    """

    if wetlab.models.LibPrepare.objects.filter(
        pk__exact=library_preparation_id
    ).exists():
        library_preparation_obj = wetlab.models.LibPrepare.objects.get(
            pk__exact=library_preparation_id
        )
        return library_preparation_obj
    else:
        return "None"


def separate_lib_prep_on_protocols(lib_prep_ids):
    """_summary_

    Parameters
    ----------
    lib_prep_ids : _type_
        _description_
    """
    protocols = {}
    lib_prep_protocols = wetlab.models.LibPrepare.objects.filter(
        pk__in=lib_prep_ids
    ).values("protocol_id__name", "pk")
    for lib_prep_protocol in lib_prep_protocols:
        if lib_prep_protocol["protocol_id__name"] not in protocols:
            protocols[lib_prep_protocol["protocol_id__name"]] = []
        protocols[lib_prep_protocol["protocol_id__name"]].append(
            lib_prep_protocol["pk"]
        )
    return protocols


def update_batch_lib_prep_sample_state(lib_prep_ids, sample_state):
    """
    Description:
        The function set the sample state having as input the list of library preparation ids
    Input:
        lib_prep_ids        # list of library preparation ids
        sample_state        # state to be set
    Return:
        None
    """
    for lib_id in lib_prep_ids:
        lib_obj = wetlab.models.LibPrepare.objects.get(pk__exact=lib_id)
        lib_obj.get_sample_obj().set_state(sample_state)

    return


def update_library_preparation_for_reuse(sample):
    """
    Description:
        The function step the reuse of the library preparation
    Input:
        sample_list        # list of samples for updating the reuse value

    Return:
        None
    """
    if wetlab.models.LibPrepare.objects.filter(
        sample_id__sample_name__exact=sample
    ).exists():
        lib_prep_obj = wetlab.models.LibPrepare.objects.filter(
            sample_id__sample_name__exact=sample
        ).last()
        lib_prep_obj.set_increase_reuse()
    return
