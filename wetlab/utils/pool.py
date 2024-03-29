# Generic imports
from datetime import date
import json

from django.db.models import Q

# Local imports
import core.utils.common
import wetlab.config
import wetlab.models


def get_lib_prep_adapter(lib_prep_ids):
    """
    Description:
        The function gets the different adapters used for the library preparations instance.
        Returns the adapter value
    Input:
        lib_prep_ids  # library preparation id list
    Return:
        returns a list of adapters in adapter_list
    """
    adapters = []
    for lib_prep_id in lib_prep_ids:
        if not wetlab.models.LibPrepare.objects.filter(pk__exact=lib_prep_id).exists():
            continue
        lib_prep_obj = wetlab.models.LibPrepare.objects.get(pk__exact=lib_prep_id)
        adapters.append(lib_prep_obj.get_adapters()[0])
    adapter_list = list(set(adapters))
    return adapter_list


def check_if_duplicated_index(lib_prep_ids):
    """
    Description:
        The function check if some of the library preparations instances has the same indexes.
        If there is not duplicated returns True else returns the list of the library preparation ids which
        have duplicated index .
    Input:
        lib_prep_ids  # library preparation id list
    Return:
        True if unique indexes are used.
        a dictionary with 'incompatible_index' as key and the list of samples names as value.
    """
    index_values = {}
    for lib_prep_id in lib_prep_ids:
        lib_prep_obj = wetlab.models.LibPrepare.objects.get(pk__exact=lib_prep_id)
        # lib_prep_count += 1
        lib_index = lib_prep_obj.get_indexes()
        combined_index = str(lib_index["i7_indexID"] + "_" + lib_index["i5_indexID"])
        """
        lib_prep_paired_end = lib_prep_obj.get_single_paired()
        if not lib_prep_paired_end in s_paired_values:
            s_paired_values[lib_prep_paired_end] = []
        s_paired_values[lib_prep_paired_end].append(lib_prep_id)
        """
        if combined_index not in index_values:
            index_values[combined_index] = []
        index_values[combined_index].append(lib_prep_id)
    if len(index_values) == len(lib_prep_ids):
        return "True"
    else:
        incompatible_samples = {}
        incompatible_index = []
        # check which index
        for key, values in index_values.items():
            if len(values) > 1:
                s_name = []
                for value in values:
                    s_name.append(
                        wetlab.models.LibPrepare.objects.get(
                            pk__exact=value
                        ).get_sample_name()
                    )
                incompatible_index.append([" and  ".join(s_name), key])
        incompatible_samples["incompatible_index"] = incompatible_index
        return incompatible_samples


def get_single_paired(lib_prep_ids):
    """
    Description:
        The function checks if at least one library preparation included in the pool
        was usin I5 index.
    Input:
        lib_prep_ids  # library preparation id list

    Return:
        PairedEnd or SingleRead
    """
    for lib_prep_id in lib_prep_ids:
        i5_index_value = wetlab.models.LibPrepare.objects.get(
            pk__exact=lib_prep_id
        ).get_i5_index()
        if i5_index_value != "":
            return "PairedEnd"
    return "SingleRead"


def define_new_pool(form_data, user_obj):
    """
    Description:
        The function performs some checks to verify the new pool can be created.

    Input:
        lib_prep_ids  # library preparation id list
    Constants:
        ERROR_LIBRARY_PREPARATION_NOT_EXISTS
        ERROR_NOT_LIBRARY_PREPARATION_SELECTED
    Functions:
        check_single_paired_compatible      # located at this file
        check_if_duplicated_index      # located at this file
        get_lib_prep_adapter      # located at this file
        generate_pool_code_id      # located at this file

    Return:
        error if some of the checks are not ok.
            ERROR if some of library preparation id does not run_exists
            incompatible_s_p_end if single reads and paired end are both in the selected library preparation ids
            incompatible_index with tuple of sample name and index values
            not_defined_lib_prep with error message in case that when reaching this steps they were not available

        protocol_parameter_list.

    """
    error = {}
    lib_prep_ids = []
    excel_json_data = json.loads(form_data["table_data"])
    heading_with_hidden = wetlab.config.HEADING_FOR_DISPLAY_SAMPLES_IN_POOL.copy()
    heading_with_hidden.append("lib_prep_id")
    c_data = core.utils.common.jspreadsheet_to_dict(
        heading_with_hidden, excel_json_data
    )
    for row in c_data:
        if row["Include in Pool"] == "True":
            lib_prep_ids.append(row["lib_prep_id"])

    if len(lib_prep_ids) == 0:
        error["ERROR"] = wetlab.config.ERROR_NOT_LIBRARY_PREPARATION_SELECTED
        return error
    duplicated_index = check_if_duplicated_index(lib_prep_ids)
    if "incompatible_index" in duplicated_index:
        error["duplicated_index"] = duplicated_index
        return error
    # check that all library preparation use the same adapter
    adapters = get_lib_prep_adapter(lib_prep_ids)
    if len(adapters) > 1:
        error["multiple_adapters"] = adapters
        return error

    pool_data = {}
    pool_data["poolName"] = form_data["poolName"]
    pool_data["platform"] = form_data["platform"]
    pool_data["poolCodeID"] = generate_pool_code_id()
    pool_data["registerUser"] = user_obj
    pool_data["adapter"] = adapters[0]
    pool_data["pairedEnd"] = get_single_paired(lib_prep_ids)
    pool_data["n_samples"] = len(lib_prep_ids)

    new_pool = wetlab.models.LibraryPool.objects.create_lib_pool(pool_data)
    # update pool_id in each library_preparation belongs the new pool
    for lib_prep in lib_prep_ids:
        lib_prep_obj = wetlab.models.LibPrepare.objects.get(pk__exact=lib_prep)
        lib_prep_obj.set_pool(new_pool)
        # return back to state completed if a library prepatarion was used for reused pool
        if lib_prep_obj.get_state() == "Reused pool":
            lib_prep_obj.set_state("Completed")

    # update the number of samples
    new_pool.set_pool_state("Selected")

    return new_pool


def generate_pool_code_id():
    """
    Description:
        The function get the date of today ( in format yyyy_mm_dd) and check if exists
        a pool_code_id with the same value.
        If yes the function increment the last digits else is set to "_1"
    Return:
        today_date + subindex
    """
    today_date = date.today().strftime("%Y_%m_%d")
    if wetlab.models.LibraryPool.objects.filter(
        pool_code_id__icontains=today_date
    ).exists():
        lib_pool = wetlab.models.LibraryPool.objects.filter(
            pool_code_id__icontains=today_date
        ).last()
        latest_code_id = lib_pool.get_pool_code_id()
        last_seq_number = int(latest_code_id.split("_")[-1])
        return str(today_date + "_" + str(last_seq_number + 1))
    else:
        return str(today_date + "_1")


def get_lib_prep_to_select_in_pool():
    """
    Description:
        The function get the library preparation instances that do not have linked to a pool and
        collect the information to display.
    Constants:
        HEADING_FOR_DISPLAY_SAMPLES_IN_POOL
    Return:
        display_list
    """
    display_list = {}
    display_list["data"] = {}

    if wetlab.models.LibPrepare.objects.filter(
        Q(lib_prep_state__lib_prep_state__exact="Completed", pools=None)
        | Q(lib_prep_state__lib_prep_state__exact="Reused pool")
    ).exists():
        lib_preparations = wetlab.models.LibPrepare.objects.filter(
            Q(lib_prep_state__lib_prep_state__exact="Completed", pools=None)
            | Q(lib_prep_state__lib_prep_state__exact="Reused pool")
        ).order_by("register_user")
        for lib_prep in lib_preparations:
            user_sample_obj = lib_prep.get_user_sample_sheet_obj()
            platform = user_sample_obj.get_sequencing_configuration_platform()
            if platform not in display_list["data"]:
                display_list["data"][platform] = []
            display_list["data"][platform].append(
                lib_prep.get_info_for_selection_in_pool()
            )

        display_list["heading"] = wetlab.config.HEADING_FOR_DISPLAY_SAMPLES_IN_POOL

    return display_list


def get_info_to_display_created_pool(pool_obj):
    information_for_created_pool = {}
    information_for_created_pool["data"] = pool_obj.get_info()
    information_for_created_pool["pool_name"] = pool_obj.get_pool_name()
    information_for_created_pool[
        "heading_pool"
    ] = wetlab.config.HEADING_FOR_DISPLAY_CREATED_POOL
    lib_prep_data = []
    if wetlab.models.LibPrepare.objects.filter(pools=pool_obj).exists():
        lib_prep_ids = wetlab.models.LibPrepare.objects.filter(pools=pool_obj)
        for lib_prep_obj in lib_prep_ids:
            lib_prep_data.append(lib_prep_obj.get_info_for_display_pool())
    information_for_created_pool["lib_prep_data"] = lib_prep_data
    information_for_created_pool[
        "heading_library_pool"
    ] = wetlab.config.HEADING_FOR_DISPLAY_LIB_PREP_IN_POOL

    return information_for_created_pool


def get_pool_obj_from_id(pool_id):
    """_summary_

    Args:
        pool_id (string): pool primary key number

    Returns:
        LibraryPool instance
    """
    return wetlab.models.LibraryPool.objects.filter(pk__exact=pool_id).last()


def get_sample_name_in_pools(pool_ids):
    """Get the list of the samples that were involved in the pools

    Args:
        pool_ids (list): pools primary key number

    Returns:
        list: which contains the list of samples in pools
    """
    return list(
        wetlab.models.LibPrepare.objects.filter(pools__pk__in=pool_ids).values_list(
            "sample_id__sample_name", flat=True
        )
    )
