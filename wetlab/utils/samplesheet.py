#!/usr/bin/env python3
# coding: utf-8
# Generic imports
import os
import re
import time

from django.conf import settings
from django.core.files.storage import FileSystemStorage

# Local imports
import wetlab.config
import wetlab.utils.common

# import wetlab.models


def samplesheet_version(samplesheet: dict) -> str:
    """
    Return samplesheet version (1 or 2) as a string

    Args:
        samplesheet (dict): samplesheet

    Returns:
        str: "1" or "2"
    """
    return samplesheet.get("Header", {}).get("FileFormatVersion", "1")


def read_file_from_path(file_path: str) -> str:
    """
    Read file from path, ensuring characters are not lost due to encoding

    :param file_path: Description
    :type file_path: str
    """
    try:
        with open(file_path, "r", encoding="utf-8-sig") as f:
            read_file = f.read()
    except Exception:
        return False
    return read_file


def write_samplesheet_to_path(samplesheet: dict, path: str) -> bool:
    """
    Write the samplesheet to a path.

    Args:
        samplesheet (dict): Sample sheet (As a dictionary)
        path (str): Path to write the file to

    Returns:
        bool: True if no issue, False if any exception is raised
    """

    string_to_write = ""
    for key, value in samplesheet.items():
        string_to_write += f"[{key}]\n"
        if (
            key in wetlab.config.TABULAR_DATA_SECTIONS_SAMPLE_SHEET.values()
            or key.endswith("_Data")
        ):
            delimiter = "\n"
            string_to_write += (
                f"{delimiter.join([','.join(row) for row in value])}{delimiter}"
            )
        else:
            for row_header, element in value.items():
                string_to_write += f"{row_header},{element}\n"
    try:
        with open(path, "w") as f:
            f.write(string_to_write)
        return True
    except Exception:
        return False


def file_read_to_dictionary(
    file_read: str,
) -> dict[str, dict[str, any] | list[list[str]]]:
    """
    Description:
        This function transforms a SampleSheet from string to dictionary,
        loading the headers as keys. If data is fully tabular (e.g. [Data]),
        load rows as list of list; if not, nested dict.
    Input:
        file_read           # content of the IEM file from user
    Constant:
        ERROR_SAMPLE_SHEET_HAS_INVALID_LINES
    Return
        ERROR if any line was not in the proper format (tabular, comma-delimited)
        samplesheet as a dictionary
    """
    samplesheet = {}
    file_read = file_read.splitlines()
    for line in file_read:
        if line.startswith("["):
            section = line.split(",")[0].lstrip("[").rstrip("]")
            samplesheet[section] = {}
            continue
        if not line or re.match("^,+$", line):
            # Empty/Filler/Artifact lines with no info
            continue
        if (
            section in wetlab.config.TABULAR_DATA_SECTIONS_SAMPLE_SHEET.values()
            or section.endswith("_Data")
        ):
            # Data is tabular; append as list of lists until next header
            if not isinstance(samplesheet[section], list):
                samplesheet[section] = []
            samplesheet[section].append(line)
        else:
            # Data is processed as key: value before adding it to dictionary
            line = line.split(",")
            try:
                # Sometimes there are empty lines (e.g. READS values)
                samplesheet[section][line[0].strip()] = (
                    ",".join(line[1:]).strip() if len(line) > 1 else ""
                )
            except IndexError:
                return {"ERROR": wetlab.config.ERROR_SAMPLE_SHEET_HAS_INVALID_LINES}
    return samplesheet


def validate_userid_in_user_iem_file(file_read, user_id_list):
    """
    Description:
        The function get if userids included in the user IEM file
    Input:
        file_read           # content of the IEM file from user
        user_id_list        # UserID list
    Constant:
        ERROR_SAMPLE_SHEET_DOES_NOT_HAVE_DESCRIPTION_FIELD
    Return
        ERROR if userid is not defined in descripion column
        userids with the ids found
    """
    users = {}
    samplesheet = file_read_to_dictionary(file_read)
    if "ERROR" in samplesheet:
        return samplesheet

    data = get_tabular_data(samplesheet)
    iskylims_user_column = wetlab.config.TABULAR_DATA_ISKYLIMS_USER_COLUMN.get(
        samplesheet_version(samplesheet)
    )
    if not data:
        users["ERROR"] = (
            wetlab.config.ERROR_SAMPLE_SHEET_DOES_NOT_HAVE_DESCRIPTION_FIELD
        )
        return users
    users_in_sample_sheet = get_user_ids_from_samplesheet(samplesheet, data)
    if not users_in_sample_sheet and iskylims_user_column:
        # Keep compatibility with old behavior for v1 sheets.
        users_in_sample_sheet = get_column_from_tabular_data(data, iskylims_user_column)

    if not users_in_sample_sheet:
        users["ERROR"] = (
            wetlab.config.ERROR_SAMPLE_SHEET_DOES_NOT_HAVE_DESCRIPTION_FIELD
        )
        return users

    userid_names = [user for user in users_in_sample_sheet if user in user_id_list]
    invalid_names = [user for user in users_in_sample_sheet if user not in user_id_list]

    if len(invalid_names) > 0:
        invalid_names = list(set(invalid_names))
        invalid_names.insert(
            0, "".join(wetlab.config.ERROR_SAMPLE_SHEET_USER_ARE_NOT_DEFINED)
        )
        users["ERROR"] = invalid_names
        return users

    users["user_ids"] = list(set(userid_names))
    return users


def delete_stored_file(input_file):
    """
    Description:
        The function delete the requested file
    Input:
        input_file  # input file to delete
    Return
        True
    """
    if os.path.exists(input_file):
        try:
            os.remove(input_file)
        except Exception:
            return False
        return True
    return False


def get_adapters(samplesheet: dict) -> tuple[str, str]:
    """
    Get the adapters from the samplesheet (v1 | v2)

    Parameters
    ----------
    samplesheet : dict
        User sample sheet read as a dictionary

    Return
    ------
    adapters : tuple
        tuple containing the values in order for adapter 1 and 2
    """
    for settings_name in wetlab.config.SETTINGS_SECTIONS_SAMPLE_SHEET:
        settings = samplesheet.get(settings_name)
        if not settings:
            continue
        for i in range(len(wetlab.config.ADAPTER_1_FIELD_NAMES)):
            try:
                adapter1 = samplesheet[wetlab.config.ADAPTER_1_FIELD_NAMES[i]]
                adapter2 = samplesheet[wetlab.config.ADAPTER_2_FIELD_NAMES[i]]
            except ValueError:
                continue
            break
        else:
            continue
        break
    else:
        return "", ""

    return adapter1, adapter2


def get_row_from_tabular_data(tabular_data: list[list], index) -> list:
    """
    Get the data headers

    Parameters
    ----------
    file_lines : list
        List of lines from sample sheet

    Return
    ------
    heading : list
        contain the values from heading row of file
    """
    try:
        return tabular_data[index]
    except IndexError:
        return []


def get_column_from_tabular_data(
    tabular_data: list[list], column_name: str, unique: bool = False
) -> list:
    """
    Docstring for get_all_values_from_tabular_data

    :param tabular_data: Description
    :type tabular_data: dict[list[list[str]]]
    """
    if not tabular_data:
        return []

    headers = tabular_data[0]
    tabular_data = [row for row in tabular_data[1:]]
    try:
        column_index = headers.index(column_name)
    except ValueError:
        return []

    values = [row[column_index] for row in tabular_data]
    return list(set(values)) if unique else values


def get_tabular_data(samplesheet: dict) -> list[list]:
    """
    Get tabular data from a samplesheet.

    Args:
        samplesheet (dict): Samplesheet with the tabular data

    Returns:
        list[list]: Tabular data, in a nested list (Matrix N*M)
    """
    data_section = wetlab.config.TABULAR_DATA_SECTIONS_SAMPLE_SHEET.get(
        samplesheet_version(samplesheet), ""
    )
    data = samplesheet.get(data_section, [])
    return [row.split(",") for row in data]


def get_tabular_data_for_section(samplesheet: dict, section_name: str) -> list[list]:
    """
    Get tabular data for a specific section in a samplesheet.

    Args:
        samplesheet (dict): Parsed samplesheet data
        section_name (str): Section name to fetch

    Returns:
        list[list]: Tabular data, in a nested list (Matrix N*M)
    """
    data = samplesheet.get(section_name, [])
    if not isinstance(data, list):
        return []
    return [row.split(",") for row in data]


def get_user_id_from_project_name(project_name: str) -> str:
    """
    Extract iSkyLIMS username from the project name suffix.
    Example: MiSeq_i100_GEN_002_20260219_testuser1 -> testuser1
    """
    if not project_name or "_" not in project_name:
        return ""
    return project_name.rsplit("_", 1)[-1].strip()


def get_user_ids_from_project_name(samplesheet: dict) -> list[str]:
    """
    Extract iSkyLIMS usernames from Cloud_Data.ProjectName values.
    """
    cloud_data = get_tabular_data_for_section(samplesheet, "Cloud_Data")
    project_names = get_column_from_tabular_data(cloud_data, "ProjectName")
    return [
        get_user_id_from_project_name(project_name) for project_name in project_names
    ]


def get_user_ids_from_samplesheet(
    samplesheet: dict, data: list[list] | None = None
) -> list:
    """
    Extract user IDs from a samplesheet.
    Priority:
    - v1: Description/custom_description column in [Data]
    - v2: custom_description/Description in [BCLConvert_Data]
    - v2 fallback: suffix from Cloud_Data.ProjectName
    """
    if data is None:
        data = get_tabular_data(samplesheet)

    if not data:
        return []

    version = samplesheet_version(samplesheet)
    iskylims_user_column = wetlab.config.TABULAR_DATA_ISKYLIMS_USER_COLUMN.get(version)

    user_id_list_db = wetlab.utils.common.get_userid_list()

    user_ids = []
    if iskylims_user_column:
        user_ids = get_column_from_tabular_data(data, iskylims_user_column)

    # Some v2 sheets still use Description rather than custom_description.
    if not user_ids and version == "2":
        user_ids = get_column_from_tabular_data(data, "Description")

    # v2 fallback: parse user id from Cloud_Data.ProjectName suffix.
    if not user_ids and version == "2":
        user_ids = get_user_ids_from_project_name(samplesheet)

    # FIXME: Due to samplesheet limitations, we need to check all possible combinations of "dash to dot"
    # FIXME: When we develop a final solution, change below.
    userid_names = []

    for user in user_ids:
        all_user_replacement_combinations = (
            wetlab.utils.common.get_all_string_replacement_combinations(
                user, old="-", new="."
            )
        )
        for user_permutation in all_user_replacement_combinations:
            if user_permutation in user_id_list_db:
                userid_names.append(user_permutation)
                break

    return userid_names


def get_projects_in_sample_sheet(samplesheet) -> list:
    """
    Description :
        get the project names defined from the file_lines
    Input:
        file_lines  # sample sheet file converted to list of lines
    Return:
        project_name_list
    """
    data = get_tabular_data(samplesheet)
    projects = get_column_from_tabular_data(data, "Sample_Project", unique=True)
    return projects


def get_reads(samplesheet: dict) -> list[str, str]:
    """
    Description :
        get the reads information from the samplesheet
    Input:
        samplesheet  # samplesheet file converted to dictionary
    Return:
        reads . List of the reads
    """
    reads_section = samplesheet.get("Reads")
    reads = []
    for key, value in reads_section.items():
        if not value:
            # Sample sheets format for read is inconsistent - V1 accepts different format for read lengths
            reads = list(reads_section.keys())
            break
        if "Read" in key:
            reads.append(value)
    return reads


def get_samples_in_sample_sheet(samplesheet: dict) -> dict:
    """
    Description :
        get the sample information from the file_lines
    Input:
        file_lines  # sample sheet file converted to list of lines
    Return:
        samples_dict contains in 'samples' key the sample name.  'sample_data'
        is a list for each sample row
    """
    samples_dict = {}
    data = get_tabular_data(samplesheet)

    samples_dict["header"] = data[0] if len(data) > 1 else ""
    samples_dict["samples"] = get_column_from_tabular_data(data, "Sample_Name")
    samples_dict["sample_data"] = [
        get_row_from_tabular_data(data, i) for i in range(1, len(data))
    ]
    return samples_dict


def get_headers(samplesheet: dict) -> list:
    data = get_tabular_data(samplesheet)
    return data[0].split(",")


def get_sample_sheet_data(file_read):
    """
    Description:
        The function reads the user sample sheet from IEM and extracts : samples, adapters, reads
        assay, index adapters, application and instrument
    Input:
        file_read    # content of the user IEM
    Constants:
        FIELDS_IN_SAMPLE_SHEET_HEADER_IEM_VERSION_5
    Functions:
        get_adapters                 # located at this file
        get_reads                    # located at this file
        get_index_adapter            # located at this file
        get_samples_in_sample_sheet  # located at this file
    Return
        sample_sheet_data dictionary with the extracted information
    """
    sample_sheet_data = {}
    samplesheet = file_read_to_dictionary(file_read)
    # initialize header data
    for item in wetlab.config.FIELDS_IN_SAMPLE_SHEET_HEADER_IEM_VERSION_5:
        sample_sheet_data[item.lower()] = samplesheet["Header"].get(item, "")

    # get adapters information
    sample_sheet_data["adapter1"], sample_sheet_data["adapter2"] = get_adapters(
        samplesheet
    )
    # get indexes adapters information
    sample_sheet_data["index_adapters"] = samplesheet.get("Headers", {}).get(
        "Index Adapters", ""
    )
    # get reads information
    sample_sheet_data["reads"] = get_reads(samplesheet)
    # get proyects in sheet_data
    sample_sheet_data["projects"] = get_projects_in_sample_sheet(samplesheet)
    # update sample sheet data
    sample_sheet_data.update(get_samples_in_sample_sheet(samplesheet))
    # include heading
    sample_sheet_data["header"] = get_headers(samplesheet)
    return sample_sheet_data


def get_sample_with_user_owner(sample_sheet_path):
    """
    Description:
        The function fetch the sample sheet and return a dictionnary with sample
        name and the user ID owner of the sample
    Input:
        sample_sheet_path    # path of the stored sample sheet
    Return
        sample_user dictionary with the extracted information
    """
    sample_user = {}
    full_path = os.path.join(settings.MEDIA_ROOT, sample_sheet_path)
    file_read = read_file_from_path(full_path)
    samplesheet = file_read_to_dictionary(file_read)
    # Retrieve tabular data with the iskylims user in the header
    data = get_tabular_data(samplesheet)
    sample_names = get_column_from_tabular_data(data, "Sample_Name")
    if not sample_names:
        # v2 sheets frequently provide Sample_ID only.
        sample_names = get_column_from_tabular_data(data, "Sample_ID")
    user_ids = get_user_ids_from_samplesheet(samplesheet, data)

    min_len = min(len(sample_names), len(user_ids))
    sample_user = {sample_names[i]: user_ids[i] for i in range(min_len)}
    return sample_user


def get_projects_in_run(in_file: str) -> dict:
    """Funtion to check if the sample sheet has a valid header. On valid file
    get project names and the user names from description column

    Args:
        in_file (str): path to the sample sheet

    Returns:
        dict: dictionary with the projects or error message
    """
    file_read = read_file_from_path(in_file)
    samplesheet = file_read_to_dictionary(file_read)
    data = get_tabular_data(samplesheet)
    sample_projects = get_column_from_tabular_data(data, "Sample_Project")
    user_ids = get_user_ids_from_samplesheet(samplesheet, data)

    # v2 fallback: project names can live in Cloud_Data.ProjectName.
    if not sample_projects and samplesheet_version(samplesheet) == "2":
        cloud_data = get_tabular_data_for_section(samplesheet, "Cloud_Data")
        sample_projects = get_column_from_tabular_data(cloud_data, "ProjectName")
        if not user_ids:
            user_ids = [
                get_user_id_from_project_name(project_name)
                for project_name in sample_projects
            ]

    min_len = min(len(sample_projects), len(user_ids))
    projects = {sample_projects[i]: user_ids[i] for i in range(min_len)}

    if not data:
        return {"ERROR": wetlab.config.ERROR_SAMPLE_SHEET_HAS_INVALID_HEADING}
    if not projects:
        return {"ERROR": wetlab.config.ERROR_SAMPLE_SHEET_DOES_NOT_HAVE_PROJECTS}

    return projects


def get_index_library_name(in_file):
    """
    Description:
        The function get the index library adapters. It searchs in the  assay
        value (used for version 4 of IEM sample sheet
        and in the Index Adapters on sample sheet version 5.
        If Index adapters is found they are used if not the assay value
    Input:
        in_file     # shample sheet file
    Output:
        library_value
    """

    file_read = read_file_from_path(in_file)
    samplesheet = file_read_to_dictionary(file_read)
    header_values = samplesheet["Header"]
    library_value = header_values.get("Assay", "") or header_values.get(
        "Index Adapters", ""
    )

    return library_value


def update_library_kit_field(library_file_name, library_kit_name, library_name):
    # FIXME This function is not used anywhere - Not going to touch it for now
    # result_directory='documents/wetlab/BaseSpaceMigrationFiles/'
    timestr = time.strftime("%Y%m%d-%H%M%S")
    tmp = re.search(r"(.*)\d{8}-\d+.*\.csv", library_file_name)
    absolute_path = str(settings.BASE_DIR + "/")
    out_file = str(
        absolute_path
        + tmp.group(1)
        + timestr
        + "_for_basespace_"
        + library_kit_name
        + ".csv"
    )
    try:
        fh_in = open(library_file_name, "r")
        fh_out = open(out_file, "w")
    except Exception:
        return "ERROR:"

    for line in fh_in:
        found_library_kit = re.search("LibraryPrepKit", line)
        # library kit found. Replace line with new library kit
        if found_library_kit:
            line = str("LibraryPrepKit," + library_name + "\n")
        fh_out.write(line)
    fh_in.close()
    fh_out.close()
    os.remove(library_file_name)
    # remove absolute path from file_name
    absolute_path = str(settings.BASE_DIR + "/")
    file_name_in_database = out_file.replace(absolute_path, "")
    return file_name_in_database


def update_sample_sheet(in_file: str, experiment_name: str):
    """
    Update experiment name in samplesheet

    Args:
        in_file (str): Path to samplesheet
        experiment_name (str): experiment name
    """
    file_read = read_file_from_path(in_file)
    samplesheet = file_read_to_dictionary(file_read)
    samplesheet["Header"]["Experiment Name"] = experiment_name
    write_samplesheet_to_path(in_file)


def create_unique_sample_id_values(in_file: str, index_file: str):
    """
    Create unique sample IDs for samples in the samplesheets by using an ongoing index.

    Args:
        in_file (str): Path to samplesheet
        index_file (str): Path to index_file
    """
    file_read = read_file_from_path(in_file)
    samplesheet = file_read_to_dictionary(file_read)

    with open(index_file, "r") as f:
        try:  # catch OSError in case of a one line file
            f.seek(-2, os.SEEK_END)
            while f.read(1) != b"\n":
                f.seek(-2, os.SEEK_CUR)
        except OSError:
            f.seek(0)
        last_line = f.readline().decode()
        index_number_str, index_letter = last_line.rstrip().split("-")
        index_number = int(index_number_str)
    data = get_tabular_data(samplesheet)
    for row in data[1:]:
        index_number += 1
        index_number = (
            index_number % 10000
        )  # Return only 4 last digits, effectively restarting at 10000
        if index_number == 0:
            # When index re-starts, we move on to the next letter
            index_letter_parts = list(index_letter)
            # Reverse order to
            for i in [1, 0]:
                ascii_index_letter = ord(index_letter_parts[i])
                ascii_index_letter += 1
                ascii_index_letter = (
                    ascii_index_letter if ascii_index_letter <= 90 else 65
                )
                index_letter_parts[i] = chr(ascii_index_letter)
                if ascii_index_letter != 65:
                    break
            # Reverse and stringify the index letter parts
            index_letter = "".join(index_letter_parts)
        # Create unique sample ID and overwrite Sample_ID
        sample_unique_id = f"{str(index_number).zfill(4)}-{index_letter}"
        row[0] = sample_unique_id

    # Dump the index value to file
    with open(index_file, "w") as f:
        f.write(sample_unique_id)

    # Dump the updated samplesheet
    write_samplesheet_to_path(samplesheet, in_file)


def set_user_names_in_sample_sheet(in_file, user_names):
    """
    Description:
        The function modifies/set the user names in the description
        column
    Input:
        in_file # sample sheet file to be updated
        user_names # dictionary having projects as key and user names
                    as their value
    Variable:
        data_line  # split line into list to set user name
        description_index # column number where is located the description
                            inside sample Sheet
        found_sample_line # flag to identify if sample heading was found
        project_index # column number where is located the project inside
                        sample Sheet

        temp_sample_sheet # temporary sample sheet to store the information
                            it will replace the in_file
    Return:
        Bool: True if successful writing, False otherwise
    """
    file_read = read_file_from_path(in_file)
    samplesheet = file_read_to_dictionary(file_read)
    data = get_tabular_data(samplesheet)
    projects = get_column_from_tabular_data("Sample_Project")
    descriptions_index = data[0].index(
        wetlab.config.TABULAR_DATA_ISKYLIMS_USER_COLUMN.get(
            samplesheet_version(samplesheet)
        )
    )
    for i in range(1, len(data)):
        data[i][descriptions_index] = user_names[projects[i - 1]]
    success_writing = write_samplesheet_to_path(samplesheet, in_file)
    return success_writing


def store_user_input_file(user_input_file):
    """
    Description:
        The function rename the file name with the present time and it stores it in
        LIBRARY_PREPARATION_SAMPLE_SHEET_DIRECTORY
    Input:
        user_input_file  # input file from user
    Constant:
        LIBRARY_PREPARATION_SAMPLE_SHEET_DIRECTORY
    Return
        stored_path_file contains the full path of the file and file_name
    """
    # create the directory if not exists
    template_dir = os.path.join(
        settings.MEDIA_ROOT, wetlab.config.LIBRARY_PREPARATION_SAMPLE_SHEET_DIRECTORY
    )
    if not os.path.exists(template_dir):
        os.makedirs(template_dir)

    file_name, file_extension = os.path.splitext(user_input_file.name)
    fs = FileSystemStorage()
    timestr = time.strftime("%Y%m%d-%H%M%S")
    # including the timestamp to the sample sheet file
    file_name = str(
        wetlab.config.LIBRARY_PREPARATION_SAMPLE_SHEET_DIRECTORY
        + file_name
        + "_"
        + timestr
        + file_extension
    )
    file_name = fs.save(file_name, user_input_file)

    # add the document directory to the input file
    stored_path_file = os.path.join(settings.MEDIA_ROOT, file_name)
    return stored_path_file, file_name


def valid_user_iem_file(file_read: str) -> bool:
    """
    Description:
        The function check if the user IEM file has a valid format by checking the headings and
        if all fields are included in the data section. in particular the description field
        where username has to be defined to assing the sample to the user.
        If there is no sample function return False
    Input:
        file_read                           # content of the input file from user
    Functions:
        get_userid_list                 # located at utils.common.py file
    Constant:
        SECTIONS_IN_IEM_SAMPLE_SHEET
    Return
        False if file cannot be read or do not have all information
    """
    samplesheet = file_read_to_dictionary(file_read)
    for section_name in samplesheet.keys():
        if (
            section_name not in wetlab.config.SECTIONS_IN_IEM_SAMPLE_SHEET
            and section_name not in wetlab.config.SECTIONS_IN_V2_SAMPLE_SHEET
        ):
            return False

    data_field_length = ""
    sample_number = 0
    data = get_tabular_data(samplesheet)

    # Check on data
    if not data:
        return False

    # Checks on data contents
    data_field_length = len(data[0])
    sample_number = len(data) - 1
    if not all([len(row) - data_field_length for row in data]):
        # All rows contain the same either the same
        return False

    if sample_number == 0:
        return False
    return True
