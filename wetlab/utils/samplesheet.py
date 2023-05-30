#!/usr/bin/env python3
# coding: utf-8
# Generic imports
import codecs
import os
import re
import time

from django.conf import settings
from django.core.files.storage import FileSystemStorage

# Local imports
import wetlab.config
import wetlab.models


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
    lines = file_read.split("\n")
    data_section_found = False
    description_index = False
    userid_names, invalid_names = [], []
    for line in lines:
        line = line.rstrip()
        if line == "":
            continue
        if "[Data]" in line:
            data_section_found = True
            continue
        if data_section_found:
            line = line.split(",")
            if not description_index:
                try:
                    description_index = line.index("Description")
                    continue
                except:
                    users[
                        "ERROR"
                    ] = wetlab.config.ERROR_SAMPLE_SHEET_DOES_NOT_HAVE_DESCRIPTION_FIELD
                    return users
            else:
                try:
                    u_name = line[description_index]
                except:
                    continue
                if u_name in user_id_list:
                    userid_names.append(u_name)
                else:
                    invalid_names.append(u_name)

    if len(invalid_names) > 0:
        invalid_names = list(set(invalid_names))
        invalid_names.insert(
            0, "".join(wetlab.config.ERROR_SAMPLE_SHEET_FOLLOWING_USER_ARE_NOT_DEFINED)
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
        except:
            return False
        return True
    return False


def get_adapters(file_lines):
    adapter1 = ""
    adapter2 = ""
    ## For accepting characters like spanish characters.
    for line in file_lines:
        if line == "":
            continue
        found_adapter = re.search("^Adapter", line)
        if found_adapter:
            adapter_code = line.split(",")[1]
            if adapter1 == "":
                adapter1 = adapter_code
                continue
            else:
                adapter2 = adapter_code
                break
        data_found = re.search("^\[Data]", line)
        if data_found:
            break

    return adapter1, adapter2


def get_index_adapter(file_lines):
    """
    Description :
        get the indexes adapter information from the file_lines
    Input:
        file_lines  # sample sheet file converted to list of lines
    Return:
        index_adapters . List of the reads
    """
    index_adapters = ""
    for line in file_lines:
        if "Index Adapters" in line:
            index_adapters = line.split(",")[1].replace('"', "")
            break
    # remove the \n or \r
    return index_adapters.rstrip()


def get_projects_in_sample_sheet(file_lines):
    """
    Description :
        get the project names defined from the file_lines
    Input:
        file_lines  # sample sheet file converted to list of lines
    Return:
        project_name_list
    """
    header_found = False
    project_list = []
    for line in file_lines:
        line = line.rstrip()
        if line == "":
            continue
        found_header = re.search("^Sample_ID,Sample_Name", line)
        if found_header:
            header_found = True
            heading = line.split(",")
            index_project_name = heading.index("Sample_Project")
            continue
        if header_found:
            try:
                project_list.append(line.split(",")[index_project_name])
            except:
                continue
    project_name_list = list(set(project_list))
    return project_name_list


def get_reads(file_lines):
    """
    Description :
        get the reads information from the file_lines
    Input:
        file_lines  # sample sheet file converted to list of lines
    Return:
        reads . List of the reads
    """
    read_found = False
    reads = []
    for line in file_lines:
        if "[Reads]" in line:
            read_found = True
            continue
        if "[Settings]" in line:
            break
        if read_found and re.search("^\w+", line):
            reads.append(line.split(",")[0])
    return reads


def get_samples_in_sample_sheet(file_lines):
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
    samples_dict["samples"] = []
    samples_dict["sample_data"] = []
    header_found = False
    for line in file_lines:
        line = line.rstrip()
        if line == "":
            continue
        found_header = re.search("^Sample_ID,Sample_Name", line)
        if found_header:
            header_found = True
            samples_dict["heading"] = line.split(",")
            index_sample_name = samples_dict["heading"].index("Sample_Name")
            continue
            ## found the index for projects
        if header_found:
            line_split = line.split(",")
            try:
                sample_name = line_split[index_sample_name].strip()
            except:
                continue
            if sample_name == "":
                continue
            data = []
            for item in line_split:
                data.append(item.strip())
            samples_dict["sample_data"].append(data)
            samples_dict["samples"].append(sample_name)
    return samples_dict


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
    # initialize data with empty values
    for item in wetlab.config.FIELDS_IN_SAMPLE_SHEET_HEADER_IEM_VERSION_5:
        sample_sheet_data[item.lower()] = ""

    file_lines = file_read.split("\n")
    for line in file_lines:
        if "IEMFileVersion" in line:
            sample_sheet_data["iem_version"] = line.split(",")[1]
            break
    # get assay information
    for line in file_lines:
        if "Assay" in line:
            sample_sheet_data["assay"] = line.split(",")[1]
            break
    # get index adapters information
    # for line in file_lines :
    #    if 'Index Adapters' in line :
    #        sample_sheet_data['index_adapters'] = line.split(',')[1]
    #        break
    # get application information
    for line in file_lines:
        if "Application" in line:
            sample_sheet_data["application"] = line.split(",")[1]
            break
    # get Instrument information
    for line in file_lines:
        if "Instrument Type" in line:
            sample_sheet_data["instrument type"] = line.split(",")[1]
            break
    # get adapters information
    sample_sheet_data["adapter1"], sample_sheet_data["adapter2"] = get_adapters(
        file_lines
    )
    # get indexes adapters information
    sample_sheet_data["index_adapters"] = get_index_adapter(file_lines)
    # get reads information
    sample_sheet_data["reads"] = get_reads(file_lines)
    # get proyects in sheet_data
    sample_sheet_data["proyects"] = get_projects_in_sample_sheet(file_lines)
    # update sample sheet data
    sample_sheet_data.update(get_samples_in_sample_sheet(file_lines))

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
    header_found = False
    full_path = os.path.join(settings.MEDIA_ROOT, sample_sheet_path)
    fh = open(full_path, "r")
    for line in fh:
        line = line.rstrip()
        if line == "":
            continue
        found_header = re.search("^Sample_ID,Sample_Name", line)
        if found_header:
            header_found = True
            sample_sheet_heading = line.split(",")
            index_sample_name = sample_sheet_heading.index("Sample_Name")
            try:
                index_description = sample_sheet_heading.index("Description")
            except:
                index_description = None
            continue
        if header_found:
            line_split = line.split(",")
            try:
                sample_name = line_split[index_sample_name].strip()
                if index_description is not None:
                    user_id = line_split[index_description].strip()
                else:
                    user_id = None
            except:
                continue
            if sample_name == "":
                continue

            sample_user[sample_name] = user_id
    fh.close()
    return sample_user


def get_projects_in_run(in_file):
    header_found = 0
    projects = {}
    fh = open(in_file, "r")
    for line in fh:
        line = line.rstrip()
        if line == "":
            continue
        found_header = re.search("^Sample_ID,Sample_Name", line)
        if found_header:
            header_found = 1
            ## found the index for projects
            p_index = line.split(",").index("Sample_Project")
            description_index = line.split(",").index("Description")
            continue
        if header_found:
            ### ignore the empty lines separated by commas
            valid_line = re.search("^\w+", line)
            if not valid_line:
                continue
            ## store the project name and the user name (Description) inside projects dict
            projects[line.split(",")[p_index]] = line.split(",")[description_index]
    fh.close()
    return projects


def get_index_library_name(in_file):
    """
    Description:
        The function get the index library adapters. It searchs in the  assay value (used for version 4 of IEM sample sheet
        and in hte Index Adapters on sample sheet version 5.
        If Index adapters is found they are used if not the assay value
    Input:
        in_file     # shample sheet file
    Output:
        library_value
    """
    library_value = ""
    ## For accepting characters like spanish characters.
    import codecs

    fh = codecs.open(in_file, "r", "utf-8")
    for line in fh:
        line = line.rstrip()
        if line == "":
            continue
        found_assay = re.search("^Assay", line)
        if found_assay:
            library_value = line.split(",")[1]
            continue
        found_adapters = re.search("^Index Adapters", line)
        if found_adapters:
            library_value = line.split(",")[1]
            break

    fh.close()
    return library_value


def update_library_kit_field(library_file_name, library_kit_name, library_name):
    # result_directory='documents/wetlab/BaseSpaceMigrationFiles/'
    timestr = time.strftime("%Y%m%d-%H%M%S")
    tmp = re.search("(.*)\d{8}-\d+.*\.csv", library_file_name)
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
    except:
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


def update_sample_sheet(in_file, experiment_name):
    out_line = str("Experiment Name," + experiment_name + "\n")
    fh_in = open(in_file, "r")
    temp = os.path.join(settings.MEDIA_ROOT, "wetlab", "tmp.txt")
    fh_out = open(temp, "w")
    experiment_line_found = False
    for line in fh_in:
        # find experiment name line
        found_experiment = re.search("^Experiment Name", line)
        if found_experiment:
            fh_out.write(out_line)
            experiment_line_found = True

        elif line == "\n" and experiment_line_found == False:
            fh_out.write(out_line)
            fh_out.write("\n")
            experiment_line_found = True
        else:
            fh_out.write(line)

    fh_in.close()
    fh_out.close()
    os.rename(temp, in_file)


def create_unique_sample_id_values(infile, index_file):
    found_sample_line = False

    fh = open(infile, "r")
    temp_sample_sheet = os.path.join(settings.MEDIA_ROOT, "wetlab", "tmp_file")
    fh_out_file = open(temp_sample_sheet, "w")
    with open(index_file) as fh_index:
        index = fh_index.readline()
        index = index.rstrip()
        index_number_str, index_letter = index.split("-")
        fh_index.close()

    for line in fh:
        if "Sample_ID" in line:
            found_sample_line = True
            fh_out_file.write(line)
            continue
        if found_sample_line:
            # discard the empty lines or the lines that contains empty lines separated by comma
            if line == "\n" or re.search("^\W", line):
                continue

            data_line = line.split(",")
            data_line[0] = str(index_number_str + "-" + index_letter)
            new_line = ",".join(data_line)
            fh_out_file.write(new_line)
            # increase the index number
            index_number = int(index_number_str) + 1
            if index_number > 9999:
                index_number = 0
                # increase a letter
                split_index_letter = list(index_letter)
                if split_index_letter[1] == "Z":
                    last_letter = chr(ord(split_index_letter[0]) + 1)
                    split_index_letter[0] = last_letter
                    split_index_letter[1] = "A"
                    index_letter = "".join(split_index_letter)
                else:
                    first_letter = chr(ord(split_index_letter[1]) + 1)
                    split_index_letter[1] = first_letter
                    index_letter = "".join(split_index_letter)

            index_number_str = str(index_number)
            index_number_str = index_number_str.zfill(4)

        else:
            fh_out_file.write(line)

    # dump the index value to file
    fh_index = open(index_file, "w")
    index_line = str(index_number_str + "-" + index_letter)
    fh_index.write(index_line)
    fh_index.close()
    fh.close()
    fh_out_file.close()
    os.rename(temp_sample_sheet, infile)


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
        True
    """
    found_sample_line = False
    temp_sample_sheet = os.path.join(settings.MEDIA_ROOT, "wetlab", "tmp_file")
    fh = open(in_file, "r")
    fh_out_file = open(temp_sample_sheet, "w")
    for line in fh:
        if "Sample_ID" in line:
            found_sample_line = True
            line = line.rstrip()
            project_index = line.split(",").index("Sample_Project")
            description_index = line.split(",").index("Description")
            fh_out_file.write(str(line + "\n"))
            continue
        if found_sample_line:
            # discard the empty lines or the lines that contains empty lines separated by comma
            if line == "\n" or re.search("^\W", line):
                continue

            data_line = line.split(",")
            project_name = data_line[project_index]
            data_line[description_index] = user_names[project_name]

            new_line = ",".join(data_line)
            fh_out_file.write(str(new_line + "\n"))

        else:
            fh_out_file.write(line)
    fh.close()
    fh_out_file.close()
    os.rename(temp_sample_sheet, in_file)
    return True


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
    # create thd directory if not exists
    template_dir = os.path.join(
        settings.MEDIA_ROOT, wetlab.config.LIBRARY_PREPARATION_SAMPLE_SHEET_DIRECTORY
    )
    if not os.path.exists(template_dir):
        os.makedirs(template_dir)

    filename, file_extension = os.path.splitext(user_input_file.name)
    fs = FileSystemStorage()
    timestr = time.strftime("%Y%m%d-%H%M%S")
    ## including the timestamp to the sample sheet file
    file_name = str(
        wetlab.config.LIBRARY_PREPARATION_SAMPLE_SHEET_DIRECTORY
        + filename
        + "_"
        + timestr
        + file_extension
    )
    filename = fs.save(file_name, user_input_file)

    ### add the document directory to the input file
    stored_path_file = os.path.join(settings.MEDIA_ROOT, file_name)
    return stored_path_file, file_name


def read_all_lines_in_sample_sheet(sample_sheet):
    """
    Description:
        The function reads the input file and return the content in a variable
    Input:
        sample_sheet    # location of sample sheet
    Return:
        read_lines
    """
    read_lines = []
    if os.path.exists(sample_sheet):
        fh = open(sample_sheet, "r")
        read_lines = fh.readlines()
        fh.close()
    return read_lines


def read_user_iem_file(in_file):
    """
    Description:
        The function reads the input file and return the content in a variable

    Return:
        file_read
    """
    fh = codecs.open(in_file, "r", "utf-8")
    try:
        file_read = fh.read()
        fh.close()
        return file_read
    except:
        fh.close()
        return False


def valid_user_iem_file(file_read):
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
    for section in wetlab.config.SECTIONS_IN_IEM_SAMPLE_SHEET:
        if section not in file_read:
            return False

    lines = file_read.split("\n")
    data_section_found = False
    data_field_length = ""
    sample_number = 0
    for line in lines:
        line = line.rstrip()
        if line == "":
            continue
        if "[Data]" in line:
            data_section_found = True
            continue
        if data_section_found:
            if data_field_length == "":
                data_field_length = len(line.split(","))
                continue
            line_split = line.split(",")
            # Allow at this step the Description field can be empth
            if (
                len(line_split) < data_field_length - 1
                or len(line_split) > data_field_length
            ):
                return False
            sample_number += 1
    if sample_number == 0:
        return False
    return True
