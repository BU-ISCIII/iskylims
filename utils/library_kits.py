#!/usr/bin/env python3
import re, os, sys, codecs
from ..wetlab_config import *



def check_index_library_file_format (input_file):
    '''
    Description:
        The function will check if the input file contains the heading
        described on the declared constant INDEX_LIBRARY_HEADING
        defined in wetlab_config
    Input:
        input_file     # contains the file name
    Variables:
        read_data # has the content of the file
        item_to_check   # interaction variable of the  INDEX_LIBRARY_HEADING
    Constants:
        INDEX_LIBRARY_HEADING # List containing the heading colunms to
                                be checked
    Return:
        False if heading was not found
        True if heading was found.

    '''
    with open (input_file , encoding="utf-8" ) as fh:
       read_data = fh.read()
    for item_to_check in INDEX_LIBRARY_HEADING :
        if item_to_check not in read_data:
            return False
    return True



def getting_index_library_name (input_file):
    '''
    Description:
        The function will read the input file and it will return
        the index library name
    Input:
        input_file     # contains the file name
    Variables:
        found_library_name # has the content of the file
        found_name   # interaction variable of the  INDEX_LIBRARY_HEADING
        library_name    # will have the name of the library fetched from file
    Return:
        library_name with the name of the library or empty '' if not found
    '''
    found_name = False
    library_name = ''
    with open (input_file , encoding="utf-8" ) as fh:
        for line in fh:
            found_library_name = re.search('^\[Name\]',line)
            if found_library_name :
                found_name = True
                continue
            if found_name:
                library_name= line.rstrip()
                break
    return library_name


def get_library_settings (input_file):
    '''
    Description:
        The function will read the input file and it will return
        library_settings with the information content file
    Input:
        input_file     # contains the file name
    Variables:
        library_settings    # dictionary to store the setting information
        found_version, found_name, found_extension and found_settings # are
                temporary created to read the next following line that
                contains the information that we are looking for
    Return:
        library_settings with the settings of the library or empty '' if not found
    '''
    library_settings ={}
    adapter_list = []
    with open (input_file, encoding="utf-8") as fh:
        for line in fh:
            if '[Version]' in line :
               found_version = True
               continue
            if found_version:
                library_settings['version'] = line.rstrip()
                found_version = False
                continue
            if '[Name]' in line :
                found_name = True
                continue
            if found_name:
                library_settings['name'] = line.rstrip()
                found_name = False
                continue
            if '[PlateExtension]' in line:
                found_extension = True
                continue
            if found_extension :
                library_settings['plate_extension'] = line.rstrip()
                found_extension = False
                continue
            if '[Settings]' in line:
                found_settings = True
                continue
            if found_settings :
                line=line.rstrip()
                line_split= line.split('\t')
                # No more adapters  have been found. copy the adapters in library_settings
                # dictionary and exit the loop
                if len (line_split) == 1:
                    library_settings['adapters'] = adapter_list
                    break
                else:
                    adapter_list.append(line_split[-1])
    return library_settings

def get_index_values (input_file):
    '''
    Description:
        The function will read the index value in the input file and it will return
        library_settings with the information content file
    Input:
        input_file     # contains the file name
    Variables:
        index_values # dictionary to store index_I7 and and index_I5
        index_7     # will store the index 7 values
        index_5     # will store the index 5 values
        found_I7 and found_I5 # are set to True to read the index
                once the I7/I5 heading is found
    Return:
        index_values with the I7 and I5 index found in the file
            Empty if not found
    '''
    index_7 , index_5 = {} , {}
    found_I7 = False
    found_I5 = False
    layout_found = False
    index_values = []
    with open (input_file , encoding="utf-8") as fh :
        for line in fh:
            if '[I7]' in line :
                found_I7 = True
                continue
            if found_I7 :
                line=line.rstrip()
                if len (line.split('\t'))> 1:
                    index_line = line.split('\t')
                    index_7[index_line[0]]= index_line[1]
                    continue
                # No more index have been found for I7. Start collecting index for I5
                else:
                    found_I7 = False
                    if '[I5]' in line:
                        found_I5 = True
                        continue
            if found_I5 :
                line=line.rstrip()
                if len (line.split('\t')) > 1:
                    index_line = line.split('\t')
                    index_5[index_line[0]]= index_line[1]
                    continue
                # No more index have been found for I5. Exit the loop
                else:
                    break
            if 'Layout' in line :
                if "SingleIndex" in line and len(I5) > 0 :
                    continue
                else:
                    layout_found = True
                    continue
            if layout_found :
                if '[' in line: # reached the end of layout. Exiting the function
                    break
                line_split = line.split('\t')
                if len(line_split) == 3 : # There are 2 index (I7 and I5)
                    index_values.append([line_split[0], line_split[1], index_7[line_split[1]], line_split[2], index_5[line_split[2]]])
                else:
                    index_values.append([line_split[0], line_split[1], index_7[line_split[1]], '', ''])

    return index_values

def remove_boom_bytes (input_file):
    s = open(input_file, mode='r', encoding='utf-8-sig').read()
    open(input_file, mode='w', encoding='utf-8').write(s)
