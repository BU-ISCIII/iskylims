import re, os, sys, codecs
import datetime, time
import itertools
from django.core.files.storage import FileSystemStorage
from iSkyLIMS_wetlab.wetlab_config import *
from iSkyLIMS_wetlab.models import *

def check_collection_index_file_format (input_file):
    '''
    Description:
        The function will check if the input file contains the heading
        described on the declared constant COLLECTION_INDEX_HEADING
        defined in wetlab_config
    Input:
        input_file     # contains the file name
    Variables:
        read_data # has the content of the file
        item_to_check   # interaction variable of the  COLLECTION_INDEX_HEADING
    Constants:
        COLLECTION_INDEX_HEADING # List containing the heading colunms to
                                be checked
    Return:
        False if heading was not found
        True if heading was found.

    '''
    with open (input_file , encoding="utf-8" ) as fh:
       read_data = fh.read()
    for item_to_check in wetlab_config.COLLECTION_INDEX_HEADING :
        if item_to_check not in read_data:
            return False
    return True

def check_collection_index_exists(collection_name):
    '''
    Description:
        Check if collection name exist on Database
    Input:
        collection_name # name to check
    Return:
        True if exists on database , else if not
    '''
    if CollectionIndexKit.objects.filter( collectionIndexName__exact = collection_name):
        return True
    return False

def get_collection_index_name (input_file):
    '''
    Description:
        The function will read the input file and it will return
        the index library name
    Input:
        input_file     # contains the file name
    Variables:
        found_library_name # has the content of the file

        library_name    # will have the name of the library fetched from file
    Return:
        library_name with the name of the library or empty '' if not found
    '''
    found_name = False
    collection_name = ''
    with open (input_file , encoding="utf-8" ) as fh:
        for line in fh:
            found_collection_name = re.search('^\[Name\]',line)
            if found_collection_name :
                found_name = True
                continue
            if found_name:
                collection_name = line.rstrip()
                break
    return collection_name


def get_collection_settings (input_file):
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
    collection_index_settings ={}
    adapter_list = []
    found_version = False
    with open (input_file, encoding="utf-8") as fh:
        for line in fh:
            line=line.rstrip()
            if line == "":
                continue
            if '[Version]' in line :
               found_version = True
               continue
            if found_version:
                collection_index_settings['version'] = line
                found_version = False
                continue
            if '[Name]' in line :
                found_name = True
                continue
            if found_name:
                collection_index_settings['name'] = line
                found_name = False
                continue
            if '[PlateExtension]' in line:
                found_extension = True
                continue
            if found_extension :
                collection_index_settings['plate_extension'] = line
                found_extension = False
                continue
            if '[Settings]' in line:
                found_settings = True
                continue
            if found_settings :
                line_split= line.split('\t')
                # No more adapters  have been found. copy the adapters in library_settings
                # dictionary and exit the loop

                if 'Adapter' in line:
                    adapter_list.append(line_split[-1])
                if (len (line_split) == 1 and len(adapter_list) > 0):

                    collection_index_settings['adapters'] = adapter_list
                    break
                else:
                    continue
    return collection_index_settings


def get_collection_index_information (collection_index_id) :

    collection_index_dict ={}
    collection_index_obj = CollectionIndexKit.objects.get(pk=collection_index_id)

    general_information = [collection_index_obj.get_collection_index_information()]
    collection_index_dict['general_information'] = general_information

    if CollectionIndexValues.objects.filter(collectionIndexKit_id = collection_index_id).exists():
        index_list = CollectionIndexValues.objects.filter(collectionIndexKit_id = collection_index_id)

        I7_indexes, I5_indexes = [], []
        index_row_info = []
        for index in index_list :
            i_values = index.get_index_value_information().split(';')
            index_row_info.append(i_values)
            # get all I7 index defined on the library
            I7_indexes.append([i_values[1],i_values[2]])
            if i_values[3] != '':
                #get all I5 index defined on the library
                I5_indexes.append([i_values[3],i_values[4]])
        I7_indexes.sort()
        collection_index_dict['I7_indexes'] = list(I7_indexes for I7_indexes,_ in itertools.groupby(I7_indexes))
        if len(I5_indexes) > 0 :
            I5_indexes.sort()
            collection_index_dict['I5_indexes'] = list(I5_indexes for I5_indexes,_ in itertools.groupby(I5_indexes))
        collection_index_dict['default_layout'] = index_row_info
        return collection_index_dict
    else:
        return False

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
            line=line.rstrip()
            if line == "":
                continue
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
                    found_I5 = False
            if 'Layout' in line :
                if "SingleIndex" in line and len(index_5) > 0 :
                    if len(index_values) == 0 :
                        continue
                    else:
                        break
                if len(index_values) > 0 :
                    break
                else:
                    layout_found = True
                    continue
            if layout_found :
                if '[' in line: # reached the end of layout. Exiting the function
                    break
                line = line.rstrip('\n')
                line_split = line.split('\t')
                if len(index_5) : # There are 2 index (I7 and I5)
                    index_values.append([line_split[0], line_split[1], index_7[line_split[1]], line_split[2], index_5[line_split[2]]])
                else:
                    index_values.append([line_split[0], line_split[1], index_7[line_split[1]], '', ''])
    return index_values

def store_collection_kits_file(collection_file):
    '''
    Description:
        The function is store the collection kit file inside COLLECTION_INDEX_KITS_DIRECTORY

    Input:
        collection_file     # contains the file to be saved
    Variables:

    Constants:
        COLLECTION_INDEX_KITS_DIRECTORY

    Functions:
        in utils.library_kits :
            -- check_index_library_file_format(saved_file)     # for checking
                            the number of index in the input file
            -- getting_index_library_name(saved_file) # gets the library name
            -- get_library_settings(saved_file) # gets the settings values
                            from the input file
            -- get_index_values(saved_file)     # gets the index value

    Return:
         Return the different information depending on the execution:
        -- Error page in case of:
            -- Uploaded file is bigger than the LIBRARY_MAXIMUM_SIZE value
            -- file uploaded does not have the right format
            -- the library already exists.
        -- library_kit_information with :
            -- ['libraries']
            ---['new_library_kit'] in case a new library kit was added
    '''

    ## fetch the file from user form and  build the file name  including
    ## the date and time on now to store in database

    split_filename=re.search('(.*)(\.\w+$)',collection_file.name)
    f_name = split_filename[1]
    f_extension = split_filename[2]
    fs_index_lib = FileSystemStorage()
    timestr = time.strftime("%Y%m%d-%H%M%S")

    ## using the MEDIA_ROOT variable defined on settings to upload the file
    file_name=os.path.join(wetlab_config.COLLECTION_INDEX_KITS_DIRECTORY ,  str(f_name + '_' +timestr + f_extension))
    filename = fs_index_lib.save(file_name,  collection_file)
    saved_file = os.path.join(settings.MEDIA_ROOT, file_name)

    return saved_file

def store_collection_settings (collection_settings, file_name) :
    # saving library settings into database
    if len(collection_settings['adapters']) == 1:
        adapter_2 = ''
    else :
        adapter_2 = collection_settings['adapters'][1]
    new_collection_settings = CollectionIndexKit(collectionIndexName = collection_settings['name'],
                                version =  collection_settings ['version'],
                                plateExtension = collection_settings['plate_extension'] ,
                                adapter1 = collection_settings['adapters'][0],
                                adapter2 = adapter_2, collectionIndexFile = file_name)
    new_collection_settings.save()
    return new_collection_settings

def remove_boom_bytes (input_file):
    s = open(input_file, mode='r', encoding='utf-8-sig').read()
    open(input_file, mode='w', encoding='utf-8').write(s)


def store_collection_indexes (collection_index, new_collection_obj):
    # saving index values into database
    for row in collection_index :
        index_to_store = CollectionIndexValues(collectionIndexKit_id = new_collection_obj,
                                defaultWell = row[0], index_7 = row[1],
                                i_7_seq = row[2], index_5 = row[3],
                                i_5_seq = row[4])
        index_to_store.save()
