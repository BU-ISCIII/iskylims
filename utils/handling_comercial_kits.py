import json
import re, os, sys, codecs
import itertools
from iSkyLIMS_core.core_config import *
from iSkyLIMS_core.models import *
from django.contrib.auth.models import User

def get_user_comercial_kits(register_user_obj, protocol_obj):
    user_kit_list = []

    if UserComercialKits.objects.filter(user = register_user_obj, basedComercial__protocol_id = protocol_obj.type).exists():
        user_kits =UserComercialKits.objects.filter(user = register_user_obj, basedComercial__protocol_id = protocol_obj.type)
        for user_kit in user_kits:
            user_kit_list.append(user_kit.get_nick_name)

    return user_kit_list

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
        COLLECTION_INDEX_HEADING # List containing the heading colunms to
                                be checked
    Return:
        False if heading was not found
        True if heading was found.

    '''
    with open (input_file , encoding="utf-8" ) as fh:
       read_data = fh.read()
    for item_to_check in COLLECTION_INDEX_HEADING :
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
    with open (input_file, encoding="utf-8") as fh:
        for line in fh:
            if '[Version]' in line :
               found_version = True
               continue
            if found_version:
                collection_index_settings['version'] = line.rstrip()
                found_version = False
                continue
            if '[Name]' in line :
                found_name = True
                continue
            if found_name:
                collection_index_settings['name'] = line.rstrip()
                found_name = False
                continue
            if '[PlateExtension]' in line:
                found_extension = True
                continue
            if found_extension :
                collection_index_settings['plate_extension'] = line.rstrip()
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
                    collection_index_settings['adapters'] = adapter_list
                    break
                else:
                    adapter_list.append(line_split[-1])
    return collection_index_settings


def index_library_information (index_library_id) :

    index_library_dict ={}
    index_library_found = IndexLibraryKit.objects.get(pk=index_library_id)
    general_information = [index_library_found.get_index_library_information().split(';')]
    index_library_dict['general_information'] = general_information

    if IndexLibraryValues.objects.filter(indexLibraryKit_id__exact = index_library_id).exists():
        index_list = IndexLibraryValues.objects.filter(indexLibraryKit_id__exact = index_library_id)

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
        index_library_dict['I7_indexes'] = list(I7_indexes for I7_indexes,_ in itertools.groupby(I7_indexes))
        if len(I5_indexes) > 0 :
            I5_indexes.sort()
            index_library_dict['I5_indexes'] = list(I5_indexes for I5_indexes,_ in itertools.groupby(I5_indexes))
        index_library_dict['default_layout'] = index_row_info
        return index_library_dict
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
        The function is store the collection kit file inisde COLLECTION_INDEX_KITS_DIRECTORY

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
    file_name=os.path.join(COLLECTION_INDEX_KITS_DIRECTORY ,  str(f_name + '_' +timestr + f_extension))
    filename = fs_index_lib.save(file_name,  collection_file)
    saved_file = os.path.join(settings.MEDIA_ROOT, file_name)

    return saved_file

def store_collection_settings (collection_settings, file_name) :
    # saving library settings into database
        if len(collection_settings['adapters']) == 1:
            adapter_2 = ''
        else :
            adapter_2 = collection_settings['adapters'][1]
        lib_settings_to_store = CollectionIndexKit(collectionIndexName = library_settings['name'],
                                    version =  library_settings ['version'],
                                    plateExtension = library_settings['plate_extension'] ,
                                    adapter1 = library_settings['adapters'][0],
                                    adapter2 = adapter_2, collectionLibraryFile = file_name)
        lib_settings_to_store.save()

def remove_boom_bytes (input_file):
    s = open(input_file, mode='r', encoding='utf-8-sig').read()
    open(input_file, mode='w', encoding='utf-8').write(s)


def store_collection_indexes (collection_index):
    # saving index values into database
    for row in collection_index :
        index_to_store = CollectionIndexValues(indexLibraryKit_id = lib_settings_to_store,
                                defaultWell = row[0], index_7 = row[1],
                                i_7_seq = row[2], index_5 = row[3],
                                i_5_seq = row[4])
        index_to_store.save()
