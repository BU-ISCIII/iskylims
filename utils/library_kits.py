#!/usr/bin/env python3
import re
'''
[Version]
1
[Name]
Nextera XT Index Kit (24 Indexes, 96 Samples)
[PlateExtension]
nexxt24
[Settings]
Adapter	CTGTCTCTTATACACATCT
[I7]
N701	TAAGGCGA
N702	CGTACTAG
N703	AGGCAGAA
N704	TCCTGAGC
N705	GGACTCCT
N706	TAGGCATG
[I5]
S502	CTCTCTAT
S503	TATCCTCT
S504	AGAGTAGA
S517	GCGTAAGA
[DefaultLayout_SingleIndex]
'''
def check_index_library_file_format (input_file):
    heading_checks = ['[Version]','[Name]', '[PlateExtension]','[Settings]', '[I7]','[I5]']
    with open (input_file ) as fh:
       read_data = fh.read()
    for check in heading_checks :
        if check not in read_data:
            return False
    
    return True



def getting_index_library_name (input_file):
    found_name = False
    library_name = ''
    with open (input_file ) as fh:
        for line in fh:
            found_name_library = re.search('^\[Name\]',line)
            if found_name_library :
                found_name = True
                continue
            if found_name:
                library_name= line.rstrip()
                break
    
    return library_name
    

def get_library_settings (input_file):
    library_settings ={}
    adapter_list = []
    with open (input_file) as fh:
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
                # No more adapters  have been found. copy the adapters in library_settings dictionary and exit the loop
                if len (line_split) == 1:
                    library_settings['adapters'] = adapter_list
                    break
                else:
                    adapter_list.append(line_split[-1])
                    
                
    return library_settings

def get_index_values (input_file):
    
    
    with open (input_file) as fh :
        index_7 , index_5 = [] , []
        found_I7 = False
        found_I5 = False
        index_values = {}
        
        for line in fh:
            if '[I7]' in line :
                found_I7 = True
                continue
            if found_I7 :
                line=line.rstrip()
                if len (line.split('\t'))> 1:
                    index_7.append(line.split('\t'))
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
                    index_5.append(line.split('\t'))
                    continue
                # No more index have been found for I5. Exit the loop
                else:
                    break
        index_values['I7'] = index_7
        index_values['I5'] = index_5
        
    return index_values
