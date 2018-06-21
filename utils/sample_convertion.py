#!/usr/bin/env python3
# coding: utf-8


import os
import re
from datetime import datetime
import time

from Bio.Seq import Seq
from django.conf import settings

from iSkyLIMS_wetlab import wetlab_config

def include_csv_header (library_kit, out_file, plate, container):
    csv_header=['FileVersion','LibraryPrepKit','ContainerType','ContainerID','Notes']

    header_settings=['1','',plate,container,'automatic generated file from iSkyLIMS']

    index_kit=csv_header.index('LibraryPrepKit')
    header_settings[index_kit]=library_kit
    out_file.write('[Header]\n')
    for i in range(len(csv_header)):
        header_line=str(csv_header[i] + ',' + header_settings[i] + '\n')
        out_file.write(header_line)
    #### adding additional line
    out_file.write('\n')

def sample_sheet_map_basespace(in_file, library_kit, library_kit_file, projects, plate):
    target_data_header_one_index = ['SampleID','Name','Species','Project','NucleicAcid',
               'Well','Index1Name','Index1Sequence']
    target_data_header_two_index = ['SampleID','Name','Species','Project','NucleicAcid',
               'Well','Index1Name','Index1Sequence','Index2Name','Index2Sequence']
    ## Note that original header does not exactly match with the real one , but it is defined like this
    ## to get an easy way to map fields in the sample sheet and the sample to import to base Space
    original_data_header_one_index = ['SampleID','Name','Plate', 'Well','Index1Name','Index1Sequence',
                       'Project','Description']
    original_data_header_two_index = ['SampleID','Name','Plate', 'Well','Index1Name','Index1Sequence',
                       'Index2Name','Index2Sequence','Project','Description']


    data_raw=[]
    well_column={}
    well_row={}
    letter_well='A'
    number_well='01'
    result_directory=wetlab_config.MIGRATION_DIRECTORY_FILES
    cwd = os.getcwd()
    data_found=0
    header_found=0
    only_one_index = False

    fh = open(in_file,'r')
    for line in fh:
        line=line.rstrip()
        if line == '':
            continue
        #import pdb; pdb.set_trace()
        date_found = re.search('^Date',line)
        if date_found :
            date_line = line.split(',')
            # check if the year contains only 2 digits
            temp_date = date_line[1].split('/')
            if len (temp_date[2]) == 2 :
                 # add the 20 century to the year
                 year_four_digits = '20' + temp_date[2]
                 temp_date[2] = year_four_digits
                 date_line[1] = '/'.join(temp_date)

            try:
                date_object = datetime.strptime(date_line[1],'%m/%d/%Y')
            except:
                date_object = datetime.strptime(date_line[1],'%d/%m/%Y')
            date_sample = date_object.strftime('%Y%m%d')
            date_found = False

        found=re.search('^\[Data\]', line)
        if found:
            data_found=1
            continue
        found_header=re.search('^Sample_ID,Sample_Name',line)
        if found_header:
            header_found=1
            continue
        if (data_found and header_found):
            dict_value_data={}
            data_split=line.split(',')
            if len(data_split) < 9 :
                # sample sheet does not include the index 2.
                only_one_index = True
                original_data_header = original_data_header_one_index
                project_index = 6
            else:
                original_data_header = original_data_header_two_index
                project_index = 8
            if (data_split[project_index] in projects):
            #if projects in line :
                for ind in range(len(original_data_header)):
                    dict_value_data[original_data_header[ind]]= data_split[ind]
               #### adding empty values of species and NucleicAccid
                dict_value_data['Species']=''
                dict_value_data['NucleicAcid']='DNA'

                if not dict_value_data['Index2Name'] in well_row:
                    well_row[dict_value_data['Index2Name']]=letter_well
                    letter_well=chr(ord(letter_well)+1)
                if not dict_value_data['Index1Name'] in well_column:
                    well_column[dict_value_data['Index1Name']]=number_well
                    number_well =str(int(number_well)+1).zfill(2)
                dict_value_data['Well']=str(well_row[dict_value_data['Index2Name']]+ well_column[dict_value_data['Index1Name']])

                data_raw.append(dict_value_data)
    fh.close()
    # containerID build on the last Sample_Plate and the date in the sample sheet
    container = str(data_split[2] + date_sample)
    data_found=0
    tmp= re.search('.*/(.*)\.csv',in_file)
    out_tmp=tmp.group(1)
    #base_dir=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    base_dir = settings.MEDIA_ROOT
    out_file = str(base_dir + result_directory +  out_tmp + '_for_basespace' + '_'+ library_kit_file + '.csv')


    #### check for validation of the shample file
    try:
        dict_value_data
    except:
        return ('Error')
    #### open file for writting the conversion file
    fh_out = open (out_file, 'w')
    #####  print csv header
    include_csv_header(library_kit,fh_out,plate,container)
    #####  print data header
    fh_out.write('[Data]\n')
    if only_one_index == True :
        target_data_header = target_data_header_one_index
    else :
        target_data_header = target_data_header_two_index
    for i in range(len(target_data_header)):
        fh_out.write(target_data_header[i])
        if i < len(target_data_header)-1:
            fh_out.write(',')
        else:
            fh_out.write('\n')

    for line in data_raw:
        #### reverse order for Index2
        if only_one_index == False :
            seq=Seq(line['Index2Sequence'])
            line['Index2Sequence']=str(seq.reverse_complement())

        for i in  range(len(target_data_header)):
            fh_out.write(line[target_data_header[i]])
            if i < len(target_data_header)-1:
                fh_out.write(',')
            else:
                fh_out.write('\n')

    fh_out.close()
    # remove the absolute path of the library file_name
    absolute_path = str(settings.BASE_DIR + '/')
    file_name_in_database = out_file.replace(absolute_path,'')
    return file_name_in_database

def get_projects_in_run(in_file):
    header_found=0
    projects={}
    fh = open(in_file,'r')
    for line in fh:
        line=line.rstrip()
        found_header=re.search('^Sample_ID,Sample_Name',line)
        if found_header:
            header_found=1
            ## found the index for projects
            p_index= line.split(',').index('Sample_Project')
            description_index=line.split(',').index('Description')
            continue
        if header_found :
            ### ignore the empty lines
            if line == '':
                continue
            ## store the project name and the user name (Description) inside projects dict
            projects[line.split(',')[p_index]]=line.split(',')[description_index]
    fh.close()
    return projects

def get_experiment_library_name (in_file):
    experiment_name = ''
    library_name = ''
    fh = open(in_file, 'r')
    for line in fh:
        line = line.rstrip()
        found_experiment = re.search('^Experiment Name',line)
        found_library = re.search('^Assay',line)
        if found_experiment :
            experiment_value = line.split(',')
            if experiment_value[1]:
                experiment_name = experiment_value[1]
                found_experiment = 0
        if found_library :
            library_value = line.split(',')
            if library_value[1]:
                library_name = library_value[1]
                found_library = 0
    fh.close()

    return experiment_name, library_name

def update_library_kit_field (library_file_name, library_kit_name, library_name):
    #result_directory='documents/wetlab/BaseSpaceMigrationFiles/'
    timestr = time.strftime("%Y%m%d-%H%M%S")
    tmp= re.search('(.*)\d{8}-\d+.*\.csv',library_file_name)
    absolute_path = str(settings.BASE_DIR + '/')
    out_file = str( absolute_path +  tmp.group(1) + timestr +'_for_basespace_'+ library_kit_name + '.csv')
    #import pdb; pdb.set_trace()
    try:
        fh_in = open (library_file_name, 'r')
        fh_out = open (out_file, 'w')
    except:
        return 'ERROR:'

    for line in fh_in:
        found_library_kit = re.search('LibraryPrepKit', line)
        # library kit found. Replace line with new library kit
        if found_library_kit:
            line = str ('LibraryPrepKit,' + library_name + '\n')
        fh_out.write(line)
    fh_in.close()
    fh_out.close()
    os.remove(library_file_name)
    # remove absolute path from file_name
    absolute_path = str(settings.BASE_DIR + '/')
    file_name_in_database = out_file.replace(absolute_path,'')
    return file_name_in_database

def update_sample_sheet (in_file, experiment_name):

    out_line = str ( 'Experiment Name,'+ experiment_name+ '\n')
    fh_in = open (in_file, 'r')
    fh_out = open ('temp.txt' ,'w')
    experiment_line_found  = False
    for line in fh_in:
        # find experiment name line
        found_experiment = re.search('^Experiment Name',line)
        if found_experiment :
            fh_out.write(out_line)
            experiment_line_found = True

        elif line == '\n' and experiment_line_found == False:
            fh_out.write(out_line)
            fh_out.write('\n')
            experiment_line_found = True
        else:
            fh_out.write(line)

    fh_in.close()
    fh_out.close()
    os.rename('temp.txt', in_file)

def create_unique_sample_id_values (infile, index_file):
    found_sample_line = False

    fh = open (infile, 'r')
    temp_sample_sheet = os.path.join(settings.MEDIA_ROOT, 'wetlab','tmp_file')
    fh_out_file = open (temp_sample_sheet, 'w')
    with open(index_file) as fh_index:
        index = fh_index.readline()
        index =index.rstrip()
        index_number_str, index_letter = index.split('-')
        fh_index.close()

    for line in fh:
        if 'Sample_ID' in line:
            found_sample_line =True
            fh_out_file.write(line)
            continue
        if found_sample_line :
            data_line = line.split(',')
            data_line [0] = str(index_number_str + '-' + index_letter)
            new_line = ','.join(data_line)
            fh_out_file.write(new_line)
            # increase the index number
            index_number = int(index_number_str) +1
            if index_number > 9999:
                index_number = 0
                #increase a letter
                split_index_letter = list(index_letter)
                if split_index_letter[1] == 'Z':
                    last_letter = chr(ord(split_index_letter[0])+1)
                    split_index_letter[0] = last_letter
                    split_index_letter[1] = 'A'
                    index_letter = ''.join(split_index_letter)
                else:
                    first_letter=chr(ord(split_index_letter[1])+1)
                    split_index_letter[1] = first_letter
                    index_letter = ''.join(split_index_letter)

            index_number_str = str(index_number)
            index_number_str = index_number_str.zfill(4)

        else:
            fh_out_file. write(line)

    #dump the index value to file
    fh_index = open(index_file, 'w')
    index_line = str(index_number_str + '-' + index_letter)
    fh_index.write(index_line)
    fh_index.close()
    fh.close()
    fh_out_file.close()
    os.rename('temp_sample_sheet', infile)

