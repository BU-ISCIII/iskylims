#!/usr/bin/env python3
# coding: utf-8


import os
import re
from datetime import datetime

from Bio.Seq import Seq



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
    
def sample_sheet_map_basespace(in_file, library_kit, projects, plate):
    target_data_header = ['SampleID','Name','Species','Project','NucleicAcid',
               'Well','Index1Name','Index1Sequence','Index2Name','Index2Sequence']
    ## Note that original header does not exactly match with the real one , but it is defined like this
    ## to get an easy way to map fields in the sample sheet and the sample to import to base Space
    original_data_header = ['SampleID','Name','Plate', 'Well','Index1Name','Index1Sequence',
                       'Index2Name','Index2Sequence','Project','Description']
    

    data_raw=[]
    well_column={}
    well_row={}
    letter_well='A'
    number_well='01'
    result_directory='documents/wetlab/BaseSpaceMigrationFiles/'
    cwd = os.getcwd()
    data_found=0
    header_found=0
  
    fh = open(in_file,'r')
    for line in fh:
        line=line.rstrip()
        if line == '':
            continue
        #import pdb; pdb.set_trace()
        date_found = re.search('^Date',line)
        if date_found :
            date_line = line.split(',')
            
            try:
                date_object = datetime.strptime(date_line[1],'%m/%d/%Y')
            except:
                date_object = datetime.strptime(date_line[1],'%d/%m/%Y')
            date_sample = date_object.strftime('%Y%m%d')
            date_found = False
        #found=re.search('Assay,(.*)',line)
        #if found:
        #    LibraryPrepKit= found.group(1)
        #    continue
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
            if (data_split[8] in projects):
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
    out_file = str(result_directory +  out_tmp + '_for_basespace' + '_'+ library_kit + '.csv')

    
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
    for i in range(len(target_data_header)):
        fh_out.write(target_data_header[i])
        if i < len(target_data_header)-1:
            fh_out.write(',')
        else:
            fh_out.write('\n')
        
    for line in data_raw:
        #### reverse order for Index2
        seq=Seq(line['Index2Sequence'])
        line['Index2Sequence']=str(seq.reverse_complement())
        
        for i in  range(len(target_data_header)):      
            fh_out.write(line[target_data_header[i]])
            if i < len(target_data_header)-1:
                fh_out.write(',')
            else:
                fh_out.write('\n')
      
    fh_out.close()
    return out_file

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
