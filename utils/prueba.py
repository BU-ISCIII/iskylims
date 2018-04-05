#!/usr/bin/env python3
#from  ..models import *

import sys, os, re
import xml.etree.ElementTree as ET
import time
import shutil
import logging
from logging.handlers import RotatingFileHandler
from interop import py_interop_run_metrics, py_interop_run, py_interop_summary, py_interop_plot
from smb.SMBConnection import SMBConnection

def open_log():
    
    LOG_FILENAME = 'testing.log'
    log_name=os.path.join('../log/', LOG_FILENAME)
    #def create_log ():
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    #create the file handler
    handler = logging.handlers.RotatingFileHandler(log_name, maxBytes=20000, backupCount=5)
    handler.setLevel(logging.DEBUG)
    
    #create a Logging format
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    #add the handlers to the logger
    logger.addHandler(handler)
    
    return logger

def open_samba_connection():
 ## open samba connection
 # There will be some mechanism to capture userID, password, client_machine_name, server_name and server_ip
 # client_machine_name can be an arbitary ASCII string
 # server_name should match the remote machine name, or else the connection will be rejected
 
    conn=SMBConnection('Luigi', 'Apple123', 'bioinfo', 'LUIGI-PC', use_ntlm_v2=True)
    conn.connect('192.168.1.3', 139)
   
    '''
    conn = SMBConnection(userid, password, client_machine_name, remote_machine_name, use_ntlm_v2 = True)
    conn.connect(server_ip, 139)
    '''
    
    return conn


runid_name='161123_NS500454_0096_AHFGV5BGXY'
run_Id_used='161123_NS500454_0096_AHFGV5BGXY'

local_dir_samba= '../tmp/processing'
logger=open_log()
logger.info('test')


demux_file='../tmp/processing/DemultiplexingStats.xml'
conversion_file='../tmp/processing/ConversionStats.xml'
run_processing_id=2



try:
    conn=open_samba_connection()
    #logger.info('Successful connection for getting size of the files' )
except:
    print('Error')
    exit(0)
remote_stats_dir= 'Data/Intensities/BaseCalls/Stats/'
samba_demux_file=os.path.join('/',run_Id_used,remote_stats_dir, 'DemultiplexingStats.xml')
logger.debug('path to fetch demultiplexingStats is %s',  samba_demux_file)

share_folder_name='Flavia'
run_folder=os.path.join('/',run_Id_used)

file_list = conn.listPath( share_folder_name, run_folder)
for sh in file_list:
    if sh.isDirectory:
        continue
    print ('file name is = ', sh.filename)
    print ('file size is = ', sh.file_size)
    
            
    # close samba connection 
conn.close()



print ('completed')

#xml_stats=parsing_statistics_xml(demux_file, conversion_file, logger)
#store_raw_xml_stats(xml_stats, run_processing_id,logger)
#process_xml_stats(xml_stats,run_processing_id, logger)
#process_binStats(local_dir_samba, run_processing_id, logger)
#create_graphics(local_dir_samba, graphic_dir, logger)


