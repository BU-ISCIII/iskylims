




def process_run_in_samplesent_state (process_list, logger):
     # prepare a dictionary with key as run_name and value the RunID
     processed_run=[]
     for run_item in process_list:
         logger.info ('Running the process sample sent state for %s', run_item)
         run_be_processed_id=RunProcess.objects.get(runName__exact=run_item).id
         logger.debug ('Run ID for the run process to be update is:  %s', run_be_processed_id)
         #run_Id_for_searching=RunningParameters.objects.get(runName_id= run_be_processed_id)
         update_run_state(run_be_processed_id, 'Process Running', logger)
         processed_run.append(run_be_processed_id)
     return processed_run


def process_run_in_processrunning_state (process_list, logger):
    processed_run=[]
    logger.debug('starting the process_run_in_processrunning_state method')
    try:
        conn=open_samba_connection()
        logger.info('check the Sucessful connection to NGS_Data before starting processing runing state method')

    except:
        return('Error')

    share_folder_name = wetlab_config.SAMBA_SHARED_FOLDER_NAME
    for run_item in process_list:
        logger.debug ('processing the run %s in process running state' , run_item)
        run_be_processed_id=RunProcess.objects.get(runName__exact=run_item).id
        run_Id_used=str(RunningParameters.objects.get(runName_id= run_be_processed_id))
        logger.debug ('found the run ID  %s' , run_Id_used )
        run_folder=os.path.join('/',run_Id_used,'Data/Intensities/BaseCalls')
        # check if runCompletion is avalilable
        logger.debug ('found the run ID  %s' , run_Id_used )
        try:  #####
            file_list = conn.listPath( share_folder_name, run_folder)
        except: #####
            logger.error('no folder found for run ID %s', run_Id_used)  #####
            continue  #####
        #import pdb; pdb.set_trace()
        found_report_directory = 0
        for sh in file_list:
            if sh.filename =='Reports' :
                logger.info('bcl2fastq has been completed for run %s', run_Id_used)
                processed_run.append(run_Id_used)
                update_run_state(run_be_processed_id, 'Bcl2Fastq Executed', logger)
                update_project_state(run_be_processed_id, 'B2FqExecuted', logger)
                # Get the time when  the Bcl2Fastq process is ending
                conversion_stats_file = os.path.join (run_Id_used,'Data/Intensities/BaseCalls/Stats/', 'ConversionStats.xml')
                conversion_attributes = conn.getAttributes(wetlab_config.SAMBA_SHARED_FOLDER_NAME ,conversion_stats_file)
                run_date = RunProcess.objects.get(pk=run_be_processed_id)
                run_date.bcl2fastq_finish_date = datetime.datetime.fromtimestamp(int(conversion_attributes.create_time)).strftime('%Y-%m-%d %H:%M:%S')
                run_date.save()
                logger.info ('Updated the Bcl2Fastq time in the run %s', run_item)
                break
            else:
                logger.debug('The directory %s has been found while looking for completion of the execution of bcl2fastq', sh.filename)
        if found_report_directory:
            logger.info('blc2fastq has been completed for the Run ID %s  it is now on Bcl2Fastq Executed state', run_Id_used)
        else:
            logger.info('blc2fastq was not finish for the Run ID %s  waiting for Bcl2Fastq to be completed', run_Id_used)


    # close samba connection
    conn.close()
    logger.info('Closing the remote connection ')
    return processed_run



def process_run_in_bcl2Fq_executed_state (process_list, logger):
    processed_run=[]
    # get the directory of samba to fetch the files
    share_folder_name = wetlab_config.SAMBA_SHARED_FOLDER_NAME
    local_dir_samba= wetlab_config.RUN_TEMP_DIRECTORY_PROCESSING
    remote_stats_dir= 'Data/Intensities/BaseCalls/Stats/'
    demux_file=os.path.join(local_dir_samba,'DemultiplexingStats.xml')
    conversion_file=os.path.join(local_dir_samba,'ConversionStats.xml')
    run_info_file=os.path.join(local_dir_samba, 'RunInfo.xml')

    logger.debug('Executing process_run_in_bcl2F_q_executed_state method')

    # check the connectivity to remote server
    try:
        conn=open_samba_connection()
        logger.info('Successful connection for updating run on bcl2F_q' )
    except:
        logger.error('ERROR:: Unable to connect to remote server')
        return 'Error'

    for run_item in process_list:

        logger.info ('Processing the process on bcl2F_q for the run %s', run_item)
        run_processing_id=RunProcess.objects.get(runName__exact=run_item).id
        run_Id_used=str(RunningParameters.objects.get(runName_id= run_processing_id))

        update_run_state(run_processing_id, 'Running Stats', logger)
        #copy the demultiplexingStats.xml file to wetlab/tmp/processing


        samba_demux_file=os.path.join('/',run_Id_used,remote_stats_dir, 'DemultiplexingStats.xml')
        logger.debug('path to fetch demultiplexingStats is %s',  samba_demux_file)
        try:
            with open(demux_file ,'wb') as demux_fp :
                conn.retrieveFile(share_folder_name, samba_demux_file, demux_fp)
        except:
            logger.error('Unable to fetch the DemultiplexingStats.xml file for RunID %s', run_Id_used)
            os.remove(demux_file)
            logger.debug('deleting DemultiplexingStats file  for RunID %s' , run_Id_used)
            continue
        logger.info('Fetched the DemultiplexingStats.xml')
        #copy the ConversionStats.xml file to wetlab/tmp/processing
        samba_conversion_file=os.path.join('/', run_Id_used,remote_stats_dir,'ConversionStats.xml')
        try:
            with open(conversion_file ,'wb') as conv_fp :
                conn.retrieveFile(share_folder_name, samba_conversion_file, conv_fp)
        except:
            logger.error('Unable to fetch the ConversionStats.xml file for RunID %s', run_Id_used)
            os.remove(conversion_file)
            os.remove(demux_file)
            logger.debug('deleting ConversionStats and DemultiplexingStats file  for RunID %s' , run_Id_used)
            continue
        logger.info('Fetched the conversionStats.xml')
        # copy RunInfo.xml  file to process the interop files
        try:
            with open(run_info_file ,'wb') as runinfo_fp :
                samba_conversion_file=os.path.join('/', run_Id_used,'RunInfo.xml')
                conn.retrieveFile(share_folder_name, samba_conversion_file, runinfo_fp)
        except:
            logger.error('Unable to fetch the RunInfo.xml file for RunID %s', run_Id_used)
            os.remove(run_info_file)
            os.remove(conversion_file)
            os.remove(demux_file)
            logger.debug('deleting RunInfo, ConversionStats and DemultiplexingStats file  for RunID %s' , run_Id_used)
            continue
        logger.info('Fetched the RunInfo.xml file')

        # copy all binary files in interop folder to local  documents/wetlab/tmp/processing/interop
        interop_local_dir_samba= os.path.join(local_dir_samba, 'InterOp')
        remote_interop_dir=os.path.join('/',run_Id_used,'InterOp')
        try:
            file_list = conn.listPath( share_folder_name, remote_interop_dir)
            logger.info('InterOp folder exists on the RunID %s', run_Id_used)
            run_parameters_file=os.path.join(local_dir_samba,'runParameters.xml')
            try:
                with open(run_parameters_file ,'wb') as runparam_fp :
                    samba_conversion_file=os.path.join('/', run_Id_used,'runParameters.xml')
                    conn.retrieveFile(share_folder_name, samba_conversion_file, runparam_fp)
                logger.info('Fetched the runParameters.xml file. Written as: '
                    +str(run_parameters_file))
            except:
                logger.error('Unable to fetch the runParameters.xml file for RunID %s', run_Id_used)
                os.remove(run_parameters_file)
                logger.debug('Deleting runParameters file  for RunID %s' , run_Id_used)

        except:
            logger.error('ERROR:: InterOP folder does not exist on RunID %s', run_Id_used)
            os.remove(run_info_file)
            os.remove(conversion_file)
            os.remove(demux_file)
            logger.debug('deleting RunInfo, ConversionStats and DemultiplexingStats file  for RunID %s' , run_Id_used)
            print('ERROR:: InterOP folder does not exist on RunID ', run_Id_used)
            continue
        error_in_interop = False
        for sh in file_list:
            if sh.isDirectory:
                continue
            else:
                interop_file_name=sh.filename
                remote_interop_file=os.path.join(remote_interop_dir, interop_file_name)
                copy_file=os.path.join(interop_local_dir_samba, interop_file_name)
                try:
                    with open(copy_file ,'wb') as cp_fp :
                        remote_file=os.path.join(remote_interop_dir,)
                        logger.debug('File %s to be copied on the local directory', interop_file_name)
                        conn.retrieveFile(share_folder_name, remote_interop_file, cp_fp)
                        logger.info('Copied %s to local Interop folder', interop_file_name)
                except:
                    logger.error("Not be able to fetch the file %s", interop_file_name)
                    os.remove(run_info_file)
                    os.remove(conversion_file)
                    os.remove(demux_file)
                    logger.debug('deleting files RunInfo, ConversionStats and DemultiplexingStats ')
                    logger.debug('because of error when fetching interop files for RunID  %s' , run_Id_used)
                    # deleting local Interop files
                    file_list_to_delete = os.list(interop_local_dir_samba)
                    logger.debug ('Deleting the local interop files ')
                    for file_name in file_list_to_delete:
                        if file_name == '.' or '..' :
                            continue
                        else:
                            os.remove(file_name)
                    error_in_interop= True
                    break
        if error_in_interop :
            continue
        else:
            # parsing the files to get the xml Stats
            logger.info('processing the XML files')
            xml_stats=parsing_statistics_xml(run_processing_id, demux_file, conversion_file, logger)
            result_of_raw_saving = store_raw_xml_stats(xml_stats,run_processing_id, logger)
            if result_of_raw_saving == 'ERROR':
                update_run_state(run_processing_id, 'ERROR-on-Raw-SavingStats', logger)
                update_project_state(run_processing_id, 'ERROR-on-Raw-SavingStats', logger)
                logger.error('Stopping process for this run an starting deleting the files')
            else:
                process_xml_stats(xml_stats,run_processing_id, logger)

                # parsing and processing the project samples
                sample_project_stats = parsing_sample_project_xml (run_processing_id,demux_file, conversion_file, logger)
                store_samples_projects (sample_project_stats, run_processing_id, logger)

                logger.info('processing interop files')
                # processing information for the interop files
                number_of_lanes = get_machine_lanes(run_processing_id)

                process_binStats(local_dir_samba, run_processing_id, logger, number_of_lanes)
                # Create graphics
                graphic_dir=os.path.join(settings.MEDIA_ROOT,wetlab_config.RUN_TEMP_DIRECTORY_PROCESSING)

                create_graphics(graphic_dir, run_processing_id, run_Id_used, logger)

                processed_run.append(run_Id_used)
                logger.info('run id %s is now on Completed state', run_Id_used)
                update_run_state(run_processing_id, 'Completed', logger)
                update_project_state(run_processing_id, 'Completed', logger)
            # clean up the used files and directories
            logger.info('starting the clean up for the copied files from remote server ')
            os.remove(demux_file)
            logger.debug('Demultiplexing file have been removed from %s', demux_file)
            os.remove(conversion_file)
            logger.debug('ConversionStats file have been removed from %s', conversion_file)
            os.remove(run_info_file)
            logger.debug('RunInfo file have been removed from %s', run_info_file)
            for file_object in os.listdir(interop_local_dir_samba):
                file_object_path = os.path.join(interop_local_dir_samba, file_object)
                if os.path.isfile(file_object_path):
                    logger.debug('Deleting file %s' , file_object_path)
                    os.remove(file_object_path)
            logger.info('xml files and binary files from InterOp folder have been removed')
            ## connect to server to get the disk space utilization of the folders
            get_run_disk_utilization (conn, run_Id_used, run_processing_id, logger)
            # Update the run with the date of the run completion
            completion_date = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            run_date_to_update = RunProcess.objects.get(pk = run_processing_id)
            run_date_to_update.process_completed_date = completion_date
            run_date_to_update.save()

    # close samba connection
    conn.close()
    logger.info('Samba connection close. Sucessful copy of files form remote server')
    return processed_run

