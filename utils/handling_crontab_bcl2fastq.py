import os
from iSkyLIMS_wetlab.wetlab_config import *
from .handling_crontab_common_functions import *


def check_demultiplexing_folder_exists(conn, run_folder, experiment_name):
    '''
    Description:
        The function check if the folder of demultiplexing files exists
    Input:
        conn                # samba connection instance
        run_folder          # run folder on the remote server
        experiment_name     # experiment name to be checked
    Constants:
        DEMULTIPLEXION_BCL2FASTQ_FOLDER
        REPORT_FOLDER
        STATS_FOLDER
        CONVERSION_STATS_FILE
    Functions:
        get_samba_application_shared_folder     # Located at utils/handling_crontab_common_functions.py
        get_samba_shared_folder                 # Located at utils/handling_crontab_common_functions.py
    Return:
        bcl2fastq_finish_date
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('%s : Starting function check_demultiplexing_folder_exists', experiment_name)
    bcl2fastq_finish_date = ''
    statistics_folder = os.path.join(get_samba_application_shared_folder(), run_folder, DEMULTIPLEXION_BCL2FASTQ_FOLDER)

    try:
        file_list = conn.listPath( get_samba_shared_folder(), statistics_folder)
    except:
        string_message = experiment_name + ' : Unable to fetch folder demultiplexing at  ' + statistics_folder
        logging_errors(string_message, True, True)
        logger.debug ('%s : End function check_demultiplexing_folder_exists with error', experiment_name)
        return {'ERROR':29}
    import pdb; pdb.set_trace()
    logger.info('%s : bcl2fastq has been completed . Collecting date when finish the bcl2fastq', experiment_name)

    s_conversion_stats = os.path.join (statistics_folder, STATS_FOLDER, CONVERSION_STATS_FILE)
    try:
        conversion_attributes = conn.getAttributes(get_samba_shared_folder() ,s_conversion_stats)
    except:
        string_message = experiment_name + ' : Unable to fetch ' + CONVERSION_STATS_FILE
        logging_errors(string_message, True, True)
        logger.debug ('%s : End function check_demultiplexing_folder_exists with error', experiment_name)
        return {'ERROR':31}
    bcl2fastq_finish_date = datetime.fromtimestamp(int(conversion_attributes.create_time)).strftime('%Y-%m-%d %H:%M:%S')
    logger.info ('%s : Collected Bcl2Fastq time ', experiment_name)

    import pdb; pdb.set_trace()
    logger.debug ('%s : End function waiting_time_expired', experiment_name)
    return bcl2fastq_finish_date
