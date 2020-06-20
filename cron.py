from datetime import datetime
from django.conf import settings


from iSkyLIMS_drylab import drylab_config
from iSkyLIMS_drylab.utils.drylab_common_functions import  open_log
from iSkyLIMS_drylab.utils.automatic_jobs import check_pending_jobs , handling_pending_jobs
import os

def check_jobs_in_queued ():

    working_path = settings.MEDIA_ROOT
    os.chdir(working_path)
    config_file = os.path.join(settings.BASE_DIR,'iSkyLIMS_drylab',  drylab_config.LOGGING_CONFIG_FILE )
    logger = open_log(config_file)
    logger.info('###########---Start Crontab-----############')
    logger.info('Start searching for pending execution preparation Pipelines')
    if check_pending_jobs():
        handling_pending_jobs()
    else:
        logger.info('No pending jobs were found')
    logger.info('###########---Start Crontab-----############')
    return
