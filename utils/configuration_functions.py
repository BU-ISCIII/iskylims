import os
from iSkyLIMS_drylab import drylab_config
from django.conf import settings


def create_email_conf_file(email_fields, application):
    '''
    Description:
        create the email configuration file . If exists the old information is deleted
    Input:
        email_fields    # Email fields settings
    Constants:
        EMAIL_CONFIGURATION_FILE_HEADING
    Functions:
        get_type_of_data # located at this file
    Return:
        True if sucessful creation
    '''
    conf_file = os.path.join(settings.BASE_DIR, application,'drylab_email_conf.py')
    #if os.path.isfile(conf_file):
    try:

        with open (conf_file, 'w') as out_fh:
            out_fh.write(drylab_config.EMAIL_CONFIGURATION_FILE_HEADING)
            for key in sorted(email_fields):
                type_of_data = get_type_of_data(email_fields[key])
                if type_of_data == 'boolean' or type_of_data == 'integer':
                    out_fh.write(key + ' = '+ email_fields[key]+ '\n')
                else:
                    out_fh.write(key + ' = \''+ email_fields[key]+ '\'\n')
            out_fh.write(drylab_config.EMAIL_CONFIGURATION_FILE_END)
    except:
        return False

    return True

def create_samba_conf_file(samba_fields, application):
    '''
    Description:
        create the samba configuration file . If exists the old information is deleted
    Input:
        samba_fields    # Samba fields settings
    Constants:
        SAMBA_CONFIGURATION_FILE_HEADING
    Functions:
        get_type_of_data # located at this file
    Return:
        True if sucessful creation
    '''
    conf_file = os.path.join(settings.BASE_DIR, application,'drylab_samba_conf.py')
    #if os.path.isfile(conf_file):
    try:

        with open (conf_file, 'w') as out_fh:
            out_fh.write(drylab_config.SAMBA_CONFIGURATION_FILE_HEADING)
            for key in sorted(samba_fields):
                type_of_data = get_type_of_data(samba_fields[key])
                if type_of_data == 'boolean' or type_of_data == 'integer':
                    out_fh.write(key + ' = '+ samba_fields[key]+ '\n')
                else:
                    out_fh.write(key + ' = \''+ samba_fields[key]+ '\'\n')
            out_fh.write(drylab_config.SAMBA_CONFIGURATION_FILE_END)
    except:
        return False

    return True

def get_email_data_from_file(application):
    '''
    Description:
        Fetch the email configuration file
    Inputs:
        application     # Application name
    Constants:
        EMAIL_CONFIGURATION_FILE_HEADING
        EMAIL_CONFIGURATION_FILE_END
    Return:
        email_data
    '''
    conf_file = os.path.join(settings.BASE_DIR, application,'drylab_email_conf.py')
    email_data = {}
    heading_found = False
    try:
        with open (conf_file, 'r') as fh:
            for line in fh.readlines():
                if not heading_found and drylab_config.EMAIL_CONFIGURATION_FILE_HEADING.split('\n')[-2] in line :
                    heading_found = True
                    continue
                if drylab_config.EMAIL_CONFIGURATION_FILE_END in line:
                    break
                if heading_found :
                    line = line.rstrip()
                    key , value = line.split(' = ')
                    email_data[key] = value.replace('\'','')
        return email_data
    except:
        email_data

def get_samba_data_from_file(application):
    '''
    Description:
        Fetch the samba configuration file
    Inputs:
        application     # Application name
    Constants:
        SAMBA_CONFIGURATION_FILE_HEADING
        SAMBA_CONFIGURATION_FILE_END
    Return:
        samba_data
    '''
    conf_file = os.path.join(settings.BASE_DIR, application,'drylab_samba_conf.py')
    samba_data = {}
    heading_found = False
    try:
        with open (conf_file, 'r') as fh:
            for line in fh.readlines():
                if not heading_found and drylab_config.SAMBA_CONFIGURATION_FILE_HEADING.split('\n')[-2] in line :
                    heading_found = True
                    continue
                if drylab_config.SAMBA_CONFIGURATION_FILE_END in line:
                    break
                if heading_found :
                    line = line.rstrip()
                    key , value = line.split(' = ')
                    samba_data[key] = value.replace('\'','')
        return samba_data
    except:
        samba_data

def open_samba_connection():
    '''
    Description:
        The function open a samba connection with the parameter settings
        defined in frylab configuration file
    Return:
        conn object for the samba connection
    '''
    try:

        conn=SMBConnection(drylab_config.SAMBA_USER_ID, drylab_config.SAMBA_USER_PASSWORD, drylab_config.SAMBA_SHARED_FOLDER_NAME,
                            drylab_config.SAMBA_REMOTE_SERVER_NAME, use_ntlm_v2=drylab_config.SAMBA_NTLM_USED, domain = drylab_config.SAMBA_DOMAIN)

        conn.connect(drylab_config.SAMBA_IP_SERVER, int(drylab_config.SAMBA_PORT_SERVER))
    except:
        return False

    return conn

def open_log(config_file):
    '''
    Description:
        The function will create the log object to write all logging information
    Input:
        logger_name    # contains the logger name that will be included
                        in the log file
    Constant:
        LOGGING_CONFIG_FILE
    Return:
        logger object
    '''
    fileConfig(config_file)
    logger = logging.getLogger(__name__)
    return logger

def get_type_of_data (data):
    '''
    Description:
        The function get always as input a string class.
        By trying to convert the input data to int or bolean it will decide the type of data
        If not possible to conver it returns string
    Return:
        type_of_data
    '''
    boolean_values = ['True', 'False', 'None']
    if data in boolean_values :
        return 'boolean'
    try:
        integer = int(data)
        return 'integer'
    except:
        return 'string'
