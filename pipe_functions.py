
import os,paramiko, fnmatch

def open_connect_sftp(host, user):
    '''
    Description: 
        Connect via ssh to the server with the user
    Input:
        host  # server to connect
        user  # user to connect to the server
    Variables:
        conn_data  # dictionary with data connection
    Return:
        ssh_client  # object SSHClient()
        sftp # object SFTPClient()
    '''
    conn_data = dict(hostname = host, username = user)
    ssh_client = paramiko.SSHClient()
    ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh_client.connect(**conn_data)
    sftp = ssh_client.open_sftp() 
    return ssh_client,sftp

def close_connect_sftp(ssh_client, sftp):
    '''
    Description:
        disconnect connection sftp
    Input:
        ssh_client  # object SSHClient()
        sftp  # object SFTPClient()
    Variables:

    Return:
    '''
    sftp.close()
    ssh_client.close()


def get_run_number(run_id):
    '''
    Description:
       get run  number from  run _id
    Input:
        sftp  # object SFTPClient()
    Variables:
    Return:
        run_number
    '''
    run_number = run_id[run_id.rfind('_')-3:run_id.rfind('_')]
    return run_number

def make_new_folder(sftp,folder):
    '''
    Description:
       make new folder UMPXXX_runxxx
    Input:
        sftp  # object SFTPClient()
        folder # new folder
    Variables:
    Return:
    '''
    sftp.mkdir(folder,modo=775)
 


def get_seq_new_folder (sftp,folder,run_id):
    '''
    Description:
       get sequential folder number XXX UMPXXX_runxxx
    Input:
        sftp  # object SFTPClient()
    Variables:
    Return:
        new_folder # sequential new folder XXX
    '''
    sftp.chdir(folder)
    directories = sftp.listdir()
    max = 0
    for dir in directories:
        if dir.startswith('UMP'):
            number = dir[3:dir.find("_")]
            if number != '':
              
            #number = int(dir[3:dir.find("_")])
                 if int(number) > max:
                      max = number
    run = run_id[run_id.rfind('_')-3:run_id.rfind('_')]
    new_folder = "UMP" + str(max) + "_run" + run
    return  new_folder

def copy_file(sftp, localfile, remotefile):
    '''
    Description:
        copy a file from localdir to remotedir
    Input:
        sftp  # object SFTPClient()
        localdir  # local  folder to be copied
        remotedir # remote folder  where local folder  will be copied
    Variables:        
        
    Return:
       num  # number files have been copied
    '''
    num = sftp.put(localfile, remotefile)
    return num

 
def copy_list_files(sftp, localdir, remotedir, extension):
    '''
    Description:
        copy fatsq files from localdir to remotedir
    Input:
        sftp  # object SFTPClient()
        localdir  # local  folder to be copied
        remotedir # remote folder  where local folder  will be copied
        extension # extension files will be copied
    Variables:        
        filename # name file 
    Return:
        count  # number files have been copied
    '''
    if not exists_path(sftp,remotedir):
        sftp.mkdir(remotedir)
    count = 0 
    
    for filename in fnmatch.filter(os.listdir(localdir),extension):
            local_file = os.path.join(localdir,filename)
            remote_file = os.path.join(remotedir,filename)
            num = copy_file(sftp,local_file,remote_file)
            count  = count + num
    return count

def exists_path(sftp, path):
    '''
    Description:
        check if path exist
    Input:
        sftp  # object SFTPClient()
        path 
    Variables:        
    Return:
        True/False
    '''
    try:
        sftp.stat(path)
    except e:
        if e.errno == errno.ENOENT:
            return False
        raise
    else:
        return True


def copy_folder_run(localdir, remotedir):
    '''
    Description:
        copy run folder from localdir to remotedir
    Input:
        
    Variables:        
    Return:
        count  # number files have been copied
    '''

