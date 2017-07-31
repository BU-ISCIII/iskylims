#!/usr/bin/env python3

def createSSHClient(server, port, user, password):
    client = paramiko.SSHClient()
    client.load_system_host_keys()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(server, port, user, password)
    return client
ssh = createSSHClient('barbarroja', 22, 'lchapado', 'chapadomaster') 
ftp = ssh.open_sftp()  
ftp.get('localfile.py', 'remotefile.py')
ftp.put('localfile', 'remotefile')
ftp.close()