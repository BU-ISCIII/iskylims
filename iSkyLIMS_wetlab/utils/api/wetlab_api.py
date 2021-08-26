from iSkyLIMS_wetlab.models import RunningParameters, Projects, SamplesInProject

def get_runfolder_from_run_name(run_name):
    '''
    Description:
        The function api get the run_name and return the folder name
    Return:
        folder_name or False if does not exist yet
    '''
    if RunningParameters.objects.filter(runName__runName_id_exact = run_name).exists():
        return RunningParameters.objects.filter(runName__runName_id_exact = run_name).get_run_folder()
    return False

def get_run_name_project_belongs(project_name):
    '''
    Description:
        The function api return the run name where the project belongs
    Return:
        run_name or False if does not exists
    '''
    return

def get_run_folder_from_user_project(project_id):
    '''
    Description:
        The function api return the run folder name for the requested project
    Input:
        project_id    # id of the project
    Return:
        run_folder or False
    '''
    if Projects.objects.filter(pk__exact = project_id).exists():
        run_obj = Projects.objects.get(pk__exact = project_id).get_run_obj()
        if RunningParameters.objects.filter(runName_id = run_obj).exists():
            return RunningParameters.objects.get(runName_id = run_obj).get_run_folder()
    return False

def get_user_projects (sharing_list):
    '''
    Description:
        The function api return the project keys and project names for the user or the sharing user list
    Input:
        sharing_list    # user list
    Return:
        project_key, project_names or False if does not exists
    '''
    project_list = []
    if Projects.objects.filter(user_id__in = sharing_list).exists():
        projects_objs = Projects.objects.filter(user_id__in = sharing_list).order_by('user_id').order_by('generatedat')
        for project in projects_objs:
            project_list.append([project.pk, project.get_project_name()])
    return project_list

def get_user_project_name (project_id_list):
    '''
    Description:
        The function api return the project names from the project_id_list
    Input:
        project_id_list     # list of project ids
    Return:
        project_names
    '''
    project_names = []
    for project_id in project_id_list:
        if Projects.objects.filter(pk__exact = project_id).exists():
            project_names.append(Projects.objects.get(pk__exact = project_id).get_project_name())
        else:
            project_names.append('')
    return project_names

def get_samples_projects (project_id_list):
    '''
    Description:
        The function api return a dictionnary having as key the project id and
        value a list of tupla sample_id, sample_name
    Input:
        project_id_list     # list of project ids
    Return:
        samples_projects
    '''
    samples_projects = {}
    for project_id in project_id_list:
        samples_projects[project_id] = []
        if Projects.objects.filter(pk__exact = project_id).exists():
            project_obj = Projects.objects.get(pk__exact = project_id)
            if SamplesInProject.objects.filter(project_id = project_obj).exists():
                samples_obj = SamplesInProject.objects.filter(project_id = project_obj)
                for sample_obj in samples_obj:
                    samples_projects[project_id].append([sample_obj.get_sample_id(), sample_obj.get_sample_name()])
    return samples_projects

def get_runs_projects_samples_and_dates(user_list_ids):
    '''
    Description:
        The function api return a dictionnary having as keys the run and the project and
        value a list of tupla sample_id, sample_name
    Input:
        user_list_ids     # user id list to get the list of samples
    Return:
        samples_data
    '''
    samples_data = []
    if SamplesInProject.objects.filter(user_id_id__in = user_list_ids).exists():
        #sample_objs = SamplesInProject.objects.filter(user_id_id__in = user_list_ids).order_by('generated_at').reverse().values_list('runProcess_id', 'runProcess_id', 'project_id','project_id','runProcess_id','sampleName','pk')
        #sample_objs = SamplesInProject.objects.filter(user_id_id__in = user_list_ids).order_by('generated_at').reverse()
        #samples_data = [[sample_obj.get_run_name(),sample_obj.get_project_name(), sample_obj.get_project_name(),sample_obj.get_project_id()]  for sample_obj in sample_objs]
        #samples_data = list(sample_objs)
        #import pdb; pdb.set_trace()
        #for pipo in SamplesInProject.objects.raw('select * from iSkyLIMS_wetlab_samplesinproject'):
        #    print (pipo.project_id.projectName)
        #import pdb; pdb.set_trace()
        query_mine =  SamplesInProject.objects.raw('select  iSkyLIMS_wetlab_samplesinproject.id , iSkyLIMS_wetlab_samplesinproject.sampleName , iSkyLIMS_wetlab_projects.projectName from iSkyLIMS_wetlab_samplesinproject inner join iSkyLIMS_wetlab_projects ON iSkyLIMS_wetlab_samplesinproject.project_id_id = iSkyLIMS_wetlab_projects.id')
        from django.db import connection
        cursor = connection.cursor()
        #cursor.execute('select  iSkyLIMS_wetlab_samplesinproject.id , iSkyLIMS_wetlab_runprocess.runName ,iSkyLIMS_wetlab_samplesinproject.sampleName , iSkyLIMS_wetlab_projects.projectName from iSkyLIMS_wetlab_samplesinproject inner join iSkyLIMS_wetlab_projects ON iSkyLIMS_wetlab_samplesinproject.project_id_id = iSkyLIMS_wetlab_projects.id inner join iSkyLIMS_wetlab_projects ON iSkyLIMS_wetlab_samplesinproject.runProcess_id = iSkyLIMS_wetlab_runprocess.id')
        #cursor.execute('select   iSkyLIMS_wetlab_runprocess.runName, iSkyLIMS_wetlab_runprocess.id ,iSkyLIMS_wetlab_projects.projectName , iSkyLIMS_wetlab_projects.id ,  iSkyLIMS_wetlab_samplesinproject.sampleName, iSkyLIMS_wetlab_samplesinproject.id  from iSkyLIMS_wetlab_samplesinproject inner join iSkyLIMS_wetlab_projects ON iSkyLIMS_wetlab_samplesinproject.project_id_id = iSkyLIMS_wetlab_projects.id inner join iSkyLIMS_wetlab_runprocess ON iSkyLIMS_wetlab_samplesinproject.runProcess_id_id = iSkyLIMS_wetlab_runprocess.id')

        #cursor.execute('select  iSkyLIMS_wetlab_samplesinproject.id , iSkyLIMS_wetlab_runprocess.runName ,iSkyLIMS_wetlab_samplesinproject.sampleName , iSkyLIMS_wetlab_projects.projectName from iSkyLIMS_wetlab_samplesinproject inner join iSkyLIMS_wetlab_projects ON iSkyLIMS_wetlab_samplesinproject.project_id_id = iSkyLIMS_wetlab_projects.id inner join iSkyLIMS_wetlab_projects ON iSkyLIMS_wetlab_samplesinproject.runProcess_id = iSkyLIMS_wetlab_runprocess.id')
        #cursor.execute('select   iSkyLIMS_wetlab_runprocess.runName, iSkyLIMS_wetlab_runprocess.id ,iSkyLIMS_wetlab_projects.projectName , iSkyLIMS_wetlab_projects.id ,  iSkyLIMS_wetlab_samplesinproject.sampleName, iSkyLIMS_wetlab_samplesinproject.id from iSkyLIMS_wetlab_samplesinproject inner join iSkyLIMS_wetlab_projects ON iSkyLIMS_wetlab_samplesinproject.project_id_id = iSkyLIMS_wetlab_projects.id inner join iSkyLIMS_wetlab_runprocess ON iSkyLIMS_wetlab_samplesinproject.runProcess_id_id = iSkyLIMS_wetlab_runprocess.id')
        user_list_str = [str(id_str) for id_str in user_list_ids]
        user_list = '( ' + ','.join(user_list_str) + ' )'
        fetch_fields = " iSkyLIMS_wetlab_runprocess.runName, iSkyLIMS_wetlab_runprocess.id ,iSkyLIMS_wetlab_projects.projectName , iSkyLIMS_wetlab_projects.id ,  iSkyLIMS_wetlab_samplesinproject.sampleName, iSkyLIMS_wetlab_samplesinproject.id, DATE_FORMAT(iSkyLIMS_wetlab_runprocess.run_completed_date,'%d/%m/%Y') AS niceDate "
        q_tables = ' FROM iSkyLIMS_wetlab_samplesinproject inner join iSkyLIMS_wetlab_projects ON iSkyLIMS_wetlab_samplesinproject.project_id_id = iSkyLIMS_wetlab_projects.id inner join iSkyLIMS_wetlab_runprocess ON iSkyLIMS_wetlab_samplesinproject.runProcess_id_id = iSkyLIMS_wetlab_runprocess.id '
        restrict_results = 'WHERE iSkyLIMS_wetlab_samplesinproject.user_id_id IN ' + user_list
        query = "SELECT " + fetch_fields + q_tables + restrict_results + 'ORDER BY (iSkyLIMS_wetlab_runprocess.run_completed_date) DESC'
        cursor.execute(query)
        #cursor.execute("select   iSkyLIMS_wetlab_runprocess.runName, iSkyLIMS_wetlab_runprocess.id ,iSkyLIMS_wetlab_projects.projectName , iSkyLIMS_wetlab_projects.id ,  iSkyLIMS_wetlab_samplesinproject.sampleName, iSkyLIMS_wetlab_samplesinproject.id, DATE_FORMAT(iSkyLIMS_wetlab_runprocess.run_completed_date,'%d/%m/%Y') AS niceDate from iSkyLIMS_wetlab_samplesinproject inner join iSkyLIMS_wetlab_projects ON iSkyLIMS_wetlab_samplesinproject.project_id_id = iSkyLIMS_wetlab_projects.id inner join iSkyLIMS_wetlab_runprocess ON iSkyLIMS_wetlab_samplesinproject.runProcess_id_id = iSkyLIMS_wetlab_runprocess.id WHERE iSkyLIMS_wetlab_samplesinproject.user_id_id IN ('1','2','3','7','8','9','16')")

        samples_data = cursor.fetchall()


        #data = ['1','2','3','4','5','6','7']
        #samples_data = [data for sample_obj in sample_objs]
            #run_obj = sample_obj.get_run_obj()
            #if run_obj.get_state() != 'Completed':
            #    continue
            #data = ['1','2','3','4','5','6','7']
            #data.append('run')
            #data.append('3')
            #data.append('project')
            #data.append('56')
            #data.append(sample_obj.get_run_name())
            #data.append(sample_obj.get_run_id())
            #data.append(sample_obj.get_project_name())
            #data.append(sample_obj.get_project_id())
            #data.append(sample_obj.get_sample_name())
            #data.append(sample_obj.get_sample_id())
            #data.append(run_obj.get_run_finish_date())
            #data.append('fecha')
            #samples_data.append(data)


    #import pdb; pdb.set_trace()
    return samples_data
