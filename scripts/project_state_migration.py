from iSkyLIMS_wetlab.models import Projects


def run ():
    
    projects_to_update = Projects.objects.all()
    p_to_update = len (projects_to_update)
    sucessful_counter = 0  
    unsucessful_projects = []
    for project in projects_to_update :
        state_of_project = project.procState
        if project.set_project_state(state_of_project):
            sucessful_counter += 1
        else:
            unsucessful_projects.append(project.get_project_name())
    if sucessful_counter != p_to_update:
        print('The following ', len(unsucessful_projects), ' projects were not updated:\n')
        for project in unsucessful_projects :
            print(project)
    else:   
        print( 'All projects are now updated ')
