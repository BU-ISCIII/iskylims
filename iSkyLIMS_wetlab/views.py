# -*- coding: utf-8 -*-
## import django
import statistics
import re, os, shutil, json
import datetime, time
from .fusioncharts.fusioncharts import FusionCharts

from django.shortcuts import get_object_or_404, render, redirect
from django.http import HttpResponse
from django.template import loader
from django.conf import settings
from django.core.files.storage import FileSystemStorage
from django.contrib.auth.models import User
from django.contrib.auth.decorators import login_required
#from django.contrib.auth.models import Group

from django_utils.models import Profile, Center
from django_utils.views import check_user_group
from .models import *

from iSkyLIMS_wetlab import wetlab_config
## import methods defined on utils.py
from .utils.sample_sheet_utils import *
from .utils.sample_functions import *
from .utils.stats_calculation import *
from .utils.stats_graphics import *
from .utils.generic_functions import *
from .utils.collection_index_functions import *
from .utils.fetching_information import *
from .utils.testing_wetlab_configuration import *
from .utils.sample_functions import *
from .utils.library_preparation import *
from .utils.pool_preparation import *
from .utils.run_preparation import  *
from .utils.additional_kits import *
from .utils.handling_statistics import *
from .utils.handling_sequencers import *
#from .utils.samplesheet_checks import *
#from .utils.wetlab_misc_utilities import normalized_data
from iSkyLIMS_core.utils.handling_samples import *
from iSkyLIMS_core.utils.handling_platforms import get_defined_platforms_and_ids
from iSkyLIMS_core.utils.generic_functions import get_inital_sample_settings_values, save_inital_sample_setting_value, send_test_email, get_email_data
from iSkyLIMS_core.utils.handling_protocols import display_protocol_list
#from iSkyLIMS_core.utils.handling_protocols import *
#from iSkyLIMS_core.utils.handling_commercial_kits import *


def index(request):
    #
    return render(request, 'iSkyLIMS_wetlab/index.html')

@login_required
def register_wetlab(request):
    #
    return render(request, 'iSkyLIMS_wetlab/index.html')


@login_required
def configuration_email(request):
    if request.user.username != 'admin':
        return redirect('/wetlab')
    email_conf_data = get_email_data()
    email_conf_data['EMAIL_ISKYLIMS'] = get_configuration_from_database('EMAIL_FOR_NOTIFICATIONS')
    if request.method == 'POST' and (request.POST['action']=='emailconfiguration'):
        result_email = send_test_email(request.POST)
        if result_email != 'OK':
            email_conf_data = get_email_data()
            email_conf_data['EMAIL_ISKYLIMS'] = request.POST['EMAIL_ISKYLIMS']
            email_conf_data['test_email'] = request.POST['test_email']
            return render(request, 'iSkyLIMS_wetlab/configurationEmail.html',{'ERROR':result_email, 'email_conf_data': email_conf_data})
        save_database_configuration_value('EMAIL_FOR_NOTIFICATIONS', request.POST['EMAIL_ISKYLIMS'])
        return render(request, 'iSkyLIMS_wetlab/configurationEmail.html',{'succesful_settings':True})
    return render(request, 'iSkyLIMS_wetlab/configurationEmail.html',{'email_conf_data': email_conf_data})

@login_required
def configuration_samba(request):
    if request.user.username != 'admin':
        return redirect('')
    samba_conf_data = get_samba_connection_data()
    if request.method == 'POST' and (request.POST['action']=='sambaconfiguration'):
        # reload configuration samba settings
        samba_user_field ={}
        for field in SAMBA_CONFIGURATION_FIELDS:
            samba_user_field[field] = request.POST[field]
        save_samba_connection_data(samba_user_field)
        try:
            open_samba_connection()
            return render(request, 'iSkyLIMS_wetlab/configurationSamba.html',{'succesful_settings':True})
        except:
            error_message = ERROR_WRONG_SAMBA_CONFIGURATION_SETTINGS
            return render(request, 'iSkyLIMS_wetlab/configurationSamba.html',{'samba_conf_data':samba_user_field, 'error_message': error_message} )
    else:
        return render(request, 'iSkyLIMS_wetlab/configurationSamba.html',{'samba_conf_data': samba_conf_data})

@login_required
def initial_settings(request):
    if not request.user.is_authenticated:
        #redirect to login webpage
        return redirect ('/accounts/login')

    if not is_wetlab_manager(request):
        return render (request,'iSkyLIMS_wetlab/error_page.html', {'content': ERROR_USER_NOT_WETLAB_MANAGER})
    initial_data = get_inital_sample_settings_values(__package__)
    form_data = {}
    if request.method == 'POST' and request.POST['action']=='defineNewSpecie':
        form_data['species'] = request.POST['specieName']
    if request.method == 'POST' and request.POST['action']=='defineNewLabRequest':
        form_data['lab_request'] = request.POST
    if request.method == 'POST' and request.POST['action']=='defineMoleculeType':
        form_data['molecule_type'] = request.POST['moleculeName']
    if request.method == 'POST' and request.POST['action']=='defineProtocolType':
        form_data['protocol_type'] = [request.POST['protocolName'], request.POST['moleculeType']]

    if form_data:
        new_inital_data = save_inital_sample_setting_value(__package__, form_data)
        if 'ERROR' in new_inital_data:
            return render(request,'iSkyLIMS_wetlab/initialSettings.html',{'initial_data': initial_data, 'ERROR': new_inital_data['ERROR']})
        return render(request,'iSkyLIMS_wetlab/initialSettings.html',{'initial_data': initial_data, 'new_setting_defined': new_inital_data})
    else:
        return render(request,'iSkyLIMS_wetlab/initialSettings.html',{'initial_data': initial_data})


@login_required
def create_nextseq_run (request):
    ## Check user == WETLAB_MANAGER: if false,  redirect to 'login' page
    if not request.user.is_authenticated:
        #redirect to login webpage
        return redirect ('/accounts/login')
    if not is_wetlab_manager(request):
        return render (request,'iSkyLIMS_wetlab/error_page.html', {'content': ERROR_USER_NOT_WETLAB_MANAGER})

    ## FIRST STEP in collecting data from the NextSeq run. Sample Sheet and experiment name are required
    if request.method == 'POST' and (request.POST['action']=='uploadFile'):
        get_user_names={}
        projects=[]
        #run_name=request.POST['runname']
        myfile = request.FILES['myfile']

        ## CHECK if file contains the extension.
        ## Error page is showed if file does not contain any extension
        split_filename=re.search('(.*)(\.\w+$)',myfile.name)
        if None == split_filename:
            return render (
                request,'iSkyLIMS_wetlab/error_page.html',
                {'content':['Uploaded file does not containt extension',
                            'Sample Sheet must have a csv extension', '', 'ADVICE:',
                            'Select the Sample file generated by Illumina Experient Manager (IEM)']})
        ext_file=split_filename.group(2)

        ## CHECK if file contains the csv extension.
        ## Error page is shown if file does not contain the csv extension
        if ext_file != '.csv':
            return render (
                request,'iSkyLIMS_wetlab/error_page.html',
                {'content':['Sample Sheet must have a csv extension', '', 'ADVICE:',
                            'Select the Sample file generated by Illumina Experient Manager (IEM)']})
        fs = FileSystemStorage()
        timestr = time.strftime("%Y%m%d-%H%M%S")
        ## including the timestamp to the sample sheet file
        # do not need to include the absolute path because django uses
        # the MEDIA_ROOT variable defined on settings to upload the file
        file_name=str(wetlab_config.RUN_SAMPLE_SHEET_DIRECTORY
                     + split_filename.group(1) + '_' + timestr + ext_file)
        filename = fs.save(file_name,  myfile)
        uploaded_file_url = fs.url(filename)

        ### add the document directory to read the csv file
        stored_file = os.path.join(settings.MEDIA_ROOT, file_name)

        ## Fetch the experiment name and the library name from the sample sheet file
        index_library_name = get_index_library_name(stored_file)
        run_name = get_experiment_name_from_file(stored_file)

        if run_name == '':
            ## define an temporary unique value for the run name
            #until the real value is get from user FORM
            run_name = timestr


        ## Check that runName is not already used in the database.
        ## Error page is showed if runName is already  defined
        if (RunProcess.objects.filter(runName__iexact = run_name)).exists():
            if RunProcess.objects.filter(runName__iexact = run_name, state__runStateName__exact ='Pre-Recorded').exists():
                ## Delete the Sample Sheet file and the row in database
                delete_run_objs = RunProcess.objects.filter(runName__iexact = run_name, state__runStateName__exact ='Pre-Recorded')
                for delete_run in delete_run_objs:
                    #sample_sheet_file = delete_run.get_sample_file()
                    ##full_path_sample_sheet_file = os.path.join(settings.MEDIA_ROOT, sample_sheet_file)
                    #os.remove(full_path_sample_sheet_file)

                    if Projects.objects.filter(runProcess = delete_run).exists():
                        project_objs = Projects.objects.filter(runProcess = delete_run)
                        for project_obj in project_objs:
                            project_obj.runProcess.remove(delete_run)

                            if project_obj.runProcess.all().count() == 0 :
                                project_obj.delete()
                    delete_run.delete()

            else:
                # delete sample sheet file
                os.remove(stored_file)
                return render (request,'iSkyLIMS_wetlab/error_page.html',
                    {'content':['Run Name is already used. ',
                        'Run Name must be unique in database.',' ',
                        'ADVICE:','Change the value in the Sample Sheet  file ']})

        ## Fetch from the Sample Sheet file the projects included in
        ## the run and the user. Error page is showed if not project/description
        ## colunms are found

        project_list = get_projects_in_run(stored_file)

        if len (project_list) == 0 :
            ## delete sample sheet file
            fs.delete(file_name)
            return render (request,'iSkyLIMS_wetlab/error_page.html',
                {'content':['Sample Sheet does not contain "Sample project" and/or "Description" fields',
                    '','ADVICE:','Check that csv file generated by Illumina Experient Manager (IEM) includes these columns']})

        ## Check if the projects are already defined on database.
        ## Error page is showed if projects are already defined on database

        project_already_defined=[]
        project_in_several_runs = get_configuration_value('PROJECTS_ALLOWED_IN_MULTIPLE_RUNS')
        for key  in project_list.keys():
            # check if project was already saved in database in Not Started State.
            # if found delete the projects, because the previous attempt to complete the run was unsuccessful
            if ( Projects.objects.filter(projectName__icontains = key).exists()):
                if project_in_several_runs != 'TRUE':
                    project_already_defined.append(key)

        if (len(project_already_defined) >0 ):
            if (len(project_already_defined)>1):
                head_text='The following projects are already defined in database:'
            else:
                head_text='The following project is already defined in database:'
            ## convert the list into string to display the user names on error page
            display_project= '  '.join(project_already_defined)
            ## delete sample sheet file before showing the error page
            fs.delete(file_name)
            return render (request,'iSkyLIMS_wetlab/error_page.html',
                {'content':[ head_text,'', display_project,'',
                    'Project names must be unique','', 'ADVICE:',
                    'Edit the Sample Sheet file to correct this error']})

        ##Once the information looks good. it will be stores in runProcess and projects table

        ## store data in runProcess table, run is in pre-recorded state
        center_requested_id = Profile.objects.get(profileUserID = request.user).profileCenter.id
        center_requested_by = Center.objects.get(pk = center_requested_id)
        new_run_obj = RunProcess(runName=run_name,sampleSheet= file_name,
                                state = RunStates.objects.get(runStateName__exact = 'Pre-Recorded'),
                                centerRequestedBy = center_requested_by)
        new_run_obj.save()
        experiment_name = '' if run_name == timestr else run_name

        ## create new project tables based on the project involved in the run and
        ## include the project information in projects variable to build the new FORM

        run_info_values ={}
        run_info_values['experiment_name'] = experiment_name
        run_info_values['index_library_name'] = index_library_name
        for key, val  in project_list.items():
            if User.objects.filter(username__exact = val).exists():
                user_id = User.objects.get(username__exact = val)
            else:
                user_id = None
            '''
            p_data = Projects(runprocess_id=RunProcess.objects.get(runName =run_name),
                            projectName=key, user_id=userid)
            p_data.save()
            '''
            if not Projects.objects.filter(projectName__iexact = key).exists():
                data = {}
                data['user_id'] = user_id
                data['projectName'] = key
                project_obj = Projects.objects.create_new_empty_project(data)
            else:
                project_obj = Projects.objects.filter(projectName__iexact = key).last()

            project_obj.add_run(new_run_obj)
            projects.append([key, val])

        run_info_values['projects_user'] = projects
        run_info_values['runname']= run_name
        ## Get the list of the library kit used (libraryKit)
        used_libraries = []
        list_libraries = LibraryKit.objects.order_by().values_list('libraryName', flat=True)
        run_info_values['used_libraryKit'] =  list_libraries

        user_names = []
        all_users = User.objects.all()
        for user in all_users :
            user_names.append(user.username)
        run_info_values['aval_users'] =  user_names
        ## displays the list of projects and the user names found on Sample Sheet
        return render(request, 'iSkyLIMS_wetlab/CreateNextSeqRun.html', {'get_user_names': run_info_values })

    ## SECOND STEP in collecting data from the NextSeq run. Confirmation /modification of data included in Sample Sheet
    elif request.method=='POST' and (request.POST['action']=='displayResult'):
        experiment_name = request.POST['experimentname']
        run_index_library_name = request.POST['runindexlibraryname']
        run_name= request.POST['runname']
        projects=request.POST.getlist('project')
        user_name=request.POST.getlist('username')
        library_kit=request.POST.getlist('libraryKit')
        project_index_kit=request.POST.getlist('projectindexlibraryname')

        ## get the sample sheet used in the run. return error if run already exists
        if not RunProcess.objects.filter (runName__exact = run_name).exists():
            return render (
                request, 'iSkyLIMS_wetlab/error_page.html',
                {'content':['You get this error page because you use the back Buttom'
                ' to return to previous page where asking for library kit name',
                'To upload again the shample sheet, use the "Upload the Run" option from the top menu']})
        run_p = RunProcess.objects.get(runName__exact = run_name)
        s_file=run_p.get_sample_file()
        ## get the different type of library kit used in the run and
        ## convert the sample sheet into Base Space. Number of converted
        ## file will be the same as the number of different lybraries use in the run
        library={}
        bs_file={}
        results=[]

        in_file = os.path.join(settings.MEDIA_ROOT,s_file)
        # Set unique Sample_ID in the sample sheet

        index_file = os.path.join(settings.MEDIA_ROOT,'wetlab', 'index_file')
        # if file does not exists create the file and assing the first value
        if not os.path.isfile(index_file) :
            with open(index_file, 'w') as index_fh:
                index_fh.write('0000-AA')
        create_unique_sample_id_values (in_file, index_file)
        # create the projects/users to update sample sheet
        user_names_in_projects ={}
        for p_index in range(len(projects)):
            user_names_in_projects[projects[p_index]] = user_name[p_index]

        set_user_names_in_sample_sheet (in_file, user_names_in_projects)
        ## build the project list for each project_library kit
        for x in range(len(project_index_kit)):
            if project_index_kit[x] in library :
                library[project_index_kit[x]].append(projects[x])
            else:
                library[project_index_kit[x]]= [projects[x]]
        ## convert the sample sheet to base space format and have different files according the library kit

        for key, value in library.items():
            lib_kit_file =key.replace(' ', '_')
            library_file = sample_sheet_map_basespace(in_file, key, lib_kit_file, value,'Plate96')
            if library_file == 'ERROR':
                # deleting the sample sheet file
                os.remove(in_file)
                # Deleting projects related to the shample sheet
                for p in range(len( projects)):
                    my_project = projects [p]
                    delete_proj=Projects.objects.get(projectName = my_project)
                    delete_proj.delete()
                # delete the run used when uploading the sample sheet
                run_p.delete()
                # show the error page
                return render (
                    request,'iSkyLIMS_wetlab/error_page.html',
                    {'content':[ 'The information on  the Library kit ', key,
                    ' For the project ', value,
                    'Does not meet the requirements to perform the conversion to import to Base Space',
                    'ADVICE', 'Check the sample sheet that was uploaded ']})
            else:
                bs_file[key] = library_file
                results.append([key, bs_file[key]])

        ## save the project information on database

        for p in range(len( projects)):
            my_project = projects [p]
            my_name = user_name[p]
            my_libkit = library_kit[p]
            library_kit_id = LibraryKit.objects.filter(libraryName__exact = library_kit[p]).last()
            update_info_proj=Projects.objects.get(projectName = my_project)
            update_info_proj.libraryKit=project_index_kit[p]
            update_info_proj.baseSpaceFile=bs_file[project_index_kit[p]]
            update_info_proj.LibraryKit_id = library_kit_id
            update_info_proj.user_id = User.objects.get(username__exact = user_name[p])
            update_info_proj.save()
        results.append(['runname', experiment_name])
        ## save the sample sheet file under tmp/recorded to be processed when run folder was created
        '''
        ########
        ## Remove the folder creation to include the sample sheet from iSkyLIMS version 2
        subfolder_name=str(run_p.id)

        temp_directory = os.path.join(settings.MEDIA_ROOT , wetlab_config.RUN_TEMP_DIRECTORY_RECORDED, subfolder_name)
        if not os.path.isdir(temp_directory):
            os.makedirs(temp_directory)
        # os.mkdir(temp_directory)
        # set group writing permission to the temporary directory
        os.chmod(temp_directory, 0o774)
        #os.mkdir(os.path.join(settings.MEDIA_ROOT, 'wetlab/tmp/recorded', subfolder_name ))
        sample_sheet_copy= os.path.join(temp_directory, 'samplesheet.csv' )
        shutil.copy(in_file,sample_sheet_copy)
        # set the group write permission to the Sample Sheet File
        os.chmod(sample_sheet_copy, 0o664)
        # update the sample sheet with the experiment name
        if run_name != experiment_name :
            update_sample_sheet (in_file, experiment_name)
        ## update the Experiment name and the state of the run to 'Recorded'
        run_p.runName = experiment_name
        run_p.index_library = run_index_library_name
        run_p.save()
        #######################
        '''
        run_p.set_run_state ('Recorded')
        sample_sheet_lines = read_all_lines_in_sample_sheet(in_file)
        sample_names_and_data = get_samples_in_sample_sheet(sample_sheet_lines)
        samples_reused = increase_reuse_if_samples_exists(sample_names_and_data['samples'])

        return render (request, 'iSkyLIMS_wetlab/CreateNextSeqRun.html', {'completed_form':results})

    return render(request, 'iSkyLIMS_wetlab/CreateNextSeqRun.html')

@login_required
def add_basespace_library (request):
    '''
    Description:
        The function is called from web, having 2 main parts:
            - User form with the information to add a new Basespace library
            - Result information as response of user submit
        Save a new basespace library name in database if it is not already defined.
    Input:
        request     # contains the request dictionary sent by django
    Variables:
        basespace_library_information ={} # returned dictionary with the information
                                to include in the web page
        basespace_library_objects # contains the object list of the basespace model
        basespace_library = [] # It is a list containing the Basespces Library names
        new_basespace_library_name # contain the new library name enter by user form
        library     # it is the new LibraryKit object
        l_kit       # is the iter variable for basespace_library_objects
    Return:
        Return the different information depending on the execution:
        -- Error page in case the library already exists.
        -- library_kit_information with :
            -- ['libraries']
            ---['new_basespace_library'] in case a new basespace library was added.
    '''

    libraries_information ={}
    basespace_library_information ={}
    basespace_library = []

    basespace_library_objects = LibraryKit.objects.all()
    if len(basespace_library_objects) >0 :
        for l_kit in basespace_library_objects :
            basespace_library.append(l_kit.libraryName)

    if request.method == 'POST' and request.POST['action'] == 'addNewBasespaceLibrary':
        new_basespace_library_name = request.POST['newBasespaceLibrary']

        ## Check that library kit is not already defined in database
        if LibraryKit.objects.filter(libraryName__icontains = new_basespace_library_name).exists():
            return render (request, 'iSkyLIMS_wetlab/error_page.html', {'content':['The Library Kit ', new_basespace_library_name, 'is already defined on the system']})

        basespace_library_information['new_basespace_library'] = new_basespace_library_name
        basespace_library.append(new_basespace_library_name)
        #save the new library on database
        b_library = LibraryKit(libraryName= new_basespace_library_name)
        b_library.save()

    basespace_library_information ['libraries'] = basespace_library
    return render(request,'iSkyLIMS_wetlab/AddBasespaceLibrary.html',{'list_of_libraries': basespace_library_information})

@login_required
def add_collection_index_kit (request):
    #get the list of the already loaded index library to be displayed
    '''
    Description:
        The function is called from web, having 2 main parts:
            - User form with the information to add a new library
            - Result information as response of user submit

    Input:
        request     # contains the request dictionary sent by django
    Variables:
        index_libraries_information     # returned dictionary with the information
                                        to include in the web page
        index_library_objects    #  contains the object list of the IndexLibraryKit model
        index_library_names     # It is a list containing the Index Library Kits names
        index_to_store      # contains the index (I7/I5) values that are stored in database
                            the same variable is used for the interaction for library_index

        library     # it is the new IndexLibraryKit object
        library_settings # settings values returned by
        l_kit       # is the iteration variable for index_library_objects
        lib_settings_to_store # is the new IndexLibraryKit object used to store
                        the information into database

        index_library_file      # contains the file provider by user in the form
        fs_index_lib    # file system index library object to store the input file
        saved_file      # contain the full path name, where the user file have been
                        stored in the server

    Constants:
        LIBRARY_KITS_DIRECTORY
        LIBRARY_MAXIMUM_SIZE
        MEDIA_ROOT
    Functions:
        in utils.library_kits :
            -- check_index_library_file_format(saved_file)     # for checking
                            the number of index in the input file
            -- getting_index_library_name(saved_file) # gets the library name
            -- get_library_settings(saved_file) # gets the settings values
                            from the input file
            -- get_index_values(saved_file)     # gets the index value

    Return:
         Return the different information depending on the execution:
        -- Error page in case of:
            -- Uploaded file is bigger than the LIBRARY_MAXIMUM_SIZE value
            -- file uploaded does not have the right format
            -- the library already exists.
        -- library_kit_information with :
            -- ['libraries']
            ---['new_library_kit'] in case a new library kit was added
    '''
    collection_index_information ={}
    collection_index_names = []

    collection_indexes = CollectionIndexKit.objects.all()
    if len(collection_indexes) > 0 :
        for c_index in collection_indexes :
            collection_index_names.append([c_index.get_id(), c_index.get_collection_index_name])

    if request.method == 'POST' and request.POST['action'] == 'addCollectionIndexKit':
        ## fetch the file from user form and  build the file name  including
        ## the date and time on now to store in database

        file_name = request.FILES['newCollectionIndexFile'].name
        saved_file = store_collection_kits_file(request.FILES['newCollectionIndexFile'])

        ## get the libary name to check if it is already defined
        if not check_collection_index_file_format(saved_file):
            os.remove(saved_file)
            return render (request, 'iSkyLIMS_wetlab/error_page.html',
                           {'content':['The Collection Index Kit file', file_name,
                                       'does not have the right format']})

        collection_name = get_collection_index_name(saved_file)
        if collection_name == '' :
            # removing the uploaded file
            os.remove(saved_file)
            return render (request, 'iSkyLIMS_wetlab/error_page.html',
                        {'content':['The Collection Index Kit file', file_name,
                                   'does not contain the  name']})

        # check if library name is already defined on database
        if check_collection_index_exists(collection_name) :
            # removing the uploaded file
            os.remove(saved_file)
            return render (request, 'iSkyLIMS_wetlab/error_page.html',
                           {'content':['The Collection Index Kit Name ', file_name,
                                       'is already defined on iSkyLIMS']})
        # Get the collection settings included in the file
        collection_settings = get_collection_settings(saved_file)
        new_collection_obj = store_collection_settings (collection_settings, file_name)
        ## get the index name and index bases for the library
        collection_index = get_index_values(saved_file)
        store_collection_indexes(collection_index, new_collection_obj)

        collection_index_information['collection_index_names'] = collection_settings['name']
        collection_index_information ['collection_index'] = collection_index_names

        return render (request, 'iSkyLIMS_wetlab/addCollectionIndexKit.html',{'collection_index_information': collection_index_information })
    else:
        collection_index_information ['collection_index'] = collection_index_names
        return render (request, 'iSkyLIMS_wetlab/addCollectionIndexKit.html',{'list_of_collection_index': collection_index_information })



@login_required
def search_run (request):
    '''
    Description:
        The function is called from web, having 2 main parts:
            - User form with the information to search runs
            - Result information can be :
                - list of the matched runs
                - run information in case that only 1 match is found
    Input:
        request     # contains the request dictionary sent by django
    Imports:
        Machines and Platform are imported from iSkyLIMS_drylab.models
            for filtering runs based on the platform
    Constants:
    ERROR_NO_MATCHES_FOR_RUN_SEARCH
    Functions:
        get_information_run() # Collects information about one run
    Variables:
        User inputs from search options
            run_name        # string characters to find in the run name
            platform_name   # platform name filter
            run_state       # state of the run
            start_date      # filter of starting date of the runs
            end_date        # filter for the end of the runs
        available_platforms # contains the list of platform defined in
                            # iSkyLIMS.models.Platform
        machine_list        # list of machines to filter on the matches runs
        platforms           # contain the object from iSkyLIMS.models.Platform
        platform_name       # has the platform get from user form

        runs_found          # runProcess object that contains the result query
                            # it is updated with the user form conditions
        r_data_display      # contains the information to display about the run
        run_list            # contains the run list that mathches te user conditions
    Return:
        Return the different information depending on the execution:
        -- Error page in case no run is founded on the matching conditions.
        -- SearchRun.html is returned with one of the following information :
            -- r_data_display   # in case that only one run is matched
            ---run_list         # in case several run matches the user conditions.

    '''
    # check user privileges
    if request.user.is_authenticated:
        try:
            groups = Group.objects.get(name=wetlab_config.WETLAB_MANAGER)
            if groups not in request.user.groups.all():
                allowed_all_runs = False
               #return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
            else:
                allowed_all_runs = True
        except:
            return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
    else:
        #redirect to login webpage
        return redirect ('/accounts/login')
    #############################################################
    ## Search for runs that fullfil the input values
    #############################################################
    run_form_data = get_run_search_fields_form()
    error_message = ERROR_NO_MATCHES_FOR_RUN_SEARCH
    if request.method == 'POST' and (request.POST['action'] == 'runsearch'):
        run_name = request.POST['runname']
        start_date = request.POST['startdate']
        end_date = request.POST['enddate']
        run_state = request.POST['runstate']
        platform_name = request.POST['platform']
        # check that some values are in the request if not return the form
        if run_name == '' and start_date == '' and end_date == '' and run_state == '' and platform_name == '' :
            return render(request, 'iSkyLIMS_wetlab/SearchRun.html')

        ### check the right format of start and end date
        if start_date != '':
            if not check_valid_date_format(start_date) :
                error_message = ERROR_INVALID_FORMAT_FOR_DATES
                return render(request, 'iSkyLIMS_wetlab/SearchRun.html', {'run_form_data': run_form_data, 'error_message' : error_message})
        if end_date != '':
            if not check_valid_date_format(end_date) :
                error_message = ERROR_INVALID_FORMAT_FOR_DATES
                return render(request, 'iSkyLIMS_wetlab/SearchRun.html', {'run_form_data': run_form_data, 'error_message' : error_message})

        ### Get all the available runs to start the filtering
        if allowed_all_runs :
            runs_found=RunProcess.objects.all().order_by('run_date').reverse()
        else:

            user_projects = Projects.objects.filter(user_id__exact = request.user.id)
            run_list =[]
            for user_project in user_projects :
                run_list.append(user_project.runprocess_id.id)
            if RunProcess.objects.filter(pk__in = run_list).exists():
                runs_found = RunProcess.objects.filter(pk__in = run_list)
            else:
                error_message = ['There are not run where ' + request.user.username + ' was involved' ]
                return render(request, 'iSkyLIMS_wetlab/SearchRun.html', {'run_form_data': run_form_data, 'error_message' : error_message})

        ### Get runs when run name is not empty
        if run_name !='':
            if (RunProcess.objects.filter(runName__iexact =run_name).exists()):
                run_name_found=RunProcess.objects.filter(runName__iexact =run_name)
                if len(run_name_found) == 1:
                    r_data_display= get_information_run(run_name_found[0])
                    return render(request, 'iSkyLIMS_wetlab/SearchRun.html', {'display_one_run': r_data_display })
            if (runs_found.filter(runName__icontains =run_name).exists()):
                runs_found=runs_found.filter(runName__icontains =run_name).order_by('runName')
            else:
                return render(request, 'iSkyLIMS_wetlab/SearchRun.html', {'run_form_data': run_form_data, 'error_message' : error_message})
        if platform_name != '' :
            sequencer_list = get_sequencer_names_from_platform(platform_name)
            if len(sequencer_list) > 0:
                runs_found = runs_found.filter(usedSequencer__sequencerName__in = sequencer_list)
            else:
                return render(request, 'iSkyLIMS_wetlab/SearchRun.html', {'run_form_data': run_form_data, 'error_message' : error_message})
        ### Check if state is not empty
        if run_state != '':
            s_state = RunStates.objects.get(runStateName__exact = run_state)
            if runs_found.filter(state__runStateName__exact = s_state).exists():
                runs_found = runs_found.filter(state__runStateName__exact = s_state).order_by('runName')
            else :
                return render(request, 'iSkyLIMS_wetlab/SearchRun.html', {'run_form_data': run_form_data, 'error_message' : error_message})
        ### Check if start_date is not empty
        if start_date !='' and end_date != '':

            if runs_found.filter(run_date__range=(start_date, end_date)).exists():
                 runs_found = runs_found.filter(run_date__range=(start_date, end_date))
            else:
                return render(request, 'iSkyLIMS_wetlab/SearchRun.html', {'run_form_data': run_form_data, 'error_message' : error_message})
        if start_date !='' and end_date == '':
            if runs_found.filter(run_date__gte = start_date).exists():
                 runs_found = runs_found.filter(run_date__gte = start_date)
            else:
                return render(request, 'iSkyLIMS_wetlab/SearchRun.html', {'run_form_data': run_form_data, 'error_message' : error_message})
        if start_date =='' and end_date != '':
            if runs_found.filter(run_date__lte = end_date).exists():
                 runs_found = runs_found.filter(run_date__lte = end_date)
            else:
                return render(request, 'iSkyLIMS_wetlab/SearchRun.html', {'run_form_data': run_form_data, 'error_message' : error_message})

        #If only 1 run mathes the user conditions, then get the project information

        if (len(runs_found)== 1) :
            return redirect ('display_run', run_id=runs_found[0].pk)

        else:
            ## collect the list of run that matches the run date
            run_list=[]
            for run_found in runs_found:
                run_list.append([run_found.get_run_id(),run_found.get_run_name(), run_found.get_run_date() ])
            return render(request, 'iSkyLIMS_wetlab/SearchRun.html', {'display_run_list': run_list })
    else:
        return render(request, 'iSkyLIMS_wetlab/SearchRun.html', {'run_form_data': run_form_data})

@login_required
def search_project (request):
    '''
    Description:
        The function is called from web, having 2 main parts:
            - User form with the information to search projects
            - Result information can be :
                - list of the matched projects
                - project information in case that only 1 match is found

    Input:
        request     # contains the request dictionary sent by django
    Imports:
        Machines and Platform are imported from iSkyLIMS_drylab.models
            for filtering runs based on the platform

    Return:
        Return the different information depending on the execution:
        -- Error page in case no run is founded on the matching conditions.
        -- SearchRun.html is returned with one of the following information :
            -- r_data_display   # in case that only one run is matched
            ---run_list         # in case several run matches the user conditions.

    '''
    project_form_data = get_project_search_fields_form()
    error_message = ERROR_NO_MATCHES_FOR_PROJECT_SEARCH
    if request.method=='POST' and (request.POST['action']=='searchproject'):
        project_name=request.POST['projectname']
        start_date=request.POST['startdate']
        end_date=request.POST['enddate']
        user_name = request.POST['username']
        sequencer_name = request.POST['sequencer']
        run_state = request.POST['runstate']
        run_process_ids = []
        # check that some values are in the request if not return the form
        if project_name == '' and start_date == '' and end_date == '' and user_name =='' and sequencer_name == '' and run_state == '':
            return render(request, 'iSkyLIMS_wetlab/SearchProject.html', {'project_form_data': project_form_data})

        if user_name !=''  and len(user_name) < 5 :
            error_message = ERROR_USER_NAME_TOO_SHORT
            return render(request, 'iSkyLIMS_wetlab/SearchProject.html', {'project_form_data': project_form_data,'error_message':error_message})

        ### check the right format of start and end date
        if start_date != '':
            if not check_valid_date_format(start_date) :
                error_message = ERROR_INVALID_FORMAT_FOR_DATES
                return render(request, 'iSkyLIMS_wetlab/SearchProject.html', {'project_form_data': project_form_data,'error_message':error_message})

        if end_date != '':
            if not check_valid_date_format(start_date) :
                error_message = ERROR_INVALID_FORMAT_FOR_DATES
                return render(request, 'iSkyLIMS_wetlab/SearchProject.html', {'project_form_data': project_form_data,'error_message':error_message})

        projects_found = Projects.objects.all()

        if project_name != '':
            projects_found = projects_found.filter(projectName__icontains = project_name)
        if sequencer_name != '':
            run_objs = RunProcess.objects.filter(usedSequencer__sequencerName__exact = sequencer_name)
            projects_found = projects_found.filter(runProcess__in = run_objs)
        if (run_state !='' ):
            run_objs = RunProcess.objects.filter(state__runStateName__exact = run_state)
            projects_found = projects_found.filter(runProcess__in = run_objs)
        if user_name != '':
            projects_found = projects_found.filter(user_id__username__icontains = user_name)
        if start_date !='' and end_date != '':
            projects_found = projects_found.filter(generatedat__range=(start_date, end_date))
        if start_date !='' and end_date == '':
            projects_found = projects_found.filter(generatedat__gte = start_date)
        if start_date =='' and end_date != '':
            projects_found = projects_found.filter(generatedat__lte = end_date)
        if len(projects_found) == 0:
            error_message = ERROR_NO_MATCHES_FOR_PROJECT_SEARCH
            return render(request, 'iSkyLIMS_wetlab/SearchProject.html', {'project_form_data': project_form_data,'error_message':error_message})

        if len (projects_found) == 1:
            return redirect ('display_project', project_id = projects_found[0].id)
        else :
            # Display a list with all projects that matches the conditions
            project_list_dict = {}
            project_list = []
            for project in projects_found :
                p_name = project.get_project_name()
                p_name_id = project.id
                project_list.append([p_name, p_name_id])
            project_list_dict ['projects'] = project_list
            return render(request, 'iSkyLIMS_wetlab/SearchProject.html', {'display_project_list': project_list_dict })

    else:
        return render(request, 'iSkyLIMS_wetlab/SearchProject.html', {'project_form_data': project_form_data})



@login_required
def retry_error_run (request):
    # check user privileges
    if request.user.is_authenticated:

        try:
            groups = Group.objects.get(name='WetlabManager')
            if groups not in request.user.groups.all():
                return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
        except:
            return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
    else:
        #redirect to login webpage
        return redirect ('/accounts/login')
    if request.method=='POST' and (request.POST['action']=='retry_correct_error'):
        run_id = request.POST['run_id']
        if (RunProcess.objects.filter(pk__exact=run_id).exists()):
            run_name_found = RunProcess.objects.get(pk__exact=run_id)
            previous_error_state = run_name_found.get_state_before_error()
            run_name_found.set_run_state(previous_error_state)
            detail_description = {}
            detail_description['information'] = SUCCESSFUL_RUN_STATE_CHANGE_FOR_RETRY
            return render (request,'iSkyLIMS_wetlab/successful_page.html', {'detail_description': detail_description , 'return_main_menu': True})
        else:
            return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['Run does not exist ']})
    else:
        #return redirect (request,'/')
        return render(request, 'iSkyLIMS_wetlab/index.html')


@login_required
def display_run (request, run_id):
    # check user privileges
    if request.user.is_authenticated:
        try:
            groups = Group.objects.get(name = wetlab_config.WETLAB_MANAGER)
            if groups not in request.user.groups.all():
                # check if user is owner of the run
                if Projects.objects.filter(runprocess_id__exact = run_id).exists():
                    projects = Projects.objects.filter(runprocess_id__exact = run_id)
                    user_list =[]
                    for project in projects:
                        user_list.append(project.user_id.id)
                    if  not request.user.id in user_list :
                        return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
                else:
                    return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['No matches have been found for the run  ']})

        except:
            return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
    else:
        #redirect to login webpage
        return redirect ('/accounts/login')
    if (RunProcess.objects.filter(pk=run_id).exists()):
        run_name_found = RunProcess.objects.get(pk=run_id)
        r_data_display  = get_information_run(run_name_found)
        return render(request, 'iSkyLIMS_wetlab/displayRun.html', {'display_one_run': r_data_display })
    else:
        return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['No matches have been found for the run  ']})

@login_required
def last_run_by_sequencer (request) :
    # check user privileges
    if request.user.is_authenticated:

        try:
            groups = Group.objects.get(name='WetlabManager')
            if groups not in request.user.groups.all():
                return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
        except:
            return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
    else:
        #redirect to login webpage
        return redirect ('/accounts/login')

    last_runs = get_last_runs_by_sequencer()
    if len(last_runs) == 0:
        return render (request, 'iSkyLIMS_wetlab/lastRunBySequencer.html', {'no_runs': 'no_runs'})
    if len(last_runs) > 1:
        return render (request, 'iSkyLIMS_wetlab/lastRunBySequencer.html', {'last_runs': last_runs})
    else:
        # if only 1 sequencer is defined, then display the information of the latest run
        return redirect ('display_run', run_id = last_runs[0][2])

@login_required
def incompleted_runs (request) :
    # check user privileges
    if request.user.is_authenticated:
        try:
            groups = Group.objects.get(name='WetlabManager')
            if groups not in request.user.groups.all():
                return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
        except:
            return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
    else:
        #redirect to login webpage
        return redirect ('/accounts/login')
    if RunProcess.objects.all().exclude(state__runStateName = 'Completed').exists() :

        display_incompleted_run = get_information_for_incompleted_run()
        return render (request, 'iSkyLIMS_wetlab/incompletedRuns.html',{'display_incompleted_run':display_incompleted_run})
        #unfinished_runs = RunProcess.objects.all().exclude(state__runStateName = 'Completed').order_by('runName')
        #for run in unfinished_runs:
        #display_incomplete_run_list[run.id] = [[run.runName, run.get_state()]]
    else:
        return render (request,'iSkyLIMS_wetlab/info_page.html', {'content':['There is no project in incompleted state' , 'All Runs are finished']})


def check_user_access (request, project_found_id ) :

    groups = Group.objects.get(name = wetlab_config.WETLAB_MANAGER)
    # check if user belongs to WetlabManager . If true allow to see the page
    if groups not in request.user.groups.all():
        #check if project belongs to the same user as the one requesting the page
        if project_found_id.user_id.id != request.user.id :
           return False
    return True



@login_required
def display_project (request, project_id):

    if (Projects.objects.filter(pk=project_id).exists()):
        project_found_id = Projects.objects.get(pk=project_id)
        if request.user.is_authenticated:
            # check that user is allow to make the change
            allowed_access = check_user_access (request, project_found_id)
            if not allowed_access :
                return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
        else:
            #redirect to login webpage
            return redirect ('/accounts/login')
        # Display the proyect information
        display_project_data  = get_information_project(project_found_id, request)
        return render(request, 'iSkyLIMS_wetlab/displayProject.html', {'display_project_data': display_project_data })
    else:
        return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['No matches have been found for the project  ' ]})
'''
@login_required
def display_sample_in_run (request, sample_run_project_id):
'''
'''
    Description:
        The function will check if the requested sample id exists, then
        it will call to get_info_sample_in_run function to collect all information
    Input:
        request     # contains the request dictionary sent by django
        sample_id   # contains the sample id to display the information
    Variables:
        sample_data_information ={} # returned dictionary with the information
                                to include in the web page
        sample_found_id  # contains the object for the sample id
    Functions:
        get_info_sample_in_run (sample_found_id)
    Return:
        Return the different information depending on the execution:
        -- Error page in case the sample id in the request does not exists.
        -- sample_data_information with the information collected by get_info_sample_in_run()
'''
'''
    if (SamplesInProject.objects.filter(pk=sample_run_project_id).exists()):
        sample_found_obj = SamplesInProject.objects.get(pk=sample_run_project_id)
        sample_data_information = get_info_sample_in_run (sample_found_obj)
        return render(request, 'iSkyLIMS_wetlab/displaySampleInRun.html',{'display_one_sample': sample_data_information })
    else:
        return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['No matches have been found for the sample  ' ]})
'''



@login_required
def display_collection_index (request, collection_index_id):

    if (CollectionIndexKit.objects.filter(pk=collection_index_id).exists()) :
        collection_index_dict = get_collection_index_information (collection_index_id)
        if collection_index_dict != False:
            return render (request, 'iSkyLIMS_wetlab/DisplayCollectionIndex.html', {'display_one_collection_index': collection_index_dict})
        else:
            return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['There are recorded information for the collection index for your request']})

    else:
        return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['No matches have been found for the Collection index ' ]})




@login_required
def search_collection_index_library (request):

    if request.method == 'POST' and (request.POST['action'] == 'searchcollectionindexkit') :
        collection_index_kit_name=request.POST['collectionindexkitname']
        adapter_1=request.POST['adapter1']
        adapter_2=request.POST['adapter2']
        index_name=request.POST['indexname']
        index_sequence=request.POST['indexbase']

        # check that some values are in the request if not return the form
        if collection_index_kit_name == '' and adapter_1 =='' and adapter_2 == '' and index_name == '' and index_sequence == '' :
            return render(request, 'iSkyLIMS_wetlab/searchCollectionIndexLibrary.html')

        if index_sequence !='' :
            if len(index_sequence) < 6 :
                return render (request,'iSkyLIMS_wetlab/searchCollectionIndexLibrary.html', {'error_message':ERROR_TOO_SHORT_INDEX_BASE_SEQUENCE})
            else:
                valid_seq_characters = ['a','A','c','C', 'g', 'G', 't','T']
                for letter in index_sequence:
                    if letter not in valid_seq_characters:
                        return render (request,'iSkyLIMS_wetlab/searchCollectionIndexLibrary.html', {'error_message':ERROR_INVALID_SEQUENCE_CHARACTERS})


        collection_indexes = CollectionIndexKit.objects.all()
        if collection_index_kit_name != '':
            if collection_indexes.filter(collectionIndexName__icontains = collection_index_kit_name).exists():
                collection_indexes = collection_indexes.filter(collectionIndexName__icontains = collection_index_kit_name)
                if len (collection_indexes) == 1:
                    return redirect ('display_collection_index', collection_index_id = collection_indexes[0].get_id())
            else:
                error_message = ERROR_NO_COLLECTION_INDEX_FOUND
                error_message.append(collection_index_kit_name)
                return render (request,'iSkyLIMS_wetlab/searchCollectionIndexLibrary.html', {'error_message':error_message})
        if adapter_1 != '':
            if collection_indexes.filter(adapter1__icontains =adapter_1).exists():
                collection_indexes = collection_indexes.filter(adapter1__icontains =adapter_1)
            else:
                return render (request,'iSkyLIMS_wetlab/searchCollectionIndexLibrary.html', {'not_found_matchs':'not_found_matchs'})
        if adapter_2 != '':
            if collection_indexes.filter(adapter1__icontains =adapter_2).exists():
                collection_indexes = collection_indexes.filter(adapter2__icontains =adapter_2)
            else:
                return render (request,'iSkyLIMS_wetlab/searchCollectionIndexLibrary.html', {'not_found_matchs':'not_found_matchs'})

        if index_name != '' or index_sequence != '' :
            collection_values = CollectionIndexValues.objects.all()
            if index_name != '':
                if collection_values.filter(index_7_contains =index_name, collectionIndexKit_id__in = collection_indexes).exists():
                    collection_values = collection_values.filter(index_7_contains =index_name, collectionIndexKit_id__in = collection_indexes)
                elif collection_values.filter(index_5_contains =index_name, collectionIndexKit_id__in = collection_indexes).exists():
                    collection_values = collection_values.filter(index_5_contains =index_name, collectionIndexKit_id__in = collection_indexes)
                else:
                    return render (request,'iSkyLIMS_wetlab/searchCollectionIndexLibrary.html', {'not_found_matchs':'not_found_matchs'})

            if index_sequence != '':
                index_found , sequence = find_index_sequence_collection_values_kit(index_sequence)
                if 'I7' in index_found:
                    collection_values = collection_values.filter(i_7_seq__icontains =sequence, collectionIndexKit_id__in = collection_indexes)
                elif 'I5' in index_found :
                    collection_values = collection_values.filter(i_5_seq__icontains =sequence, collectionIndexKit_id__in = collection_indexes)
                else:
                    return render (request,'iSkyLIMS_wetlab/searchCollectionIndexLibrary.html', {'not_found_matchs':'not_found_matchs'})

            if len(collection_values) == 1:
                return redirect ('display_collection_index', collection_index_id = collection_values[0].get_collection_index_id)
            else:
                matched_collection_index = []
                collection_index_id_list = []

                for collection_value in collection_values :
                    if collection_value.get_collection_index_id() not in collection_index_id_list:
                        collection_index_id_list.append(collection_value.get_collection_index_id())
                        matched_collection_index.append([collection_value.get_collection_index_id(), collection_value.get_collection_index_name()])
                if len (matched_collection_index) == 1:
                    return redirect ('display_collection_index', collection_index_id = matched_collection_index[0][0])
                else:
                    return render (request,'iSkyLIMS_wetlab/searchCollectionIndexLibrary.html', {'matched_collection_index': matched_collection_index})
        else:
            if len(collection_indexes) == 1:
                return redirect ('display_collection_index', collection_index_id = collection_indexes[0].get_id())
            else:
                matched_collection_index = []
                for collection_index in collection_indexes:
                    matched_collection_index.append([collection_index.get_id(), collection_index.get_collection_index_name()])
                return render (request,'iSkyLIMS_wetlab/searchCollectionIndexLibrary.html', {'matched_collection_index': matched_collection_index})

    else:
        return render (request, 'iSkyLIMS_wetlab/searchCollectionIndexLibrary.html')


@login_required
def change_run_name (request, run_id):
    if RunProcess.objects.filter(pk=run_id).exists():
        run = RunProcess.objects.get(pk = run_id)
        if not request.user.is_authenticated :
            return redirect ('/accounts/login')
        # check if user is allow to make the change
        groups = Group.objects.get(name='WetlabManager')
        # check if user belongs to WetlabManager . If true allow to see the page
        if groups not in request.user.groups.all():
            return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
        if request.method == 'POST' and request.POST['action'] == 'change_run_name':
            new_run_name = request.POST['runName']
            if new_run_name == '':
                return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['Empty value is not allowed for the Run Name ']})
            if RunProcess.objects.filter(runName__exact = new_run_name).exists():
                return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['The given Run Name is already in use', 'Go back to the previous page and change the run name']})
            changed_run_name ={}
            old_run_name = run.runName
            run.runName = new_run_name
            run.save()
            changed_run_name ['new_run_name'] = [[new_run_name, run_id]]
            changed_run_name ['old_run_name'] = old_run_name

            return render (request, 'iSkyLIMS_wetlab/ChangeRunName.html', {'changed_run_name': changed_run_name})
        else:
            form_change_run_name ={}
            form_change_run_name['run_name'] = run.runName
            return render (request, 'iSkyLIMS_wetlab/ChangeRunName.html', {'form_change_run_name':form_change_run_name})
    else:
        return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['There is no Run for your query  ' ]})

@login_required
def change_project_libKit (request, project_id) :
    # check if project exists
    if Projects.objects.filter(pk = project_id).exists():
        project = Projects.objects.get(pk = project_id)
        if not request.user.is_authenticated:
            #redirect to login webpage
            return redirect ('/accounts/login')
        # check that user is allow to make the change
        allowed_access = check_user_access (request, project)
        if not allowed_access :
            return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})


        if request.method == 'POST' and request.POST['action'] == 'change_project_libKit':
            new_library_name = request.POST['projectlibkit']
            old_library_name = project.get_index_library_name()
            if old_library_name == new_library_name :
                return render (request, 'iSkyLIMS_wetlab/info_page.html', {'content': ['The library kit from the input text is the same to the existing defined for this project', 'No change is done']})
            # check if there is no other project in the same Run with the same Library Kit
            # if the library is shared with other project then error message is displayed

            if not check_user_group (request, 'WetlabManager') :
                project_run_id = project.runprocess_id.id
                project_lib_kit = project.libraryKit
                if Projects.objects.filter(runprocess_id = project_run_id).exclude(pk = project_id).exists():
                    all_project_with_same_run_id = Projects.objects.filter(runprocess_id = project_run_id).exclude(pk = project_id)
                    # there are more than 1 project on the same Run.
                    # check if these projects have in common the same library kit
                    other_lib_kits = []
                    for other_project in all_project_with_same_run_id:
                        other_lib_kits.append(other_project.libraryKit)
                    if  project_lib_kit in other_lib_kits:
                        message = str('The library Kit ' + old_library_name + 'is shared with other projects in the same Run ')

                        return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':[message, '', 'Contact with your administrator .']})
            old_lib_kit_file = project.baseSpaceFile
            new_file_name = new_library_name.replace(' ' , '_')
            #
            new_file = update_library_kit_field(old_lib_kit_file,new_file_name,new_library_name)
            if new_file == 'ERROR':
                return render (request, 'iSkyLIMS_wetlab/error_page.html', {'content':['']})
            # update the database with new file
            project.baseSpaceFile = new_file
            project.libraryKit = new_library_name
            project.save()
            # Preparing the data to show in the web page
            new_file = project.baseSpaceFile
            change_library_kit_dict ={}
            change_library_kit_dict['project']= project.projectName
            change_library_kit_dict['library_name'] = new_library_name
            change_library_kit_dict['file_to_download'] = new_file

            #
            return render (request, 'iSkyLIMS_wetlab/ChangeProjectLibraryKit.html',{'changed_lib_kit':change_library_kit_dict})
        else:
            form_change_lib_kit ={}
            project_data =[]
            project_name = project.projectName
            form_change_lib_kit['project_name'] = project_name

            project_info_text = ['Run Name', 'Project Name', 'Project date', 'User Name', 'Library Kit']
            project_values = project.get_p_info_change_library().split(';')

            for item in range (len(project_info_text)):
                project_data.append([project_info_text[item], project_values[item]])
                form_change_lib_kit['project_data'] = project_data
            return render (request, 'iSkyLIMS_wetlab/ChangeProjectLibraryKit.html',{'form_change_lib_kit': form_change_lib_kit})
    else:
        return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['No project has been found for changing the library Kit ' ]})

'''
@login_required
def change_run_libKit (request, run_id):
    #check if run exist
    if RunProcess.objects.filter(pk = run_id).exists():
        run_obj = RunProcess.objects.get(pk = run_id)
        if not request.user.is_authenticated :
            return redirect ('/accounts/login')
        # check if user is allow to make the change
        groups = Group.objects.get(name='WetlabManager')
        # check if user belongs to WetlabManager . If true allow to see the page
        if groups not in request.user.groups.all():
            return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
        if request.method == 'POST' and request.POST['action'] == 'change_run_libKit':
            new_library_kit = request.POST.getlist('runlibraryKit')
            projects_name = request.POST.getlist('projectInRun')
            changed_lib_kit_dict = {}

            # check if there is only one library kit associated to the run
            if len(new_library_kit) == 1 :
                # Check if new library kit was set in the form
                project = Projects.objects.get(projectName__exact = projects_name[0])
                old_library_kit = project.get_index_library_name()
                if new_library_kit[0] == old_library_kit :
                    return render (request, 'iSkyLIMS_wetlab/info_page.html', {'content': ['The library kit from the input text is the same to the existing defined for this project', 'No change is done']})
                # change the library name
                old_lib_kit_file = project.baseSpaceFile
                new_file_name = new_library_kit[0].replace(' ' , '_')
                #

                new_file = update_library_kit_field(old_lib_kit_file,new_file_name,new_library_kit[0])
                if new_file == 'ERROR':
                    return render (request, 'iSkyLIMS_wetlab/error_page.html', {'content':['']})
                # update the database with new file
                project.baseSpaceFile = new_file
                project.libraryKit = new_library_kit[0]
                project.save()

            else:
                old_files_to_be_deleted = []
                # check if any of the library has change
                need_to_be_updated = False
                for item in range(len(projects_name)):
                    project = Projects.objects.get(projectName__exact = projects_name[item])
                    old_library_kit = project.libraryKit
                    # get the library kit file name to delete it later
                    old_file = project.baseSpaceFile
                    # build the list to delete later the old library kit files
                    if old_file not in old_files_to_be_deleted :
                        old_files_to_be_deleted.append(old_file)
                    if new_library_kit[item] != old_library_kit:
                        need_to_be_updated = True
                if not need_to_be_updated:
                    return render (request, 'iSkyLIMS_wetlab/info_page.html', {'content': ['The library kits from the input text are the same to the existing ones defined for the Run', 'No change is done']})
                # get the Sample Sheet file related to the run
                sample_file = run.get_sample_file()

                lib_kit_dict = {}
                in_file=str('documents/' + sample_file)
                #
                ## build the project list for each library kit
                for x in range(len(new_library_kit)):
                    if new_library_kit[x] in lib_kit_dict :
                        lib_kit_dict[new_library_kit[x]].append(projects_name[x])
                    else:
                        lib_kit_dict[new_library_kit[x]]= [projects_name[x]]

                ## convert the sample sheet to base space format and have different files according the library kit
                #

                for key, value in lib_kit_dict.items():
                    lib_kit_file =key.replace(' ', '_')
                    library_file = sample_sheet_map_basespace(in_file, key, lib_kit_file, value,'Plate96')
                    if library_file == 'ERROR':
                        return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':[ 'The information on  the Library kit ', key,' For the project ', value,
                          'could not be changed to the new value of the library kit ','ADVICE', 'Contact your administrator']})
                    for updated_project_name in value :
                        p_updated = Projects.objects.get(projectName__exact = updated_project_name)
                        p_updated.libraryKit = key
                        p_updated.baseSpaceFile = library_file
                        p_updated.save()
                # delete old library kit files
                for delete_file in old_files_to_be_deleted :
                    os.remove(delete_file)
                # prepare the information to be displayed
            changed_lib_kit_dict = {}
            run_data = {}
            for item in range(len(projects_name)) :
                project = Projects.objects.get(projectName__exact = projects_name[item])
                p_name =project.projectName
                lib_name = project.libraryKit
                lib_file = project.baseSpaceFile
                run_data[p_name] = ([[lib_name, lib_file]])

            changed_lib_kit_dict['run_name'] = run.get_run_name()
            changed_lib_kit_dict['run_data'] = run_data

            return render (request, 'iSkyLIMS_wetlab/ChangeRunLibraryKit.html',{'changed_lib_kit': changed_lib_kit_dict})

        else:
            form_change_lib_kit = {}
            run_library_data = []
            project_list = []
            #library_kit_list = []
            library_kit_dict = {}
            form_change_lib_kit['run_name'] = run_obj.get_run_name()
            # get the library Kits used in run
            import pdb; pdb.set_trace()
            run_obj.get_index_library()
                pass
            project_list = Projects.objects.filter(runprocess_id = run_id)
            for project in project_list :
                lib_kit = project.libraryKit
                run_library_data.append([project.projectName, lib_kit])
                if lib_kit in library_kit_dict :
                    library_kit_dict[project.libraryKit].append(project.projectName)
                else :
                    library_kit_dict[project.libraryKit] = [project.projectName]

            form_change_lib_kit['run_data'] = library_kit_dict
            form_change_lib_kit['run_library_data'] = run_library_data

            return render (request, 'iSkyLIMS_wetlab/ChangeRunLibraryKit.html',{'form_change_lib_kit': form_change_lib_kit})
    else:
        return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['No run has been found for changing the library Kit ' ]})
'''

@login_required
def stats_experiment (request):
    return render (request, 'iSkyLIMS_wetlab/StatsPerExperiment.html', {})

@login_required
def stats_per_sequencer (request):
    sequencer_names = get_sequencer_installed_names()
    if request.method == 'POST':
        sequencer = request.POST['sequencer']
        start_date=request.POST['startdate']
        end_date=request.POST['enddate']

        if start_date != '':
            if not check_valid_date_format(start_date):
                error_message = ERROR_INVALID_FORMAT_FOR_DATES
                return render (request, 'iSkyLIMS_wetlab/StatsPerSequencer.html', {'sequencer_names':sequencer_names, 'error_message': error_message})
        if end_date !='' :
            if not check_valid_date_format(end_date):
                error_message = ERROR_INVALID_FORMAT_FOR_DATES
                return render (request, 'iSkyLIMS_wetlab/StatsPerSequencer.html', {'sequencer_names':sequencer_names, 'error_message': error_message})
        runs_using_sequencer = get_sequencers_run_from_time_interval(sequencer, start_date, end_date )
        if len(runs_using_sequencer) == 0:
            error_message = ERROR_NO_MATCHES_FOR_SEQUENCER_STATS
            return render (request, 'iSkyLIMS_wetlab/StatsPerSequencer.html', {'sequencer_names':sequencer_names, 'error_message': error_message})
        sequencer_data = get_stats_sequencer_data_from_selected_runs (runs_using_sequencer,sequencer, start_date, end_date )



        '''


                    # Get data from researcher projects




                    sequencer_statistics = {}

                    sequencer_statistics ['run_heading'] = ['Run name' , 'Run Date' , 'Sequencer Model' , 'State' , 'Space Img(Mb)' , 'Space Fasta(Mb)' , 'Space Other(Mb)']
                    runs_data = []
                    sequencer_run = {}
                    runs_data_list = []
                    runs_time_list = []
                    for run_machine in run_by_machine:
                        ru_seq_model = run_machine.get_run_used_sequencer()
                        ru_data=run_machine.get_info_process().split(';')
                        #runName, state, requested_center, generated_date, rundate,completed_date, bcl2fastq_date, finish_date, useSpaceImgMb, useSpaceFastaMb,useSpaceOtherMb
                        runs_data_list.append([ru_data[0], ru_data[4] , ru_seq_model , ru_data[3] , ru_data[8] , ru_data[9] , ru_data[10]])
                        # time the run took from run_date to completed_date
                        runs_time_list.append([ru_data[0], ru_data[5]-ru_data[4]])

                    sequencer_run[ru_seq_model] = runs_data_list
                    runs_data.append(sequencer_run)
                    sequencer_statistics['run_data'] = runs_data
                    # Get data from machine projects

                    sequencer_statistics ['projects_heading'] = ['Project name', 'Date', 'Libraty Kit','Samples', 'Cluster PF', 'Yield Mb', '% Q> 30', 'Mean','Sequencer ID']
                    projects_data =[]
                    p_machine_date , p_machine_num_sample = {} , {}
                    p_machine_lib_kit, p_machine_sequencer = {}, {}
                    p_machine_q30_dict, p_machine_mean_dict = {} , {}
                    p_machine_yield_mb_dict, p_machine_cluster_pf_dict ={} , {}
                    projects_name_dict , projects_id_list = {} , {}



                    for project_machine in projects_by_machine:

                        q_30_list , mean_q_list = [] , []
                        yield_mb_list,  cluster_pf_list = [], []
                        p_name = project_machine.get_project_name()

                        sequencer_in_project = project_machine.runprocess_id.get_run_used_sequencer()
                        if not sequencer_in_project in projects_name_dict :
                            p_machine_num_sample[sequencer_in_project] ={}
                            p_machine_sequencer[sequencer_in_project] ={}
                            p_machine_date[sequencer_in_project] ={}
                            p_machine_lib_kit[sequencer_in_project] ={}
                            p_machine_q30_dict[sequencer_in_project] ={}
                            p_machine_mean_dict[sequencer_in_project] ={}
                            p_machine_yield_mb_dict[sequencer_in_project] ={}
                            p_machine_cluster_pf_dict[sequencer_in_project] ={}
                            projects_name_dict[sequencer_in_project] = []
                            projects_id_list[sequencer_in_project] =  []
                        projects_name_dict[sequencer_in_project].append(p_name)
                        m_project_id = project_machine.id
                        projects_id_list[sequencer_in_project].append(m_project_id)
                        p_machine_num_sample[sequencer_in_project][p_name] = StatsFlSummary.objects.get(project_id__exact = m_project_id).sampleNumber
                        p_machine_date [sequencer_in_project][p_name] = project_machine.get_date()
                        p_machine_lib_kit[sequencer_in_project][p_name]= project_machine.get_index_library_name()
                        p_machine_sequencer[sequencer_in_project][p_name] = str(project_machine.runprocess_id.sequencerModel)
                        lanes_in_project = StatsLaneSummary.objects.filter( project_id__exact = m_project_id)
                        for lane in lanes_in_project :
                            q_30_value, mean_q_value , yield_mb_value , cluster_pf_value = lane.get_stats_info()
                            q_30_list.append(float(q_30_value))
                            mean_q_list.append(float(mean_q_value))
                            yield_mb_list.append(float(yield_mb_value.replace(',','')))
                            cluster_pf_list.append(float(cluster_pf_value.replace(',','')))
                        p_machine_q30_dict[sequencer_in_project] [p_name]= format(statistics.mean(q_30_list), '.2f')
                        p_machine_mean_dict[sequencer_in_project][p_name] = format(statistics.mean(mean_q_list), '.2f')
                        p_machine_yield_mb_dict[sequencer_in_project][p_name] = round(sum(yield_mb_list))
                        p_machine_cluster_pf_dict[sequencer_in_project][p_name] = round(sum(cluster_pf_list))

                    # Create the table with projects executed by the machine
                    machine_seq_graphs, machine_graphs = [], []
                    for sequencer, projects_name_list in projects_name_dict.items() :
                        sequencer_proj = {}
                        proj_data =[]
                        for project_name in projects_name_list :
                            proj_data.append([project_name, p_machine_date[sequencer][project_name], p_machine_lib_kit[sequencer][project_name],
                                    p_machine_num_sample[sequencer][project_name], '{0:,}'.format(int(p_machine_cluster_pf_dict[sequencer][project_name])),
                                    '{0:,}'.format(int(p_machine_yield_mb_dict[sequencer][project_name])), p_machine_q30_dict[sequencer][project_name],
                                    p_machine_mean_dict[sequencer][project_name], p_machine_sequencer[sequencer][project_name]])
                        sequencer_proj[sequencer] = proj_data
                        projects_data.append(sequencer_proj)
                        # create the graphic for q30 quality
                        theme = 'ocean'
                        heading = 'Graphics for Q > 30 for machine ' + m_name
                        sub_caption = 'Sequencer ' + sequencer
                        x_axis_name = 'Projects'
                        y_axis_name = 'Q 30 (in %)'

                        data_source = column_graphic_simple (heading, sub_caption, x_axis_name, y_axis_name, theme, p_machine_q30_dict[sequencer])
                        seq_chart = sequencer + 'q30_chart'
                        seq_graph = sequencer + 'q30_graph'
                        q30_machine_seq_graph = FusionCharts("column3d", seq_graph , "500", "350",seq_chart , "json", data_source).render()

                        machine_seq_graphs.append([seq_chart, q30_machine_seq_graph])
                        # create the graphic for mean quality
                        theme = 'carbon'
                        heading = 'Graphics for Mean quality for machine ' + m_name
                        sub_caption = 'Sequencer ' + sequencer
                        x_axis_name = 'Projects'
                        y_axis_name = 'Mean Quality'
                        data_source = column_graphic_simple (heading, sub_caption, x_axis_name, y_axis_name, theme, p_machine_mean_dict[sequencer])
                        seq_chart = sequencer + 'mean_q_chart'
                        seq_graph = sequencer + 'mean_q_graph'
                        mean_q_machine_seq_graph = FusionCharts("column3d", seq_graph , "500", "350", seq_chart, "json", data_source).render()
                        machine_seq_graphs.append([seq_chart, mean_q_machine_seq_graph])

                        # create the graphic for yield Mb
                        theme = 'zune'
                        heading = 'Graphics for Yield Mb for machine ' + m_name
                        sub_caption = 'Sequencer ' + sequencer
                        x_axis_name = 'Projects'
                        y_axis_name = 'Yield Mb'
                        data_source = column_graphic_simple (heading, sub_caption, x_axis_name, y_axis_name, theme, p_machine_yield_mb_dict[sequencer])
                        seq_chart = sequencer + 'yield_mb_chart'
                        seq_graph = sequencer + 'yield_mb_graph'
                        yield_mb_machine_graph = FusionCharts("column3d", seq_graph , "500", "350", seq_chart, "json", data_source).render()
                        machine_seq_graphs.append([seq_chart, yield_mb_machine_graph])
                        # create the graphic for cluster Pf
                        theme = 'ocean'
                        heading = 'Graphics for Cluster Pf for machine ' + m_name
                        sub_caption = 'Sequencer ' + sequencer
                        x_axis_name = 'Projects'
                        y_axis_name = 'Cluster Pf'
                        data_source = column_graphic_simple (heading, sub_caption, x_axis_name, y_axis_name, theme, p_machine_cluster_pf_dict[sequencer])
                        seq_chart = sequencer + 'cluster_pf_chart'
                        seq_graph = sequencer + 'cluster_pf_graph'
                        cluster_pf_machine_graph = FusionCharts("column3d", seq_graph , "500", "350", seq_chart, "json", data_source).render()
                        machine_seq_graphs.append([seq_chart, cluster_pf_machine_graph])


                        machine_graphs.append(machine_seq_graphs)

                    sequencer_statistics ['machine_graph'] = machine_graphs
                    sequencer_statistics ['machine_name'] = m_name
                    sequencer_statistics['projects_data'] = projects_data


                    #collecting data for comparation graphics

                    # Calculating the mean for all projects performed by researcher
                    comp_q30_dict, comp_mean_q_dict = {} , {}
                    comp_yield_mb_dict, comp_cluster_pf_dict = {} , {}

                    q30_val = p_machine_q30_dict.values()
                    mean_q_val = p_machine_mean_dict.values()
                    yield_mb_val = p_machine_yield_mb_dict.values()
                    cluster_pf = p_machine_cluster_pf_dict.values()

                    for sequencer in projects_name_dict.keys() :
                        #sequencer_proj = {}
                        #proj_data =[]
                        comp_q30_dict[sequencer] , comp_mean_q_dict [sequencer]= {} , {}
                        comp_yield_mb_dict[sequencer], comp_cluster_pf_dict [sequencer] = {}, {}
                        comp_q30_dict [sequencer][m_name] = format(statistics.mean( [float(x) for x in list(p_machine_q30_dict[sequencer].values())]),'.2f')
                        comp_mean_q_dict[sequencer] [m_name] = format(statistics.mean( [float(x) for x in list(p_machine_mean_dict[sequencer].values())]),'.2f')

                        comp_yield_mb_dict[sequencer] [m_name] = sum(list(p_machine_yield_mb_dict[sequencer].values()))
                        comp_cluster_pf_dict[sequencer] [m_name] = sum(list(p_machine_cluster_pf_dict[sequencer].values()))

                    total_q_30_list, total_mean_q_list = [] , []
                    total_yield_mb_list, total_cluster_pf_list = [] , []
                    total_lanes_summary = {}

                    for sequencer in projects_name_dict.keys() :
                        runs_sequencer = RunProcess.objects.filter(sequencerModel__machineName__exact = sequencer)
                        run_sequencer_id_list = []
                        for run in runs_sequencer :
                            run_sequencer_id_list.append(run.pk)

                        if StatsLaneSummary.objects.filter(runprocess_id__in  = run_sequencer_id_list).exclude(defaultAll__isnull = False).exclude(project_id__in = projects_id_list[sequencer]).exists():
                            total_lanes_summary[sequencer] = StatsLaneSummary.objects.filter(runprocess_id__in  = run_sequencer_id_list).exclude(defaultAll__isnull = False).exclude(project_id__in = projects_id_list[sequencer])
                        else:
                            total_lanes_summary[sequencer] = ''


                    comp_graphs, comp_seq_graphs = [] , []
                    for sequencer in projects_name_dict.keys() :
                        for lane_summary in total_lanes_summary[sequencer] :
                            q_30_value, mean_q_value , yield_mb_value , cluster_pf_value = lane_summary.get_stats_info()
                            total_q_30_list.append(float(q_30_value))
                            total_mean_q_list.append(float(mean_q_value))
                            total_yield_mb_list.append(int(yield_mb_value.replace(',','')))
                            total_cluster_pf_list.append(int(cluster_pf_value.replace(',','')))


                    sequencer_statistics ['comp_graphs'] = comp_graphs

                    # Sequencer graphic utilization
                    sequencer_used = {}
                    for sequencer in projects_name_dict.keys() :
                        sequencer_used[sequencer] = len( projects_name_dict[sequencer])

                    theme = 'ocean'
                    heading = 'Sequencer utilization for machine ' + m_name
                    sub_caption = ''
                    data_source = pie_graphic_standard (heading, sub_caption, theme, sequencer_used)
                    sequencer_pie_graph = FusionCharts("pie3d", "sequencer_pie_graph" , "500", "400", "sequencer_pie_chart", "json", data_source).render()
                    sequencer_statistics ['sequencer_pie_graph'] = sequencer_pie_graph

                    return  render(request, 'iSkyLIMS_wetlab/StatsPerMachine.html', {'sequencer_statistics' : sequencer_statistics})
        '''
        return render (request, 'iSkyLIMS_wetlab/StatsPerSequencer.html', {'sequencer_data':sequencer_data})

    else:
        return render (request, 'iSkyLIMS_wetlab/StatsPerSequencer.html', {'sequencer_names':sequencer_names})


@login_required
def stats_per_researcher (request):
    if request.method == 'POST':
        r_name = request.POST['researchername']
        start_date=request.POST['startdate']
        end_date=request.POST['enddate']

        researcher_statistics = get_researcher_statistics(r_name, start_date, end_date)
        if 'ERROR' in researcher_statistics:
            error_message = researcher_statistics
            return render (request,'iSkyLIMS_wetlab/StatsPerResearcher.html', {'researcher_statistics':error_message})

        return  render(request, 'iSkyLIMS_wetlab/StatsPerResearcher.html', {'researcher_statistics' : researcher_statistics})

    else:
        return render (request, 'iSkyLIMS_wetlab/StatsPerResearcher.html', {})


@login_required
def stats_per_time (request):
    if request.method=='POST':
        start_date=request.POST['startdate']
        end_date=request.POST['enddate']
         ### check the right format of start and end date
        if start_date != '' and not check_valid_date_format(start_date):
            errror_message = wetlab_config.ERROR_INVALID_FORMAT_FOR_DATES
            return render (request,'iSkyLIMS_wetlab/StatsPerTime.html', {'ERROR':error_message})
        if end_date != '' and not check_valid_date_format(start_date):
            errror_message = wetlab_config.ERROR_INVALID_FORMAT_FOR_DATES
            return render (request,'iSkyLIMS_wetlab/StatsPerTime.html', {'ERROR':error_message})
        #############################################################
        #### searching for runs were match the state and start and end date
        #############################################################
        if (start_date != '' and end_date != ''):
            stat_per_time ={}
            if (RunProcess.objects.filter( state__runStateName='Completed', run_date__range=(start_date, end_date)).exists()):
                run_stats_list=RunProcess.objects.filter(state__runStateName='Completed', run_date__range=(start_date, end_date)).order_by('run_date')

                run_list={}
                run_date_name ={}
                ## get the run names that matches de conditions
                for run in run_stats_list:
                    #run_list.append([run.get_run_name(),run.id])
                    run_date = str(run.run_date)
                    if run_date in run_date_name:
                        run_date_name [run_date] +=1
                    else:
                        run_date_name[run_date] = 1
                    run_list [run.id] = [[run.get_run_name(), run_date]]
                    #
                stat_per_time ['run_names'] = run_list
                if len (run_stats_list) == 1:
                    number_of_runs = '1 Run'
                else:
                    number_of_runs = str(len (run_stats_list)) + '  Runs'
                stat_per_time ['number_of_runs'] = number_of_runs
                stat_per_time ['dates'] = start_date + ' and  ' + end_date
                #
                ############################################################
                ### define the graphics for found run in the period
                heading = 'Runs found during the period ' + str(start_date) + ' and ' + str(end_date)
                sub_caption = ''
                x_axis_name = 'Date'
                y_axis_name = 'Number of runs'
                run_period_chart_number = 'run_period_chart-1'
                run_period_index_graph = 'exq1'

                data_source = researcher_project_column_graphic (heading, sub_caption, x_axis_name, y_axis_name, 'ocean', run_date_name)
                stat_per_time['run_period_graphic'] = FusionCharts("column3d", run_period_index_graph , "550", "350", run_period_chart_number, "json", data_source).render()

                ####### end creation run preparation graphics

                #############################################################
                ### collect statistics for Projects
                if RunProcess.objects.filter(state__runStateName__exact = 'Completed', run_completed_date__range=(start_date, end_date)).exists():
                    run_objs =  RunProcess.objects.filter(state__runStateName__exact = 'Completed', run_completed_date__range=(start_date, end_date))
                    if Projects.objects.filter(runProcess__in = run_objs).exists():
                        project_found_list = Projects.objects.filter(runProcess__in = run_objs)

                        project_list={}
                        project_date_name ={}
                        ## get the project names that matches de conditions
                        for project in project_found_list:
                            project_run_date = str(project.project_run_date)
                            if project_run_date in project_date_name:
                                project_date_name [project_run_date] +=1
                            else:
                                project_date_name[project_run_date] = 1
                            project_list [project.id] = [[project.get_project_name(), project_run_date]]
                            #
                        stat_per_time ['project_names'] = project_list
                        if len (project_found_list) == 1:
                            number_of_projects = '1 Project'
                        else:
                            number_of_projects = str(len (project_found_list)) + '  Projects'
                        stat_per_time ['number_of_projects'] = number_of_projects
                        stat_per_time ['dates'] = start_date + ' and  ' + end_date
                        #
                        ############################################################
                        ### define the graphics for found run in the period
                        heading = 'Projects found during the period ' + str(start_date) + ' and ' + str(end_date)
                        sub_caption = ''
                        x_axis_name = 'Date'
                        y_axis_name = 'Number of Projects'
                        run_period_chart_number = 'project_period_chart-1'
                        run_period_index_graph = 'project_period-1'

                        data_source = researcher_project_column_graphic (heading, sub_caption, x_axis_name, y_axis_name, 'carbon', project_date_name)
                        stat_per_time['project_period_graphic'] = FusionCharts("column3d", run_period_index_graph , "550", "350", run_period_chart_number, "json", data_source).render()

                ####### end creation run preparation graphics
                #############################################################
                ### collect statistics for unkow Barcodes
                #top_unbarcode_list = []
                count_unbarcode  = {}

                for run in run_stats_list:
                    run_obj = run.get_run_id()
                    run_param_obj = RunningParameters.objects.get(runName_id = run_obj)
                    lanes_in_sequencer = int(run_param_obj.get_number_of_lanes())
                    top_unbarcode_all_runs  = {}
                    for lane_number in range (1, lanes_in_sequencer +1):
                        lane_unbarcodes = RawTopUnknowBarcodes.objects.filter(runprocess_id =run, lane_number__exact = lane_number)
                        for lane_unbarcode in lane_unbarcodes :
                            if not lane_number in count_unbarcode :
                                count_unbarcode[lane_number] = {}
                            unbarcode_num , unknown_barcode,  = lane_unbarcode.get_unknow_barcodes().split(';')
                            value_unbarcode = int(unbarcode_num.replace(',',''))
                            if not unknown_barcode in count_unbarcode[lane_number] :
                                count_unbarcode[lane_number][unknown_barcode] = value_unbarcode
                            else:
                                count_unbarcode[lane_number][unknown_barcode] += value_unbarcode
                            if not unknown_barcode in top_unbarcode_all_runs :
                                top_unbarcode_all_runs[unknown_barcode] = value_unbarcode
                            else:
                                top_unbarcode_all_runs[unknown_barcode] += value_unbarcode

                themes = ['','ocean','fint','carbon','zune', '']
                # prepare the column graphic for nunber of top Unknow Barcode
                unbar_lane_chart = []
                for lane_number in range (1, lanes_in_sequencer +1):
                    heading = 'Number of undetermined barcode sequence in lane ' + str(lane_number)
                    chart_number = 'chart-' + str(lane_number)
                    render_number = 'ex'+ str(lane_number)
                    lane_chart = 'lane_chart'+ str(lane_number)
                    data_source = graphic_for_unbarcodes(heading , themes[lane_number] , count_unbarcode[lane_number])
                    lane_graphic = FusionCharts("column3d", render_number , "500", "400", chart_number, "json", data_source)
                    unbar_lane_chart.append([chart_number, str(lane_number),lane_graphic.render()])
                stat_per_time ['unbar_lane_chart'] = unbar_lane_chart


                # prepare the pie graphic for the number of top Unknow Barcode per sequence
                data_source = pie_graphic ('Number of count for the Undetermined Sequences', 'fint',top_unbarcode_all_runs)
                unknow_pie3d = FusionCharts("pie3d", "ex5" , "500", "400", "chart-5", "json", data_source)
                stat_per_time ['unknow_pie3d'] = unknow_pie3d.render()

                #########################
                ### Insert information for disk space utilization
                run_disk_utilization ={}
                for run_disk_stats in run_stats_list :
                    run_name_disk = run_disk_stats.runName
                    run_disk_utilization[run_name_disk] = run_disk_stats.get_disk_space_utilization()


                heading = 'Disk space used for each Run found during the period ' + str(start_date) + ' and ' + str(end_date)
                sub_caption = ''
                x_axis_name = 'Date'
                y_axis_name = 'Disk space used (MB)'
                disk_space_period_chart_number = 'disk_usage_chart-1'
                disk_space_period_index_graph = 'diskusage1'

                data_source = researcher_project_column_graphic (heading, sub_caption, x_axis_name, y_axis_name, 'carbon', run_disk_utilization)
                stat_per_time['disk_space_period_graphic'] = FusionCharts("column3d", disk_space_period_index_graph , "950", "350", disk_space_period_chart_number, "json", data_source).render()

                #
                return render(request, 'iSkyLIMS_wetlab/StatsPerTime.html', {'display_stats_per_time': stat_per_time })

            else:
                return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['No matches have been found for Runs created between', start_date, ' and the ',  end_date ]})
        else:
            return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':'Start date and End Date cannot be empty '})


    return render (request,'iSkyLIMS_wetlab/StatsPerTime.html')

def get_list_of_libraries_values (library_found, q30_comparations, mean_comparations , n_bases_comparations) :

    for project_to_compare in library_found :
        library_to_compare_name = project_to_compare.get_index_library_name()
        #project_to_compare_id = project_to_compare.id
        q30_compare_lib, mean_compare_lib, yield_mb_compare_lib = [], [] , []

        ### This line must changed to handle project name is reused in several runs
        run_used_in_project = project_to_compare.runProcess.all().last()

        run_param_obj = RunningParameters.objects.get(runName_id = run_used_in_project)
        # get the number of lanes by quering the SequencerModel in the RunProcess
        #number_of_lanes = project_to_compare.runprocess_id.get_sequencing_lanes()
        number_of_lanes = int(run_param_obj.get_number_of_lanes())
        for lane_number in range (1,number_of_lanes + 1):
            try:
                lane_in_project = StatsLaneSummary.objects.get(project_id = project_to_compare, lane__exact = lane_number)
            except:
                continue
            q_30_value, mean_q_value, yield_mb , cluster_pf = lane_in_project.get_stats_info()
            q30_compare_lib.append(float(q_30_value))
            mean_compare_lib.append(float(mean_q_value))
            yield_mb_compare_lib.append(float(yield_mb.replace(',','')))
        if library_to_compare_name in q30_comparations:
            q30_tmp_list =[float(q30_comparations [library_to_compare_name]), statistics.mean (q30_compare_lib)]
            q30_comparations [library_to_compare_name] = format(statistics.mean (q30_tmp_list), '.2f')
            mean_tmp_list = [float(mean_comparations [library_to_compare_name]), statistics.mean (mean_compare_lib)]
            mean_comparations [library_to_compare_name] = format(statistics.mean (mean_tmp_list), '.2f')
            n_bases_list =[float(n_bases_comparations [library_to_compare_name]), sum (yield_mb_compare_lib)]

            n_bases_comparations [library_to_compare_name] = format(statistics.mean (n_bases_list), '.2f')
        else:
            q30_comparations [library_to_compare_name] = format(statistics.mean (q30_compare_lib), '.2f')
            mean_comparations [library_to_compare_name] = format(statistics.mean (mean_compare_lib), '.2f')
            n_bases_comparations [library_to_compare_name] = format(statistics.mean (yield_mb_compare_lib), '.2f')


@login_required
def stats_per_library (request):
    if request.method=='POST' :
        library_kit_name=request.POST['libraryKitName']
        start_date=request.POST['startdate']
        end_date=request.POST['enddate']
        # check that some values are in the request if not return the form
        if library_kit_name == '' and start_date == '' and end_date == '' :
            return render(request, 'iSkyLIMS_wetlab/StatsPerLibrary.html')

        if library_kit_name !=''  and len(library_kit_name) < 5 :
            error_message = ERROR_TOO_SHORT_INDEX_LIBRAY_NAME
            return render(request, 'iSkyLIMS_wetlab/StatsPerLibrary.html', {'error_message':error_message})

        ### check the right format of start and end date
        if start_date != '':
            if not check_valid_date_format(start_date) :
                error_message = ERROR_INVALID_FORMAT_FOR_DATES
                return render(request, 'iSkyLIMS_wetlab/StatsPerLibrary.html', {'error_message':error_message})
        if end_date != '':
            if not check_valid_date_format(end_date) :
                error_message = ERROR_INVALID_FORMAT_FOR_DATES
                return render(request, 'iSkyLIMS_wetlab/StatsPerLibrary.html', {'error_message':error_message})

        if library_kit_name != '':
            if Projects.objects.filter(libraryKit__icontains = library_kit_name, runProcess__state__runStateName__exact = 'Completed').exists():
                library_found = Projects.objects.filter(libraryKit__icontains = library_kit_name, runProcess__state__runStateName__exact = 'Completed')
            else:
                error_message = ERROR_NO_MATCHES_FOR_LIBRARY_STATISTICS
                return render(request, 'iSkyLIMS_wetlab/StatsPerLibrary.html' , {'error_message':error_message})
        else:
            library_found = Projects.objects.filter(runProcess__state__runStateName__exact = 'Completed')
        if (start_date != '' and end_date != ''):
            if library_found.filter(project_run_date__range=(start_date, end_date)).exists():
                 library_found = library_found.filter(project_run_date__range=(start_date, end_date))
            else:
                error_message = ERROR_NO_MATCHES_FOR_LIBRARY_STATISTICS
                return render(request, 'iSkyLIMS_wetlab/StatsPerLibrary.html', {'error_message':error_message})
        if start_date !='' and end_date == '':
            if library_found.filter(project_run_date__gte = start_date).exists():
                 library_found = library_found.filter(project_run_date__gte = start_date)
                 #
            else:
                error_message = ERROR_NO_MATCHES_FOR_LIBRARY_STATISTICS
                return render(request, 'iSkyLIMS_wetlab/StatsPerLibrary.html', {'error_message':error_message})
        if start_date =='' and end_date != '':
            if library_found.filter(project_run_date__lte = end_date).exists():
                #
                library_found = library_found.filter(project_run_date__lte = end_date)
            else:
                error_message = ERROR_NO_MATCHES_FOR_LIBRARY_STATISTICS
                return render(request, 'iSkyLIMS_wetlab/StatsPerLibrary.html', {'error_message':error_message})

        #Collecting the statistics for the selected library
        # Get the projects which are using the library kit
        #
        library_stats ={}
        projects_name_in_library =[]
        q_30_list , mean_q_list , yield_mb_list = [] , [] ,[]
        # Getting 1 library. Library could be in several projects. Information is collected per lane and by project
        #check if only 1 library kit matches the query
        library_names ={}
        for library in library_found :
            library_names [library.get_index_library_name()] = 1
        #
        if len(library_names) == 1:
            # There is only 1 library in the query. Results displays all projects data which have this library kit
            mean_lane_graphic ={}
            for project in library_found :
                projects_name_in_library.append(project.get_project_name())
            q30_in_lib, mean_in_lib, yield_mb_in_lib = [], [] , []
            for lane_number in range (1,5):
                q_30_lane , mean_q_lane , yield_mb_lane = {} , {} ,{}
                for project in library_found :
                    project_id = project.get_project_id()
                    # Get quality information for each Lane summary of the project id
                    #
                    lane_in_project = StatsLaneSummary.objects.get(project_id__exact = project_id, lane__exact = lane_number)
                    q_30_value, mean_q_value, yield_mb , cluster_pf = lane_in_project.get_stats_info()
                    project_name = project.get_project_name()
                    q_30_lane[project_name] = q_30_value
                    q30_in_lib.append(float(q_30_value))
                    mean_q_lane[project_name] = mean_q_value
                    mean_in_lib.append(float(mean_q_value))
                    yield_mb_lane[project_name] = yield_mb.replace(',','')
                    yield_mb_in_lib.append(float(yield_mb.replace(',','')))
                    #
                # creating the Yield MBases graphics
                chart_number = 'chart-' + str(lane_number)
                render_number = 'ex'+ str(lane_number)
                heading = 'Number of MBases in the projects for Lane ' + str(lane_number)
                data_source = graphic_for_library_kit (heading, 'projects in lane ' ,'Project Names', 'Number of M bases', 'ocean', yield_mb_lane)
                yield_mb_lane_graphic = FusionCharts("column3d", render_number , "500", "300", chart_number, "json", data_source)
                #
                yield_graphic = 'yield_mb_graphic' + str(lane_number)
                library_stats [yield_graphic] = yield_mb_lane_graphic.render()

                # creating the Q30 graphics
                chart_number = 'q30-chart-' + str(lane_number)
                render_number = 'q30-ex'+ str(lane_number)
                heading = 'Percent of bases > Q30 in the projects for Lane ' + str(lane_number)
                data_source = graphic_for_library_kit (heading, 'projects in lane ' ,'Project Names', 'Percent of Q 30', 'zune', q_30_lane)
                q30_lane_graphic = FusionCharts("column3d", render_number , "400", "300", chart_number, "json", data_source)
                #
                q30_graphic = 'q30_graphic' + str(lane_number)
                library_stats [q30_graphic] = q30_lane_graphic.render()

                # creating the Mean graphics
                chart_number = 'mean-chart-' + str(lane_number)
                render_number = 'mean-ex'+ str(lane_number)
                heading = 'Mean Quality Score in the projects for Lane ' + str(lane_number)
                data_source = graphic_for_library_kit (heading, 'projects in lane ' ,'Project Names', 'Percent of Q 30', 'carbon', mean_q_lane)
                mean_lane_graphic = FusionCharts("column3d", render_number , "400", "300", chart_number, "json", data_source)
                #
                mean_graphic = 'mean_graphic' + str(lane_number)
                library_stats [mean_graphic] = mean_lane_graphic.render()


            library_name = project.get_index_library_name()
            library_stats['library_name'] = library_name
            library_stats['project_names'] = projects_name_in_library
            #
            ########################################################################
            # set the data for the library under study
            ########################################################################
            q30_comparations , mean_comparations , n_bases_comparations = {}, {} , {}
            q30_comparations [library_name] = format(statistics.mean (q30_in_lib), '.2f')
            mean_comparations [library_name] = format(statistics.mean (mean_in_lib), '.2f')
            n_bases_comparations [library_name] = format(statistics.mean (yield_mb_in_lib), '.2f')
            error_in_library_to_compare = ''
            # get the data for the libraries to compare with
            if start_date == '' and end_date == '':
                if Projects.objects.filter(runprocess_id__state__runStateName__exact = 'Completed').exclude(libraryKit__exact = library_name).exists():
                    libraries_to_compare = Projects.objects.filter(runprocess_id__state__runStateName__exact = 'Completed').exclude(libraryKit__exact = library_name)
                else :
                    error_in_library_to_compare ='No other library have been found for doing the comparison. '

            if start_date != '' and end_date == '':
                if Projects.objects.filter(runprocess_id__state__runStateName__exact = 'Completed', generatedat__gte = start_date).exclude(libraryKit__exact = library_name).exists():
                    libraries_to_compare = Projects.objects.filter(runprocess_id__state__runStateName__exact = 'Completed', generatedat__gte = start_date).exclude(libraryKit__exact = library_name)
                else :
                    error_in_library_to_compare ='No other library have been found for doing the comparison, with the starting date  ' + start_date

            if start_date == '' and end_date != '':
                if Projects.objects.filter(runprocess_id__state__runStateName__exact = 'Completed', generatedat__lte = end_date).exclude(libraryKit__exact = library_name).exists():
                    libraries_to_compare = Projects.objects.filter(runprocess_id__state__runStateName__exact = 'Completed', generatedat__lte = end_date).exclude(libraryKit__exact = library_name)
                else :
                    error_in_library_to_compare ='No other library have been found for doing the comparison ending with  ' + end_date

            if start_date != '' and end_date != '':
                if Projects.objects.filter(runprocess_id__state__runStateName__exact = 'Completed', generatedat__range =(start_date, end_date)).exclude(libraryKit__exact = library_name).exists():
                    libraries_to_compare = Projects.objects.filter(runprocess_id__state__runStateName__exact = 'Completed', generatedat__range =(start_date, end_date)).exclude(libraryKit__exact = library_name)
                else :
                    error_in_library_to_compare ='No other library have been found for doing the comparison for the start date  ' + start_date + '  and with the ending date  ' + end_date

            if error_in_library_to_compare == '':
                for project_to_compare in libraries_to_compare :
                    library_to_compare_name = project_to_compare.get_index_library_name()
                    project_to_compare_id = project_to_compare.get_project_id()
                    #q_30_lane , mean_q_lane , yield_mb_lane = {} , {} ,{}
                    q30_compare_lib, mean_compare_lib, yield_mb_compare_lib = [], [] , []

                    run_obj = project_to_compare.get_run_obj()
                    run_param_obj = RunningParameters.objects.get(run_id = run_obj)
                    lanes_in_sequencer = int(run_param_obj.get_number_of_lanes())
                    for lane_number in range (1,lanes_in_sequencer+ 1):

                        lane_in_project = StatsLaneSummary.objects.get(project_id__exact = project_to_compare_id, lane__exact = lane_number)
                        q_30_value, mean_q_value, yield_mb , cluster_pf = lane_in_project.get_stats_info()
                        q30_compare_lib.append(float(q_30_value))
                        mean_compare_lib.append(float(mean_q_value))
                        yield_mb_compare_lib.append(float(yield_mb.replace(',','')))
                    if library_to_compare_name in q30_comparations:
                        q30_tmp_list =[float(q30_comparations [library_to_compare_name]), statistics.mean (q30_compare_lib)]
                        q30_comparations [library_to_compare_name] = format(statistics.mean (q30_tmp_list), '.2f')
                        mean_tmp_list = [float(mean_comparations [library_to_compare_name]), statistics.mean (mean_compare_lib)]
                        mean_comparations [library_to_compare_name] = format(statistics.mean (mean_tmp_list), '.2f')
                        n_bases_list =[float(n_bases_comparations [library_to_compare_name]), statistics.mean (yield_mb_compare_lib)]
                        n_bases_comparations [library_to_compare_name] = format(statistics.mean (n_bases_list), '.2f')
                    else:
                        q30_comparations [library_to_compare_name] = format(statistics.mean (q30_compare_lib), '.2f')
                        mean_comparations [library_to_compare_name] = format(statistics.mean (mean_compare_lib), '.2f')
                        n_bases_comparations [library_to_compare_name] = format(statistics.mean (yield_mb_compare_lib), '.2f')

            else:
                library_stats ['error_library'] = error_in_library_to_compare


            heading = 'Comparison of Percent of bases > Q30  '
            data_source = graphic_for_library_kit (heading, 'Q30 comparison ' ,'Library Names', 'Percent of Q 30', '', q30_comparations)
            comp_q30_lib_graphic = FusionCharts("column3d", 'comp-q30-1' , "500", "300", 'comp-q30-chart-1', "json", data_source)
            library_stats ['comp_q30_graphic'] = comp_q30_lib_graphic.render()

            heading = 'Comparison of Mean Quality Score '
            data_source = graphic_for_library_kit (heading, 'Mean Quality Score comparison ' ,'Library Names', 'Mean Quality Score', '', mean_comparations)
            comp_mean_lib_graphic = FusionCharts("column3d", 'comp-mean-1' , "500", "300", 'comp-mean-chart-1', "json", data_source)
            library_stats ['comp_mean_graphic'] = comp_mean_lib_graphic.render()

            heading = 'Number of Bases comparison'
            data_source = graphic_for_library_kit (heading, 'Number of Bases comparison ' ,'Library Names', 'Number of Bases ', '', n_bases_comparations)
            comp_mean_lib_graphic = FusionCharts("column3d", 'comp-n_bases-1' , "500", "300", 'comp-n_bases-chart-1', "json", data_source)
            library_stats ['comp_n_bases_graphic'] = comp_mean_lib_graphic.render()

            return render (request,'iSkyLIMS_wetlab/StatsPerLibrary.html', {'display_library_stats': library_stats })
        else:
            library_list_stats ={}
            libraries_found_name =[]
            # get the library names that match with the searching criteria
            for library in library_found :
                lib_name =library.get_index_library_name ()
                if not lib_name in libraries_found_name :
                    libraries_found_name.append(lib_name)
            #
            library_list_stats['library_names'] = libraries_found_name
            q30_comparations , mean_comparations , n_bases_comparations = {}, {} , {}
            ###
            # get the data for displaying the libraries found in the form request
            ###
            get_list_of_libraries_values (library_found, q30_comparations, mean_comparations , n_bases_comparations)

            heading = 'Comparison of Percent of bases > Q30  '
            data_source = graphic_for_library_kit (heading, 'Q30 comparison ' ,'Library Names', 'Percent of Q 30', '', q30_comparations)
            comp_q30_lib_graphic = FusionCharts("column3d", 'comp-q30-1' , "500", "300", 'comp-q30-chart-1', "json", data_source)
            #
            library_list_stats ['comp_q30_graphic'] = comp_q30_lib_graphic.render()

            heading = 'Comparison of Mean Quality Score '
            data_source = graphic_for_library_kit (heading, 'Mean Quality Score comparison ' ,'Library Names', 'Mean Quality Score', '', mean_comparations)
            comp_mean_lib_graphic = FusionCharts("column3d", 'comp-mean-1' , "500", "300", 'comp-mean-chart-1', "json", data_source)
            #
            library_list_stats ['comp_mean_graphic'] = comp_mean_lib_graphic.render()

            heading = 'Number of Bases comparison'
            data_source = graphic_for_library_kit (heading, 'Number of Bases comparison ' ,'Library Names', 'Number of Bases ', '', n_bases_comparations)
            comp_mean_lib_graphic = FusionCharts("column3d", 'comp-n_bases-1' , "500", "300", 'comp-n_bases-chart-1', "json", data_source)
            #
            library_list_stats ['comp_n_bases_graphic'] = comp_mean_lib_graphic.render()
            ###
            # get the data for displaying the libraries found in the form request
            ###
            all_libraries = Projects.objects.filter(runProcess__state__runStateName__exact = 'Completed')
            if (start_date != '' and end_date != ''):
                if all_libraries.filter(generatedat__range=(start_date, end_date)).exists():
                     library_found = library_found.filter(generatedat__range=(start_date, end_date))
            if start_date !='' and end_date == '':
                if all_libraries.filter(generatedat__gte = start_date).exists():
                     all_libraries = library_found.filter(generatedat__gte = start_date)
            if start_date =='' and end_date != '':
                if all_libraries.filter(generatedat__lte = end_date).exists():
                    #
                    all_libraries = library_found.filter(generatedat__lte = end_date)

            q30_comparations , mean_comparations , n_bases_comparations = {}, {} , {}
            get_list_of_libraries_values (all_libraries, q30_comparations, mean_comparations , n_bases_comparations)
            #
            heading = 'Library kits of Percent of bases > Q30  '
            data_source = graphic_for_library_kit (heading, 'Q30 library kits ' ,'Library Names', 'Percent of Q 30', '', q30_comparations)
            lib_q30_lib_graphic = FusionCharts("column3d", 'lib-q30-lib' , "500", "300", 'lib-q30-chart-1', "json", data_source)
            #
            library_list_stats ['lib_q30_graphic'] = lib_q30_lib_graphic.render()

            heading = 'Library kits of Mean Quality Score '
            data_source = graphic_for_library_kit (heading, 'Mean Quality Score Library kits ' ,'Library Names', 'Mean Quality Score', '', mean_comparations)
            lib_mean_lib_graphic = FusionCharts("column3d", 'lib-mean-lib' , "500", "300", 'lib-mean-chart-1', "json", data_source)
            #
            library_list_stats ['lib_mean_graphic'] = lib_mean_lib_graphic.render()

            heading = 'Number of Bases per Library kits'
            data_source = graphic_for_library_kit (heading, 'Number of Bases per Library kits ' ,'Library Names', 'Number of Bases ', '', n_bases_comparations)
            lib_mean_lib_graphic = FusionCharts("column3d", 'lib-n_bases-lib' , "500", "300", 'lib-n_bases-chart-1', "json", data_source)
            #
            library_list_stats ['lib_n_bases_graphic'] = lib_mean_lib_graphic.render()

            # Create the graphic for number of time that library has been used

            count_libraries = {}
            for library_used in all_libraries :
                lib_name = library_used.get_index_library_name()
                if not lib_name in count_libraries :
                    count_libraries[lib_name] = 1
                else:
                    count_libraries[lib_name] +=1
            data_source = pie_graphic ('Library utilization in projects', 'fint',count_libraries)
            libraries_kit_utilization = FusionCharts("pie3d", "lib_kit_utilization_graph-1" , "500", "400", "lib_kit_utilization_chart-1", "json", data_source)
            library_list_stats ['libraries_kit_utilization'] = libraries_kit_utilization.render()

            return render (request,'iSkyLIMS_wetlab/StatsPerLibrary.html', {'display_list_of_library_stats': library_list_stats })

    else:
        return render (request,'iSkyLIMS_wetlab/StatsPerLibrary.html')

@login_required
def annual_report (request) :
    # check user privileges
    if request.user.is_authenticated:
        try:
            groups = Group.objects.get(name='WetlabManager')
            if groups not in request.user.groups.all():
                return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
        except:
            return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
    else:
        #redirect to login webpage
        return redirect ('/accounts/login')

    if request.method=='POST' :
        year_selected = int(request.POST['yearselected'])
        # get the current year to compare with the input
        present_year = datetime.datetime.now().year
        if year_selected > present_year:
            return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['Annual Report cannot be done on the future  ',
                            'the input year in the Form  ',year_selected , 'is not allowed']})

        completed_run_in_year = RunProcess.objects.filter(run_date__year = year_selected, state__runStateName__exact = 'Completed')
        #
        uncompleted_run_in_year = RunProcess.objects.filter(run_date__year = year_selected).exclude(state__runStateName__exact = 'Completed')
        if len (completed_run_in_year)  == 0 and len (uncompleted_run_in_year) == 0:
            return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['Annual Report cannot be generated because there is no runs performed the year ', year_selected ]})

        annual_report_information = {}
        annual_report_information['year'] = year_selected
        number_of_runs = {}
        number_of_runs['Completed Runs'] = 0
        number_of_runs['Not Finish Runs'] = 0
        if len ( completed_run_in_year) > 0 :
            completed_run = []
            for run in completed_run_in_year :
                completed_run.append(run.get_run_name)
            annual_report_information['completed_run'] = completed_run
            number_of_runs['Completed Runs'] = len ( completed_run_in_year)
        if len ( uncompleted_run_in_year) > 0 :
            uncompleted_run = []
            for run_uncompleted in uncompleted_run_in_year :
                uncompleted_run.append(run_uncompleted.get_run_name)
            annual_report_information['uncompleted_run'] = uncompleted_run
            number_of_runs['Not Finish Runs'] = len ( uncompleted_run_in_year)
        # prepare the pie graphic for the number of completed/ unfinished runs
        data_source = pie_graphic_standard('Number of Runs performed on the year', "",'ocean',number_of_runs)
        graphic_completed_run = FusionCharts("pie3d", "ex1" , "400", "300", "chart-1", "json", data_source)
        annual_report_information ['graphic_completed_run'] = graphic_completed_run.render()

        #
        ### Collecting information from StatsRunSummary
        run_found_bin_summary_year = StatsRunSummary.objects.filter(stats_summary_run_date__year = year_selected, level__exact = 'Total')
        q30_year, aligned_year, error_rate_year  = {} , {} , {}
        for run_bin_summary in run_found_bin_summary_year :
            bin_summary_data = run_bin_summary.get_bin_run_summary().split(';')
            run_name = run_bin_summary.runprocess_id.get_run_name()
            aligned_year[run_name]= bin_summary_data[2]
            error_rate_year[run_name]= bin_summary_data[3]
            q30_year[run_name]= bin_summary_data[5]
        annual_report_information ['aligned_data'] = aligned_year
        annual_report_information ['error_rate_data'] = error_rate_year
        annual_report_information ['q30_data'] = q30_year
        # graphics for StatsRunSummary
        heading = 'Aligned % for the runs done on year '+ str(year_selected )
        data_source = column_graphic_for_year_report (heading, 'Aligned  ' , 'Run names ', 'Aligned (in %)', 'ocean', aligned_year)
        aligned_year_graphic = FusionCharts("column3d", 'aligned_year' , "600", "300", 'aligned_chart-3', "json", data_source)
        annual_report_information ['aligned_graphic'] = aligned_year_graphic.render()

        heading = 'Error Rate for the runs done on year '+ str(year_selected )
        data_source = column_graphic_for_year_report (heading, 'Error rate ' , 'Run names ', 'Error rate', 'carbon', error_rate_year)
        error_rate_year_graphic = FusionCharts("column3d", 'error_rate_year' , "600", "300", 'error_rate_chart-4', "json", data_source)
        annual_report_information ['error_rate_graphic'] = error_rate_year_graphic.render()

        heading = '>Q30 for the runs done on year '+ str(year_selected )
        data_source = column_graphic_for_year_report (heading, 'Q30  ' , 'Run names ', '>Q 30 (in %)', 'fint', q30_year)
        q30_year_graphic = FusionCharts("column3d", 'q30_year' , "600", "300", 'q30_chart-2', "json", data_source)
        #
        annual_report_information ['q30_graphic'] = q30_year_graphic.render()
        #

        # Get the information for investigator name and the projects done
        # number_proyects_investigator contains a dict with 3 ranges 1-5, 6-10, more than 11
        investigator_projects = Projects.objects.filter(project_run_date__year = year_selected).order_by('user_id')
        project_by_user = {}
        investigator_5_project, investigator_10_project, investigator_more_10_project = {}, {} , {}
        #
        for investigator in investigator_projects:
            user_name = investigator.get_user_name()
            if user_name in project_by_user:
                project_by_user [user_name].append(investigator.get_project_name())
            else:
                project_by_user [user_name]=([investigator.get_project_name()])
        for key, value in project_by_user.items():
            if len(value) <= 5 :
                investigator_5_project[key]= value
            elif len (value) <=10:
                investigator_10_project[key]= value
            else:
                investigator_more_10_project[key]= value
        annual_report_information['user_5_projects'] = investigator_5_project
        annual_report_information['user_10_projects'] = investigator_10_project
        annual_report_information['user_more_10_projects'] = investigator_more_10_project

        # Create the bar graphic for user projects
        p_user_year ={}
        p_user_year['1 - 5']= len(investigator_5_project)
        p_user_year['6 - 10']= len(investigator_10_project)
        p_user_year['more than 10']= len(investigator_more_10_project)
        heading = 'Projects done per investigator on year '+ str(year_selected )
        data_source = column_graphic_for_year_report (heading, '  ' , 'Projects ', 'number of users', 'ocean', p_user_year)
        p_user_year_graphic = FusionCharts("column3d", 'bar_project_user_year' , "400", "300", 'p_user_chart-1', "json", data_source)
        annual_report_information ['p_user_year_graphic'] = p_user_year_graphic.render()

        data_source = pie_graphic_standard (heading, 'Percentage' ,'carbon', p_user_year)
        pie_p_user_year_graphic = FusionCharts("pie3d", "pie_project_user_year" , "400", "300", "p_user_chart-2", "json", data_source)
        annual_report_information ['pie_p_user_year_graphic'] = pie_p_user_year_graphic.render()
        #
        return render (request, 'iSkyLIMS_wetlab/AnnualReport.html',{'display_annual_report': annual_report_information})
    else:
        return render (request, 'iSkyLIMS_wetlab/AnnualReport.html')

@login_required
def monthly_report (request) :
    # check user privileges
    if request.user.is_authenticated:
        try:
            groups = Group.objects.get(name='WetlabManager')
            if groups not in request.user.groups.all():
                return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
        except:
            return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
    else:
        #redirect to login webpage
        return redirect ('/accounts/login')

    if request.method=='POST' :

        input_value = request.POST['month_year_selected']
        browser_used = request.META['HTTP_USER_AGENT']
        if 'Firefox' in browser_used :
            try:
                datetime.datetime.strptime(input_value, '%m-%Y')
                month_selected, year_selected = input_value.split('-')
            except:
                return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['Input field does not have the right format  ',
                            'the right input format is MM-YYYY   the entry ', input_value , ' is not allowed']})

        else:
            try:
                datetime.datetime.strptime(input_value, '%Y-%m')
                year_selected , month_selected = input_value.split('-')
            except:
                return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['Monthly Report input field does not have the right format  ',
                            'the right input format is MM-YYYY  ' ,' the entry ', input_value , ' is not allowed']})

        # get the current year to compare with the input
        present_year = datetime.datetime.now().year
        #
        if (int(year_selected) > present_year) :
            return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['Monthly Report cannot be done on the future  ',
                            'the input year in the Form  ', year_selected , 'is not allowed']})

        present_month = datetime.datetime.now().month
        if (int(year_selected) == present_year) and (int(month_selected) > present_month) :
            return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['Monthly Report cannot be done on the future  ',
                            'the input month in the Form  ', month_selected , 'is not allowed']})

        completed_run_in_year_month = RunProcess.objects.filter(run_date__year = year_selected,  run_date__month = month_selected ,state__runStateName__exact = 'Completed')
        #
        uncompleted_run_in_year_month = RunProcess.objects.filter(run_date__year = year_selected, run_date__month = month_selected).exclude(state__runStateName__exact = 'Completed')
        if len (completed_run_in_year_month)  == 0 and len (uncompleted_run_in_year_month) == 0:
            return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['Montly Report cannot be generated because there is no runs performed the year ', year_selected ]})

        monthly_report_information = {}
        monthly_report_information['month_year'] = str(month_selected + '  ' +  year_selected )
        number_of_runs = {}
        number_of_runs['Completed Runs'] = 0
        number_of_runs['Not Finish Runs'] = 0
        if len ( completed_run_in_year_month) > 0 :
            completed_run = []
            for run in completed_run_in_year_month :
                completed_run.append(run.get_run_name)
            monthly_report_information['completed_run'] = completed_run
            number_of_runs['Completed Runs'] = len ( completed_run_in_year_month)
        if len ( uncompleted_run_in_year_month) > 0 :
            uncompleted_run = []
            for run_uncompleted in uncompleted_run_in_year_month :
                uncompleted_run.append(run_uncompleted.get_run_name)
            monthly_report_information['uncompleted_run'] = uncompleted_run
            number_of_runs['Not Finish Runs'] = len ( uncompleted_run_in_year_month)
        # prepare the pie graphic for the number of completed/ unfinished runs
        heading = str ('Graphics of the Runs performed on the ' + month_selected + ' - ' + year_selected)
        data_source = pie_graphic_standard(heading, "",'ocean',number_of_runs)
        graphic_completed_run = FusionCharts("pie3d", "ex1" , "400", "300", "chart-1", "json", data_source)
        monthly_report_information ['graphic_completed_run'] = graphic_completed_run.render()

        # Get the information for investigator name and the projects done
        # number_proyects_investigator contains a dict with 3 ranges 1, 2, more than 2
        investigator_projects = Projects.objects.filter(project_run_date__year = year_selected, project_run_date__month = month_selected).order_by('user_id')
        project_by_user = {}
        investigator_1_project, investigator_2_projects, investigator_more_2_projects = {}, {} , {}
        #
        for investigator in investigator_projects:
            user_name = investigator.get_user_name()
            if user_name in project_by_user:
                project_by_user [user_name].append(investigator.get_project_name())
            else:
                project_by_user [user_name]=([investigator.get_project_name()])
        for key, value in project_by_user.items():

            if len(value) == 1 :
                investigator_1_project[key]= value
            elif len (value) == 2:
                investigator_2_projects[key]= value
            else:
                investigator_more_2_projects[key]= value
        monthly_report_information['user_1_project'] = investigator_1_project
        monthly_report_information['user_2_projects'] = investigator_2_projects
        monthly_report_information['user_more_2_projects'] = investigator_more_2_projects

        # Create the bar graphic for user projects
        p_user_month ={}
        p_user_month['1 project']= len(investigator_1_project)
        p_user_month['2 projects']= len(investigator_2_projects)
        p_user_month['more than 2']= len(investigator_more_2_projects)

        heading = 'Projects done per investigator on '+ str(month_selected + ' - ' + year_selected )
        data_source = column_graphic_for_year_report (heading, '  ' , 'Projects ', 'number of users', 'ocean', p_user_month)
        p_user_monthly_graphic = FusionCharts("column3d", 'bar_project_user_month' , "400", "300", 'p_user_chart-1', "json", data_source)
        monthly_report_information ['p_user_monthly_graphic'] = p_user_monthly_graphic.render()

        data_source = pie_graphic_standard (heading, 'Percentage' ,'carbon', p_user_month)
        pie_p_user_monthly_graphic = FusionCharts("pie3d", "pie_project_user_month" , "400", "300", "p_user_chart-2", "json", data_source)
        monthly_report_information ['pie_p_user_monthly_graphic'] = pie_p_user_monthly_graphic.render()

        ### Collecting information from StatsRunSummary
        run_found_bin_summary_month = StatsRunSummary.objects.filter(stats_summary_run_date__year = year_selected, stats_summary_run_date__month = month_selected, level__exact = 'Total')
        q30_month, aligned_month, error_rate_month  = {} , {} , {}
        for run_bin_summary in run_found_bin_summary_month :
            bin_summary_data = run_bin_summary.get_bin_run_summary().split(';')
            run_name = run_bin_summary.runprocess_id.get_run_name()
            aligned_month[run_name]= bin_summary_data[2]
            error_rate_month[run_name]= bin_summary_data[3]
            q30_month[run_name]= bin_summary_data[5]
        monthly_report_information ['aligned_data'] = aligned_month
        monthly_report_information ['error_rate_data'] = error_rate_month
        monthly_report_information ['q30_data'] = q30_month
        # graphics for StatsRunSummary
        heading = 'Aligned % for the runs done on '+ str(month_selected + ' - ' + year_selected)
        data_source = column_graphic_for_year_report (heading, 'Aligned  ' , 'Run names ', 'Aligned (in %)', 'ocean', aligned_month)
        aligned_month_graphic = FusionCharts("column3d", 'aligned_year' , "600", "300", 'aligned_chart-3', "json", data_source)
        monthly_report_information ['aligned_graphic'] = aligned_month_graphic.render()

        heading = 'Error Rate for the runs done on  '+ str(month_selected + ' - ' + year_selected)
        data_source = column_graphic_for_year_report (heading, 'Error rate ' , 'Run names ', 'Error rate', 'carbon', error_rate_month)
        error_rate_month_graphic = FusionCharts("column3d", 'error_rate_year' , "600", "300", 'error_rate_chart-4', "json", data_source)
        monthly_report_information ['error_rate_graphic'] = error_rate_month_graphic.render()

        heading = '>Q30 for the runs done on  '+ str(month_selected + ' - ' + year_selected)
        data_source = column_graphic_for_year_report (heading, 'Q30  ' , 'Run names ', '>Q 30 (in %)', 'fint', q30_month)
        q30_month_graphic = FusionCharts("column3d", 'q30_year' , "600", "300", 'q30_chart-2', "json", data_source)

        monthly_report_information ['q30_graphic'] = q30_month_graphic.render()

        return render (request, 'iSkyLIMS_wetlab/MonthlyReport.html',{'display_monthly_report': monthly_report_information})
    else:
        return render (request, 'iSkyLIMS_wetlab/MonthlyReport.html')

@login_required
def quarter_report (request) :
    # check user privileges
    if request.user.is_authenticated:
        try:
            groups = Group.objects.get(name='WetlabManager')
            if groups not in request.user.groups.all():
                return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
        except:
            return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
    else:
        #redirect to login webpage
        return redirect ('/accounts/login')

    if request.method=='POST' :
        year_selected = request.POST['yearselected']
        quarter_selected = int(request.POST['quarter'])
        quarter_string = ['', 'First Quarter (January -- March) ', 'Second Quarter (April -- June) ',
                            'Third Quarter (July -- September) ', 'Fourth Quarter (October -- Decemmber) ' ]
        days_in_end_quarter = ['0','31','30','30','31']
        start_quarter = str(quarter_selected *3 -2)
        end_quarter = str(quarter_selected *3)
        start_date = str(year_selected + '-' + start_quarter + '-01')
        end_date = str(year_selected + '-' + end_quarter + '-' + days_in_end_quarter[quarter_selected])
        # get the current year to compare with the input
        present_year = datetime.datetime.now().year
        if int (year_selected) > present_year:
            return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['Quarter Report cannot be done on the future  ',
                            'the input year in the Form  ',year_selected , 'is not allowed']})

        present_month = datetime.datetime.now().month
        if (int(year_selected) == present_year) and (int(end_quarter) > present_month) :
            return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['Quater Report cannot be done on the future  ',
                            'the selected Quarter ', quarter_string [quarter_selected] + str(year_selected) , 'is not allowed']})

        #
        completed_run_in_quarter = RunProcess.objects.filter( run_date__range =(start_date, end_date) , state__runStateName = 'Completed')
        #
        uncompleted_run_in_quarter = RunProcess.objects.filter(run_date__range =(start_date, end_date)).exclude(state__runStateName = 'Completed')
        if len (completed_run_in_quarter)  == 0 and len (uncompleted_run_in_quarter) == 0:
            return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['Quater Report cannot be generated because there is no runs performed the Quarter ',
                            quarter_string [quarter_selected] + str(year_selected) ]})

        quarter_report_information = {}
        quarter_report_information['quarter_year'] = quarter_string [quarter_selected] + str(year_selected)
        number_of_runs = {}
        number_of_runs['Completed Runs'] = 0
        number_of_runs['Not Finish Runs'] = 0
        if len ( completed_run_in_quarter) > 0 :
            completed_run = []
            for run in completed_run_in_quarter :
                completed_run.append(run.get_run_name)
            quarter_report_information['completed_run'] = completed_run
            number_of_runs['Completed Runs'] = len ( completed_run_in_quarter)
        if len ( uncompleted_run_in_quarter) > 0 :
            uncompleted_run = []
            for run_uncompleted in uncompleted_run_in_quarter :
                uncompleted_run.append(run_uncompleted.get_run_name)
            quarter_report_information['uncompleted_run'] = uncompleted_run
            number_of_runs['Not Finish Runs'] = len ( uncompleted_run_in_quarter)
        # prepare the pie graphic for the number of completed/ unfinished runs
        data_source = pie_graphic_standard('Number of Runs performed on the year', "",'ocean',number_of_runs)
        graphic_completed_run = FusionCharts("pie3d", "ex1" , "400", "300", "chart-1", "json", data_source)
        quarter_report_information ['graphic_completed_run'] = graphic_completed_run.render()

        #
        ### Collecting information from StatsRunSummary
        run_found_bin_summary_quarter = StatsRunSummary.objects.filter(stats_summary_run_date__range = (start_date, end_date), level__exact = 'Total')
        q30_quarter, aligned_quarter, error_rate_quarter  = {} , {} , {}
        for run_bin_summary in run_found_bin_summary_quarter :
            bin_summary_data = run_bin_summary.get_bin_run_summary().split(';')
            run_name = run_bin_summary.runprocess_id.get_run_name()
            aligned_quarter[run_name]= bin_summary_data[2]
            error_rate_quarter[run_name]= bin_summary_data[3]
            q30_quarter[run_name]= bin_summary_data[5]
        quarter_report_information ['aligned_data'] = aligned_quarter
        quarter_report_information ['error_rate_data'] = error_rate_quarter
        quarter_report_information ['q30_data'] = q30_quarter
        # graphics for StatsRunSummary
        heading = 'Aligned % for the runs done on  ' + quarter_string [quarter_selected] + str(year_selected)
        data_source = column_graphic_for_year_report (heading, 'Aligned  ' , 'Run names ', 'Aligned (in %)', 'ocean', aligned_quarter)
        aligned_quarter_graphic = FusionCharts("column3d", 'aligned_year' , "600", "300", 'aligned_chart-3', "json", data_source)
        quarter_report_information ['aligned_graphic'] = aligned_quarter_graphic.render()

        heading = 'Error Rate for the runs done on year '+ quarter_string [quarter_selected] + str(year_selected)
        data_source = column_graphic_for_year_report (heading, 'Error rate ' , 'Run names ', 'Error rate', 'carbon', error_rate_quarter)
        error_rate_quarter_graphic = FusionCharts("column3d", 'error_rate_year' , "600", "300", 'error_rate_chart-4', "json", data_source)
        quarter_report_information ['error_rate_graphic'] = error_rate_quarter_graphic.render()

        heading = '>Q30 for the runs done on year '+ quarter_string [quarter_selected] + str(year_selected)
        data_source = column_graphic_for_year_report (heading, 'Q30  ' , 'Run names ', '>Q 30 (in %)', 'fint', q30_quarter)
        q30_quarter_graphic = FusionCharts("column3d", 'q30_year' , "600", "300", 'q30_chart-2', "json", data_source)
        #
        quarter_report_information ['q30_graphic'] = q30_quarter_graphic.render()
        #

        # Get the information for investigator name and the projects done
        # number_proyects_investigator contains a dict with 3 ranges 1-5, 6-10, more than 11
        investigator_projects = Projects.objects.filter(project_run_date__range = (start_date, end_date)).order_by('user_id')
        project_by_user = {}
        investigator_5_project, investigator_10_project, investigator_more_10_project = {}, {} , {}
        #
        for investigator in investigator_projects:
            user_name = investigator.get_user_name()
            if user_name in project_by_user:
                project_by_user [user_name].append(investigator.get_project_name())
            else:
                project_by_user [user_name]=([investigator.get_project_name()])
        for key, value in project_by_user.items():
            if len(value) <= 5 :
                investigator_5_project[key]= value
            elif len (value) <=10:
                investigator_10_project[key]= value
            else:
                investigator_more_10_project[key]= value
        quarter_report_information['user_5_projects'] = investigator_5_project
        quarter_report_information['user_10_projects'] = investigator_10_project
        quarter_report_information['user_more_10_projects'] = investigator_more_10_project

        # Create the bar graphic for user projects
        p_user_quarter ={}
        p_user_quarter['1 - 5']= len(investigator_5_project)
        p_user_quarter['6 - 10']= len(investigator_10_project)
        p_user_quarter['more than 10']= len(investigator_more_10_project)
        heading = 'Projects done per investigator on year '+ str(year_selected )
        data_source = column_graphic_for_year_report (heading, '  ' , 'Projects ', 'number of users', 'ocean', p_user_quarter)
        p_user_quarter_graphic = FusionCharts("column3d", 'bar_project_user_year' , "400", "300", 'p_user_chart-1', "json", data_source)
        quarter_report_information ['p_user_year_graphic'] = p_user_quarter_graphic.render()

        data_source = pie_graphic_standard (heading, 'Percentage' ,'carbon', p_user_quarter)
        pie_p_user_quarter_graphic = FusionCharts("pie3d", "pie_project_user_year" , "400", "300", "p_user_chart-2", "json", data_source)
        quarter_report_information ['pie_p_user_year_graphic'] = pie_p_user_quarter_graphic.render()
        #
        return render (request, 'iSkyLIMS_wetlab/QuarterReport.html',{'display_quarter_report': quarter_report_information})
    else:
        return render (request, 'iSkyLIMS_wetlab/QuarterReport.html')

'''
def open_samba_connection ():

    from smb.SMBConnection import SMBConnection
    ##conn=SMBConnection('bioinfocifs', 'fCdEg979I-W.gUx-teDr', 'NGS_Data', 'quibitka', use_ntlm_v2=True)
    conn=SMBConnection(wetlab_config.SAMBA_USER_ID, wetlab_config.SAMBA_USER_PASSWORD, wetlab_config.SAMBA_SHARED_FOLDER_NAME,wetlab_config.SAMBA_REMOTE_SERVER_NAME, use_ntlm_v2=True)
    ##conn.connect('172.21.7.11', 445)
    conn.connect(wetlab_config.SAMBA_IP_SERVER, 445)
    return conn
'''
def get_size_dir (directory, conn, ):
    count_file_size = 0
    file_list = conn.listPath(wetlab_config.SAMBA_SHARED_FOLDER_NAME, directory)
    for sh_file in file_list:
        if sh_file.isDirectory:
            if (sh_file.filename == '.' or sh_file.filename == '..'):
                continue

            sub_directory = os.path.join (directory,sh_file.filename)
            count_file_size += get_size_dir (sub_directory, conn)
        else:
            count_file_size += sh_file.file_size

    return count_file_size




def update_tables (request):
    #### Update the run date for the projects. StatsBinRunRead and  StatsRunSummary
    #### tables when they were not updated. It takes the run date from the run date
    '''
    run_founds = RunProcess.objects.all()
    for run in run_founds :
        run_id = run.id
        run_date = run.run_date
        projects_to_update = Projects.objects.filter(runprocess_id__exact = run_id)
        for project in projects_to_update :
            project.project_run_date = run_date
            project.save()
        stats_run_to_update = StatsRunSummary.objects.filter(runprocess_id__exact = run_id)
        for stats_run in stats_run_to_update :
            stats_run.stats_summary_run_date = run_date
            stats_run.save()
        stats_read_to_update = StatsBinRunRead.objects.filter(runprocess_id__exact = run_id)
        for  stats_read in stats_read_to_update :
            stats_read.stats_read_run_date = run_date
            stats_read.save()
    return render(request, 'iSkyLIMS_wetlab/info_page.html', {'content':['The tables have been updated']})
    '''
    ### Update the disc space used of each run
    ### It will connect to quibitka to get the size of the file for each run


    #conn=open_samba_connection()
    if RunProcess.objects.filter(state__runStateName ='Completed', useSpaceImgMb = 0).exists():
        conn = open_samba_connection()
        run_list_be_updated = RunProcess.objects.filter(state__runStateName = 'Completed' , useSpaceImgMb =0 )
        for run_be_updated in run_list_be_updated:
            run_id = run_be_updated.id
            run_parameter_id=RunningParameters.objects.get(pk=run_id)

            runID_value = run_parameter_id.RunID
            get_full_list = conn.listPath(wetlab_config.SAMBA_SHARED_FOLDER_NAME ,runID_value)
            rest_of_dir_size = 0
            data_dir_size = 0
            images_dir_size = 0
            in_mega_bytes = 1024*1024
            #
            for item_list in get_full_list:
                if item_list.filename == '.' or item_list.filename == '..':
                    continue
                if item_list.filename == 'Data':
                    dir_data = os.path.join(runID_value,'Data')
                    data_dir_size = get_size_dir(dir_data , conn)
                    continue

                elif item_list.filename == 'Images':
                    dir_images = os.path.join(runID_value, 'Images')
                    images_dir_size = get_size_dir(dir_images , conn)
                    continue

                if item_list.isDirectory:
                    item_dir = os.path.join(runID_value, item_list.filename)
                    rest_of_dir_size += get_size_dir(item_dir, conn)
                else:
                    rest_of_dir_size += item_list.file_size
            #
            # format file space and save it into database
            data_dir_size_formated = '{0:,}'.format(round(data_dir_size/in_mega_bytes))
            images_dir_size_formated = '{0:,}'.format(round(images_dir_size/in_mega_bytes))
            rest_of_dir_size_formated = '{0:,}'.format(round(rest_of_dir_size/in_mega_bytes))
            run_be_updated.useSpaceImgMb= images_dir_size_formated
            run_be_updated.useSpaceFastaMb= data_dir_size_formated
            run_be_updated.useSpaceOtherMb= rest_of_dir_size_formated
            #
            run_be_updated.save()

        '''

        get_full_list = conn.listPath('NGS_Data' ,run_name)
        rest_of_dir_size = 0
        data_dir_size = 0
        images_dir_size = 0


        conn.close()

        '''
        return render(request, 'iSkyLIMS_wetlab/info_page.html', {'content':['The Disk space usage have been updated']})
    else:
        return render(request, 'iSkyLIMS_wetlab/error_page.html', {'content':['There is no tables which requiered to update with Disk space usage information']})

def update_tables_date (request):
    if RunProcess.objects.filter(state__runStateName ='Completed', run_finish_date = None).exists():
        #
        conn = open_samba_connection()
        run_list_be_updated = RunProcess.objects.filter(state__runStateName = 'Completed' , run_finish_date = None )
        for run_be_updated in run_list_be_updated:
            run_id = run_be_updated.id
            run_parameter_id=RunningParameters.objects.get(pk=run_id)
            runID_value = run_parameter_id.RunID
            completion_file = os.path.join(runID_value, 'RunCompletionStatus.xml')
            try:
                completion_attributes = conn.getAttributes(wetlab_config.SAMBA_SHARED_FOLDER_NAME ,completion_file)
                # fetching the time creation on the RunCompletionStatus.xml for Run finish datetime
                run_be_updated.run_finish_date = datetime.datetime.fromtimestamp(int(completion_attributes.create_time)).strftime('%Y-%m-%d %H:%M:%S')
            except:
                pass
            conversion_stats_file = os.path.join (runID_value,'Data/Intensities/BaseCalls/Stats/', 'ConversionStats.xml')
            try:
                conversion_attributes = conn.getAttributes(wetlab_config.SAMBA_SHARED_FOLDER_NAME ,conversion_stats_file)
                #
                run_be_updated.bcl2fastq_finish_date = datetime.datetime.fromtimestamp(int(conversion_attributes.create_time)).strftime('%Y-%m-%d %H:%M:%S')
            except:
                pass
            finish_process_date = StatsRunSummary.objects.filter(runprocess_id__exact = run_id)
            run_be_updated.process_completed_date = finish_process_date[0].generatedat

            run_be_updated.save()

        return render(request, 'iSkyLIMS_wetlab/info_page.html', {'content':['The dates for the Runs have been updated']})
    else:
        return render(request, 'iSkyLIMS_wetlab/error_page.html', {'content':['There is no tables which requiered to update with date information']})


@login_required
def configuration_test (request):
    # check user privileges
    if request.user.is_authenticated:
        if not request.user.is_staff or not request.user.is_superuser:
            return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['You do have the enough privileges to see this page ','Contact with your administrator .']})
    else:
        #redirect to login webpage
        return redirect ('/accounts/login')
    if request.method=='POST' and request.POST['action'] == 'basicTest':
        test_results = {}
        wetlab_config_file = os.path.join(settings.BASE_DIR, 'iSkyLIMS_wetlab', 'wetlab_config.py')
        test_results['iSkyLIMS_settings'] = get_iSkyLIMS_settings()
        test_results['config_file'] = get_config_file(wetlab_config_file)
        test_results['attr_files'] = get_files_attribute(os.path.join(settings.MEDIA_ROOT, 'wetlab'))
        test_results['database_access'] = check_access_database()
        test_results['samba_connection'] = check_samba_connection()

        test_results['basic_checks_ok'] = 'OK'
        #if test_results['config_file']  and test_results['attr_files']  and test_results['database_access'] and test_results['samba_connection']:
        for result in test_results :
            if test_results[result] == 'NOK':
                test_results['basic_checks_ok'] = 'NOK'
                break
        available_run_test = []
        if RunConfigurationTest.objects.all().exists():
            run_test_objs = RunConfigurationTest.objects.all()
            for run_test_obj in run_test_objs:
                available_run_test.append([run_test_obj.get_run_test_name(),run_test_obj.get_run_test_id()])
        return render (request,'iSkyLIMS_wetlab/ConfigurationTest.html', {'test_results': test_results, 'available_run_test':available_run_test})


    elif request.method=='POST' and request.POST['action'] == 'executeRunTest':
        if RunConfigurationTest.objects.filter(pk__exact = request.POST['runTest']).exists():
            run_test_obj =  RunConfigurationTest.objects.filter(pk__exact = request.POST['runTest']).last()
            run_test_folder = run_test_obj.get_run_test_folder()
            run_test_name = run_test_obj.get_run_test_name()
        if not folder_test_exists(run_test_folder):
            return render(request,'iSkyLIMS_wetlab/ConfigurationTest.html',{'error':wetlab_config.ERROR_NOT_FOLDER_RUN_TEST_WAS_FOUND})
        run_test_result = execute_test_for_testing_run(run_test_name,run_test_folder)
        run_test_result['run_test_name'] = run_test_name
        logger = logging.getLogger(__name__)
        if 'ERROR' in run_test_result :
            log_trace = []
            with open (logging.getLoggerClass().root.handlers[0].baseFilename, 'r') as fh :
                for line in fh :
                    if run_test_name in line:
                        line = line.replace('\n', '')
                        log_trace.append(line)

            return render (request,'iSkyLIMS_wetlab/ConfigurationTest.html', {'run_test_result': run_test_result, 'log_trace': log_trace})
        else:
            return render (request,'iSkyLIMS_wetlab/ConfigurationTest.html', {'run_test_result': run_test_result})
    elif request.method=='POST' and request.POST['action'] == 'deleteTestRun':
        if RunConfigurationTest.objects.filter(runTestName__exact = request.POST['deleteRun']).exists():
            run_test_objs = RunConfigurationTest.objects.filter(runTestName__exact = request.POST['deleteRun'])
            for run_test_obj in run_test_objs:
                delete_test_run (run_test_obj)
            return render(request,'iSkyLIMS_wetlab/ConfigurationTest.html')
    else:
        return render(request,'iSkyLIMS_wetlab/ConfigurationTest.html')

@login_required
def create_protocol (request):
    ## Check user == WETLAB_MANAGER: if false,  redirect to 'login' page
    if request.user.is_authenticated:
        if not is_wetlab_manager(request):
            return render (
                request,'iSkyLIMS_wetlab/error_page.html',
                {'content':['You do not have enough privileges to see this page ',
                            'Contact with your administrator .']})
    else:
        #redirect to login webpage
        return redirect ('/accounts/login')
    # get the list of defined protocols
    defined_protocols, other_protocol_list = display_available_protocols (__package__)
    additional_kits = get_additional_kits_list (__package__)
    defined_protocol_types = display_protocol_types (__package__)


    if request.method == 'POST' and request.POST['action'] == 'addNewProtocol':
        new_protocol = request.POST['newProtocolName']
        protocol_type = request.POST['protocolType']
        description = request.POST['description']

        if check_if_protocol_exists (new_protocol, __package__):
            return render ( request,'iSkyLIMS_wetlab/error_page.html',{'content':['Protocol Name ', new_protocol,
                            'Already exists.']})
        new_protocol_id = create_new_protocol(new_protocol, protocol_type, description, __package__)

        return render(request, 'iSkyLIMS_wetlab/createProtocol.html',{'defined_protocols': defined_protocols,
                            'defined_protocol_types':defined_protocol_types, 'new_defined_protocol': new_protocol,
                            'new_protocol_id':new_protocol_id,  'other_protocol_list' :other_protocol_list})

    return render(request, 'iSkyLIMS_wetlab/createProtocol.html',{'defined_protocols': defined_protocols,
                        'defined_protocol_types':defined_protocol_types, 'other_protocol_list' :other_protocol_list,
                        'additional_kits': additional_kits})


@login_required
def define_sample_projects (request):
    ## Check user == WETLAB_MANAGER: if false,  redirect to 'login' page
    if request.user.is_authenticated:
        if not is_wetlab_manager(request):
            return render (
                request,'iSkyLIMS_wetlab/error_page.html',
                {'content':['You do not have enough privileges to see this page ',
                            'Contact with your administrator .']})
    else:
        #redirect to login webpage
        return redirect ('/accounts/login')
    # get the information of defined sample Projects
    defined_samples_projects = get_info_for_defined_sample_projects (__package__)

    if request.method == 'POST' and request.POST['action'] == 'addNewSampleProject':
        sample_project_name = request.POST['sampleProyectName']
        #description = request.POST['description']

        if check_if_sample_project_exists (sample_project_name, __package__):
            error_message = ERROR_SAMPLE_PROJECT_ALREADY_EXISTS
            return render ( request,'iSkyLIMS_wetlab/createSampleProjects.html',{'defined_samples_projects': defined_samples_projects,
                                    'error_message' : error_message})
        new_sample_project_id = create_new_sample_project (request.POST, __package__)
        new_defined_sample_project = sample_project_name
        return render(request, 'iSkyLIMS_wetlab/createSampleProjects.html',{'defined_samples_projects': defined_samples_projects,
                             'new_sample_project_id': new_sample_project_id, 'new_defined_sample_project' : new_defined_sample_project})

    return render(request, 'iSkyLIMS_wetlab/createSampleProjects.html',{'defined_samples_projects': defined_samples_projects})

def define_additional_kits(request, protocol_id):
    ## Check user == WETLAB_MANAGER: if false,  redirect to 'login' page
    if request.user.is_authenticated:
        if not is_wetlab_manager(request):
            return render ( request,'iSkyLIMS_wetlab/error_page.html',
                {'content':['You do not have enough privileges to see this page ',
                            'Contact with your administrator .']})
    else:
        #redirect to login webpage
        return redirect ('/accounts/login')

    additional_kits = define_table_for_additional_kits(protocol_id)
    if request.method == 'POST' and request.POST['action'] == 'defineAdditionalKits':

        recorded_additional_kits = set_additional_kits(request.POST, request.user)
        if len(recorded_additional_kits) == 0:
            return render(request, 'iSkyLIMS_wetlab/defineAdditionalKits.html', {'additional_kits':additional_kits})

        return render(request, 'iSkyLIMS_wetlab/defineAdditionalKits.html', {'recorded_additional_kits':recorded_additional_kits})

    else:
        if not check_if_protocol_exists(protocol_id, __package__):
            return render ( request,'iSkyLIMS_wetlab/error_page.html',
                        {'content':['The requested Protocol does not exist',
                            'Create the protocol name before assigning additional kits for protocol.']})

        return render(request, 'iSkyLIMS_wetlab/defineAdditionalKits.html', {'additional_kits':additional_kits})


@login_required
def display_sample_project(request,sample_project_id):

    samples_project_data = get_info_to_display_sample_project (sample_project_id)
    if 'ERROR' in samples_project_data :
        error_message = ERROR_SAMPLE_PROJECT_DOES_NOT_EXISTS
        return render (request,'iSkyLIMS_wetlab/error_page.html', {'content': error_message })
    return render(request, 'iSkyLIMS_wetlab/displaySampleProject.html',{'samples_project_data': samples_project_data})

@login_required
def display_protocol (request, protocol_id):
    if not is_wetlab_manager(request):
        return render (request,'iSkyLIMS_wetlab/error_page.html',
            {'content':['You do not have enough privileges to see this page ',
                        'Contact with your administrator .']})
    if not check_if_protocol_exists(protocol_id, __package__):
        return render (request,'iSkyLIMS_wetlab/error_page.html',
            {'content':['The protocol that you are trying to get ',
                        'DOES NOT exists .']})
    protocol_data = get_all_protocol_info (protocol_id)
    kit_data = get_all_additional_kit_info(protocol_id)


    return render(request, 'iSkyLIMS_wetlab/displayProtocol.html', {'protocol_data': protocol_data, 'kit_data': kit_data})


@login_required
def define_protocol_parameters (request, protocol_id):

    ## Check user == WETLAB_MANAGER: if false,  redirect to 'login' page
    if request.user.is_authenticated:
        if not is_wetlab_manager(request):
            return render ( request,'iSkyLIMS_wetlab/error_page.html',
                {'content':['You do not have enough privileges to see this page ',
                            'Contact with your administrator .']})
    else:
        #redirect to login webpage
        return redirect ('/accounts/login')

    if request.method == 'POST' and request.POST['action'] == 'define_protocol_parameters':

        recorded_prot_parameters = set_protocol_parameters(request)

        return render(request, 'iSkyLIMS_wetlab/defineProtocolParameters.html', {'recorded_prot_parameters':recorded_prot_parameters})

    else:
        if not check_if_protocol_exists(protocol_id, __package__):
            return render ( request,'iSkyLIMS_wetlab/error_page.html',
                        {'content':['The requested Protocol does not exist',
                            'Create the protocol name before assigning custom protocol parameters.']})


        prot_parameters = define_table_for_prot_parameters(protocol_id)
        return render(request, 'iSkyLIMS_wetlab/defineProtocolParameters.html', {'prot_parameters':prot_parameters})

@login_required
def add_commercial_kit (request):
    app_name = __package__.split('.')[0]
    defined_protocols = get_defined_protocols(app_name, False)
    defined_platforms = get_defined_platforms_and_ids('NGS')
    commercial_kits_data = get_data_for_commercial_kits('NGS')

    if request.method == 'POST' and request.POST['action'] == 'addCommercialKit':
        if get_commercial_kit_id (request.POST['kitName']) :

            return render(request, 'iSkyLIMS_wetlab/addCommercialKit.html',{'defined_protocols': defined_protocols, 'invalid_name': request.POST['kitName']})
        new_kit = store_commercial_kit(request.POST)
        new_kit_data = get_commercial_kit_basic_data(new_kit)
        return render(request, 'iSkyLIMS_wetlab/addCommercialKit.html',{'new_kit_data': new_kit_data})
    else:
        return render(request, 'iSkyLIMS_wetlab/addCommercialKit.html',{'defined_protocols': defined_protocols, 'defined_platforms' : defined_platforms,
                                'commercial_kits_data': commercial_kits_data})

@login_required
def add_user_lot_commercial_kit (request):
    defined_kits = get_defined_commercial_kits()
    if request.method == 'POST' and request.POST['action'] == 'addUserLotKit':
        if get_lot_user_commercial_kit_id (request.POST['barCode']) :
            return render(request, 'iSkyLIMS_wetlab/addUserLotCommercialKit.html',{'defined_kits': defined_kits, 'invalid_name': request.POST['nickName']})
        new_lot_kit = store_lot_user_commercial_kit(request.POST, request.user)
        new_lot_kit_data = get_lot_user_commercial_kit_basic_data(new_lot_kit)
        return render(request, 'iSkyLIMS_wetlab/addUserLotCommercialKit.html',{'new_lot_kit_data':new_lot_kit_data})
    else:
        return render(request, 'iSkyLIMS_wetlab/addUserLotCommercialKit.html',{'defined_kits':defined_kits})



@login_required
def pending_to_update(request):
    pending = {}
    # get the samples in defined state

    pending['defined'] = get_samples_in_defined_state('')
    pending['extract_molecule'] = get_samples_in_extracted_molecule_state(request.user)
    pending ['graphic_pending_samples'] = pending_samples_for_grafic(pending).render()

    return render(request, 'iSkyLIMS_wetlab/pendingToUpdate.html', {'pending':pending})


@login_required
def record_samples(request):
    '''
    Functions :
        analyze_input_samples  : located at iSkyLIMS_core/handling_samples.py
        analyze_input_sample_project_fields  : located at iSkyLIMS_core/handling_samples.py
        prepare_sample_input_table : located at iSkyLIMS_core/utils/handling_samples.py
        get_codeID_for_resequencing : located at iSkyLIMS_wetlab/utils/sample_functions.py
        prepare_sample_project_input_table :  located at iSkyLIMS_core/utils/handling_samples.py
        analyze_reprocess_data  : located at iSkyLIMS_wetlab/utils/sample_functions.py
        get_info_for_reprocess_samples : located at iSkyLIMS_core/utils/handling_samples.py
    '''
    ## Record new samples
    if request.method == 'POST' and request.POST['action'] == 'recordsample':
        sample_recorded = analyze_input_samples (request, __package__)
        # if no samples are in any of the options, displays the inital page

        if (not 'defined_samples' in sample_recorded and not 'pre_defined_samples' in sample_recorded and not 'invalid_samples' in sample_recorded and not 'incomplete_samples' in sample_recorded) :
            sample_information = prepare_sample_input_table(__package__)
            return render(request, 'iSkyLIMS_wetlab/recordSample.html',{'sample_information':sample_information})

        if 'sample_id_for_action' in sample_recorded :
            sample_recorded.update(get_codeID_for_resequencing(sample_recorded))
        if 'incomplete_samples' in sample_recorded :
            sample_recorded.update(prepare_sample_input_table(__package__))
            sample_recorded['number_of_samples'] = len(sample_recorded['incomplete_samples'])
        if 'pre_defined_samples_id' in sample_recorded:
            sample_recorded.update(prepare_sample_project_input_table(sample_recorded['pre_defined_samples_id']))
        return render(request, 'iSkyLIMS_wetlab/recordSample.html',{'sample_recorded':sample_recorded})


    ## Request to reprocess the samples
    elif request.method == 'POST' and request.POST['action'] == 'reprocessSamples':
        samples = {}

        to_be_reprocessed_ids = request.POST['invalidSamplesID'].split(',')
        reprocess_id = request.POST['sampleIDforAction']
        json_data = json.loads(request.POST['reprocess_data'])

        result = analyze_reprocess_data(json_data[0], reprocess_id, request.user)
        if result == 'Invalid options':
            to_be_reprocessed_ids.insert(0,reprocess_id)
            sample_recorded = get_info_for_reprocess_samples(to_be_reprocessed_ids, reprocess_id)
            sample_recorded['invalid_samples_id'] = request.POST['invalidSamplesID']
            sample_recorded['sample_id_for_action'] = reprocess_id
            sample_recorded.update(get_codeID_for_resequencing(sample_recorded))
            sample_recorded['reprocess_result'] = 'False'
        else:
            if to_be_reprocessed_ids[0] == '':

                return render(request, 'iSkyLIMS_wetlab/recordSample.html',{'all_sucessful_reprocess':True})
            else:
                next_to_be_process_id = str(to_be_reprocessed_ids[0])
                sample_recorded = get_info_for_reprocess_samples(to_be_reprocessed_ids, next_to_be_process_id)
                sample_processed_id = to_be_reprocessed_ids.pop()
                sample_recorded['invalid_samples_id'] = ','.join(to_be_reprocessed_ids)
                sample_recorded['sample_id_for_action'] = next_to_be_process_id
                sample_recorded.update(get_codeID_for_resequencing(sample_recorded))
                sample_recorded['reprocess_result'] = 'True'
                return render(request, 'iSkyLIMS_wetlab/recordSample.html',{'sample_recorded':sample_recorded})

        if len(require_to_update) > 0 :
            for key, value in require_to_update.items():
                if value == 'newLibPreparation':
                    if not libraryPreparation.objects.filter(sample_id__pk__exact = key):
                        continue
                    lib_prep_obj = libraryPreparation.objects.get(sample_id__pk__exact = key)
                    #lib_prep_obj.set_state('Reused')
                    lib_prep_obj.set_increase_reuse()
                elif value == 'newPool':
                    pass
                else:
                    continue

        return render(request, 'iSkyLIMS_wetlab/recordSample.html',{'reprocess_result':reprocess_result})

    ## display the form to show the samples in pre-defined state that user requested to complete
    elif request.method == 'POST' and request.POST['action'] == 'select_samples_pre_defined':
        if 'samples_in_list' in request.POST :
            pre_defined_samples_id = request.POST.getlist('samples')
        sample_recorded = prepare_sample_project_input_table(pre_defined_samples_id)
        return render(request, 'iSkyLIMS_wetlab/recordSample.html',{'sample_recorded':sample_recorded})

    ## Add the additional information related to the project
    elif request.method == 'POST' and request.POST['action'] == 'sampleprojectdata':
        sample_recorded = analyze_input_sample_project_fields(request.POST)

        if request.POST['pending_pre_defined'] != '':
            sample_recorded.update(prepare_sample_project_input_table(request.POST['pending_pre_defined'].split(',')))
            return render(request, 'iSkyLIMS_wetlab/recordSample.html',{'sample_recorded':sample_recorded})
        else:
            return render(request, 'iSkyLIMS_wetlab/recordSample.html',{'sample_recorded':sample_recorded})

    ## Form to get the new samples
    else:
        sample_information = prepare_sample_input_table(__package__)
        return render(request, 'iSkyLIMS_wetlab/recordSample.html',{'sample_information':sample_information})

@login_required
def define_sample_projects_fields (request, sample_project_id):
    ## Check user == WETLAB_MANAGER: if false,  redirect to 'login' page
    if request.user.is_authenticated:
        if not is_wetlab_manager(request):
            return render (
                request,'iSkyLIMS_wetlab/error_page.html',
                {'content':['You do not have enough privileges to see this page ',
                            'Contact with your administrator .']})
    else:
        #redirect to login webpage
        return redirect ('/accounts/login')
    # get the list of defined sample Projects

    if request.method == 'POST' and request.POST['action'] == 'defineSampleProjectFields':

        sample_project_field_data = set_sample_project_fields(request.POST)

        return render(request, 'iSkyLIMS_wetlab/defineSampleProjectFields.html', {'sample_project_field_data':sample_project_field_data})

    else:
        if not check_if_sample_project_id_exists(sample_project_id):
            return render ( request,'iSkyLIMS_wetlab/error_page.html',
                        {'content':['The requested Protocol does not exist',
                            'Create the protocol name before assigning custom protocol parameters.']})


        sample_project_data = define_table_for_sample_project_fields(sample_project_id)
        return render(request, 'iSkyLIMS_wetlab/defineSampleProjectFields.html', {'sample_project_data':sample_project_data})


@login_required
def modify_additional_kits(request, protocol_id):
    ## Check user == WETLAB_MANAGER: if false,  redirect to 'login' page
    if request.user.is_authenticated:
        if not is_wetlab_manager(request):
            return render (
                request,'iSkyLIMS_wetlab/error_page.html',
                {'content':['You do not have enough privileges to see this page ',
                            'Contact with your administrator .']})
    else:
        #redirect to login webpage
        return redirect ('/accounts/login')

    if request.method == 'POST' and request.POST['action'] == 'modifyAdditionalKits':
        additional_kits_data_saved = modify_fields_in_additional_kits(request.POST, request.user)
        return render(request, 'iSkyLIMS_wetlab/modifyAdditionalKits.html', {'additional_kits_data_saved':additional_kits_data_saved})
    else:
        if not check_if_protocol_exists(protocol_id, __package__):
            return render ( request,'iSkyLIMS_wetlab/error_page.html',
                        {'content':['The requested additional kits do not exist',
                            'Create the addtional kits before.']})
        additional_kits_data = get_additional_kits_data_to_modify(protocol_id)
        return render(request, 'iSkyLIMS_wetlab/modifyAdditionalKits.html', {'additional_kits_data':additional_kits_data})




@login_required
def modify_protocol_fields(request, protocol_id):
    ## Check user == WETLAB_MANAGER: if false,  redirect to 'login' page
    if request.user.is_authenticated:
        if not is_wetlab_manager(request):
            return render (
                request,'iSkyLIMS_wetlab/error_page.html',
                {'content':['You do not have enough privileges to see this page ',
                            'Contact with your administrator .']})
    else:
        #redirect to login webpage
        return redirect ('/accounts/login')

    if request.method == 'POST' and request.POST['action'] == 'modifyProtocolFields':
        protocol_field_saved = modify_fields_in_protocol(request.POST)
        return render(request, 'iSkyLIMS_wetlab/modifyProtocolFields.html', {'protocol_field_saved':protocol_field_saved})
    else:
        if not check_if_protocol_exists(protocol_id, __package__):
            return render ( request,'iSkyLIMS_wetlab/error_page.html',
                        {'content':['The requested Protocol does not exist',
                            'Create the protocol name before assigning custom parameters.']})
        protocol_field = get_protocol_fields(protocol_id)
        return render(request, 'iSkyLIMS_wetlab/modifyProtocolFields.html', {'protocol_field':protocol_field})

@login_required
def modify_sample_project_fields(request, sample_project_id):
    ## Check user == WETLAB_MANAGER: if false,  redirect to 'login' page

    if request.user.is_authenticated:
        if not is_wetlab_manager(request):
            return render (
                request,'iSkyLIMS_wetlab/error_page.html',
                {'content':['You do not have enough privileges to see this page ',
                            'Contact with your administrator .']})
    else:
        #redirect to login webpage
        return redirect ('/accounts/login')

    if request.method == 'POST' and request.POST['action'] == 'modifySampleProjectFields':
        sample_project_field_saved = modify_fields_in_sample_project(request.POST)
        return render(request, 'iSkyLIMS_wetlab/modifySampleProjectFields.html', {'sample_project_field_saved':sample_project_field_saved})

    else:
        if not check_if_sample_project_id_exists(sample_project_id):
            return render ( request,'iSkyLIMS_wetlab/error_page.html',
                        {'content':['The requested Sample project does not exist',
                            'Create the sample project name before assigning custom sample project parameters.']})
        sample_project_field = get_parameters_sample_project(sample_project_id)
        return render(request, 'iSkyLIMS_wetlab/modifySampleProjectFields.html', {'sample_project_field':sample_project_field})


@login_required
def define_molecule_uses (request):
    '''
    Functions:
        display_molecule_use : located at iSkyLIMS_core/utils/handling_samples.py
        record_molecule_use : located at iSkyLIMS_core/utils/handling_samples.py
    '''
    molecule_use_data = display_molecule_use(__package__)
    if request.method == 'POST' and request.POST['action'] == 'record_molecule_use':
        molecule_use_data.update(record_molecule_use (request.POST, __package__))

    return render(request, 'iSkyLIMS_wetlab/defineMoleculeUses.html', {'molecule_use_data':molecule_use_data})

@login_required
def define_type_of_samples (request):
    '''
    Functions:
        display_sample_types : located at iSkyLIMS_core/utils/handling_samples.py
        save_type_of_sample : located at iSkyLIMS_core/utils/handling_samples.py
    '''
    sample_types = display_sample_types (__package__)
    if request.method == 'POST' and request.POST['action'] == 'addNewSampleType':
        sample_types.update(save_type_of_sample(request.POST, __package__))

    return render(request, 'iSkyLIMS_wetlab/defineTypeOfSamples.html', {'sample_types':sample_types})


@login_required
def display_sample (request, sample_id):
    '''
    Functions:
        get_all_sample_information : located at iSkyLIMS_core/utils/handling_samples.py
        get_all_library_information  located at iSkyLIMS_wetlab/utils/library_preparation.py
        get_additional_kits_used_in_sample   located at iSkyLIMS_wetlab/utils/additional_kits.py
        get_sample_in_project_obj_from_sample_name  # located at iSkyLIMS_wetlab/utils/sample_functions.py
    '''
    sample_information = get_all_sample_information(sample_id, True)
    if not 'Error' in sample_information:
        sample_information.update(get_molecule_lot_kit_in_sample(sample_id))
        sample_information.update(get_all_library_information(sample_id))
        sample_information.update(get_additional_kits_used_in_sample(sample_id))
    else:
        sample_information = {}
    sample_obj =get_sample_obj_from_id(sample_id)
    if sample_obj:
        sample_name = sample_obj.get_sample_name()
        run_sample_obj = get_sample_in_project_obj_from_sample_name(sample_name)
        if run_sample_obj:
            sample_information.update(get_info_sample_in_run(run_sample_obj))

    if len(sample_information) == 0  :
        return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['No Sample was found']})
    else:
        return render(request, 'iSkyLIMS_wetlab/displaySample.html',{'sample_information':sample_information})

@login_required
def display_sample_in_run (request, sample_run_id):
    '''
    Functions:
        get_info_sample_in_run
    '''
    sample_run_obj = get_sample_in_project_obj_from_id(sample_run_id)
    if not sample_run_obj :
        return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['No Sample was found']})
    sample_information = get_info_sample_in_run(sample_run_obj)
    return render(request, 'iSkyLIMS_wetlab/displaySample.html',{'sample_information':sample_information})



@login_required
def display_type_of_sample(request, sample_type_id):
    '''
    Functions:
        get_type_of_sample_information : located at iSkyLIMS_core/utils/handling_samples.py
    '''
    type_of_sample_data = get_type_of_sample_information(sample_type_id)
    return render(request, 'iSkyLIMS_wetlab/displayTypeOfSample.html',{'type_of_sample_data':type_of_sample_data})

@login_required
def handling_library_preparations(request):
    '''
    Functions:
        analyze_and_store_input_additional_kits : located at utils/additional_kits.py
        get_samples_for_library_preparation : located at utils/library_preparation.py
        check_users_exists
        extract_user_sample_sheet_data   : located at utils/library_preparation.py
        get_additional_kits_from_lib_prep   : located at utils/additional_kits.py
        get_data_for_library_preparation_in_defined : located at iSkyLIMS_core/utils/handling_samples.py
        get_type_of_sample_information : located at iSkyLIMS_core/utils/handling_samples.py
        get_library_preparation_heading_for_samples : located at utils/library_preparation.py
        get_protocols_for_library_preparation : located at utils/library_preparation.py
        get_samples_in_lib_prep_state :  located at utils/library_preparation.py
        validate_sample_sheet_data  :  located at utils/library_preparation.py
        get_list_of_collection_kits     : located at utils/collection_index_function.py
    '''
    # get the information for returning the uploaded file in case errors in the sample sheet
    samples_in_lib_prep = get_samples_for_library_preparation()

    if request.method == 'POST' and request.POST['action'] == 'assignProtocol':
        samples_in_lib_prep_protocol = extract_protocol_library_preparation_form(request.POST)
        if len(samples_in_lib_prep_protocol) == 0 :
            return render (request, 'iSkyLIMS_wetlab/handlingLibraryPreparations.html', {'stored_lib_prep':stored_lib_prep})
        library_preparation_objs = create_library_preparation_instance(samples_in_lib_prep_protocol, request.user)
        lib_prep_protocol_parameters = get_protocol_parameters_for_library_preparation(library_preparation_objs)
        return render (request, 'iSkyLIMS_wetlab/handlingLibraryPreparations.html', {'lib_prep_protocol_parameters':lib_prep_protocol_parameters})

    # add protocol parameters for the user selected library preparation on defined state
    if request.method == 'POST' and request.POST['action'] == 'addProtocolParameter':
        lib_prep_ids = request.POST.getlist('libpreparation')
        library_preparation_objs = []
        for lib_prep_id in lib_prep_ids:
            library_preparation_objs.append(get_lib_prep_obj_from_id(lib_prep_id))
        lib_prep_protocol_parameters = get_protocol_parameters_for_library_preparation(library_preparation_objs)
        return render (request, 'iSkyLIMS_wetlab/handlingLibraryPreparations.html', {'lib_prep_protocol_parameters':lib_prep_protocol_parameters})

    # store the parameter library preparation protocol
    if request.method == 'POST' and request.POST['action'] == 'recordProtocolParamters':

        stored_params = analyze_and_store_input_param_values (request.POST)
        if 'ERROR' in stored_params:
            error_message = stored_params['ERROR']
            lib_prep_ids = request.POST['lib_prep_ids'].split(',')
            library_preparation_objs = []
            for lib_prep_id in lib_prep_ids:
                library_preparation_objs.append(get_lib_prep_obj_from_id(lib_prep_id))
            lib_prep_protocol_parameters = get_protocol_parameters_for_library_preparation(library_preparation_objs)
            # restore the user data
            lib_prep_protocol_parameters['data'] = json.loads(request.POST['protocol_data'])
            return render (request, 'iSkyLIMS_wetlab/handlingLibraryPreparations.html', {'ERROR': error_message, 'lib_prep_protocol_parameters':lib_prep_protocol_parameters})
        return render (request, 'iSkyLIMS_wetlab/handlingLibraryPreparations.html', {'stored_params':stored_params})

    if request.method == 'POST' and request.POST['action'] == 'importsamplesheet':
        data = {}
        data['full_path_file'], data['file_name'] = store_user_input_file (request.FILES['uploadfile'])
        file_read = read_user_iem_file(data['full_path_file'])
        if not valid_user_iem_file(file_read):
            # Error found when extracting data from sample sheet
            data['ERROR'] = wetlab_config.ERROR_INVALID_FILE_FORMAT
            if not delete_stored_file(data['full_path_file']):
                data['ERROR'].append(ERROR_UNABLE_TO_DELETE_USER_FILE)
            return render (request, 'iSkyLIMS_wetlab/handlingLibraryPreparations.html', {'ERROR':data['ERROR'], 'samples_in_lib_prep':samples_in_lib_prep})
        user_in_description = get_configuration_value('DESCRIPTION_IN_SAMPLE_SHEET_MUST_HAVE_USERNAME')
        if user_in_description == 'TRUE':
            user_id_in_s_sheet = extract_userids_from_sample_sheet_data(file_read)
            if 'ERROR' in user_id_in_s_sheet:
                if not delete_stored_file(data['full_path_file']):
                    user_id_in_s_sheet['ERROR'].append(ERROR_UNABLE_TO_DELETE_USER_FILE)
                return render (request, 'iSkyLIMS_wetlab/handlingLibraryPreparations.html', {'ERROR':user_id_list['ERROR'], 'samples_in_lib_prep':samples_in_lib_prep})
        else:
            user_id_in_s_sheet = []
        sample_sheet_data = get_sample_sheet_data(file_read)

        valid_data = validate_sample_sheet_data(sample_sheet_data)
        if 'ERROR' in valid_data:
            if not delete_stored_file(data['full_path_file']):
                valid_data['ERROR'].append(ERROR_UNABLE_TO_DELETE_USER_FILE)
            return render (request, 'iSkyLIMS_wetlab/handlingLibraryPreparations.html', {'ERROR':valid_data['ERROR'], 'samples_in_lib_prep':samples_in_lib_prep})

        platform = request.POST['platform']
        configuration = request.POST[request.POST['platform']]
        sample_sheet_data['file_name'] = data['file_name']
        sample_sheet_data['userid_names'] = user_id_in_s_sheet
        lib_prep_sample_sheet_obj = store_library_preparation_sample_sheet(sample_sheet_data, request.user, platform, configuration)
        #stored_lib_prep_sample = store_library_preparation_samples(sample_sheet_data,  request.user, request.POST['lib_protocols'], lib_prep_sample_sheet_obj)
        display_sample_sheet = format_sample_sheet_to_display_in_form(sample_sheet_data)
        display_sample_sheet['lib_prep_user_sample_sheet'] = lib_prep_sample_sheet_obj.get_user_sample_sheet_id()
        display_sample_sheet['platform'] = platform
        display_sample_sheet['iem_version'] = sample_sheet_data['iem_version']
        if user_in_description == 'TRUE':
            display_sample_sheet['user_list'] = get_userid_list()
        return render (request, 'iSkyLIMS_wetlab/handlingLibraryPreparations.html', {'display_sample_sheet':display_sample_sheet})


    if request.method == 'POST' and request.POST['action'] == 'storeIndexSample':
        store_data_result = store_confirmation_library_preparation_index(request.POST)
        if 'ERROR' in store_data_result:
            return render (request, 'iSkyLIMS_wetlab/handlingLibraryPreparations.html', {'ERROR':valid_data['ERROR'], 'samples_in_lib_prep':samples_in_lib_prep})
        stored_index = 'True'
        return render (request, 'iSkyLIMS_wetlab/handlingLibraryPreparations.html', {'stored_index':stored_index})

    if request.method == 'POST' and request.POST['action'] == 'libpreparationdefined':
        lib_prep_defined = request.POST.getlist('libpreparation')
        lib_protocols = get_protocol_from_library_id(lib_prep_defined[0])
        stored_lib_prep = get_library_preparation_heading_for_samples(lib_prep_defined, lib_protocols)
        return render (request, 'iSkyLIMS_wetlab/handlingLibraryPreparations.html', {'stored_lib_prep':stored_lib_prep})




        # return render (request, 'iSkyLIMS_wetlab/handlingLibraryPreparations.html', {'stored_params':stored_params})
    if request.method == 'POST' and request.POST['action'] == 'assignAdditionalKits':
        lib_prep_ids = request.POST.getlist('libpreparation')
        additional_kits = get_additional_kits_from_lib_prep(lib_prep_ids)
        return render (request, 'iSkyLIMS_wetlab/handlingLibraryPreparations.html', {'additional_kits':additional_kits})

    if request.method =='POST' and request.POST['action'] == 'storeAdditionalKits':
        stored_additional_kits = analyze_and_store_input_additional_kits (request.POST)
        if 'ERROR' in stored_additional_kits:
            error_message = stored_additional_kits['ERROR']
            lib_prep_ids = request.POST['lib_prep_ids'].split(',')
            additional_kits = get_additional_kits_from_lib_prep(lib_prep_ids)
            additional_kits['data'] = json.loads(request.POST['protocol_data'])
            return render (request, 'iSkyLIMS_wetlab/handlingLibraryPreparations.html', {'ERROR': error_message,'additional_kits':additional_kits})

        return render (request, 'iSkyLIMS_wetlab/handlingLibraryPreparations.html', {'stored_additional_kits':stored_additional_kits})
    else:

        return render (request, 'iSkyLIMS_wetlab/handlingLibraryPreparations.html', {'samples_in_lib_prep':samples_in_lib_prep})


def handling_molecules(request):
    '''
    Functions:
        get_samples_in_state : located at iSkyLIMS_core/utils/handling_samples.py
        create_table_to_select_molecules : located at iSkyLIMS_core/utils/handling_samples.py
        display_molecule_protocol_parameters  : located at iSkyLIMS_core/utils/handling_samples.py
    '''

    if request.method == 'POST' and request.POST['action'] == 'selectedMolecules':
        # If no samples are selected , call again this function to display again the sample list

        samples = get_selected_recorded_samples (request.POST['selected_samples'])
        if len(samples) == 0:
            return redirect ('handling_molecules')
        molecule_protocol = get_table_record_molecule (samples, __package__)
        if 'ERROR' in molecule_protocol :
            return render (request, 'iSkyLIMS_wetlab/error_page.html',
                {'content':['There was no valid sample selected ']})
        molecule_protocol['samples'] = ','.join(samples)

        return render(request, 'iSkyLIMS_wetlab/handlingMolecules.html',{'molecule_protocol':molecule_protocol})

    elif request.method == 'POST' and request.POST['action'] == 'updateMoleculeProtocol':
        molecule_recorded = record_molecules (request.POST, request.user, __package__)

        if 'molecule_code_ids' in request.POST and request.POST['molecule_code_ids'] != '' and 'molecule_code_ids' in molecule_recorded :
            # Add the already recorded molecules to the new ones
            molecule_recorded['molecule_code_ids'] += ',' + request.POST['molecule_code_ids']
            molecule_recorded['molecule_ids'] += ',' + request.POST['molecule_ids']
        if 'incomplete_sample_ids' in molecule_recorded:
            ## collect the information to select in the option fields
            molecule_recorded.update(get_table_record_molecule (molecule_recorded['incomplete_sample_ids'], __package__))

            return render(request, 'iSkyLIMS_wetlab/handlingMolecules.html',{'molecule_recorded':molecule_recorded})

        show_molecule_parameters = display_molecule_protocol_parameters(molecule_recorded['molecule_ids'].split(','),request.user)
        return render(request, 'iSkyLIMS_wetlab/handlingMolecules.html',{'molecule_recorded':molecule_recorded, 'show_molecule_parameters': show_molecule_parameters})
        #     added_molecule_protocol_parameters, sample_updated_list = add_molecule_protocol_parameters(request)
        #if 'pending' in request.POST :
        #molecules = request.POST['pending'].split(',')
        #show_molecule_parameters = display_molecule_protocol_parameters(molecules,request.user)
        #return render(request, 'iSkyLIMS_wetlab/setMoleculeValues.html',{'added_molecule_protocol_parameters':added_molecule_protocol_parameters, 'show_molecule_parameters':show_molecule_parameters})
        #else:
        #return render(request, 'iSkyLIMS_wetlab/setMoleculeValues.html',{'added_molecule_protocol_parameters':added_molecule_protocol_parameters})


    elif request.method == 'POST' and request.POST['action'] == 'selectedOwnerMolecules':
        # If no samples are selected , call again this function to display again the sample list
        if not 'molecules' in request.POST :
            return redirect ('handling_molecules')
        molecules = request.POST.getlist('molecules')
        # Set to true to reuse the html Code
        molecule_recorded = True
        show_molecule_parameters = display_molecule_protocol_parameters(molecules, request.user)
        return render(request, 'iSkyLIMS_wetlab/handlingMolecules.html',{'molecule_recorded':molecule_recorded, 'show_molecule_parameters': show_molecule_parameters})

    elif request.method == 'POST' and request.POST['action'] == 'addMoleculeParameters':
        molecule_parameters_updated = add_molecule_protocol_parameters(request.POST)
        if 'pending' in request.POST :
            molecules = request.POST['pending'].split(',')
            show_molecule_parameters = display_molecule_protocol_parameters(molecules,request.user)
            return render(request, 'iSkyLIMS_wetlab/handlingMolecules.html',{'molecule_parameters_updated':molecule_parameters_updated, 'show_molecule_parameters': show_molecule_parameters})
        else:
            return render(request, 'iSkyLIMS_wetlab/handlingMolecules.html',{'molecule_parameters_updated':molecule_parameters_updated})

    elif request.method == 'POST' and request.POST['action'] == 'requestMoleculeUse':
        molecule_use = set_molecule_use(request.POST, __package__)
        return render(request, 'iSkyLIMS_wetlab/handlingMolecules.html',{'molecule_use':molecule_use})

    else:
        sample_availables , user_molecules, request_molecule_use = '' ,'' ,''
        samples_list = get_samples_in_state ('Defined')
        if samples_list :
            sample_availables = create_table_to_select_molecules (samples_list)

        user_owner_molecules = get_molecule_in_state ('Defined' , request.user)
        if len(user_owner_molecules) > 0 :
            user_molecules = create_table_user_molecules (user_owner_molecules)

        samples_pending_use = get_samples_in_state ('Pending for use')
        if samples_pending_use:
            request_molecule_use = create_table_pending_use(samples_pending_use, __package__)
        # check if there are defined the type
        molecule_use_defined = check_if_molecule_use_defined(__package__)

        return render(request, 'iSkyLIMS_wetlab/handlingMolecules.html',{'sample_availables': sample_availables, 'user_molecules': user_molecules ,
                                'molecule_use_defined' : molecule_use_defined,'request_molecule_use':request_molecule_use})

    return

@login_required
def repeat_library_preparation(request):
    '''
    Functions:
    analyze_reprocess_data  : located at utils/sample_functions.py
    '''

    if  request.method == 'POST' and request.POST['action'] == 'repeat_library_preparation':
        lib_prep_id = request.POST['lib_prep_id']
        molecule_code_id = request.POST['molecule_code_id']
        sample_id = request.POST['sample_id']
        result = analyze_reprocess_data([molecule_code_id, 'New Library Preparation'], sample_id, request.user)
        detail_description = {}
        if result == 'Invalid options':
            detail_description['heading'] = ERROR_UNABLE_SAVE_REQUEST
            detail_description['information'] = [ERROR_INVALID_PARAMETERS_WHEN_REUSING_LIB_PREP]
            return render(request, 'iSkyLIMS_wetlab/error_page.html',{'detail_description':detail_description})
        detail_description['information'] = SUCCESSFUL_REUSE_MOLECULE_EXTRACTION
        return render(request, 'iSkyLIMS_wetlab/successful_page.html',{'detail_description':detail_description})
    # return to the main page because the page was not requested for the right page
    return redirect('')

@login_required
def repeat_molecule_extraction(request):
    '''
    Functions:
    analyze_reprocess_data  : located at utils/sample_functions.py
    get_table_record_molecule : located at iSkyLIMS_core/utils/handling_samples.py
    '''

    if  request.method == 'POST' and request.POST['action'] == 'repeat_extraction':
        sample_id = request.POST['sample_id']
        if analyze_reprocess_data(['New Extraction'], sample_id, request.user):
            molecule_protocol = get_table_record_molecule ([sample_id], __package__)
            molecule_protocol['samples'] = sample_id
            # create a copy of the request, to allow to modify it
            #request.POST = request.POST.copy()
            #request.POST['action'] = 'selectedMolecules'
            #request.POST['samples'] = sample_id

            return render(request, 'iSkyLIMS_wetlab/handlingMolecules.html',{'molecule_protocol':molecule_protocol})
    # return to the main page because the page was not requested for the right page
    return redirect('')


@login_required
def repeat_pool (request):
    '''
    Functions:
    analyze_reprocess_data  : located at utils/sample_functions.py
    '''
    if  request.method == 'POST' and request.POST['action'] == 'repeat_pool':

        lib_prep_obj = get_lib_prep_obj_from_id (request.POST['lib_prep_id'])
        lib_prep_code_id =  lib_prep_obj.get_lib_prep_code()
        molecule_code_id = lib_prep_obj.get_molecule_code_id()
        sample_id = lib_prep_obj.get_sample_id()

        result = analyze_reprocess_data([molecule_code_id, lib_prep_code_id, 'New Pool'], sample_id, request.user)
        detail_description = {}
        if result == 'Invalid options':
            detail_description['heading'] = ERROR_UNABLE_SAVE_REQUEST
            detail_description['information'] = [ERROR_INVALID_PARAMETERS_WHEN_REUSING_LIB_PREP]
            return render(request, 'iSkyLIMS_wetlab/error_page.html',{'detail_description':detail_description})
        detail_description['information'] = SUCCESSFUL_REUSE_LIB_PREP
        return render(request, 'iSkyLIMS_wetlab/successful_page.html',{'detail_description':detail_description})
    # return to the main page because the page was not requested for the right page
    return redirect('')

@login_required
def search_sample (request):
    '''
    Functions:
        get_sample_states  : located at iSkyLIMS_core/utils/handling_samples.py
        check_valid_date_format : located at utils/generic_functions.py
        search_samples          : located at iSkyLIMS_core/utils/handling_samples.py
        search_run_samples      : located at utils/sample_functions.py
    '''
    search_data = {}
    search_data['s_state'] = get_sample_states()

    if  request.method == 'POST' and request.POST['action'] == 'searchsample':
        sample_name=request.POST['samplename']
        start_date=request.POST['startdate']
        end_date=request.POST['enddate']
        user_name = request.POST['username']
        sample_state =request.POST['sampleState']

        # check that some values are in the request if not return the form
        if user_name == '' and start_date == '' and end_date == '' and sample_name =='' and sample_state == '':
            return render(request, 'iSkyLIMS_wetlab/searchSample.html',{'search_data':search_data})

        if user_name !=''  and len(user_name) <5 :
            return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['The user name must contains at least 5 caracters ',
                    'ADVICE:', 'write the full user name to get a better match']})
        ### check the right format of start and end date

        if start_date != '' and not check_valid_date_format(start_date):
            return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['The format for the "Start Date Search" Field is incorrect ',
                    'ADVICE:', 'Use the format  (DD-MM-YYYY)']})
        if end_date != '' and not check_valid_date_format(end_date):
            return render (request,'iSkyLIMS_wetlab/error_page.html', {'content':['The format for the "End Date Search" Field is incorrect ',
                     'ADVICE:', 'Use the format  (DD-MM-YYYY)']})
        ### Get projects when sample name is not empty
        sample_list = search_samples(sample_name, user_name, sample_state, start_date, end_date )
        if sample_state == '':
            run_sample_list = search_run_samples(sample_name, user_name, start_date, end_date)
        else:
            run_sample_list = ''

        if len(sample_list) == 0 and len(run_sample_list) == 0:
            return render (request,'iSkyLIMS_wetlab/searchSample.html', {'no_samples':ERROR_NO_SAMPLE_FOUND, 'sample_list':sample_list, 'search_data':search_data})
        elif len(sample_list) == 1 and len(run_sample_list) == 0:
            return redirect ('display_sample' , sample_id = sample_list[0])
        elif len(sample_list) == 0 and len(run_sample_list) == 1:
            #get_info_sample_in_run (run_sample_obj)
            return redirect ('display_sample_in_run' , sample_run_id = run_sample_list[0])
        else:
            # get the sample information to select it , because there are also matches on run_sample
            if len(sample_list) == 1:
                sample_obj = get_sample_obj_from_id(sample_list[0])
                sample_list = [sample_obj.get_info_for_searching()]
            if len(run_sample_list) == 1:
                run_sample_obj = get_sample_in_project_obj_from_id(run_sample_list[0])
                run_sample_list = [run_sample_obj.get_info_for_searching()]
            return render(request, 'iSkyLIMS_wetlab/searchSample.html',{'sample_list':sample_list , 'run_sample_list':run_sample_list})

    else:
        return render(request, 'iSkyLIMS_wetlab/searchSample.html',{'search_data':search_data})

@login_required
def set_molecule_values(request):
    if request.method == 'POST' and request.POST['action'] == 'continueWithMolecule':
        if request.POST['samples'] == '':
            return render (request,'iSkyLIMS_wetlab/error_page.html',
                {'content':['There was no sample selected ']})
        if  'samples_in_list' in request.POST:
            samples = request.POST.getlist('samples')
        else:
            samples = request.POST['samples'].split(',')

        molecule_protocol = get_table_record_molecule (samples, __package__)
        if 'ERROR' in molecule_protocol :
            return render (request, 'iSkyLIMS_wetlab/error_page.html',
                {'content':['There was no valid sample selected ']})

        molecule_protocol['samples'] = ','.join(samples)

        return render(request, 'iSkyLIMS_wetlab/setMoleculeValues.html',{'molecule_protocol':molecule_protocol})

    elif request.method == 'POST' and request.POST['action'] == 'updateMoleculeProtocol':

        molecule_recorded = record_molecules (request)

        if not 'heading' in molecule_recorded:
            samples = request.POST['samples'].split(',')
            molecule_protocol = get_table_record_molecule (samples, __package__)
            molecule_protocol['data'] = molecule_recorded['incomplete_molecules']
            molecule_protocol['samples'] = ','.join(samples)
            return render(request, 'iSkyLIMS_wetlab/setMoleculeValues.html',{'molecule_protocol':molecule_protocol})
        else:
            if 'incomplete_molecules' in molecule_recorded:
                samples = molecule_recorded['incomplete_molecules_ids'].split(',')
                molecule_recorded.update(get_table_record_molecule (samples))

            return render(request, 'iSkyLIMS_wetlab/setMoleculeValues.html',{'molecule_recorded':molecule_recorded})

    elif request.method == 'POST' and request.POST['action'] == 'displayMoleculeParameters':
        if  'samples_in_list' in request.POST:
            molecules = request.POST.getlist('molecules')
        else:
            molecules = request.POST['molecules'].split(',')
        show_molecule_parameters = display_molecule_protocol_parameters(molecules, request.user)
        return render(request, 'iSkyLIMS_wetlab/setMoleculeValues.html',{'show_molecule_parameters':show_molecule_parameters})

    elif request.method == 'POST' and request.POST['action'] == 'addMoleculeParameters':
        added_molecule_protocol_parameters, sample_updated_list = add_molecule_protocol_parameters(request)
        if 'pending' in request.POST :
            molecules = request.POST['pending'].split(',')
            show_molecule_parameters = display_molecule_protocol_parameters(molecules,request.user)
            return render(request, 'iSkyLIMS_wetlab/setMoleculeValues.html',{'added_molecule_protocol_parameters':added_molecule_protocol_parameters, 'show_molecule_parameters':show_molecule_parameters})
        else:
            return render(request, 'iSkyLIMS_wetlab/setMoleculeValues.html',{'added_molecule_protocol_parameters':added_molecule_protocol_parameters})

    else:
        register_user = request.user.username
        display_list = get_defined_samples (register_user)

        return render(request, 'iSkyLIMS_wetlab/setMoleculeValues.html',{'display_list': display_list})
    return render(request, 'iSkyLIMS_wetlab/setMoleculeValues.html',{})

'''   Replaced by handling_library_preparations
@login_required
def set_library_preparation(request):

    if request.method == 'POST' and request.POST['action'] == 'importsamplesheet':
        #protocol = request.POST['lib_protocols']
        #single_paired = request.POST['singlePairedEnd']
        #read_length = request.POST['readlength']
        extension_file = '.csv'
        reg_user = request.user.username
        stored_file , file_name = store_user_input_file(request.FILES['importsamplesheet'], extension_file)
        samples_in_s_sheet = get_samples_in_sample_sheet(stored_file)
        not_defined_samples = []
        not_allowed_lib_prep = []
        # Check if samples inside sample sheet file are valid to add the library preparation data
        for sample in samples_in_s_sheet['sample_names']:
            sample_obj = get_sample_instance (sample, reg_user)
            if not sample_obj :
                not_defined_samples.append(sample)
                continue

            if LibraryPreparation.objects.filter(sample_id = sample_obj , libPrepState__libPrepState__exact = 'Defined').exists():
                not_allowed_lib_prep.append(sample)

        if (not_defined_samples or not_allowed_lib_prep):
            not_defined_samples_str, not_allowed_lib_prep_str = '' , ''
            if not_defined_samples:
                not_defined_samples_str = 'Samples not defined: ' + ', '.join(not_defined_samples)
            if not_allowed_lib_prep:
                not_allowed_lib_prep_str = 'Samples in wrong state : ' + ', '.join(not_allowed_lib_prep)
            return render ( request,'iSkyLIMS_wetlab/error_page.html',
                {'content':[not_defined_samples_str , not_allowed_lib_prep_str ]})



        #index_adapters = get_indexes_adapters (stored_file)
        if not CollectionIndexKit.objects.filter(collectionIndexName__iexact = index_adapters).exists():
            return render ( request,'iSkyLIMS_wetlab/error_page.html',
                {'content':['Collection Index Kit ' , index_adapters , 'is not defined' ]})

        library_prep_workflow = get_index_library_name(stored_file)
        adapter1, adapter2 = get_adapters(stored_file)
        assay = get_assay_from_file (stored_file)
        # store user sample sheet in database
        #user_sample_sheet_data = {}
        #stored_lib_prep = {}
        #stored_lib_prep['data'] = []


        extracted_data_list = extract_sample_data (samples_in_s_sheet)
        # Check if exists duplicate index in sample sheet
        duplicate_index = find_duplicate_index(extracted_data_list)
        if duplicate_index:
            detail_description = {}
            detail_description ['heading'] = ['Index', 'Samples']
            detail_description['information'] =  duplicate_index
            return render ( request,'iSkyLIMS_wetlab/error_page.html',
                {'content':['Found some samples, has the same index' ],
                'detail_description': detail_description })


        stored_lib_prep = prepare_lib_prep_table_new_run (index_adapters, request, extracted_data_list,
                                    file_name, assay, adapter1, adapter2)

        return render (request, 'iSkyLIMS_wetlab/setLibraryPreparation.html', {'stored_lib_prep':stored_lib_prep})

    elif request.method == 'POST' and request.POST['action'] == 'addLibPrepParam':

        if  'lib_prep_in_list' in request.POST:
            lib_prep_ids = request.POST.getlist('lib_prep_id')
            if len('lib_prep_in_list') == 1:
                lib_prep_ids = list(request.POST['lib_prep_id'])
        else:
            lib_prep_ids = request.POST['lib_prep_id'].split(',')
        stored_lib_prep ={}
        stored_lib_prep['data'] =[]

        protocol_obj = LibraryPreparation.objects.get(pk = lib_prep_ids[0]).get_protocol_obj()
        parameter_heading = get_protocol_parameters(protocol_obj)
        length_heading = len(HEADING_FIX_FOR_ADDING_LIB_PARAMETERS) + len (parameter_heading)
        for lib_id in lib_prep_ids:
            library_preparation_obj = LibraryPreparation.objects.get(pk = lib_id)
            data = ['']*length_heading
            data[0] = library_preparation_obj.get_sample_name()
            data[1] = library_preparation_obj.get_lib_prep_code()
            #data[2] = library_preparation_obj.get_collection_index_kit()

            stored_lib_prep['data'].append(data)
        stored_lib_prep['lib_prep_id'] = ','.join(lib_prep_ids)
        stored_lib_prep['heading'] = HEADING_FIX_FOR_ADDING_LIB_PARAMETERS
        stored_lib_prep['par_heading'] = parameter_heading
        stored_lib_prep['heading_in_excel'] = ','.join(HEADING_FIX_FOR_ADDING_LIB_PARAMETERS + parameter_heading)
        register_user_obj = User.objects.get(username__exact = request.user.username)
        stored_lib_prep['reagents_kits'] = get_lot_commercial_kits(register_user_obj, protocol_obj)

        return render (request, 'iSkyLIMS_wetlab/setLibraryPreparation.html', {'stored_lib_prep':stored_lib_prep})

    elif request.method == 'POST' and request.POST['action'] == 'recordProtocolParamters':

        stored_params = analyze_input_param_values (request)

        return render (request, 'iSkyLIMS_wetlab/setLibraryPreparation.html', {'stored_params':stored_params})

    elif request.method == 'POST' and request.POST['action'] == 'addLibraryProtocol':
        if  'molecules_in_list' in request.POST:
            molecules = request.POST.getlist('molecules')
        else:
            molecules = request.POST['molecules'].split(',')
        prot_lib_json_data = json.loads(request.POST['protocol_data'])
        heading_in_excel = request.POST['heading_in_excel'].split(',')


        for i in range(len(molecules)) :
            valid_molecules = []
            protocol_name = prot_lib_json_data[i][heading_in_excel.index('Protocol used')]
            if protocol_name == '':
                continue
            if not Protocols.objects.filter(name__exact = protocol_name).exists():
                continue
            protocol_obj = Protocols.objects.get(name__exact = protocol_name)
            if not ProtocolLibraryParameters.objects.filter(protocol_id = protocol_obj).exists():
                continue
            valid_molecules.append(molecules[i])
            prot_lib_parameters = ProtocolLibraryParameters.objects.filter(protocol_id = protocol_obj)

'''

@login_required
def set_library_values (request):
    fix_headings = ['DNA Code ID', 'Protocol', 'Extraction Kit']

    if request.method == 'POST' and request.POST['action'] == 'continueWithDNA':
        lib_preparation_data = {}
    else:
        register_user = request.user.username
        grouped_samples_obj = get_samples_for_library_definition (register_user)
        s_list = []
        all_sample_list = []
        display_list = {}

        for key in grouped_samples_obj.keys():
            s_list = []
            for sample in grouped_samples_obj[key]:
                s_info = sample.get_sample_definition_information().split(';')
                s_info.append(str(sample.pk))
                s_list.append(s_info)
            all_sample_list.append(s_list)

        display_list['lib_kits'] = get_available_lib_kit(register_user)

        display_list['list_of_samples'] = all_sample_list
        display_list['heading'] = ['Registered date ','Sample Code ID', 'Type', 'DNA/RNA', 'Protocol', 'Library Kit']

        return render(request, 'iSkyLIMS_wetlab/setLibraryValues.html',{'display_list': display_list})


@login_required
def create_pool (request):
    ## Check user == WETLAB_MANAGER: if false,  redirect to 'login' page
    if request.user.is_authenticated:
        if not is_wetlab_manager(request):
            return render ( request,'iSkyLIMS_wetlab/error_page.html',
                {'content':['You do not have enough privileges to see this page ',
                            'Contact with your administrator .']})
    else:
        #redirect to login webpage
        return redirect ('/accounts/login')
    # collect the information for collecting
    display_list = get_lib_prep_to_select_in_pool()

    if request.method == 'POST' and request.POST['action'] == 'createPool':
        new_pool = define_new_pool(request.POST,  request.user)

        if not isinstance(new_pool, LibraryPool) :
            display_list.update(new_pool)
            return  render(request, 'iSkyLIMS_wetlab/createPool.html',{'display_list': display_list})
        information_for_created_pool = get_info_to_display_created_pool(new_pool )
        return  render(request, 'iSkyLIMS_wetlab/createPool.html',{'information_for_created_pool': information_for_created_pool})

    else:
        return  render(request, 'iSkyLIMS_wetlab/createPool.html',{'display_list': display_list})


@login_required
def create_new_run (request):
    if request.user.is_authenticated:
        if not is_wetlab_manager(request):
            return render ( request,'iSkyLIMS_wetlab/error_page.html',
                {'content':['You do not have enough privileges to see this page ',
                            'Contact with your administrator .']})
    else:
        #redirect to login webpage
        return redirect ('/accounts/login')

    if request.method == 'POST' and request.POST['action'] == 'createNewRun':
        display_pools_for_run = display_available_pools()
        if not 'poolID' in request.POST :
            error_message = wetlab_config.ERROR_NO_POOL_WAS_SELECTED_IN_FORM
            return  render(request, 'iSkyLIMS_wetlab/CreateNewRun.html',{'display_pools_for_run': display_pools_for_run, 'ERROR': error_message})
        compatibility = check_valid_data_for_creation_run (request.POST, request.user)
        if 'ERROR' in compatibility :
            return render (request,'iSkyLIMS_wetlab/CreateNewRun.html',{'display_pools_for_run':display_pools_for_run, 'ERROR': compatibility['ERROR']})
        #compatible_index = check_index_compatible(lib_prep_ids)
        display_sample_information = create_run_in_pre_recorded_and_get_data_for_confirmation(request.POST, request.user)

        return  render(request, 'iSkyLIMS_wetlab/CreateNewRun.html',{'display_sample_information': display_sample_information})

    elif request.method == 'POST' and request.POST['action'] == 'continueWithRun':
        run_id = request.POST['run_ids']
        experiment_name = get_experiment_name(run_id)
        pool_objs = LibraryPool.objects.filter(runProcess_id__exact = run_id)
        pool_ids = []
        for pool in pool_objs :
            pool_ids.append(pool.get_id())
        lib_prep_ids = get_library_prep_in_pools (pool_ids)

        display_sample_information = get_library_preparation_data_in_run(lib_prep_ids, pool_ids)
        display_sample_information.update(get_stored_user_sample_sheet(lib_prep_ids))
        display_sample_information['experiment_name'] = experiment_name
        display_sample_information['run_process_id'] = run_id
        return  render(request, 'iSkyLIMS_wetlab/CreateNewRun.html',{'display_sample_information': display_sample_information})

    elif request.method == 'POST' and request.POST['action'] ==  'storeDataNewRun':
        run_obj = get_run_obj_from_id(request.POST['run_process_id'])
        if run_obj.get_state() != 'Pre-Recorded':
            exp_name = run_obj.get_run_name()
            error_message = ERROR_RUN_NAME_CREATED_ALREADY.copy()
            error_message.insert(1,exp_name)
            display_pools_for_run = display_available_pools()
            return  render(request, 'iSkyLIMS_wetlab/CreateNewRun.html',{'display_pools_for_run': display_pools_for_run, 'ERROR': error_message})
        run_data = collect_data_and_update_library_preparation_samples_for_run(request.POST, request.user)

        projects_objs = create_new_projects_added_to_run(run_data['projects'], run_data['run_obj'],request.user)
        if 'ERROR' in projects_objs:
            display_pools_for_run = display_available_pools()
            return  render(request, 'iSkyLIMS_wetlab/CreateNewRun.html',{'display_pools_for_run': display_pools_for_run, 'ERROR': projects_objs['ERROR']})
        ## store data in runProcess table, run is in pre-recorded state
        #center_requested_id = Profile.objects.get(profileUserID = request.user).profileCenter.id
        #center_requested_by = Center.objects.get(pk = center_requested_id)

        #update_run_with_sample_sheet(request.POST['run_process_id'], run_data['sample_sheet'])

        #run_obj = run_data['run_obj']

        run_obj.set_run_state('Recorded')

        sample_sheet_name = store_confirmation_sample_sheet(run_data)
        # update the sample state for each one in the run
        pools_obj = LibraryPool.objects.filter(runProcess_id = run_obj)


        for pool_obj in pools_obj:
            pool_obj.set_pool_state('Used')
        update_batch_lib_prep_sample_state(run_data['lib_prep_ids'],  'Sequencing')
        created_new_run = {}
        created_new_run['exp_name'] = run_data['exp_name']
        created_new_run['run_process_id'] = request.POST['run_process_id']
        created_new_run['sample_sheet'] = sample_sheet_name

        return  render(request, 'iSkyLIMS_wetlab/CreateNewRun.html',{'created_new_run': created_new_run})
    else:
        display_pools_for_run = display_available_pools()
        return  render(request, 'iSkyLIMS_wetlab/CreateNewRun.html',{'display_pools_for_run': display_pools_for_run})


@login_required
def pending_sample_preparations(request):
    pending = {}
    # get the samples in defined state
    pending['defined'] = get_samples_in_defined_state('')
    pending['extract_molecule'] = get_samples_in_extracted_molecule_state('')
    pending['create_library_preparation'] = get_samples_in_lib_prep_state()
    #pending['lib_prep_protocols'] = get_protocol_lib()
    # get the library preparation in defined state
    pending['add_lib_prep_parameters'] = get_lib_prep_to_add_parameters()
    pending ['graphic_pending_samples'] = pending_samples_for_grafic(pending).render()
    return render(request, 'iSkyLIMS_wetlab/pendingSamplePreparations.html',{'pending': pending})

@login_required
def compare_samples(request):
    user_is_wetlab_manager = is_wetlab_manager(request)
    samples_data = get_list_of_samples_in_projects(request.user, user_is_wetlab_manager)
    samples_data['user'] = request.user.username
    if request.method == 'POST' and request.POST['action'] == 'compareSamples':
        selected_sample_objs = analyze_compare_samples_form(request.POST['table_data'])
        if len(selected_sample_objs) == 0:
            error_message = ERROR_NO_SAMPLES_SELECTED
            return render(request, 'iSkyLIMS_wetlab/compareSamples.html',{'ERROR':error_message, 'samples_data': samples_data})
        compared_data = get_comparation_sample_information(selected_sample_objs)

        return render(request, 'iSkyLIMS_wetlab/compareSamples.html',{'compared_data': compared_data})
    else:
        return render(request, 'iSkyLIMS_wetlab/compareSamples.html',{'samples_data': samples_data})

@login_required
def user_commercial_kit_inventory(request):
    expired_kit = get_expired_lot_user_kit(request.user)
    valid_kit = get_valid_lot_user_kit(request.user)
    if request.method == 'POST' and request.POST['action'] == 'runOutUserLotKit':
        selected_user_kits = request.POST.getlist('userKit')
        if len(selected_user_kits) == 0:
            return render(request, 'iSkyLIMS_wetlab/userCommercialKitInventory.html',{'expired_kit': expired_kit,
                                    'valid_kit': valid_kit, 'user_name': request.user.username})
        run_out_kits = set_user_lot_kit_to_run_out(selected_user_kits)
        return render(request, 'iSkyLIMS_wetlab/userCommercialKitInventory.html',{'run_out_kits': run_out_kits})

    else:
        return render(request, 'iSkyLIMS_wetlab/userCommercialKitInventory.html',{'expired_kit': expired_kit,
                                'valid_kit': valid_kit, 'user_name': request.user.username})


@login_required
def search_user_lot_kit(request):
    protocol_list = display_protocol_list()
    platform_list = get_defined_platforms_and_ids('NGS')
    if request.method == 'POST' and request.POST['action'] ==  'searchuserkit':
        if (request.POST['expired'] == '' and  request.POST['lotNumber'] == '' and request.POST['commercial'] == '' and request.POST['protocol'] == ''  and request.POST['platform'] == '') and 'exclude_runout' not in request.POST:
            return render(request, 'iSkyLIMS_wetlab/searchUserLotKit.html',{'protocol_list': protocol_list, 'platform_list' :platform_list, })

        if request.POST['expired'] != '' and not check_valid_date_format(request.POST['expired']):
            error_message = ERROR_INVALID_FORMAT_FOR_DATES
            return render (request,'iSkyLIMS_wetlab/searchUserLotKit.html', {'protocol_list': protocol_list, 'platform_list' :platform_list, 'ERROR': error_message})

        user_kits_objs = search_user_lot_kit_from_user_form(request.POST)
        if user_kits_objs == 'No defined':
            error_message = ERROR_NO_USER_LOT_KIT_DEFINED
            return render (request,'iSkyLIMS_wetlab/searchUserLotKit.html', {'protocol_list': protocol_list, 'platform_list' :platform_list, 'ERROR': error_message})
        if len(user_kits_objs) > 1:
            display_user_kit_list = display_user_lot_kit_information_from_query_list(user_kits_objs)
            return render(request, 'iSkyLIMS_wetlab/searchUserLotKit.html',{'display_user_kit_list': display_user_kit_list})
        elif len(user_kits_objs) == 0 :
            error_message = ERROR_NO_MATCHES_FOR_USER_LOT_KIT
            return render (request,'iSkyLIMS_wetlab/searchUserLotKit.html', {'protocol_list': protocol_list, 'platform_list' :platform_list, 'ERROR': error_message})
        else:
            display_one_user_kit = user_kits_objs[0].get_user_lot_kit_id()
        return redirect ('display_user_lot_kit', user_kit_id = display_one_user_kit)
    else:
        return render(request, 'iSkyLIMS_wetlab/searchUserLotKit.html',{'protocol_list': protocol_list , 'platform_list' :platform_list})

@login_required
def display_user_lot_kit(request, user_kit_id):
    user_kit_obj = get_user_lot_commercial_kit_obj_from_id(user_kit_id)
    if user_kit_obj == None:
        return render(request, 'iSkyLIMS_wetlab/error_page.html',{'content': ['Invalid User Lot Commercial Kit']})
    user_lot_kit_data = get_user_lot_kit_data_to_display(user_kit_obj)
    return render(request, 'iSkyLIMS_wetlab/displayUserLotKit.html',{'user_lot_kit_data': user_lot_kit_data})


@login_required
def sequencer_configuration(request):
    if not request.user.is_authenticated:
        #redirect to login webpage
        return redirect ('/accounts/login')

    if not is_wetlab_manager(request):
        return render (request,'iSkyLIMS_wetlab/error_page.html', {'content': ERROR_USER_NOT_WETLAB_MANAGER})
    sequencer_info = get_list_sequencer_configuration()
    sequencer_info['platforms'] = get_platform_data()
    sequencer_info['sequencer_names']= get_defined_sequencers()

    if request.method =='POST' and request.POST['action'] == 'addNewSequencer':
        new_sequencer = define_new_sequencer(request.POST)
        if 'ERROR' in new_sequencer :
            return render(request, 'iSkyLIMS_wetlab/sequencerConfiguration.html', {'sequencer_info': sequencer_info, 'ERROR': new_sequencer})
        return render(request, 'iSkyLIMS_wetlab/sequencerConfiguration.html', {'sequencer_info': sequencer_info, 'new_defined_sequencer': new_sequencer})
    if request.method =='POST' and request.POST['action'] == 'addNewConfiguration':
        new_defined_configuration = define_new_seq_configuration(request.POST)
        if 'ERROR' in new_defined_configuration :
            return render(request, 'iSkyLIMS_wetlab/sequencerConfiguration.html', {'sequencer_info': sequencer_info, 'ERROR': new_defined_configuration})
        return render(request, 'iSkyLIMS_wetlab/sequencerConfiguration.html', {'sequencer_info': sequencer_info, 'new_defined_configuration': new_defined_configuration})
    else:
        return render(request, 'iSkyLIMS_wetlab/sequencerConfiguration.html', {'sequencer_info':sequencer_info} )

@login_required
def sequencer_inventory(request):
    if not request.user.is_authenticated:
        #redirect to login webpage
        return redirect ('/accounts/login')
    if request.method =='POST' and request.POST['action'] == 'setRunOutDate':
        pass
    else:
        sequencer_data = get_sequencer_inventory_data()
        return render(request, 'iSkyLIMS_wetlab/sequencerInventory.html', {'sequencer_data':sequencer_data} )

    return
