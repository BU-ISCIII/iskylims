# Generic imports
import datetime
import json
import os
import re
import statistics
import time

from django.conf import settings
from django.contrib.auth.decorators import login_required
from django.contrib.auth.models import User
from django.core.files.storage import FileSystemStorage
from django.shortcuts import render, redirect
from django_utils.models import Center, Profile
from django_utils.views import check_user_group
from core.utils.common import (
    get_email_data,
    get_inital_sample_settings_values,
    save_inital_sample_setting_value,
    send_test_email,
)

# Local imports
import core.utils.load_batch
import core.utils.platforms
import core.utils.protocols
import core.utils.samples
import core.fusioncharts.fusioncharts
import wetlab.config
import wetlab.fusioncharts.fusioncharts
from .utils.collection_index import *
# updated --
import wetlab.utils.common
import logging
# end updated --

from .utils.fetch_info import *
from .utils.library import *
from .utils.pool import *
from .utils.run import *
from .utils.sample import *
from .utils.samplesheet import *
from .utils.sequencers import *
from .utils.statistics import *
from .utils.stats_graphs import *
from .utils.test_conf import *
import wetlab.models

import wetlab.utils.additional_kits
import wetlab.utils.collection_index



def index(request):
    org_name = wetlab.utils.common.get_configuration_from_database("ORGANIZATION_NAME")
    return render(
        request, "wetlab/index.html", {"organization_name": org_name}
    )


@login_required
def configuration_email(request):
    if request.user.username != "admin":
        return redirect("/wetlab")
    email_conf_data = get_email_data()
    email_conf_data["EMAIL_ISKYLIMS"] = wetlab.utils.common.get_configuration_from_database(
        "EMAIL_FOR_NOTIFICATIONS"
    )
    if request.method == "POST" and (request.POST["action"] == "emailconfiguration"):
        result_email = send_test_email(request.POST)
        if result_email != "OK":
            email_conf_data = wetlab.utils.common.get_email_data()
            email_conf_data["EMAIL_ISKYLIMS"] = request.POST["EMAIL_ISKYLIMS"]
            email_conf_data["test_email"] = request.POST["test_email"]
            return render(
                request,
                "wetlab/configuration_email.html",
                {"ERROR": result_email, "email_conf_data": email_conf_data},
            )
        wetlab.utils.common.save_database_configuration_value(
            "EMAIL_FOR_NOTIFICATIONS", request.POST["EMAIL_ISKYLIMS"]
        )
        return render(
            request,
            "wetlab/configuration_email.html",
            {"succesful_settings": True},
        )
    return render(
        request,
        "wetlab/configuration_email.html",
        {"email_conf_data": email_conf_data},
    )


@login_required
def configuration_samba(request):
    if request.user.username != "admin":
        return redirect("/wetlab")
    samba_conf_data = wetlab.utils.common.get_samba_connection_data()
    if request.method == "POST" and (request.POST["action"] == "sambaconfiguration"):
        # reload configuration samba settings
        samba_user_field = {}
        for field in SAMBA_CONFIGURATION_FIELDS:
            # set values for switch fields
            try:
                if request.POST[field] == "on":
                    samba_user_field[field] = "True"
                else:
                    samba_user_field[field] = request.POST[field]
            except KeyError:
                samba_user_field[field] = "False"
        wetlab.utils.common.save_samba_connection_data(samba_user_field)
        try:
            wetlab.utils.common.open_samba_connection()
            return render(
                request,
                "wetlab/configuration_samba.html",
                {"succesful_settings": True},
            )
        except Exception:
            error_message = wetlab.config.ERROR_WRONG_SAMBA_CONFIGURATION_SETTINGS
            return render(
                request,
                "wetlab/configuration_samba.html",
                {"samba_conf_data": samba_user_field, "error_message": error_message},
            )
    else:
        return render(
            request,
            "wetlab/configuration_samba.html",
            {"samba_conf_data": samba_conf_data},
        )


@login_required
def configuration_test(request):
    # check user privileges
    if request.user.is_authenticated:
        if request.user.username != "admin":
            return redirect("/wetlab")
    if request.method == "POST" and request.POST["action"] == "basicTest":
        test_results = {}
        config_file = os.path.join(
            settings.BASE_DIR, "wetlab", "config.py"
        )
        test_results["iSkyLIMS_settings"] = get_iSkyLIMS_settings()
        test_results["config_file"] = get_config_file(config_file)
        test_results["attr_files"] = get_files_attribute(
            os.path.join(settings.MEDIA_ROOT, "wetlab")
        )
        test_results["database_access"] = check_access_database()
        test_results["samba_connection"] = check_samba_connection()

        test_results["basic_checks_ok"] = "OK"
        for result in test_results:
            if test_results[result] == "NOK":
                test_results["basic_checks_ok"] = "NOK"
                break
        available_run_test = []
        if RunConfigurationTest.objects.all().exists():
            run_test_objs = RunConfigurationTest.objects.all()
            for run_test_obj in run_test_objs:
                available_run_test.append(
                    [run_test_obj.get_run_test_name(), run_test_obj.get_run_test_id()]
                )
        return render(
            request,
            "wetlab/configuration_test.html",
            {"test_results": test_results, "available_run_test": available_run_test},
        )

    elif request.method == "POST" and request.POST["action"] == "executeRunTest":
        if RunConfigurationTest.objects.filter(
            pk__exact=request.POST["runTest"]
        ).exists():
            run_test_obj = RunConfigurationTest.objects.filter(
                pk__exact=request.POST["runTest"]
            ).last()
            run_test_folder = run_test_obj.get_run_test_folder()
            run_test_name = run_test_obj.get_run_test_name()
        if not folder_test_exists(run_test_folder):
            return render(
                request,
                "wetlab/ConfigurationTest.html",
                {"error": config.ERROR_NOT_FOLDER_RUN_TEST_WAS_FOUND},
            )
        run_test_result = execute_test_for_testing_run(run_test_name)
        run_test_result["run_test_name"] = run_test_name
        if "ERROR" in run_test_result:
            log_trace = []
            with open(
                logging.getLoggerClass().root.handlers[0].baseFilename, "r"
            ) as fh:
                for line in fh:
                    if run_test_name in line:
                        line = line.replace("\n", "")
                        log_trace.append(line)

            return render(
                request,
                "wetlab/configuration_test.html",
                {"run_test_result": run_test_result, "log_trace": log_trace},
            )
        else:
            return render(
                request,
                "wetlab/configuration_test.html",
                {"run_test_result": run_test_result},
            )
    elif request.method == "POST" and request.POST["action"] == "deleteTestRun":
        if RunConfigurationTest.objects.filter(
            run_test_name__exact=request.POST["deleteRun"]
        ).exists():
            if RunProcess.objects.filter(
                run_name__exact=request.POST["deleteRun"]
            ).exists():
                run_test_objs = RunProcess.objects.filter(
                    run_name__exact=request.POST["deleteRun"]
                )
                for run_test_obj in run_test_objs:
                    delete_test_run(run_test_obj)
                    delete_successful = {"run_name": request.POST["deleteRun"]}
                return render(
                    request,
                    "wetlab/ConfigurationTest.html",
                    {"delete_successful": delete_successful},
                )
            return render(request, "wetlab/configuration_test.html")
    else:
        return render(request, "wetlab/configuration_test.html")



@login_required
def initial_settings(request):
    if not request.user.is_authenticated:
        # redirect to login webpage
        return redirect("/accounts/login")

    if not wetlab.utils.common.is_wetlab_manager(request):
        return render(
            request,
            "wetlab/error_page.html",
            {"content": wetlab.config.ERROR_USER_NOT_WETLAB_MANAGER},
        )
    initial_data = get_inital_sample_settings_values(__package__)
    form_data = {}
    if request.method == "POST" and request.POST["action"] == "defineNewSpecie":
        form_data["species"] = request.POST["specieName"]
    if request.method == "POST" and request.POST["action"] == "defineNewState":
        form_data["state"] = request.POST["stateName"]
    if request.method == "POST" and request.POST["action"] == "defineNewCity":
        form_data["city"] = request.POST
    if request.method == "POST" and request.POST["action"] == "defineNewLabRequest":
        form_data["lab_request"] = request.POST
    if request.method == "POST" and request.POST["action"] == "defineMoleculeType":
        form_data["molecule_type"] = request.POST["moleculeName"]
    if request.method == "POST" and request.POST["action"] == "defineProtocolType":
        form_data["protocol_type"] = [
            request.POST["protocolName"],
            request.POST["moleculeType"],
        ]

    if form_data:
        new_inital_data = save_inital_sample_setting_value(__package__, form_data)
        if "ERROR" in new_inital_data:
            return render(
                request,
                "wetlab/initial_settings.html",
                {"initial_data": initial_data, "ERROR": new_inital_data["ERROR"]},
            )
        return render(
            request,
            "wetlab/initial_settings.html",
            {"initial_data": initial_data, "new_setting_defined": new_inital_data},
        )
    else:
        return render(
            request,
            "wetlab/initial_settings.html",
            {"initial_data": initial_data},
        )


@login_required
def create_nextseq_run(request):
    # Check user == WETLAB_MANAGER: if false,  redirect to 'login' page
    if not wetlab.utils.common.is_wetlab_manager(request):
        return render(
            request,
            "wetlab/error_page.html",
            {"content": wetlab.config.ERROR_USER_NOT_WETLAB_MANAGER},
        )

    # FIRST STEP in collecting data from the NextSeq run. Sample Sheet and experiment name are required
    if request.method == "POST" and (request.POST["action"] == "uploadFile"):
        projects = []
        # run_name=request.POST['runname']
        myfile = request.FILES["myfile"]

        # CHECK if file contains the extension.
        # Error page is showed if file does not contain any extension
        split_filename = re.search("(.*)(\.\w+$)", myfile.name)
        if split_filename is None:
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "Uploaded file does not containt extension",
                        "Sample Sheet must have a csv extension",
                        "",
                        "ADVICE:",
                        "Select the Sample file generated by Illumina Experient Manager (IEM)",
                    ]
                },
            )
        ext_file = split_filename.group(2)

        # CHECK if file contains the csv extension.
        # Error page is shown if file does not contain the csv extension
        if ext_file != ".csv":
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "Sample Sheet must have a csv extension",
                        "",
                        "ADVICE:",
                        "Select the Sample file generated by Illumina Experient Manager (IEM)",
                    ]
                },
            )
        fs = FileSystemStorage()
        timestr = time.strftime("%Y%m%d-%H%M%S")
        # including the timestamp to the sample sheet file
        # do not need to include the absolute path because django uses
        # the MEDIA_ROOT variable defined on settings to upload the file
        file_name = str(
            wetlab.config.RUN_SAMPLE_SHEET_DIRECTORY
            + split_filename.group(1)
            + "_"
            + timestr
            + ext_file
        )
        fs.save(file_name, myfile)

        # add the document directory to read the csv file
        stored_file = os.path.join(settings.MEDIA_ROOT, file_name)

        # Fetch the experiment name and the library name from the sample sheet file
        index_library_name = get_index_library_name(stored_file)
        run_name = wetlab.utils.common.get_experiment_name_from_file(stored_file)

        if run_name == "":
            # define an temporary unique value for the run name
            # until the real value is get from user FORM
            run_name = timestr

        # Check that runName is not already used in the database.
        # Error page is showed if runName is already  defined

        if (RunProcess.objects.filter(run_name__iexact=run_name)).exists():
            if RunProcess.objects.filter(
                run_name__iexact=run_name, state__run_state_name__exact="Pre-Recorded"
            ).exists():
                # Delete the Sample Sheet file and the row in database
                delete_run_objs = RunProcess.objects.filter(
                    run_name__iexact=run_name, state__run_state_name__exact="Pre-Recorded"
                )
                for delete_run in delete_run_objs:
                    # sample_sheet_file = delete_run.get_sample_file()
                    # full_path_sample_sheet_file = os.path.join(settings.MEDIA_ROOT, sample_sheet_file)
                    # os.remove(full_path_sample_sheet_file)

                    if Projects.objects.filter(run_process=delete_run).exists():
                        project_objs = Projects.objects.filter(run_process=delete_run)
                        for project_obj in project_objs:
                            project_obj.run_process.remove(delete_run)

                            if project_obj.run_process.all().count() == 0:
                                project_obj.delete()
                    delete_run.delete()

            else:
                # delete sample sheet file
                os.remove(stored_file)
                return render(
                    request,
                    "wetlab/error_page.html",
                    {
                        "content": [
                            "Run Name is already used. ",
                            "Run Name must be unique in database.",
                            " ",
                            "ADVICE:",
                            "Change the value in the Sample Sheet  file ",
                        ]
                    },
                )

        # Fetch from the Sample Sheet file the projects included in
        # the run and the user. Error page is showed if not project/description
        # colunms are found

        project_list = get_projects_in_run(stored_file)

        if len(project_list) == 0:
            # delete sample sheet file
            fs.delete(file_name)
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        'Sample Sheet does not contain "Sample project" and/or "Description" fields',
                        "",
                        "ADVICE:",
                        "Check that csv file generated by Illumina Experient Manager (IEM) includes these columns",
                    ]
                },
            )

        # Check if the projects are already defined on database.
        # Error page is showed if projects are already defined on database

        project_already_defined = []
        project_in_several_runs = wetlab.utils.common.get_configuration_value(
            "PROJECTS_ALLOWED_IN_MULTIPLE_RUNS"
        )
        for key in project_list.keys():
            # check if project was already saved in database in Not Started State.
            # if found delete the projects, because the previous attempt to complete the run was unsuccessful
            if Projects.objects.filter(project_name__icontains=key).exists():
                   project_already_defined.append(key)

        if len(project_already_defined) > 0 and project_in_several_runs != "TRUE":
            if len(project_already_defined) > 1:
                head_text = "The following projects are already defined in database:"
            else:
                head_text = "The following project is already defined in database:"
            # convert the list into string to display the user names on error page
            display_project = "  ".join(project_already_defined)
            # delete sample sheet file before showing the error page
            fs.delete(file_name)
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        head_text,
                        "",
                        display_project,
                        "",
                        "Project names must be unique",
                        "",
                        "ADVICE:",
                        "Edit the Sample Sheet file to correct this error",
                    ]
                },
            )

        # Once the information looks good. it will be stores in runProcess and projects table

        # store data in runProcess table, run is in pre-recorded state
        center_requested_id = Profile.objects.get(
            profile_user_id=request.user
        ).profile_center.id
        center_requested_by = Center.objects.get(pk=center_requested_id)
        new_run_obj = RunProcess(
            run_name=run_name,
            sample_sheet=file_name,
            state=RunStates.objects.get(run_state_name__exact="Pre-Recorded"),
            center_requested_by=center_requested_by,
        )
        new_run_obj.save()
        experiment_name = "" if run_name == timestr else run_name

        # create new project tables based on the project involved in the run and
        # include the project information in projects variable to build the new FORM

        run_info_values = {}
        run_info_values["experiment_name"] = experiment_name
        run_info_values["index_library_name"] = index_library_name
        for key, val in project_list.items():
            if User.objects.filter(username__exact=val).exists():
                user_id = User.objects.get(username__exact=val)
            else:
                user_id = None

            if not Projects.objects.filter(project_name__iexact=key).exists():
                data = {}
                data["user_id"] = user_id
                data["projectName"] = key
                project_obj = Projects.objects.create_new_empty_project(data)
            else:
                project_obj = Projects.objects.filter(project_name__iexact=key).last()

            project_obj.add_run(new_run_obj)
            projects.append([key, val])

        run_info_values["projects_user"] = projects
        run_info_values["runname"] = run_name
        # Get the list of the library kit used (libraryKit)
        list_libraries = LibraryKit.objects.order_by().values_list(
            "library_name", flat=True
        )
        run_info_values["used_libraryKit"] = list_libraries

        user_names = []
        all_users = User.objects.all()
        for user in all_users:
            user_names.append(user.username)
        run_info_values["aval_users"] = user_names
        # displays the list of projects and the user names found on Sample Sheet
        return render(
            request,
            "wetlab/create_next_seq_run.html",
            {"get_user_names": run_info_values},
        )

    # SECOND STEP in collecting data from the NextSeq run. Confirmation /modification of data included in Sample Sheet
    elif request.method == "POST" and (request.POST["action"] == "displayResult"):
        experiment_name = request.POST["experimentname"]
        run_name = request.POST["runname"]
        projects = request.POST.getlist("project")
        user_name = request.POST.getlist("username")
        library_kit = request.POST.getlist("libraryKit")
        project_index_kit = request.POST.getlist("projectindexlibraryname")

        # get the sample sheet used in the run. return error if run already exists
        if not RunProcess.objects.filter(run_name__exact=run_name).exists():
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "You get this error page because you use the back Buttom"
                        " to return to previous page where asking for library kit name",
                        'To upload again the shample sheet, use the "Upload the Run" option from the top menu',
                    ]
                },
            )
        run_p = RunProcess.objects.get(run_name__exact=run_name)
        s_file = run_p.get_sample_file()
        # get the different type of library kit used in the run and
        # convert the sample sheet into Base Space. Number of converted
        # file will be the same as the number of different lybraries use in the run
        library = {}
        results = []

        in_file = os.path.join(settings.MEDIA_ROOT, s_file)
        # Set unique Sample_ID in the sample sheet

        index_file = os.path.join(settings.MEDIA_ROOT, "wetlab", "index_file")
        # if file does not exists create the file and assing the first value
        if not os.path.isfile(index_file):
            with open(index_file, "w") as index_fh:
                index_fh.write("0000-AA")
        create_unique_sample_id_values(in_file, index_file)
        # create the projects/users to update sample sheet
        user_names_in_projects = {}
        for p_index in range(len(projects)):
            user_names_in_projects[projects[p_index]] = user_name[p_index]

        set_user_names_in_sample_sheet(in_file, user_names_in_projects)
        # build the project list for each project_library kit
        for x in range(len(project_index_kit)):
            if project_index_kit[x] in library:
                library[project_index_kit[x]].append(projects[x])
            else:
                library[project_index_kit[x]] = [projects[x]]
        # save the project information on database
        for p in range(len(projects)):
            my_project = projects[p]
            library_kit_id = LibraryKit.objects.filter(
                library_name__exact=library_kit[p]
            ).last()
            update_info_proj = Projects.objects.get(project_name=my_project)
            update_info_proj.libraryKit = project_index_kit[p]
            # removed the link to base space file
            # update_info_proj.baseSpaceFile=bs_file[project_index_kit[p]]
            update_info_proj.baseSpaceFile = None
            update_info_proj.LibraryKit_id = library_kit_id
            update_info_proj.user_id = User.objects.get(username__exact=user_name[p])
            update_info_proj.save()
        results.append(["runname", experiment_name])
        run_p.set_run_state("Recorded")
        sample_sheet_lines = read_all_lines_in_sample_sheet(in_file)
        sample_names_and_data = get_samples_in_sample_sheet(sample_sheet_lines)
        increase_reuse_if_samples_exists(sample_names_and_data["samples"])

        return render(
            request,
            "wetlab/create_next_seq_run.html",
            {"completed_form": results},
        )

    return render(request, "wetlab/create_next_seq_run.html")


@login_required
def add_basespace_library(request):
    """
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
    """

    basespace_library_information = {}
    basespace_library = []

    basespace_library_objects = LibraryKit.objects.all()
    if len(basespace_library_objects) > 0:
        for l_kit in basespace_library_objects:
            basespace_library.append(l_kit.library_name)

    if request.method == "POST" and request.POST["action"] == "addNewBasespaceLibrary":
        new_basespace_library_name = request.POST["newBasespaceLibrary"]

        # Check that library kit is not already defined in database
        if LibraryKit.objects.filter(
            library_name__icontains=new_basespace_library_name
        ).exists():
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "The Library Kit ",
                        new_basespace_library_name,
                        "is already defined on the system",
                    ]
                },
            )

        basespace_library_information[
            "new_basespace_library"
        ] = new_basespace_library_name
        basespace_library.append(new_basespace_library_name)
        # save the new library on database
        b_library = LibraryKit(libraryName=new_basespace_library_name)
        b_library.save()

    basespace_library_information["libraries"] = basespace_library
    return render(
        request,
        "wetlab/AddBasespaceLibrary.html",
        {"list_of_libraries": basespace_library_information},
    )


@login_required
def add_collection_index_kit(request):
    # get the list of the already loaded index library to be displayed
    """
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
    """
    collection_index_information = {}
    collection_index_names = []

    collection_indexes = CollectionIndexKit.objects.all()
    if len(collection_indexes) > 0:
        for c_index in collection_indexes:
            collection_index_names.append(
                [c_index.get_id(), c_index.get_collection_index_name]
            )

    if request.method == "POST" and request.POST["action"] == "addCollectionIndexKit":
        # fetch the file from user form and  build the file name  including
        # the date and time on now to store in database

        file_name = request.FILES["newCollectionIndexFile"].name
        saved_file = wetlab.utils.collection_index.store_collection_kits_file(request.FILES["newCollectionIndexFile"])

        # get the libary name to check if it is already defined
        if not wetlab.utils.collection_index.check_collection_index_file_format(os.path.join(settings.MEDIA_ROOT,saved_file)):
            os.remove(saved_file)
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "The Collection Index Kit file",
                        file_name,
                        "does not have the right format",
                    ]
                },
            )

        collection_name = wetlab.utils.collection_index.get_collection_index_name(os.path.join(settings.MEDIA_ROOT,saved_file))
        if collection_name == "":
            # removing the uploaded file
            os.remove(saved_file)
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "The Collection Index Kit file",
                        file_name,
                        "does not contain the  name",
                    ]
                },
            )

        # check if library name is already defined on database
        if wetlab.utils.collection_index.check_collection_index_exists(collection_name):
            # removing the uploaded file
            os.remove(saved_file)
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "The Collection Index Kit Name ",
                        file_name,
                        "is already defined on iSkyLIMS",
                    ]
                },
            )
        # Get the collection settings included in the file
        collection_settings = wetlab.utils.collection_index.get_collection_settings(os.path.join(settings.MEDIA_ROOT,saved_file))
        new_collection_obj = wetlab.utils.collection_index.store_collection_settings(collection_settings, saved_file)
        # get the index name and index bases for the library
        collection_index = wetlab.utils.collection_index.get_index_values(os.path.join(settings.MEDIA_ROOT,saved_file))
        wetlab.utils.collection_index.store_collection_indexes(collection_index, new_collection_obj)

        collection_index_information["collection_index_names"] = collection_settings[
            "name"
        ]
        collection_index_information["collection_index"] = collection_index_names

        return render(
            request,
            "wetlab/addCollectionIndexKit.html",
            {"collection_index_information": collection_index_information},
        )
    else:
        collection_index_information["collection_index"] = collection_index_names
        return render(
            request,
            "wetlab/addCollectionIndexKit.html",
            {"list_of_collection_index": collection_index_information},
        )


@login_required
def search_run(request):
    """
    Description:
        The function is called from web, having 2 main parts:
            - User form with the information to search runs
            - Result information can be :
                - list of the matched runs
                - run information in case that only 1 match is found
    Input:
        request     # contains the request dictionary sent by django
    Imports:
        Machines and Platform are imported from drylab.models
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
        -- search_run.html is returned with one of the following information :
            -- r_data_display   # in case that only one run is matched
            ---run_list         # in case several run matches the user conditions.

    """
    # check user privileges
    groups = Group.objects.filter(name=wetlab.config.WETLAB_MANAGER).last()
    if groups not in request.user.groups.all():
        allowed_all_runs = False
    else:
        allowed_all_runs = True

    # Search for runs that fullfil the input values
    ########
    run_form_data = get_run_search_fields_form()
    error_message = wetlab.config.ERROR_NO_MATCHES_FOR_RUN_SEARCH
    if request.method == "POST" and (request.POST["action"] == "runsearch"):
        run_name = request.POST["runname"]
        start_date = request.POST["startdate"]
        end_date = request.POST["enddate"]
        run_state = request.POST["runstate"]
        platform_name = request.POST["platform"]
        # check that some values are in the request if not return the form
        if (
            run_name == ""
            and start_date == ""
            and end_date == ""
            and run_state == ""
            and platform_name == ""
        ):
            return render(request, "wetlab/search_run.html")

        # check the right format of start and end date
        if start_date != "":
            if not check_valid_date_format(start_date):
                error_message = wetlab.config.ERROR_INVALID_FORMAT_FOR_DATES
                return render(
                    request,
                    "wetlab/search_run.html",
                    {"run_form_data": run_form_data, "error_message": error_message},
                )
        if end_date != "":
            if not check_valid_date_format(end_date):
                error_message = wetlab.config.ERROR_INVALID_FORMAT_FOR_DATES
                return render(
                    request,
                    "wetlab/search_run.html",
                    {"run_form_data": run_form_data, "error_message": error_message},
                )

        # Get all the available runs to start the filtering
        if allowed_all_runs:
            runs_found = RunProcess.objects.all().order_by("run_date").reverse()
        else:
            user_ids = wetlab.utils.common.get_allowed_user_for_sharing(request.user)
            user_projects = Projects.objects.filter(user_id__in=user_ids)
            run_list = []
            for user_project in user_projects:
                # run_list.append(user_project.runprocess_id.id)
                run_objs = user_project.runProcess.all()
                for run_obj in run_objs:
                    run_list.append(run_obj.get_run_id())
            if RunProcess.objects.filter(pk__in=run_list).exists():
                runs_found = RunProcess.objects.filter(pk__in=run_list)
            else:
                error_message = [
                    "There are not run where " + request.user.username + " was involved"
                ]
                return render(
                    request,
                    "wetlab/search_run.html",
                    {"run_form_data": run_form_data, "error_message": error_message},
                )

        # Get runs when run name is not empty
        if run_name != "":
            if RunProcess.objects.filter(run_name__iexact=run_name).exists():
                run_name_found = RunProcess.objects.filter(run_name__iexact=run_name)
                if len(run_name_found) == 1:
                    return redirect("display_run", run_id=run_name_found[0].pk)
            if runs_found.filter(run_name__icontains=run_name).exists():
                runs_found = runs_found.filter(run_name__icontains=run_name).order_by(
                    "run_name"
                )
            else:
                return render(
                    request,
                    "wetlab/search_run.html",
                    {"run_form_data": run_form_data, "error_message": error_message},
                )
        if platform_name != "":
            sequencer_list = get_sequencer_names_from_platform(platform_name)
            if len(sequencer_list) > 0:
                runs_found = runs_found.filter(
                    used_sequencer__sequencer_name__in=sequencer_list
                )
            else:
                return render(
                    request,
                    "wetlab/search_run.html",
                    {"run_form_data": run_form_data, "error_message": error_message},
                )
        # Check if state is not empty
        if run_state != "":
            s_state = RunStates.objects.filter(run_state_name__exact=run_state).last()
            if runs_found.filter(state__run_state_name__exact=s_state).exists():
                runs_found = runs_found.filter(
                    state__run_state_name__exact=s_state
                ).order_by("run_name")
            else:
                return render(
                    request,
                    "wetlab/search_run.html",
                    {"run_form_data": run_form_data, "error_message": error_message},
                )
        # Check if start_date is not empty
        if start_date != "" and end_date != "":
            if runs_found.filter(run_date__range=(start_date, end_date)).exists():
                runs_found = runs_found.filter(run_date__range=(start_date, end_date))
            else:
                return render(
                    request,
                    "wetlab/search_run.html",
                    {"run_form_data": run_form_data, "error_message": error_message},
                )
        if start_date != "" and end_date == "":
            if runs_found.filter(run_date__gte=start_date).exists():
                runs_found = runs_found.filter(run_date__gte=start_date)
            else:
                return render(
                    request,
                    "wetlab/search_run.html",
                    {"run_form_data": run_form_data, "error_message": error_message},
                )
        if start_date == "" and end_date != "":
            if runs_found.filter(run_date__lte=end_date).exists():
                runs_found = runs_found.filter(run_date__lte=end_date)
            else:
                return render(
                    request,
                    "wetlab/search_run.html",
                    {"run_form_data": run_form_data, "error_message": error_message},
                )

        # If only 1 run mathes the user conditions, then get the project information

        if len(runs_found) == 1:
            return redirect("display_run", run_id=runs_found[0].pk)

        else:
            # collect the list of run that matches the run date
            run_list = []
            for run_found in runs_found:
                run_list.append(
                    [
                        run_found.get_run_id(),
                        run_found.get_run_name(),
                        run_found.get_run_date(),
                    ]
                )
            return render(
                request,
                "wetlab/search_run.html",
                {"display_run_list": run_list},
            )
    else:
        return render(
            request, "wetlab/search_run.html", {"run_form_data": run_form_data}
        )


@login_required
def search_project(request):
    """
    Description:
        The function is called from web, having 2 main parts:
            - User form with the information to search projects
            - Result information can be :
                - list of the matched projects
                - project information in case that only 1 match is found

    Input:
        request     # contains the request dictionary sent by django
    Imports:
        Machines and Platform are imported from drylab.models
            for filtering runs based on the platform

    Return:
        Return the different information depending on the execution:
        -- Error page in case no run is founded on the matching conditions.
        -- search_run.html is returned with one of the following information :
            -- r_data_display   # in case that only one run is matched
            ---run_list         # in case several run matches the user conditions.

    """
    project_form_data = get_project_search_fields_form()
    error_message = wetlab.config.ERROR_NO_MATCHES_FOR_PROJECT_SEARCH
    if request.method == "POST" and (request.POST["action"] == "searchproject"):
        project_name = request.POST["projectname"]
        start_date = request.POST["startdate"]
        end_date = request.POST["enddate"]
        user_name = request.POST["username"]
        sequencer_name = request.POST["sequencer"]
        run_state = request.POST["runstate"]
        # check that some values are in the request if not return the form
        if (
            project_name == ""
            and start_date == ""
            and end_date == ""
            and user_name == ""
            and sequencer_name == ""
            and run_state == ""
        ):
            return render(
                request,
                "wetlab/search_project.html",
                {"project_form_data": project_form_data},
            )

        if user_name != "" and len(user_name) < 5:
            error_message = wetlab.config.ERROR_USER_NAME_TOO_SHORT
            return render(
                request,
                "wetlab/search_project.html",
                {
                    "project_form_data": project_form_data,
                    "error_message": error_message,
                },
            )

        # check the right format of start and end date
        if start_date != "":
            if not check_valid_date_format(start_date):
                error_message = wetlab.config.ERROR_INVALID_FORMAT_FOR_DATES
                return render(
                    request,
                    "wetlab/search_project.html",
                    {
                        "project_form_data": project_form_data,
                        "error_message": error_message,
                    },
                )

        if end_date != "":
            if not check_valid_date_format(start_date):
                error_message = wetlab.config.ERROR_INVALID_FORMAT_FOR_DATES
                return render(
                    request,
                    "wetlab/search_project.html",
                    {
                        "project_form_data": project_form_data,
                        "error_message": error_message,
                    },
                )

        projects_found = Projects.objects.all()

        if project_name != "":
            projects_found = projects_found.filter(project_name__icontains=project_name)
        if sequencer_name != "":
            run_objs = RunProcess.objects.filter(
                used_sequencer__sequencer_name__exact=sequencer_name
            )
            projects_found = projects_found.filter(run_process__in=run_objs)
        if run_state != "":
            run_objs = RunProcess.objects.filter(state__run_state_name__exact=run_state)
            projects_found = projects_found.filter(run_process__in=run_objs)
        if user_name != "":
            # check if user has a shared user
            if user_name == request.user.username:
                p_shared_list = wetlab.utils.common.get_allowed_user_for_sharing(request.user)
                projects_found = projects_found.filter(user_id__in=p_shared_list)
            else:
                projects_found = projects_found.filter(
                    user_id__username__icontains=user_name
                )
        if start_date != "" and end_date != "":
            projects_found = projects_found.filter(
                generated_at__range=(start_date, end_date)
            )
        if start_date != "" and end_date == "":
            projects_found = projects_found.filter(generated_at__gte=start_date)
        if start_date == "" and end_date != "":
            projects_found = projects_found.filter(generated_at__lte=end_date)
        if len(projects_found) == 0:
            error_message = wetlab.config.ERROR_NO_MATCHES_FOR_PROJECT_SEARCH
            return render(
                request,
                "wetlab/search_project.html",
                {
                    "project_form_data": project_form_data,
                    "error_message": error_message,
                },
            )

        if len(projects_found) == 1:
            return redirect("display_project", project_id=projects_found[0].id)
        else:
            # Display a list with all projects that matches the conditions
            project_list_dict = {}
            project_list = []
            for project in projects_found:
                p_name = project.get_project_name()
                p_name_id = project.id
                project_list.append([p_name, p_name_id])
            project_list_dict["projects"] = project_list
            return render(
                request,
                "wetlab/search_project.html",
                {"display_project_list": project_list_dict},
            )

    else:
        return render(
            request,
            "wetlab/search_project.html",
            {"project_form_data": project_form_data},
        )


@login_required
def retry_error_run(request):
    # check user privileges
    if request.user.is_authenticated:
        try:
            groups = Group.objects.get(name="WetlabManager")
            if groups not in request.user.groups.all():
                return render(
                    request,
                    "wetlab/error_page.html",
                    {
                        "content": [
                            "You do have the enough privileges to see this page ",
                            "Contact with your administrator .",
                        ]
                    },
                )
        except Exception:
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "You do have the enough privileges to see this page ",
                        "Contact with your administrator .",
                    ]
                },
            )
    else:
        # redirect to login webpage
        return redirect("/accounts/login")
    if request.method == "POST" and (request.POST["action"] == "retry_correct_error"):
        run_id = request.POST["run_id"]
        if RunProcess.objects.filter(pk__exact=run_id).exists():
            run_name_found = RunProcess.objects.get(pk__exact=run_id)
            previous_error_state = run_name_found.get_state_before_error()
            run_name_found.set_run_state(previous_error_state)
            detail_description = {}
            detail_description["information"] = wetlab.config.SUCCESSFUL_RUN_STATE_CHANGE_FOR_RETRY
            return render(
                request,
                "wetlab/successful_page.html",
                {"detail_description": detail_description, "return_main_menu": True},
            )
        else:
            return render(
                request,
                "wetlab/error_page.html",
                {"content": ["Run does not exist "]},
            )
    else:
        # return redirect (request,'/')
        return render(request, "wetlab/index.html")


@login_required
def skip_cancel_situation(request):
    # check user privileges
    if request.user.is_authenticated:
        try:
            groups = Group.objects.get(name="WetlabManager")
            if groups not in request.user.groups.all():
                return render(
                    request,
                    "wetlab/error_page.html",
                    {
                        "content": [
                            "You do have the enough privileges to see this page ",
                            "Contact with your administrator .",
                        ]
                    },
                )
        except Exception:
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "You do have the enough privileges to see this page ",
                        "Contact with your administrator .",
                    ]
                },
            )
    else:
        # redirect to login webpage
        return redirect("/accounts/login")
    if request.method == "POST" and (request.POST["action"] == "skip_cancel_situation"):
        run_id = request.POST["run_id"]
        if RunProcess.objects.filter(pk__exact=run_id).exists():
            run_name_found = RunProcess.objects.get(pk__exact=run_id)
            run_name_found.set_run_state("Sample Sent")
            run_name_found.set_forced_continue_on_error()
            detail_description = {}
            detail_description["information"] = wetlab.config.SUCCESSFUL_RUN_STATE_CHANGE_FOR_RETRY
            return render(
                request,
                "wetlab/successful_page.html",
                {"detail_description": detail_description, "return_main_menu": True},
            )
        else:
            return render(
                request,
                "wetlab/error_page.html",
                {"content": ["Run does not exist "]},
            )
    else:
        # return redirect (request,'/')
        return render(request, "wetlab/index.html")


@login_required
def display_run(request, run_id):
    # check user privileges
    if not request.user.is_authenticated:
        # redirect to login webpage
        return redirect("/accounts/login")
    groups = Group.objects.get(name=wetlab.config.WETLAB_MANAGER)
    if groups not in request.user.groups.all():
        # check if user is owner of the run or belongs to the shared user
        shared_user_ids = wetlab.utils.common.get_allowed_user_for_sharing(request.user)
        if Projects.objects.filter(run_process__exact=run_id).exists():
            projects = Projects.objects.filter(run_process__exact=run_id)
            allowed = False
            for project in projects:
                if int(project.get_user_center_name()) in shared_user_ids:
                    allowed = True
                    break
            if not allowed:
                return render(
                    request,
                    "wetlab/error_page.html",
                    {
                        "content": [
                            "You do have the enough privileges to see this page ",
                            "Contact with your administrator .",
                        ]
                    },
                )
        else:
            return render(
                request,
                "wetlab/error_page.html",
                {"content": ["No matches have been found for the run "]},
            )

    if RunProcess.objects.filter(pk=run_id).exists():
        run_name_found = RunProcess.objects.get(pk=run_id)
        r_data_display = get_information_run(run_name_found)
        return render(
            request,
            "wetlab/displayRun.html",
            {"display_one_run": r_data_display},
        )
    else:
        return render(
            request,
            "wetlab/error_page.html",
            {"content": ["No matches have been found for the run  "]},
        )


@login_required
def last_run_by_sequencer(request):
    # check user privileges
    if request.user.is_authenticated:
        try:
            groups = Group.objects.get(name="WetlabManager")
            if groups not in request.user.groups.all():
                return render(
                    request,
                    "wetlab/error_page.html",
                    {
                        "content": [
                            "You do have the enough privileges to see this page ",
                            "Contact with your administrator .",
                        ]
                    },
                )
        except Exception:
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "You do have the enough privileges to see this page ",
                        "Contact with your administrator .",
                    ]
                },
            )
    else:
        # redirect to login webpage
        return redirect("/accounts/login")

    last_runs = get_last_runs_by_sequencer()
    if len(last_runs) == 0:
        return render(
            request, "wetlab/last_run_by_sequencer.html", {"no_runs": "no_runs"}
        )
    if len(last_runs) > 1:
        return render(
            request, "wetlab/last_run_by_sequencer.html", {"last_runs": last_runs}
        )
    else:
        # if only 1 sequencer is defined, then display the information of the latest run
        return redirect("display_run", run_id=last_runs[0][2])


@login_required
def incompleted_runs(request):
    # check user privileges
    if request.user.is_authenticated:
        try:
            groups = Group.objects.get(name="WetlabManager")
            if groups not in request.user.groups.all():
                return render(
                    request,
                    "wetlab/error_page.html",
                    {
                        "content": [
                            "You do have the enough privileges to see this page ",
                            "Contact with your administrator .",
                        ]
                    },
                )
        except Exception:
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "You do have the enough privileges to see this page ",
                        "Contact with your administrator .",
                    ]
                },
            )
    else:
        # redirect to login webpage
        return redirect("/accounts/login")
    if RunProcess.objects.all().exclude(state__run_state_name="Completed").exists():
        display_incompleted_run = get_information_for_incompleted_run()
        return render(
            request,
            "wetlab/incompleted_runs.html",
            {"display_incompleted_run": display_incompleted_run},
        )
    else:
        return render(request,"wetlab/incompleted_runs.html")


def check_user_access(request, project_found_id):
    groups = Group.objects.get(name=wetlab.config.WETLAB_MANAGER)
    # check if user belongs to WetlabManager . If true allow to see the page
    if groups not in request.user.groups.all():
        # check if project belongs to the same user as the one requesting the page
        if project_found_id.user_id.id != request.user.id:
            return False
    return True


@login_required
def display_project(request, project_id):
    if not request.user.is_authenticated:
        # redirect to login webpage
        return redirect("/accounts/login")
    if Projects.objects.filter(pk=project_id).exists():
        project_obj = Projects.objects.filter(pk=project_id).last()

        # check that user is allow to see the project
        groups = Group.objects.get(name=wetlab.config.WETLAB_MANAGER)
        if groups not in request.user.groups.all():
            p_shared_list = wetlab.utils.common.get_allowed_user_for_sharing(request.user)
            if int(project_obj.get_user_center_name()) not in p_shared_list:
                return render(
                    request,
                    "wetlab/error_page.html",
                    {
                        "content": [
                            "You do have the enough privileges to see this page ",
                            "Contact with your administrator .",
                        ]
                    },
                )
        # Display the proyect information
        display_project_data = get_information_project(project_obj, request)
        return render(
            request,
            "wetlab/displayProject.html",
            {"display_project_data": display_project_data},
        )
    else:
        return render(
            request,
            "wetlab/error_page.html",
            {"content": ["No matches have been found for the project  "]},
        )


@login_required
def display_collection_index(request, collection_index_id):
    if wetlab.models.CollectionIndexKit.objects.filter(pk=collection_index_id).exists():
        collection_index_dict = wetlab.utils.collection_index.get_collection_index_information(collection_index_id)
        if collection_index_dict is not False:
            return render(
                request,
                "wetlab/display_collection_index.html",
                {"display_one_collection_index": collection_index_dict},
            )
        else:
            return render(
                request,
                "wetlab/display_collection_index.html",
                {
                    "error_message":
                        "There are recorded information for the collection index for your request"
                },
            )
    else:
        return render(
            request,
            "wetlab/display_collection_index.html",
            {"error_message": "There is no information for the Collection index "},
        )


@login_required
def search_collection_index_library(request):
    if request.method == "POST" and (
        request.POST["action"] == "searchcollectionindexkit"
    ):
        collection_index_kit_name = request.POST["collectionindexkitname"]
        adapter_1 = request.POST["adapter1"]
        adapter_2 = request.POST["adapter2"]
        index_name = request.POST["indexname"]
        index_sequence = request.POST["indexbase"]

        # check that some values are in the request if not return the form
        if (
            collection_index_kit_name == ""
            and adapter_1 == ""
            and adapter_2 == ""
            and index_name == ""
            and index_sequence == ""
        ):
            return render(request, "wetlab/search_collection_index_library.html")

        if index_sequence != "":
            if len(index_sequence) < 6:
                return render(
                    request,
                    "wetlab/search_collection_index_library.html",
                    {"error_message": wetlab.config.ERROR_TOO_SHORT_INDEX_BASE_SEQUENCE},
                )
            else:
                valid_seq_characters = ["a", "A", "c", "C", "g", "G", "t", "T"]
                for letter in index_sequence:
                    if letter not in valid_seq_characters:
                        return render(
                            request,
                            "wetlab/search_collection_index_library.html",
                            {"error_message": wetlab.config.ERROR_INVALID_SEQUENCE_CHARACTERS},
                        )

        collection_indexes = CollectionIndexKit.objects.all()
        if collection_index_kit_name != "":
            if collection_indexes.filter(
                collection_index_name__icontains=collection_index_kit_name
            ).exists():
                collection_indexes = collection_indexes.filter(
                    collection_index_name__icontains=collection_index_kit_name
                )
                if len(collection_indexes) == 1:
                    return redirect(
                        "display_collection_index",
                        collection_index_id=collection_indexes[0].get_id(),
                    )
            else:
                return render(
                    request,
                    "wetlab/search_collection_index_library.html",
                    {"not_found_matchs": "not_found_matchs"},
                )
        if adapter_1 != "":
            if collection_indexes.filter(adapter_1__icontains=adapter_1).exists():
                collection_indexes = collection_indexes.filter(
                    adapter_1__icontains=adapter_1
                )
            else:
                return render(
                    request,
                    "wetlab/search_collection_index_library.html",
                    {"not_found_matchs": "not_found_matchs"},
                )
        if adapter_2 != "":
            if collection_indexes.filter(adapter_1__icontains=adapter_2).exists():
                collection_indexes = collection_indexes.filter(
                    adapter_2__icontains=adapter_2
                )
            else:
                return render(
                    request,
                    "wetlab/search_collection_index_library.html",
                    {"not_found_matchs": "not_found_matchs"},
                )

        if index_name != "" or index_sequence != "":
            collection_values = CollectionIndexValues.objects.all()
            if index_name != "":
                if collection_values.filter(
                    index_7_contains=index_name,
                    collection_index_kit_id__in=collection_indexes,
                ).exists():
                    collection_values = collection_values.filter(
                        index_7_contains=index_name,
                        collection_index_kit_id__in=collection_indexes,
                    )
                elif collection_values.filter(
                    index_5_contains=index_name,
                    collection_index_kit_id__in=collection_indexes,
                ).exists():
                    collection_values = collection_values.filter(
                        index_5_contains=index_name,
                        collection_index_kit_id__in=collection_indexes,
                    )
                else:
                    return render(
                        request,
                        "wetlab/search_collection_index_library.html",
                        {"not_found_matchs": "not_found_matchs"},
                    )

            if index_sequence != "":
                index_found, sequence = find_index_sequence_collection_values_kit(
                    index_sequence
                )
                if "I7" in index_found:
                    collection_values = collection_values.filter(
                        i_7_seq__icontains=sequence,
                        collection_index_kit_id__in=collection_indexes,
                    )
                elif "I5" in index_found:
                    collection_values = collection_values.filter(
                        i_5_seq__icontains=sequence,
                        collection_index_kit_id__in=collection_indexes,
                    )
                else:
                    return render(
                        request,
                        "wetlab/search_collection_index_library.html",
                        {"not_found_matchs": "not_found_matchs"},
                    )

            if len(collection_values) == 1:
                return redirect(
                    "display_collection_index",
                    collection_index_id=collection_values[0].get_collection_index_id,
                )
            else:
                matched_collection_index = []
                collection_index_id_list = []

                for collection_value in collection_values:
                    if (
                        collection_value.get_collection_index_id()
                        not in collection_index_id_list
                    ):
                        collection_index_id_list.append(
                            collection_value.get_collection_index_id()
                        )
                        matched_collection_index.append(
                            [
                                collection_value.get_collection_index_id(),
                                collection_value.get_collection_index_name(),
                            ]
                        )
                if len(matched_collection_index) == 1:
                    return redirect(
                        "display_collection_index",
                        collection_index_id=matched_collection_index[0][0],
                    )
                else:
                    return render(
                        request,
                        "wetlab/search_collection_index_library.html",
                        {"matched_collection_index": matched_collection_index},
                    )
        else:
            if len(collection_indexes) == 1:
                return redirect(
                    "display_collection_index",
                    collection_index_id=collection_indexes[0].get_id(),
                )
            else:
                matched_collection_index = []
                for collection_index in collection_indexes:
                    matched_collection_index.append(
                        [
                            collection_index.get_id(),
                            collection_index.get_collection_index_name(),
                        ]
                    )
                return render(
                    request,
                    "wetlab/search_collection_index_library.html",
                    {"matched_collection_index": matched_collection_index},
                )

    else:
        return render(request, "wetlab/search_collection_index_library.html")


@login_required
def change_run_name(request, run_id):
    if RunProcess.objects.filter(pk=run_id).exists():
        run = RunProcess.objects.get(pk=run_id)
        if not request.user.is_authenticated:
            return redirect("/accounts/login")
        # check if user is allow to make the change
        groups = Group.objects.get(name="WetlabManager")
        # check if user belongs to WetlabManager . If true allow to see the page
        if groups not in request.user.groups.all():
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "You do have the enough privileges to see this page ",
                        "Contact with your administrator .",
                    ]
                },
            )
        if request.method == "POST" and request.POST["action"] == "change_run_name":
            new_run_name = request.POST["runName"]
            if new_run_name == "":
                return render(
                    request,
                    "wetlab/error_page.html",
                    {"content": ["Empty value is not allowed for the Run Name "]},
                )
            if RunProcess.objects.filter(run_name__exact=new_run_name).exists():
                return render(
                    request,
                    "wetlab/error_page.html",
                    {
                        "content": [
                            "The given Run Name is already in use",
                            "Go back to the previous page and change the run name",
                        ]
                    },
                )
            changed_run_name = {}
            old_run_name = run.run_name
            run.run_name = new_run_name
            run.save()
            changed_run_name["new_run_name"] = [[new_run_name, run_id]]
            changed_run_name["old_run_name"] = old_run_name

            return render(
                request,
                "wetlab/ChangeRunName.html",
                {"changed_run_name": changed_run_name},
            )
        else:
            form_change_run_name = {}
            form_change_run_name["run_name"] = run.run_name
            return render(
                request,
                "wetlab/ChangeRunName.html",
                {"form_change_run_name": form_change_run_name},
            )
    else:
        return render(
            request,
            "wetlab/error_page.html",
            {"content": ["There is no Run for your query  "]},
        )


@login_required
def change_project_libKit(request, project_id):
    # check if project exists
    if Projects.objects.filter(pk=project_id).exists():
        project = Projects.objects.get(pk=project_id)
        if not request.user.is_authenticated:
            # redirect to login webpage
            return redirect("/accounts/login")
        # check that user is allow to make the change
        allowed_access = check_user_access(request, project)
        if not allowed_access:
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "You do have the enough privileges to see this page ",
                        "Contact with your administrator .",
                    ]
                },
            )

        if (
            request.method == "POST"
            and request.POST["action"] == "change_project_libKit"
        ):
            new_library_name = request.POST["projectlibkit"]
            old_library_name = project.get_index_library_name()
            if old_library_name == new_library_name:
                return render(
                    request,
                    "wetlab/info_page.html",
                    {
                        "content": [
                            "The library kit from the input text is the same to the existing defined for this project",
                            "No change is done",
                        ]
                    },
                )
            # check if there is no other project in the same Run with the same Library Kit
            # if the library is shared with other project then error message is displayed

            if not check_user_group(request, "WetlabManager"):
                project_run_id = project.runprocess_id.id
                project_lib_kit = project.libraryKit
                if (
                    Projects.objects.filter(runprocess_id=project_run_id)
                    .exclude(pk=project_id)
                    .exists()
                ):
                    all_project_with_same_run_id = Projects.objects.filter(
                        runprocess_id=project_run_id
                    ).exclude(pk=project_id)
                    # there are more than 1 project on the same Run.
                    # check if these projects have in common the same library kit
                    other_lib_kits = []
                    for other_project in all_project_with_same_run_id:
                        other_lib_kits.append(other_project.libraryKit)
                    if project_lib_kit in other_lib_kits:
                        message = str(
                            "The library Kit "
                            + old_library_name
                            + "is shared with other projects in the same Run "
                        )

                        return render(
                            request,
                            "wetlab/error_page.html",
                            {
                                "content": [
                                    message,
                                    "",
                                    "Contact with your administrator .",
                                ]
                            },
                        )
            old_lib_kit_file = project.baseSpaceFile
            new_file_name = new_library_name.replace(" ", "_")
            #
            new_file = update_library_kit_field(
                old_lib_kit_file, new_file_name, new_library_name
            )
            if new_file == "ERROR":
                return render(
                    request, "wetlab/error_page.html", {"content": [""]}
                )
            # update the database with new file
            project.baseSpaceFile = new_file
            project.libraryKit = new_library_name
            project.save()
            # Preparing the data to show in the web page
            new_file = project.baseSpaceFile
            change_library_kit_dict = {}
            change_library_kit_dict["project"] = project.projectName
            change_library_kit_dict["library_name"] = new_library_name
            change_library_kit_dict["file_to_download"] = new_file

            #
            return render(
                request,
                "wetlab/ChangeProjectLibraryKit.html",
                {"changed_lib_kit": change_library_kit_dict},
            )
        else:
            form_change_lib_kit = {}
            project_data = []
            project_name = project.projectName
            form_change_lib_kit["project_name"] = project_name

            project_info_text = [
                "Run Name",
                "Project Name",
                "Project date",
                "User Name",
                "Library Kit",
            ]
            project_values = project.get_p_info_change_library().split(";")

            for item in range(len(project_info_text)):
                project_data.append([project_info_text[item], project_values[item]])
                form_change_lib_kit["project_data"] = project_data
            return render(
                request,
                "wetlab/ChangeProjectLibraryKit.html",
                {"form_change_lib_kit": form_change_lib_kit},
            )
    else:
        return render(
            request,
            "wetlab/error_page.html",
            {"content": ["No project has been found for changing the library Kit "]},
        )


@login_required
def stats_experiment(request):
    return render(request, "wetlab/StatsPerExperiment.html", {})


@login_required
def stats_per_sequencer(request):
    sequencer_names = get_sequencer_installed_names()
    if request.method == "POST":
        sequencer = request.POST["sequencer"]
        start_date = request.POST["startdate"]
        end_date = request.POST["enddate"]

        if start_date != "":
            if not check_valid_date_format(start_date):
                error_message = wetlab.config.ERROR_INVALID_FORMAT_FOR_DATES
                return render(
                    request,
                    "wetlab/StatsPerSequencer.html",
                    {
                        "sequencer_names": sequencer_names,
                        "error_message": error_message,
                    },
                )
        if end_date != "":
            if not check_valid_date_format(end_date):
                error_message = wetlab.config.ERROR_INVALID_FORMAT_FOR_DATES
                return render(
                    request,
                    "wetlab/StatsPerSequencer.html",
                    {
                        "sequencer_names": sequencer_names,
                        "error_message": error_message,
                    },
                )
        runs_using_sequencer = get_sequencers_run_from_time_interval(
            sequencer, start_date, end_date
        )
        if len(runs_using_sequencer) == 0:
            error_message = wetlab.config.ERROR_NO_MATCHES_FOR_SEQUENCER_STATS
            return render(
                request,
                "wetlab/StatsPerSequencer.html",
                {"sequencer_names": sequencer_names, "error_message": error_message},
            )
        sequencer_data = get_stats_sequencer_data_from_selected_runs(
            runs_using_sequencer, sequencer, start_date, end_date
        )

        return render(
            request,
            "wetlab/StatsPerSequencer.html",
            {"sequencer_data": sequencer_data},
        )

    else:
        return render(
            request,
            "wetlab/StatsPerSequencer.html",
            {"sequencer_names": sequencer_names},
        )


@login_required
def stats_per_researcher(request):
    if request.method == "POST":
        r_name = request.POST["researchername"]
        start_date = request.POST["startdate"]
        end_date = request.POST["enddate"]

        researcher_statistics = get_researcher_statistics(r_name, start_date, end_date)
        if "ERROR" in researcher_statistics:
            error_message = researcher_statistics["ERROR"]
            return render(
                request,
                "wetlab/StatsPerResearcher.html",
                {"researcher_statistics": error_message},
            )

        return render(
            request,
            "wetlab/StatsPerResearcher.html",
            {"researcher_statistics": researcher_statistics},
        )

    else:
        return render(request, "wetlab/StatsPerResearcher.html", {})


@login_required
def stats_per_time(request):
    if request.method == "POST":
        start_date = request.POST["startdate"]
        end_date = request.POST["enddate"]
        # check the right format of start and end date
        if start_date != "" and not check_valid_date_format(start_date):
            error_message = wetlab.config.ERROR_INVALID_FORMAT_FOR_DATES
            return render(
                request, "wetlab/StatsPerTime.html", {"ERROR": error_message}
            )
        if end_date != "" and not check_valid_date_format(start_date):
            error_message = wetlab.config.ERROR_INVALID_FORMAT_FOR_DATES
            return render(
                request, "wetlab/StatsPerTime.html", {"ERROR": error_message}
            )
        ########
        # searching for runs were match the state and start and end date
        ########
        if start_date != "" and end_date != "":
            stat_per_time = {}
            if RunProcess.objects.filter(
                state__run_state_name="Completed", run_date__range=(start_date, end_date)
            ).exists():
                run_stats_list = RunProcess.objects.filter(
                    state__run_state_name="Completed",
                    run_date__range=(start_date, end_date),
                ).order_by("run_date")

                run_list = {}
                run_date_name = {}
                # get the run names that matches de conditions
                for run in run_stats_list:
                    # run_list.append([run.get_run_name(),run.id])
                    run_date = str(run.run_date)
                    if run_date in run_date_name:
                        run_date_name[run_date] += 1
                    else:
                        run_date_name[run_date] = 1
                    run_list[run.id] = [[run.get_run_name(), run_date]]
                    #
                stat_per_time["run_names"] = run_list
                if len(run_stats_list) == 1:
                    number_of_runs = "1 Run"
                else:
                    number_of_runs = str(len(run_stats_list)) + "  Runs"
                stat_per_time["number_of_runs"] = number_of_runs
                stat_per_time["dates"] = start_date + " and  " + end_date
                #
                ########
                # define the graphics for found run in the period
                heading = (
                    "Runs found during the period "
                    + str(start_date)
                    + " and "
                    + str(end_date)
                )
                sub_caption = ""
                x_axis_name = "Date"
                y_axis_name = "Number of runs"
                run_period_chart_number = "run_period_chart-1"
                run_period_index_graph = "exq1"

                data_source = researcher_project_column_graphic(
                    heading,
                    sub_caption,
                    x_axis_name,
                    y_axis_name,
                    "ocean",
                    run_date_name,
                )
                stat_per_time["run_period_graphic"] = core.fusioncharts.fusioncharts.FusionCharts(
                    "column3d",
                    run_period_index_graph,
                    "550",
                    "350",
                    run_period_chart_number,
                    "json",
                    data_source,
                ).render()

                # end creation run preparation graphics

                ########
                # collect statistics for Projects
                if RunProcess.objects.filter(
                    state__run_state_name__exact="Completed",
                    run_completed_date__range=(start_date, end_date),
                ).exists():
                    run_objs = RunProcess.objects.filter(
                        state__run_state_name__exact="Completed",
                        run_completed_date__range=(start_date, end_date),
                    )
                    if Projects.objects.filter(run_process__in=run_objs).exists():
                        project_found_list = Projects.objects.filter(
                            run_process__in=run_objs
                        )

                        project_list = {}
                        project_date_name = {}
                        # get the project names that matches de conditions
                        for project in project_found_list:
                            project_run_date = str(project.project_run_date)
                            if project_run_date in project_date_name:
                                project_date_name[project_run_date] += 1
                            else:
                                project_date_name[project_run_date] = 1
                            project_list[project.id] = [
                                [project.get_project_name(), project_run_date]
                            ]
                            #
                        stat_per_time["project_names"] = project_list
                        if len(project_found_list) == 1:
                            number_of_projects = "1 Project"
                        else:
                            number_of_projects = (
                                str(len(project_found_list)) + "  Projects"
                            )
                        stat_per_time["number_of_projects"] = number_of_projects
                        stat_per_time["dates"] = start_date + " and  " + end_date
                        #
                        ########
                        # define the graphics for found run in the period
                        heading = (
                            "Projects found during the period "
                            + str(start_date)
                            + " and "
                            + str(end_date)
                        )
                        sub_caption = ""
                        x_axis_name = "Date"
                        y_axis_name = "Number of Projects"
                        run_period_chart_number = "project_period_chart-1"
                        run_period_index_graph = "project_period-1"

                        data_source = researcher_project_column_graphic(
                            heading,
                            sub_caption,
                            x_axis_name,
                            y_axis_name,
                            "carbon",
                            project_date_name,
                        )
                        stat_per_time["project_period_graphic"] = core.fusioncharts.fusioncharts.FusionCharts(
                            "column3d",
                            run_period_index_graph,
                            "550",
                            "350",
                            run_period_chart_number,
                            "json",
                            data_source,
                        ).render()

                # end creation run preparation graphics
                ########
                # collect statistics for unkow Barcodes
                # top_unbarcode_list = []
                count_unbarcode = {}

                for run in run_stats_list:
                    run_obj = run.get_run_id()
                    run_param_obj = RunningParameters.objects.get(run_name_id=run_obj)
                    lanes_in_sequencer = int(run_param_obj.get_number_of_lanes())
                    top_unbarcode_all_runs = {}
                    for lane_number in range(1, lanes_in_sequencer + 1):
                        lane_unbarcodes = RawTopUnknowBarcodes.objects.filter(
                            runprocess_id=run, lane_number__exact=lane_number
                        )
                        for lane_unbarcode in lane_unbarcodes:
                            if lane_number not in count_unbarcode:
                                count_unbarcode[lane_number] = {}
                            (
                                unbarcode_num,
                                unknown_barcode,
                            ) = lane_unbarcode.get_unknow_barcodes().split(";")
                            value_unbarcode = int(unbarcode_num.replace(",", ""))
                            if unknown_barcode not in count_unbarcode[lane_number]:
                                count_unbarcode[lane_number][
                                    unknown_barcode
                                ] = value_unbarcode
                            else:
                                count_unbarcode[lane_number][
                                    unknown_barcode
                                ] += value_unbarcode
                            if unknown_barcode not in top_unbarcode_all_runs:
                                top_unbarcode_all_runs[
                                    unknown_barcode
                                ] = value_unbarcode
                            else:
                                top_unbarcode_all_runs[
                                    unknown_barcode
                                ] += value_unbarcode

                themes = ["", "ocean", "fint", "carbon", "zune", ""]
                # prepare the column graphic for nunber of top Unknow Barcode
                unbar_lane_chart = []
                for lane_number in range(1, lanes_in_sequencer + 1):
                    heading = "Number of undetermined barcode sequence in lane " + str(
                        lane_number
                    )
                    chart_number = "chart-" + str(lane_number)
                    render_number = "ex" + str(lane_number)
                    data_source = graphic_for_unbarcodes(
                        heading, themes[lane_number], count_unbarcode[lane_number]
                    )
                    lane_graphic = core.fusioncharts.fusioncharts.FusionCharts(
                        "column3d",
                        render_number,
                        "500",
                        "400",
                        chart_number,
                        "json",
                        data_source,
                    )
                    unbar_lane_chart.append(
                        [chart_number, str(lane_number), lane_graphic.render()]
                    )
                stat_per_time["unbar_lane_chart"] = unbar_lane_chart

                # prepare the pie graphic for the number of top Unknow Barcode per sequence
                data_source = pie_graphic(
                    "Number of count for the Undetermined Sequences",
                    "fint",
                    top_unbarcode_all_runs,
                )
                unknow_pie3d = core.fusioncharts.fusioncharts.FusionCharts(
                    "pie3d", "ex5", "500", "400", "chart-5", "json", data_source
                )
                stat_per_time["unknow_pie3d"] = unknow_pie3d.render()

                ####
                # Insert information for disk space utilization
                run_disk_utilization = {}
                for run_disk_stats in run_stats_list:
                    run_name_disk = run_disk_stats.run_name
                    run_disk_utilization[
                        run_name_disk
                    ] = run_disk_stats.get_disk_space_utilization()

                heading = (
                    "Disk space used for each Run found during the period "
                    + str(start_date)
                    + " and "
                    + str(end_date)
                )
                sub_caption = ""
                x_axis_name = "Date"
                y_axis_name = "Disk space used (MB)"
                disk_space_period_chart_number = "disk_usage_chart-1"
                disk_space_period_index_graph = "diskusage1"

                data_source = researcher_project_column_graphic(
                    heading,
                    sub_caption,
                    x_axis_name,
                    y_axis_name,
                    "carbon",
                    run_disk_utilization,
                )
                stat_per_time["disk_space_period_graphic"] = core.fusioncharts.fusioncharts.FusionCharts(
                    "column3d",
                    disk_space_period_index_graph,
                    "950",
                    "350",
                    disk_space_period_chart_number,
                    "json",
                    data_source,
                ).render()

                #
                return render(
                    request,
                    "wetlab/StatsPerTime.html",
                    {"display_stats_per_time": stat_per_time},
                )

            else:
                return render(
                    request,
                    "wetlab/error_page.html",
                    {
                        "content": [
                            "No matches have been found for Runs created between",
                            start_date,
                            " and the ",
                            end_date,
                        ]
                    },
                )
        else:
            return render(
                request,
                "wetlab/error_page.html",
                {"content": "Start date and End Date cannot be empty "},
            )

    return render(request, "wetlab/StatsPerTime.html")


def get_list_of_libraries_values(
    library_found, q30_comparations, mean_comparations, n_bases_comparations
):
    for project_to_compare in library_found:
        library_to_compare_name = project_to_compare.get_index_library_name()
        # project_to_compare_id = project_to_compare.id
        q30_compare_lib, mean_compare_lib, yield_mb_compare_lib = [], [], []

        # This line must changed to handle project name is reused in several runs
        run_used_in_project = project_to_compare.runProcess.all().last()

        run_param_obj = RunningParameters.objects.get(run_name_id=run_used_in_project)
        # get the number of lanes by quering the SequencerModel in the RunProcess
        # number_of_lanes = project_to_compare.runprocess_id.get_sequencing_lanes()
        number_of_lanes = int(run_param_obj.get_number_of_lanes())
        for lane_number in range(1, number_of_lanes + 1):
            try:
                lane_in_project = StatsLaneSummary.objects.get(
                    project_id=project_to_compare, lane__exact=lane_number
                )
            except Exception:
                continue
            (
                q_30_value,
                mean_q_value,
                yield_mb,
                cluster_pf,
            ) = lane_in_project.get_stats_info()
            q30_compare_lib.append(float(q_30_value))
            mean_compare_lib.append(float(mean_q_value))
            yield_mb_compare_lib.append(float(yield_mb.replace(",", "")))
        if library_to_compare_name in q30_comparations:
            q30_tmp_list = [
                float(q30_comparations[library_to_compare_name]),
                statistics.mean(q30_compare_lib),
            ]
            q30_comparations[library_to_compare_name] = format(
                statistics.mean(q30_tmp_list), ".2f"
            )
            mean_tmp_list = [
                float(mean_comparations[library_to_compare_name]),
                statistics.mean(mean_compare_lib),
            ]
            mean_comparations[library_to_compare_name] = format(
                statistics.mean(mean_tmp_list), ".2f"
            )
            n_bases_list = [
                float(n_bases_comparations[library_to_compare_name]),
                sum(yield_mb_compare_lib),
            ]

            n_bases_comparations[library_to_compare_name] = format(
                statistics.mean(n_bases_list), ".2f"
            )
        else:
            q30_comparations[library_to_compare_name] = format(
                statistics.mean(q30_compare_lib), ".2f"
            )
            mean_comparations[library_to_compare_name] = format(
                statistics.mean(mean_compare_lib), ".2f"
            )
            n_bases_comparations[library_to_compare_name] = format(
                statistics.mean(yield_mb_compare_lib), ".2f"
            )


@login_required
def stats_per_library(request):
    if request.method == "POST":
        library_kit_name = request.POST["libraryKitName"]
        start_date = request.POST["startdate"]
        end_date = request.POST["enddate"]
        # check that some values are in the request if not return the form
        if library_kit_name == "" and start_date == "" and end_date == "":
            return render(request, "wetlab/StatsPerLibrary.html")

        if library_kit_name != "" and len(library_kit_name) < 5:
            error_message = wetlab.config.ERROR_TOO_SHORT_INDEX_LIBRAY_NAME
            return render(
                request,
                "wetlab/StatsPerLibrary.html",
                {"error_message": error_message},
            )

        # check the right format of start and end date
        if start_date != "":
            if not check_valid_date_format(start_date):
                error_message = wetlab.config.ERROR_INVALID_FORMAT_FOR_DATES
                return render(
                    request,
                    "wetlab/StatsPerLibrary.html",
                    {"error_message": error_message},
                )
        if end_date != "":
            if not check_valid_date_format(end_date):
                error_message = wetlab.config.ERROR_INVALID_FORMAT_FOR_DATES
                return render(
                    request,
                    "wetlab/StatsPerLibrary.html",
                    {"error_message": error_message},
                )

        if library_kit_name != "":
            if Projects.objects.filter(
                library_kit__icontains=library_kit_name,
                run_process__state__run_state_name__exact="Completed",
            ).exists():
                library_found = Projects.objects.filter(
                    library_kit__icontains=library_kit_name,
                    run_process__state__run_state_name__exact="Completed",
                )
            else:
                error_message = wetlab.config.ERROR_NO_MATCHES_FOR_LIBRARY_STATISTICS
                return render(
                    request,
                    "wetlab/StatsPerLibrary.html",
                    {"error_message": error_message},
                )
        else:
            library_found = Projects.objects.filter(
                run_process__state__run_state_name__exact="Completed"
            )
        if start_date != "" and end_date != "":
            if library_found.filter(
                project_run_date__range=(start_date, end_date)
            ).exists():
                library_found = library_found.filter(
                    project_run_date__range=(start_date, end_date)
                )
            else:
                error_message = wetlab.config.ERROR_NO_MATCHES_FOR_LIBRARY_STATISTICS
                return render(
                    request,
                    "wetlab/StatsPerLibrary.html",
                    {"error_message": error_message},
                )
        if start_date != "" and end_date == "":
            if library_found.filter(project_run_date__gte=start_date).exists():
                library_found = library_found.filter(project_run_date__gte=start_date)
                #
            else:
                error_message = wetlab.config.ERROR_NO_MATCHES_FOR_LIBRARY_STATISTICS
                return render(
                    request,
                    "wetlab/StatsPerLibrary.html",
                    {"error_message": error_message},
                )
        if start_date == "" and end_date != "":
            if library_found.filter(project_run_date__lte=end_date).exists():
                #
                library_found = library_found.filter(project_run_date__lte=end_date)
            else:
                error_message = wetlab.config.ERROR_NO_MATCHES_FOR_LIBRARY_STATISTICS
                return render(
                    request,
                    "wetlab/StatsPerLibrary.html",
                    {"error_message": error_message},
                )

        # Collecting the statistics for the selected library
        # Get the projects which are using the library kit
        #
        library_stats = {}
        projects_name_in_library = []
        # Getting 1 library. Library could be in several projects. Information is collected per lane and by project
        # check if only 1 library kit matches the query
        library_names = {}
        for library in library_found:
            library_names[library.get_index_library_name()] = 1
        #
        if len(library_names) == 1:
            # There is only 1 library in the query. Results displays all projects data which have this library kit
            mean_lane_graphic = {}
            for project in library_found:
                projects_name_in_library.append(project.get_project_name())
            q30_in_lib, mean_in_lib, yield_mb_in_lib = [], [], []
            for lane_number in range(1, 5):
                q_30_lane, mean_q_lane, yield_mb_lane = {}, {}, {}
                for project in library_found:
                    project_id = project.get_project_id()
                    # Get quality information for each Lane summary of the project id
                    #
                    lane_in_project = StatsLaneSummary.objects.get(
                        project_id__exact=project_id, lane__exact=lane_number
                    )
                    (
                        q_30_value,
                        mean_q_value,
                        yield_mb,
                        cluster_pf,
                    ) = lane_in_project.get_stats_info()
                    project_name = project.get_project_name()
                    q_30_lane[project_name] = q_30_value
                    q30_in_lib.append(float(q_30_value))
                    mean_q_lane[project_name] = mean_q_value
                    mean_in_lib.append(float(mean_q_value))
                    yield_mb_lane[project_name] = yield_mb.replace(",", "")
                    yield_mb_in_lib.append(float(yield_mb.replace(",", "")))
                    #
                # creating the Yield MBases graphics
                chart_number = "chart-" + str(lane_number)
                render_number = "ex" + str(lane_number)
                heading = "Number of MBases in the projects for Lane " + str(
                    lane_number
                )
                data_source = graphic_for_library_kit(
                    heading,
                    "projects in lane ",
                    "Project Names",
                    "Number of M bases",
                    "ocean",
                    yield_mb_lane,
                )
                yield_mb_lane_graphic = core.fusioncharts.fusioncharts.FusionCharts(
                    "column3d",
                    render_number,
                    "500",
                    "300",
                    chart_number,
                    "json",
                    data_source,
                )
                #
                yield_graphic = "yield_mb_graphic" + str(lane_number)
                library_stats[yield_graphic] = yield_mb_lane_graphic.render()

                # creating the Q30 graphics
                chart_number = "q30-chart-" + str(lane_number)
                render_number = "q30-ex" + str(lane_number)
                heading = "Percent of bases > Q30 in the projects for Lane " + str(
                    lane_number
                )
                data_source = graphic_for_library_kit(
                    heading,
                    "projects in lane ",
                    "Project Names",
                    "Percent of Q 30",
                    "zune",
                    q_30_lane,
                )
                q30_lane_graphic = core.fusioncharts.fusioncharts.FusionCharts(
                    "column3d",
                    render_number,
                    "400",
                    "300",
                    chart_number,
                    "json",
                    data_source,
                )
                #
                q30_graphic = "q30_graphic" + str(lane_number)
                library_stats[q30_graphic] = q30_lane_graphic.render()

                # creating the Mean graphics
                chart_number = "mean-chart-" + str(lane_number)
                render_number = "mean-ex" + str(lane_number)
                heading = "Mean Quality Score in the projects for Lane " + str(
                    lane_number
                )
                data_source = graphic_for_library_kit(
                    heading,
                    "projects in lane ",
                    "Project Names",
                    "Percent of Q 30",
                    "carbon",
                    mean_q_lane,
                )
                mean_lane_graphic = core.fusioncharts.fusioncharts.FusionCharts(
                    "column3d",
                    render_number,
                    "400",
                    "300",
                    chart_number,
                    "json",
                    data_source,
                )
                #
                mean_graphic = "mean_graphic" + str(lane_number)
                library_stats[mean_graphic] = mean_lane_graphic.render()

            library_name = project.get_index_library_name()
            library_stats["library_name"] = library_name
            library_stats["project_names"] = projects_name_in_library
            #
            #########
            # set the data for the library under study
            #########
            q30_comparations, mean_comparations, n_bases_comparations = {}, {}, {}
            q30_comparations[library_name] = format(statistics.mean(q30_in_lib), ".2f")
            mean_comparations[library_name] = format(
                statistics.mean(mean_in_lib), ".2f"
            )
            n_bases_comparations[library_name] = format(
                statistics.mean(yield_mb_in_lib), ".2f"
            )
            error_in_library_to_compare = ""
            # get the data for the libraries to compare with
            if start_date == "" and end_date == "":
                if (
                    Projects.objects.filter(
                        runprocess_id__state__run_state_name__exact="Completed"
                    )
                    .exclude(libraryKit__exact=library_name)
                    .exists()
                ):
                    libraries_to_compare = Projects.objects.filter(
                        runprocess_id__state__run_state_name__exact="Completed"
                    ).exclude(library_kit__exact=library_name)
                else:
                    error_in_library_to_compare = (
                        "No other library have been found for doing the comparison. "
                    )

            if start_date != "" and end_date == "":
                if (
                    Projects.objects.filter(
                        runprocess_id__state__run_state_name__exact="Completed",
                        generate_dat__gte=start_date,
                    )
                    .exclude(library_kit__exact=library_name)
                    .exists()
                ):
                    libraries_to_compare = Projects.objects.filter(
                        runprocess_id__state__run_state_name__exact="Completed",
                        generate_dat__gte=start_date,
                    ).exclude(library_kit__exact=library_name)
                else:
                    error_in_library_to_compare = (
                        "No other library have been found for doing the comparison, with the starting date  "
                        + start_date
                    )

            if start_date == "" and end_date != "":
                if (
                    Projects.objects.filter(
                        runprocess_id__state__run_state_name__exact="Completed",
                        generate_dat__lte=end_date,
                    )
                    .exclude(library_kit__exact=library_name)
                    .exists()
                ):
                    libraries_to_compare = Projects.objects.filter(
                        runprocess_id__state__run_state_name__exact="Completed",
                        generate_dat__lte=end_date,
                    ).exclude(library_kit__exact=library_name)
                else:
                    error_in_library_to_compare = (
                        "No other library have been found for doing the comparison ending with  "
                        + end_date
                    )

            if start_date != "" and end_date != "":
                if (
                    Projects.objects.filter(
                        runprocess_id__state__run_state_name__exact="Completed",
                        generatedat__range=(start_date, end_date),
                    )
                    .exclude(library_kit__exact=library_name)
                    .exists()
                ):
                    libraries_to_compare = Projects.objects.filter(
                        runprocess_id__state__run_state_name__exact="Completed",
                        generate_dat__range=(start_date, end_date),
                    ).exclude(library_kit__exact=library_name)
                else:
                    error_in_library_to_compare = (
                        "No other library have been found for doing the comparison for the start date  "
                        + start_date
                        + "  and with the ending date  "
                        + end_date
                    )

            if error_in_library_to_compare == "":
                for project_to_compare in libraries_to_compare:
                    library_to_compare_name = (
                        project_to_compare.get_index_library_name()
                    )
                    project_to_compare_id = project_to_compare.get_project_id()
                    # q_30_lane , mean_q_lane , yield_mb_lane = {} , {} ,{}
                    q30_compare_lib, mean_compare_lib, yield_mb_compare_lib = [], [], []

                    run_obj = project_to_compare.get_run_obj()
                    run_param_obj = RunningParameters.objects.get(run_id=run_obj)
                    lanes_in_sequencer = int(run_param_obj.get_number_of_lanes())
                    for lane_number in range(1, lanes_in_sequencer + 1):
                        lane_in_project = StatsLaneSummary.objects.get(
                            project_id__exact=project_to_compare_id,
                            lane__exact=lane_number,
                        )
                        (
                            q_30_value,
                            mean_q_value,
                            yield_mb,
                            cluster_pf,
                        ) = lane_in_project.get_stats_info()
                        q30_compare_lib.append(float(q_30_value))
                        mean_compare_lib.append(float(mean_q_value))
                        yield_mb_compare_lib.append(float(yield_mb.replace(",", "")))
                    if library_to_compare_name in q30_comparations:
                        q30_tmp_list = [
                            float(q30_comparations[library_to_compare_name]),
                            statistics.mean(q30_compare_lib),
                        ]
                        q30_comparations[library_to_compare_name] = format(
                            statistics.mean(q30_tmp_list), ".2f"
                        )
                        mean_tmp_list = [
                            float(mean_comparations[library_to_compare_name]),
                            statistics.mean(mean_compare_lib),
                        ]
                        mean_comparations[library_to_compare_name] = format(
                            statistics.mean(mean_tmp_list), ".2f"
                        )
                        n_bases_list = [
                            float(n_bases_comparations[library_to_compare_name]),
                            statistics.mean(yield_mb_compare_lib),
                        ]
                        n_bases_comparations[library_to_compare_name] = format(
                            statistics.mean(n_bases_list), ".2f"
                        )
                    else:
                        q30_comparations[library_to_compare_name] = format(
                            statistics.mean(q30_compare_lib), ".2f"
                        )
                        mean_comparations[library_to_compare_name] = format(
                            statistics.mean(mean_compare_lib), ".2f"
                        )
                        n_bases_comparations[library_to_compare_name] = format(
                            statistics.mean(yield_mb_compare_lib), ".2f"
                        )

            else:
                library_stats["error_library"] = error_in_library_to_compare

            heading = "Comparison of Percent of bases > Q30  "
            data_source = graphic_for_library_kit(
                heading,
                "Q30 comparison ",
                "Library Names",
                "Percent of Q 30",
                "",
                q30_comparations,
            )
            comp_q30_lib_graphic = core.fusioncharts.fusioncharts.FusionCharts(
                "column3d",
                "comp-q30-1",
                "500",
                "300",
                "comp-q30-chart-1",
                "json",
                data_source,
            )
            library_stats["comp_q30_graphic"] = comp_q30_lib_graphic.render()

            heading = "Comparison of Mean Quality Score "
            data_source = graphic_for_library_kit(
                heading,
                "Mean Quality Score comparison ",
                "Library Names",
                "Mean Quality Score",
                "",
                mean_comparations,
            )
            comp_mean_lib_graphic = core.fusioncharts.fusioncharts.FusionCharts(
                "column3d",
                "comp-mean-1",
                "500",
                "300",
                "comp-mean-chart-1",
                "json",
                data_source,
            )
            library_stats["comp_mean_graphic"] = comp_mean_lib_graphic.render()

            heading = "Number of Bases comparison"
            data_source = graphic_for_library_kit(
                heading,
                "Number of Bases comparison ",
                "Library Names",
                "Number of Bases ",
                "",
                n_bases_comparations,
            )
            comp_mean_lib_graphic = core.fusioncharts.fusioncharts.FusionCharts(
                "column3d",
                "comp-n_bases-1",
                "500",
                "300",
                "comp-n_bases-chart-1",
                "json",
                data_source,
            )
            library_stats["comp_n_bases_graphic"] = comp_mean_lib_graphic.render()

            return render(
                request,
                "wetlab/StatsPerLibrary.html",
                {"display_library_stats": library_stats},
            )
        else:
            library_list_stats = {}
            libraries_found_name = []
            # get the library names that match with the searching criteria
            for library in library_found:
                lib_name = library.get_index_library_name()
                if lib_name not in libraries_found_name:
                    libraries_found_name.append(lib_name)
            #
            library_list_stats["library_names"] = libraries_found_name
            q30_comparations, mean_comparations, n_bases_comparations = {}, {}, {}
            #
            # get the data for displaying the libraries found in the form request
            #
            get_list_of_libraries_values(
                library_found, q30_comparations, mean_comparations, n_bases_comparations
            )

            heading = "Comparison of Percent of bases > Q30  "
            data_source = graphic_for_library_kit(
                heading,
                "Q30 comparison ",
                "Library Names",
                "Percent of Q 30",
                "",
                q30_comparations,
            )
            comp_q30_lib_graphic = core.fusioncharts.fusioncharts.FusionCharts(
                "column3d",
                "comp-q30-1",
                "500",
                "300",
                "comp-q30-chart-1",
                "json",
                data_source,
            )
            #
            library_list_stats["comp_q30_graphic"] = comp_q30_lib_graphic.render()

            heading = "Comparison of Mean Quality Score "
            data_source = graphic_for_library_kit(
                heading,
                "Mean Quality Score comparison ",
                "Library Names",
                "Mean Quality Score",
                "",
                mean_comparations,
            )
            comp_mean_lib_graphic = core.fusioncharts.fusioncharts.FusionCharts(
                "column3d",
                "comp-mean-1",
                "500",
                "300",
                "comp-mean-chart-1",
                "json",
                data_source,
            )
            #
            library_list_stats["comp_mean_graphic"] = comp_mean_lib_graphic.render()

            heading = "Number of Bases comparison"
            data_source = graphic_for_library_kit(
                heading,
                "Number of Bases comparison ",
                "Library Names",
                "Number of Bases ",
                "",
                n_bases_comparations,
            )
            comp_mean_lib_graphic = core.fusioncharts.fusioncharts.FusionCharts(
                "column3d",
                "comp-n_bases-1",
                "500",
                "300",
                "comp-n_bases-chart-1",
                "json",
                data_source,
            )
            #
            library_list_stats["comp_n_bases_graphic"] = comp_mean_lib_graphic.render()
            #
            # get data for displaying the libraries found in the request form
            #
            all_libraries = Projects.objects.filter(
                run_process__state__run_state_name__exact="Completed"
            )
            if start_date != "" and end_date != "":
                if all_libraries.filter(
                    generate_dat__range=(start_date, end_date)
                ).exists():
                    library_found = library_found.filter(
                        generate_dat__range=(start_date, end_date)
                    )
            if start_date != "" and end_date == "":
                if all_libraries.filter(generate_dat__gte=start_date).exists():
                    all_libraries = library_found.filter(generate_dat__gte=start_date)
            if start_date == "" and end_date != "":
                if all_libraries.filter(generate_dat__lte=end_date).exists():
                    all_libraries = library_found.filter(generate_dat__lte=end_date)

            q30_comparations, mean_comparations, n_bases_comparations = {}, {}, {}
            get_list_of_libraries_values(
                all_libraries, q30_comparations, mean_comparations, n_bases_comparations
            )
            #
            heading = "Library kits of Percent of bases > Q30  "
            data_source = graphic_for_library_kit(
                heading,
                "Q30 library kits ",
                "Library Names",
                "Percent of Q 30",
                "",
                q30_comparations,
            )
            lib_q30_lib_graphic = core.fusioncharts.fusioncharts.FusionCharts(
                "column3d",
                "lib-q30-lib",
                "500",
                "300",
                "lib-q30-chart-1",
                "json",
                data_source,
            )
            #
            library_list_stats["lib_q30_graphic"] = lib_q30_lib_graphic.render()

            heading = "Library kits of Mean Quality Score "
            data_source = graphic_for_library_kit(
                heading,
                "Mean Quality Score Library kits ",
                "Library Names",
                "Mean Quality Score",
                "",
                mean_comparations,
            )
            lib_mean_lib_graphic = core.fusioncharts.fusioncharts.FusionCharts(
                "column3d",
                "lib-mean-lib",
                "500",
                "300",
                "lib-mean-chart-1",
                "json",
                data_source,
            )
            #
            library_list_stats["lib_mean_graphic"] = lib_mean_lib_graphic.render()

            heading = "Number of Bases per Library kits"
            data_source = graphic_for_library_kit(
                heading,
                "Number of Bases per Library kits ",
                "Library Names",
                "Number of Bases ",
                "",
                n_bases_comparations,
            )
            lib_mean_lib_graphic = core.fusioncharts.fusioncharts.FusionCharts(
                "column3d",
                "lib-n_bases-lib",
                "500",
                "300",
                "lib-n_bases-chart-1",
                "json",
                data_source,
            )
            #
            library_list_stats["lib_n_bases_graphic"] = lib_mean_lib_graphic.render()

            # Create the graphic for number of time that library has been used

            count_libraries = {}
            for library_used in all_libraries:
                lib_name = library_used.get_index_library_name()
                if lib_name not in count_libraries:
                    count_libraries[lib_name] = 1
                else:
                    count_libraries[lib_name] += 1
            data_source = pie_graphic(
                "Library utilization in projects", "fint", count_libraries
            )
            libraries_kit_utilization = core.fusioncharts.fusioncharts.FusionCharts(
                "pie3d",
                "lib_kit_utilization_graph-1",
                "500",
                "400",
                "lib_kit_utilization_chart-1",
                "json",
                data_source,
            )
            library_list_stats[
                "libraries_kit_utilization"
            ] = libraries_kit_utilization.render()

            return render(
                request,
                "wetlab/StatsPerLibrary.html",
                {"display_list_of_library_stats": library_list_stats},
            )

    else:
        return render(request, "wetlab/StatsPerLibrary.html")


@login_required
def annual_report(request):
    # check user privileges
    if request.user.is_authenticated:
        try:
            groups = Group.objects.get(name="WetlabManager")
            if groups not in request.user.groups.all():
                return render(
                    request,
                    "wetlab/error_page.html",
                    {
                        "content": [
                            "You do have the enough privileges to see this page ",
                            "Contact with your administrator .",
                        ]
                    },
                )
        except Exception:
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "You do have the enough privileges to see this page ",
                        "Contact with your administrator .",
                    ]
                },
            )
    else:
        # redirect to login webpage
        return redirect("/accounts/login")

    if request.method == "POST":
        year_selected = int(request.POST["yearselected"])
        # get the current year to compare with the input
        present_year = datetime.now().year
        if year_selected > present_year:
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "Annual Report cannot be done on the future  ",
                        "the input year in the Form  ",
                        year_selected,
                        "is not allowed",
                    ]
                },
            )

        completed_run_in_year = RunProcess.objects.filter(
            run_date__year=year_selected, state__run_state_name__exact="Completed"
        )
        #
        uncompleted_run_in_year = RunProcess.objects.filter(
            run_date__year=year_selected
        ).exclude(state__run_state_name__exact="Completed")
        if len(completed_run_in_year) == 0 and len(uncompleted_run_in_year) == 0:
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "Annual Report cannot be generated because there is no runs performed the year ",
                        year_selected,
                    ]
                },
            )

        annual_report_information = {}
        annual_report_information["year"] = year_selected
        number_of_runs = {}
        number_of_runs["Completed Runs"] = 0
        number_of_runs["Not Finish Runs"] = 0
        if len(completed_run_in_year) > 0:
            completed_run = []
            for run in completed_run_in_year:
                completed_run.append(run.get_run_name)
            annual_report_information["completed_run"] = completed_run
            number_of_runs["Completed Runs"] = len(completed_run_in_year)
        if len(uncompleted_run_in_year) > 0:
            uncompleted_run = []
            for run_uncompleted in uncompleted_run_in_year:
                uncompleted_run.append(run_uncompleted.get_run_name)
            annual_report_information["uncompleted_run"] = uncompleted_run
            number_of_runs["Not Finish Runs"] = len(uncompleted_run_in_year)
        # prepare the pie graphic for the number of completed/ unfinished runs
        data_source = pie_graphic_standard(
            "Number of Runs performed on the year", "", "ocean", number_of_runs
        )
        graphic_completed_run = core.fusioncharts.fusioncharts.FusionCharts(
            "pie3d", "ex1", "400", "300", "chart-1", "json", data_source
        )
        annual_report_information[
            "graphic_completed_run"
        ] = graphic_completed_run.render()

        #
        # Collecting information from StatsRunSummary
        run_found_bin_summary_year = StatsRunSummary.objects.filter(
            stats_summary_run_date__year=year_selected, level__exact="Total"
        )
        q30_year, aligned_year, error_rate_year = {}, {}, {}
        for run_bin_summary in run_found_bin_summary_year:
            bin_summary_data = run_bin_summary.get_bin_run_summary().split(";")
            run_name = run_bin_summary.runprocess_id.get_run_name()
            aligned_year[run_name] = bin_summary_data[2]
            error_rate_year[run_name] = bin_summary_data[3]
            q30_year[run_name] = bin_summary_data[5]
        annual_report_information["aligned_data"] = aligned_year
        annual_report_information["error_rate_data"] = error_rate_year
        annual_report_information["q30_data"] = q30_year
        # graphics for StatsRunSummary
        heading = "Aligned % for the runs done on year " + str(year_selected)
        data_source = column_graphic_for_year_report(
            heading, "Aligned  ", "Run names ", "Aligned (in %)", "ocean", aligned_year
        )
        aligned_year_graphic = core.fusioncharts.fusioncharts.FusionCharts(
            "column3d",
            "aligned_year",
            "600",
            "300",
            "aligned_chart-3",
            "json",
            data_source,
        )
        annual_report_information["aligned_graphic"] = aligned_year_graphic.render()

        heading = "Error Rate for the runs done on year " + str(year_selected)
        data_source = column_graphic_for_year_report(
            heading,
            "Error rate ",
            "Run names ",
            "Error rate",
            "carbon",
            error_rate_year,
        )
        error_rate_year_graphic = core.fusioncharts.fusioncharts.FusionCharts(
            "column3d",
            "error_rate_year",
            "600",
            "300",
            "error_rate_chart-4",
            "json",
            data_source,
        )
        annual_report_information[
            "error_rate_graphic"
        ] = error_rate_year_graphic.render()

        heading = ">Q30 for the runs done on year " + str(year_selected)
        data_source = column_graphic_for_year_report(
            heading, "Q30  ", "Run names ", ">Q 30 (in %)", "fint", q30_year
        )
        q30_year_graphic = core.fusioncharts.fusioncharts.FusionCharts(
            "column3d", "q30_year", "600", "300", "q30_chart-2", "json", data_source
        )
        #
        annual_report_information["q30_graphic"] = q30_year_graphic.render()
        #

        # Get the information for investigator name and the projects done
        # number_proyects_investigator contains a dict with 3 ranges 1-5, 6-10, more than 11
        investigator_projects = Projects.objects.filter(
            project_run_date__year=year_selected
        ).order_by("user_id")
        project_by_user = {}
        (
            investigator_5_project,
            investigator_10_project,
            investigator_more_10_project,
        ) = ({}, {}, {})
        #
        for investigator in investigator_projects:
            user_name = investigator.get_user_name()
            if user_name in project_by_user:
                project_by_user[user_name].append(investigator.get_project_name())
            else:
                project_by_user[user_name] = [investigator.get_project_name()]
        for key, value in project_by_user.items():
            if len(value) <= 5:
                investigator_5_project[key] = value
            elif len(value) <= 10:
                investigator_10_project[key] = value
            else:
                investigator_more_10_project[key] = value
        annual_report_information["user_5_projects"] = investigator_5_project
        annual_report_information["user_10_projects"] = investigator_10_project
        annual_report_information[
            "user_more_10_projects"
        ] = investigator_more_10_project

        # Create the bar graphic for user projects
        p_user_year = {}
        p_user_year["1 - 5"] = len(investigator_5_project)
        p_user_year["6 - 10"] = len(investigator_10_project)
        p_user_year["more than 10"] = len(investigator_more_10_project)
        heading = "Projects done per investigator on year " + str(year_selected)
        data_source = column_graphic_for_year_report(
            heading, "  ", "Projects ", "number of users", "ocean", p_user_year
        )
        p_user_year_graphic = core.fusioncharts.fusioncharts.FusionCharts(
            "column3d",
            "bar_project_user_year",
            "400",
            "300",
            "p_user_chart-1",
            "json",
            data_source,
        )
        annual_report_information["p_user_year_graphic"] = p_user_year_graphic.render()

        data_source = pie_graphic_standard(heading, "Percentage", "carbon", p_user_year)
        pie_p_user_year_graphic = core.fusioncharts.fusioncharts.FusionCharts(
            "pie3d",
            "pie_project_user_year",
            "400",
            "300",
            "p_user_chart-2",
            "json",
            data_source,
        )
        annual_report_information[
            "pie_p_user_year_graphic"
        ] = pie_p_user_year_graphic.render()
        #
        return render(
            request,
            "wetlab/AnnualReport.html",
            {"display_annual_report": annual_report_information},
        )
    else:
        return render(request, "wetlab/AnnualReport.html")


@login_required
def monthly_report(request):
    # check user privileges
    if request.user.is_authenticated:
        try:
            groups = Group.objects.get(name="WetlabManager")
            if groups not in request.user.groups.all():
                return render(
                    request,
                    "wetlab/error_page.html",
                    {
                        "content": [
                            "You do have the enough privileges to see this page ",
                            "Contact with your administrator .",
                        ]
                    },
                )
        except Exception:
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "You do have the enough privileges to see this page ",
                        "Contact with your administrator .",
                    ]
                },
            )
    else:
        # redirect to login webpage
        return redirect("/accounts/login")

    if request.method == "POST":
        input_value = request.POST["month_year_selected"]
        browser_used = request.META["HTTP_USER_AGENT"]
        if "Firefox" in browser_used:
            try:
                datetime.strptime(input_value, "%m-%Y")
                month_selected, year_selected = input_value.split("-")
            except Exception:
                return render(
                    request,
                    "wetlab/error_page.html",
                    {
                        "content": [
                            "Input field does not have the right format  ",
                            "the right input format is MM-YYYY   the entry ",
                            input_value,
                            " is not allowed",
                        ]
                    },
                )

        else:
            try:
                datetime.strptime(input_value, "%Y-%m")
                year_selected, month_selected = input_value.split("-")
            except Exception:
                return render(
                    request,
                    "wetlab/error_page.html",
                    {
                        "content": [
                            "Monthly Report input field does not have the right format  ",
                            "the right input format is MM-YYYY  ",
                            " the entry ",
                            input_value,
                            " is not allowed",
                        ]
                    },
                )

        # get the current year to compare with the input
        present_year = datetime.now().year
        #
        if int(year_selected) > present_year:
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "Monthly Report cannot be done on the future  ",
                        "the input year in the Form  ",
                        year_selected,
                        "is not allowed",
                    ]
                },
            )

        present_month = datetime.now().month
        if (int(year_selected) == present_year) and (
            int(month_selected) > present_month
        ):
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "Monthly Report cannot be done on the future  ",
                        "the input month in the Form  ",
                        month_selected,
                        "is not allowed",
                    ]
                },
            )

        completed_run_in_year_month = RunProcess.objects.filter(
            run_date__year=year_selected,
            run_date__month=month_selected,
            state__run_state_name__exact="Completed",
        )
        #
        uncompleted_run_in_year_month = RunProcess.objects.filter(
            run_date__year=year_selected, run_date__month=month_selected
        ).exclude(state__run_state_name__exact="Completed")
        if (
            len(completed_run_in_year_month) == 0
            and len(uncompleted_run_in_year_month) == 0
        ):
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "Montly Report cannot be generated because there is no runs performed the year ",
                        year_selected,
                    ]
                },
            )

        monthly_report_information = {}
        monthly_report_information["month_year"] = str(
            month_selected + "  " + year_selected
        )
        number_of_runs = {}
        number_of_runs["Completed Runs"] = 0
        number_of_runs["Not Finish Runs"] = 0
        if len(completed_run_in_year_month) > 0:
            completed_run = []
            for run in completed_run_in_year_month:
                completed_run.append(run.get_run_name)
            monthly_report_information["completed_run"] = completed_run
            number_of_runs["Completed Runs"] = len(completed_run_in_year_month)
        if len(uncompleted_run_in_year_month) > 0:
            uncompleted_run = []
            for run_uncompleted in uncompleted_run_in_year_month:
                uncompleted_run.append(run_uncompleted.get_run_name)
            monthly_report_information["uncompleted_run"] = uncompleted_run
            number_of_runs["Not Finish Runs"] = len(uncompleted_run_in_year_month)
        # prepare the pie graphic for the number of completed/ unfinished runs
        heading = str(
            "Graphics of the Runs performed on the "
            + month_selected
            + " - "
            + year_selected
        )
        data_source = pie_graphic_standard(heading, "", "ocean", number_of_runs)
        graphic_completed_run = core.fusioncharts.fusioncharts.FusionCharts(
            "pie3d", "ex1", "400", "300", "chart-1", "json", data_source
        )
        monthly_report_information[
            "graphic_completed_run"
        ] = graphic_completed_run.render()

        # Get the information for investigator name and the projects done
        # number_proyects_investigator contains a dict with 3 ranges 1, 2, more than 2
        investigator_projects = Projects.objects.filter(
            project_run_date__year=year_selected, project_run_date__month=month_selected
        ).order_by("user_id")
        project_by_user = {}
        (
            investigator_1_project,
            investigator_2_projects,
            investigator_more_2_projects,
        ) = ({}, {}, {})
        #
        for investigator in investigator_projects:
            user_name = investigator.get_user_name()
            if user_name in project_by_user:
                project_by_user[user_name].append(investigator.get_project_name())
            else:
                project_by_user[user_name] = [investigator.get_project_name()]
        for key, value in project_by_user.items():
            if len(value) == 1:
                investigator_1_project[key] = value
            elif len(value) == 2:
                investigator_2_projects[key] = value
            else:
                investigator_more_2_projects[key] = value
        monthly_report_information["user_1_project"] = investigator_1_project
        monthly_report_information["user_2_projects"] = investigator_2_projects
        monthly_report_information[
            "user_more_2_projects"
        ] = investigator_more_2_projects

        # Create the bar graphic for user projects
        p_user_month = {}
        p_user_month["1 project"] = len(investigator_1_project)
        p_user_month["2 projects"] = len(investigator_2_projects)
        p_user_month["more than 2"] = len(investigator_more_2_projects)

        heading = "Projects done per investigator on " + str(
            month_selected + " - " + year_selected
        )
        data_source = column_graphic_for_year_report(
            heading, "  ", "Projects ", "number of users", "ocean", p_user_month
        )
        p_user_monthly_graphic = core.fusioncharts.fusioncharts.FusionCharts(
            "column3d",
            "bar_project_user_month",
            "400",
            "300",
            "p_user_chart-1",
            "json",
            data_source,
        )
        monthly_report_information[
            "p_user_monthly_graphic"
        ] = p_user_monthly_graphic.render()

        data_source = pie_graphic_standard(
            heading, "Percentage", "carbon", p_user_month
        )
        pie_p_user_monthly_graphic = core.fusioncharts.fusioncharts.FusionCharts(
            "pie3d",
            "pie_project_user_month",
            "400",
            "300",
            "p_user_chart-2",
            "json",
            data_source,
        )
        monthly_report_information[
            "pie_p_user_monthly_graphic"
        ] = pie_p_user_monthly_graphic.render()

        # Collecting information from StatsRunSummary
        run_found_bin_summary_month = StatsRunSummary.objects.filter(
            stats_summary_run_date__year=year_selected,
            stats_summary_run_date__month=month_selected,
            level__exact="Total",
        )
        q30_month, aligned_month, error_rate_month = {}, {}, {}
        for run_bin_summary in run_found_bin_summary_month:
            bin_summary_data = run_bin_summary.get_bin_run_summary().split(";")
            run_name = run_bin_summary.runprocess_id.get_run_name()
            aligned_month[run_name] = bin_summary_data[2]
            error_rate_month[run_name] = bin_summary_data[3]
            q30_month[run_name] = bin_summary_data[5]
        monthly_report_information["aligned_data"] = aligned_month
        monthly_report_information["error_rate_data"] = error_rate_month
        monthly_report_information["q30_data"] = q30_month
        # graphics for StatsRunSummary
        heading = "Aligned % for the runs done on " + str(
            month_selected + " - " + year_selected
        )
        data_source = column_graphic_for_year_report(
            heading, "Aligned  ", "Run names ", "Aligned (in %)", "ocean", aligned_month
        )
        aligned_month_graphic = core.fusioncharts.fusioncharts.FusionCharts(
            "column3d",
            "aligned_year",
            "600",
            "300",
            "aligned_chart-3",
            "json",
            data_source,
        )
        monthly_report_information["aligned_graphic"] = aligned_month_graphic.render()

        heading = "Error Rate for the runs done on  " + str(
            month_selected + " - " + year_selected
        )
        data_source = column_graphic_for_year_report(
            heading,
            "Error rate ",
            "Run names ",
            "Error rate",
            "carbon",
            error_rate_month,
        )
        error_rate_month_graphic = core.fusioncharts.fusioncharts.FusionCharts(
            "column3d",
            "error_rate_year",
            "600",
            "300",
            "error_rate_chart-4",
            "json",
            data_source,
        )
        monthly_report_information[
            "error_rate_graphic"
        ] = error_rate_month_graphic.render()

        heading = ">Q30 for the runs done on  " + str(
            month_selected + " - " + year_selected
        )
        data_source = column_graphic_for_year_report(
            heading, "Q30  ", "Run names ", ">Q 30 (in %)", "fint", q30_month
        )
        q30_month_graphic = core.fusioncharts.fusioncharts.FusionCharts(
            "column3d", "q30_year", "600", "300", "q30_chart-2", "json", data_source
        )

        monthly_report_information["q30_graphic"] = q30_month_graphic.render()

        return render(
            request,
            "wetlab/MonthlyReport.html",
            {"display_monthly_report": monthly_report_information},
        )
    else:
        return render(request, "wetlab/MonthlyReport.html")


@login_required
def quarter_report(request):
    # check user privileges
    if request.user.is_authenticated:
        try:
            groups = Group.objects.get(name="WetlabManager")
            if groups not in request.user.groups.all():
                return render(
                    request,
                    "wetlab/error_page.html",
                    {
                        "content": [
                            "You do have the enough privileges to see this page ",
                            "Contact with your administrator .",
                        ]
                    },
                )
        except Exception:
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "You do have the enough privileges to see this page ",
                        "Contact with your administrator .",
                    ]
                },
            )
    else:
        # redirect to login webpage
        return redirect("/accounts/login")

    if request.method == "POST":
        year_selected = request.POST["yearselected"]
        quarter_selected = int(request.POST["quarter"])
        quarter_string = [
            "",
            "First Quarter (January -- March) ",
            "Second Quarter (April -- June) ",
            "Third Quarter (July -- September) ",
            "Fourth Quarter (October -- Decemmber) ",
        ]
        days_in_end_quarter = ["0", "31", "30", "30", "31"]
        start_quarter = str(quarter_selected * 3 - 2)
        end_quarter = str(quarter_selected * 3)
        start_date = str(year_selected + "-" + start_quarter + "-01")
        end_date = str(
            year_selected
            + "-"
            + end_quarter
            + "-"
            + days_in_end_quarter[quarter_selected]
        )
        # get the current year to compare with the input
        present_year = datetime.now().year
        if int(year_selected) > present_year:
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "Quarter Report cannot be done on the future  ",
                        "the input year in the Form  ",
                        year_selected,
                        "is not allowed",
                    ]
                },
            )

        present_month = datetime.now().month
        if (int(year_selected) == present_year) and (int(end_quarter) > present_month):
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "Quater Report cannot be done on the future  ",
                        "the selected Quarter ",
                        quarter_string[quarter_selected] + str(year_selected),
                        "is not allowed",
                    ]
                },
            )

        #
        completed_run_in_quarter = RunProcess.objects.filter(
            run_date__range=(start_date, end_date), state__run_state_name="Completed"
        )
        #
        uncompleted_run_in_quarter = RunProcess.objects.filter(
            run_date__range=(start_date, end_date)
        ).exclude(state__run_state_name="Completed")
        if len(completed_run_in_quarter) == 0 and len(uncompleted_run_in_quarter) == 0:
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "Quater Report cannot be generated because there is no runs performed the Quarter ",
                        quarter_string[quarter_selected] + str(year_selected),
                    ]
                },
            )

        quarter_report_information = {}
        quarter_report_information["quarter_year"] = quarter_string[
            quarter_selected
        ] + str(year_selected)
        number_of_runs = {}
        number_of_runs["Completed Runs"] = 0
        number_of_runs["Not Finish Runs"] = 0
        if len(completed_run_in_quarter) > 0:
            completed_run = []
            for run in completed_run_in_quarter:
                completed_run.append(run.get_run_name)
            quarter_report_information["completed_run"] = completed_run
            number_of_runs["Completed Runs"] = len(completed_run_in_quarter)
        if len(uncompleted_run_in_quarter) > 0:
            uncompleted_run = []
            for run_uncompleted in uncompleted_run_in_quarter:
                uncompleted_run.append(run_uncompleted.get_run_name)
            quarter_report_information["uncompleted_run"] = uncompleted_run
            number_of_runs["Not Finish Runs"] = len(uncompleted_run_in_quarter)
        # prepare the pie graphic for the number of completed/ unfinished runs
        data_source = pie_graphic_standard(
            "Number of Runs performed on the year", "", "ocean", number_of_runs
        )
        graphic_completed_run = core.fusioncharts.fusioncharts.FusionCharts(
            "pie3d", "ex1", "400", "300", "chart-1", "json", data_source
        )
        quarter_report_information[
            "graphic_completed_run"
        ] = graphic_completed_run.render()

        #
        # Collecting information from StatsRunSummary
        run_found_bin_summary_quarter = StatsRunSummary.objects.filter(
            stats_summary_run_date__range=(start_date, end_date), level__exact="Total"
        )
        q30_quarter, aligned_quarter, error_rate_quarter = {}, {}, {}
        for run_bin_summary in run_found_bin_summary_quarter:
            bin_summary_data = run_bin_summary.get_bin_run_summary().split(";")
            run_name = run_bin_summary.runprocess_id.get_run_name()
            aligned_quarter[run_name] = bin_summary_data[2]
            error_rate_quarter[run_name] = bin_summary_data[3]
            q30_quarter[run_name] = bin_summary_data[5]
        quarter_report_information["aligned_data"] = aligned_quarter
        quarter_report_information["error_rate_data"] = error_rate_quarter
        quarter_report_information["q30_data"] = q30_quarter
        # graphics for StatsRunSummary
        heading = (
            "Aligned % for the runs done on  "
            + quarter_string[quarter_selected]
            + str(year_selected)
        )
        data_source = column_graphic_for_year_report(
            heading,
            "Aligned  ",
            "Run names ",
            "Aligned (in %)",
            "ocean",
            aligned_quarter,
        )
        aligned_quarter_graphic = core.fusioncharts.fusioncharts.FusionCharts(
            "column3d",
            "aligned_year",
            "600",
            "300",
            "aligned_chart-3",
            "json",
            data_source,
        )
        quarter_report_information["aligned_graphic"] = aligned_quarter_graphic.render()

        heading = (
            "Error Rate for the runs done on year "
            + quarter_string[quarter_selected]
            + str(year_selected)
        )
        data_source = column_graphic_for_year_report(
            heading,
            "Error rate ",
            "Run names ",
            "Error rate",
            "carbon",
            error_rate_quarter,
        )
        error_rate_quarter_graphic = core.fusioncharts.fusioncharts.FusionCharts(
            "column3d",
            "error_rate_year",
            "600",
            "300",
            "error_rate_chart-4",
            "json",
            data_source,
        )
        quarter_report_information[
            "error_rate_graphic"
        ] = error_rate_quarter_graphic.render()

        heading = (
            ">Q30 for the runs done on year "
            + quarter_string[quarter_selected]
            + str(year_selected)
        )
        data_source = column_graphic_for_year_report(
            heading, "Q30  ", "Run names ", ">Q 30 (in %)", "fint", q30_quarter
        )
        q30_quarter_graphic = core.fusioncharts.fusioncharts.FusionCharts(
            "column3d", "q30_year", "600", "300", "q30_chart-2", "json", data_source
        )
        #
        quarter_report_information["q30_graphic"] = q30_quarter_graphic.render()
        #

        # Get the information for investigator name and the projects done
        # number_proyects_investigator contains a dict with 3 ranges 1-5, 6-10, more than 11
        investigator_projects = Projects.objects.filter(
            project_run_date__range=(start_date, end_date)
        ).order_by("user_id")
        project_by_user = {}
        (
            investigator_5_project,
            investigator_10_project,
            investigator_more_10_project,
        ) = ({}, {}, {})
        #
        for investigator in investigator_projects:
            user_name = investigator.get_user_name()
            if user_name in project_by_user:
                project_by_user[user_name].append(investigator.get_project_name())
            else:
                project_by_user[user_name] = [investigator.get_project_name()]
        for key, value in project_by_user.items():
            if len(value) <= 5:
                investigator_5_project[key] = value
            elif len(value) <= 10:
                investigator_10_project[key] = value
            else:
                investigator_more_10_project[key] = value
        quarter_report_information["user_5_projects"] = investigator_5_project
        quarter_report_information["user_10_projects"] = investigator_10_project
        quarter_report_information[
            "user_more_10_projects"
        ] = investigator_more_10_project

        # Create the bar graphic for user projects
        p_user_quarter = {}
        p_user_quarter["1 - 5"] = len(investigator_5_project)
        p_user_quarter["6 - 10"] = len(investigator_10_project)
        p_user_quarter["more than 10"] = len(investigator_more_10_project)
        heading = "Projects done per investigator on year " + str(year_selected)
        data_source = column_graphic_for_year_report(
            heading, "  ", "Projects ", "number of users", "ocean", p_user_quarter
        )
        p_user_quarter_graphic = core.fusioncharts.fusioncharts.FusionCharts(
            "column3d",
            "bar_project_user_year",
            "400",
            "300",
            "p_user_chart-1",
            "json",
            data_source,
        )
        quarter_report_information[
            "p_user_year_graphic"
        ] = p_user_quarter_graphic.render()

        data_source = pie_graphic_standard(
            heading, "Percentage", "carbon", p_user_quarter
        )
        pie_p_user_quarter_graphic = core.fusioncharts.fusioncharts.FusionCharts(
            "pie3d",
            "pie_project_user_year",
            "400",
            "300",
            "p_user_chart-2",
            "json",
            data_source,
        )
        quarter_report_information[
            "pie_p_user_year_graphic"
        ] = pie_p_user_quarter_graphic.render()
        #
        return render(
            request,
            "wetlab/QuarterReport.html",
            {"display_quarter_report": quarter_report_information},
        )
    else:
        return render(request, "wetlab/QuarterReport.html")


def get_size_dir(
    directory,
    conn,
):
    count_file_size = 0
    file_list = conn.listPath(wetlab.config.SAMBA_SHARED_FOLDER_NAME, directory)
    for sh_file in file_list:
        if sh_file.isDirectory:
            if sh_file.filename == "." or sh_file.filename == "..":
                continue

            sub_directory = os.path.join(directory, sh_file.filename)
            count_file_size += get_size_dir(sub_directory, conn)
        else:
            count_file_size += sh_file.file_size

    return count_file_size



@login_required
def create_protocol(request):
    # Check user == WETLAB_MANAGER: if false,  redirect to 'login' page
    if request.user.is_authenticated:
        if not wetlab.utils.common.is_wetlab_manager(request):
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "You do not have enough privileges to see this page ",
                        "Contact with your administrator .",
                    ]
                },
            )
    # get the list of defined protocols
    defined_protocols, other_protocol_list = display_available_protocols(__package__)
    additional_kits = wetlab.utils.additional_kits.get_additional_kits_list(__package__)
    defined_protocol_types = display_protocol_types(__package__)

    if request.method == "POST" and request.POST["action"] == "addNewProtocol":
        new_protocol = request.POST["newProtocolName"]
        protocol_type = request.POST["protocolType"]
        description = request.POST["description"]

        if check_if_protocol_exists(new_protocol, __package__):
            return render(
                request,
                "wetlab/create_protocol.html",
                {"ERROR": "Protocol Name " + new_protocol + "Already exists."}
            )
        new_protocol_id = create_new_protocol(
            new_protocol, protocol_type, description, __package__
        )

        return render(
            request,
            "wetlab/create_protocol.html",
            {
                "defined_protocols": defined_protocols,
                "defined_protocol_types": defined_protocol_types,
                "new_defined_protocol": new_protocol,
                "new_protocol_id": new_protocol_id,
                "other_protocol_list": other_protocol_list,
            },
        )

    return render(
        request,
        "wetlab/create_protocol.html",
        {
            "defined_protocols": defined_protocols,
            "defined_protocol_types": defined_protocol_types,
            "other_protocol_list": other_protocol_list,
            "additional_kits": additional_kits,
        },
    )


@login_required
def define_sample_projects(request):
    # Check user == WETLAB_MANAGER: if false,  redirect to 'login' page
    if request.user.is_authenticated:
        if not wetlab.utils.common.is_wetlab_manager(request):
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "You do not have enough privileges to see this page ",
                        "Contact with your administrator .",
                    ]
                },
            )
    # get the information of defined sample Projects
    defined_samples_projects = core.utils.samples.get_info_for_defined_sample_projects(__package__)

    if request.method == "POST" and request.POST["action"] == "addNewSampleProject":
        sample_project_name = request.POST["sampleProyectName"]
        # description = request.POST['description']

        if core.utils.samples.check_if_sample_project_exists(sample_project_name, __package__):
            error_message = wetlab.config.ERROR_SAMPLE_PROJECT_ALREADY_EXISTS
            return render(
                request,
                "wetlab/createSampleProjects.html",
                {
                    "defined_samples_projects": defined_samples_projects,
                    "error_message": error_message,
                },
            )
        new_sample_project_id = core.utils.samples.create_new_sample_project(request.POST, __package__)
        new_defined_sample_project = sample_project_name
        return render(
            request,
            "wetlab/createSampleProjects.html",
            {
                "defined_samples_projects": defined_samples_projects,
                "new_sample_project_id": new_sample_project_id,
                "new_defined_sample_project": new_defined_sample_project,
            },
        )

    return render(
        request,
        "wetlab/createSampleProjects.html",
        {"defined_samples_projects": defined_samples_projects},
    )


def define_additional_kits(request, protocol_id):
    # Check user == WETLAB_MANAGER: if false,  redirect to 'login' page
    if request.user.is_authenticated:
        if not wetlab.utils.common.is_wetlab_manager(request):
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "You do not have enough privileges to see this page ",
                        "Contact with your administrator .",
                    ]
                },
            )
    else:
        # redirect to login webpage
        return redirect("/accounts/login")

    additional_kits = wetlab.utils.additional_kits.define_table_for_additional_kits(protocol_id)
    if request.method == "POST" and request.POST["action"] == "defineAdditionalKits":
        recorded_additional_kits = wetlab.utils.additional_kits.set_additional_kits(request.POST, request.user)
        if len(recorded_additional_kits) == 0:
            return render(
                request,
                "wetlab/defineAdditionalKits.html",
                {"additional_kits": additional_kits},
            )

        return render(
            request,
            "wetlab/defineAdditionalKits.html",
            {"recorded_additional_kits": recorded_additional_kits},
        )

    else:
        if not check_if_protocol_exists(protocol_id, __package__):
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "The requested Protocol does not exist",
                        "Create the protocol name before assigning additional kits for protocol.",
                    ]
                },
            )

        return render(
            request,
            "wetlab/defineAdditionalKits.html",
            {"additional_kits": additional_kits},
        )


@login_required
def display_sample_project(request, sample_project_id):
    samples_project_data = core.utils.samples.get_info_to_display_sample_project(sample_project_id)
    if "ERROR" in samples_project_data:
        error_message = samples_project_data["ERROR"]
        return render(
            request, "wetlab/error_page.html", {"content": error_message}
        )
    return render(
        request,
        "wetlab/displaySampleProject.html",
        {"samples_project_data": samples_project_data},
    )


@login_required
def display_protocol(request, protocol_id):
    if not wetlab.utils.common.is_wetlab_manager(request):
        return render(
            request,
            "wetlab/error_page.html",
            {
                "content": [
                    "You do not have enough privileges to see this page ",
                    "Contact with your administrator .",
                ]
            },
        )
    if not check_if_protocol_exists(protocol_id, __package__):
        return render(
            request,
            "wetlab/error_page.html",
            {
                "content": [
                    "The protocol that you are trying to get ",
                    "DOES NOT exists .",
                ]
            },
        )
    protocol_data = get_all_protocol_info(protocol_id)
    kit_data = wetlab.utils.additional_kits.get_all_additional_kit_info(protocol_id)

    return render(
        request,
        "wetlab/displayProtocol.html",
        {"protocol_data": protocol_data, "kit_data": kit_data},
    )


@login_required
def define_protocol_parameters(request, protocol_id):
    # Check user == WETLAB_MANAGER: if false,  redirect to 'login' page
    if request.user.is_authenticated:
        if not wetlab.utils.common.is_wetlab_manager(request):
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "You do not have enough privileges to see this page ",
                        "Contact with your administrator .",
                    ]
                },
            )
    else:
        # redirect to login webpage
        return redirect("/accounts/login")

    if (
        request.method == "POST"
        and request.POST["action"] == "define_protocol_parameters"
    ):
        recorded_prot_parameters = set_protocol_parameters(request)

        return render(
            request,
            "wetlab/defineProtocolParameters.html",
            {"recorded_prot_parameters": recorded_prot_parameters},
        )

    else:
        if not check_if_protocol_exists(protocol_id, __package__):
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "The requested Protocol does not exist",
                        "Create the protocol name before assigning custom protocol parameters.",
                    ]
                },
            )

        prot_parameters = define_table_for_prot_parameters(protocol_id)
        return render(
            request,
            "wetlab/defineProtocolParameters.html",
            {"prot_parameters": prot_parameters},
        )


@login_required
def add_commercial_kit(request):
    app_name = __package__.split(".")[0]
    defined_protocols = get_defined_protocols(app_name, False)
    defined_platforms = core.utils.platforms.get_defined_platforms_and_ids("NGS")
    commercial_kits_data = get_data_for_commercial_kits("NGS")

    if request.method == "POST" and request.POST["action"] == "addCommercialKit":
        if get_commercial_kit_id(request.POST["kitName"]):
            return render(
                request,
                "wetlab/addCommercialKit.html",
                {
                    "defined_protocols": defined_protocols,
                    "invalid_name": request.POST["kitName"],
                },
            )
        new_kit = store_commercial_kit(request.POST)
        new_kit_data = get_commercial_kit_basic_data(new_kit)
        return render(
            request,
            "wetlab/addCommercialKit.html",
            {"new_kit_data": new_kit_data},
        )
    else:
        return render(
            request,
            "wetlab/addCommercialKit.html",
            {
                "defined_protocols": defined_protocols,
                "defined_platforms": defined_platforms,
                "commercial_kits_data": commercial_kits_data,
            },
        )


@login_required
def add_user_lot_commercial_kit(request):
    defined_kits = get_defined_commercial_kits()
    if request.method == "POST" and request.POST["action"] == "addUserLotKit":
        if get_lot_user_commercial_kit_id(request.POST["barCode"]):
            return render(
                request,
                "wetlab/addUserLotCommercialKit.html",
                {
                    "defined_kits": defined_kits,
                    "invalid_name": request.POST["nickName"],
                },
            )
        new_lot_kit = store_lot_user_commercial_kit(request.POST, request.user)
        new_lot_kit_data = get_lot_user_commercial_kit_basic_data(new_lot_kit)
        return render(
            request,
            "wetlab/addUserLotCommercialKit.html",
            {"new_lot_kit_data": new_lot_kit_data},
        )
    else:
        return render(
            request,
            "wetlab/addUserLotCommercialKit.html",
            {"defined_kits": defined_kits},
        )


@login_required
def pending_to_update(request):
    pending = {}
    # get the samples in defined state

    pending["defined"] = core.utils.samples.get_samples_in_defined_state("")
    pending["extract_molecule"] = core.utils.samples.get_samples_in_extracted_molecule_state(request.user)
    pending["graphic_pending_samples"] = pending_samples_for_grafic(pending).render()

    return render(request, "wetlab/pendingToUpdate.html", {"pending": pending})


@login_required
def record_samples(request):
    """
    Functions :
        analyze_input_samples
        analyze_input_sample_project_fields
        prepare_sample_input_table
        get_codeID_for_resequencing
        prepare_sample_project_input_table
        analyze_reprocess_data
        get_info_for_reprocess_samples
    """
    # Record new samples
    if request.method == "POST" and request.POST["action"] == "recordsample":
        sample_recorded = core.utils.samples.analyze_input_samples(request, __package__)
        # if no samples are in any of the options, displays the inital page

        if (
            "defined_samples" not in sample_recorded
            and "pre_defined_samples" not in sample_recorded
            and "invalid_samples" not in sample_recorded
            and "incomplete_samples" not in sample_recorded
        ):
            sample_information = core.utils.samples.prepare_sample_input_table(__package__)
            return render(
                request,
                "wetlab/record_sample.html",
                {"sample_information": sample_information},
            )

        if "sample_id_for_action" in sample_recorded:
            sample_recorded.update(get_codeID_for_resequencing(sample_recorded))
        if "incomplete_samples" in sample_recorded:
            sample_recorded.update(core.utils.samples.prepare_sample_input_table(__package__))
            sample_recorded["number_of_samples"] = len(
                sample_recorded["incomplete_samples"]
            )
        if "pre_defined_samples_id" in sample_recorded:
            sample_recorded.update(
                core.utils.samples.prepare_sample_project_input_table(
                    sample_recorded["pre_defined_samples_id"]
                )
            )
        return render(
            request,
            "wetlab/record_sample.html",
            {"sample_recorded": sample_recorded},
        )

    # Request to reprocess the samples
    elif request.method == "POST" and request.POST["action"] == "reprocessSamples":
        to_be_reprocessed_ids = request.POST["invalidSamplesID"].split(",")
        reprocess_id = request.POST["sampleIDforAction"]
        json_data = json.loads(request.POST["reprocess_data"])

        result = analyze_reprocess_data(json_data[0], reprocess_id, request.user)
        if result == "Invalid options":
            to_be_reprocessed_ids.insert(0, reprocess_id)
            sample_recorded = core.utils.samples.get_info_for_reprocess_samples(
                to_be_reprocessed_ids, reprocess_id
            )
            sample_recorded["invalid_samples_id"] = request.POST["invalidSamplesID"]
            sample_recorded["sample_id_for_action"] = reprocess_id
            sample_recorded.update(get_codeID_for_resequencing(sample_recorded))
            sample_recorded["reprocess_result"] = "False"
        else:
            if to_be_reprocessed_ids[0] == "":
                return render(
                    request,
                    "wetlab/record_sample.html",
                    {"all_sucessful_reprocess": True},
                )
            else:
                next_to_be_process_id = str(to_be_reprocessed_ids[0])
                sample_recorded = core.utils.samples.get_info_for_reprocess_samples(
                    to_be_reprocessed_ids, next_to_be_process_id
                )
                sample_recorded["invalid_samples_id"] = ",".join(to_be_reprocessed_ids)
                sample_recorded["sample_id_for_action"] = next_to_be_process_id
                sample_recorded.update(get_codeID_for_resequencing(sample_recorded))
                sample_recorded["reprocess_result"] = "True"
                return render(
                    request,
                    "wetlab/record_sample.html",
                    {"sample_recorded": sample_recorded},
                )

        # NOT WORKING.This is not working so commented for now. require_to_update not defined anywhere.
        """
        if len(require_to_update) > 0:
            for key, value in require_to_update.items():
                if value == "newLibPreparation":
                    if not libraryPreparation.objects.filter(sample_id__pk__exact=key):
                        continue
                    lib_prep_obj = libraryPreparation.objects.get(
                        sample_id__pk__exact=key
                    )
                    # lib_prep_obj.set_state('Reused')
                    lib_prep_obj.set_increase_reuse()
                elif value == "newPool":
                    pass
                else:
                    continue
        """

        return render(
            request,
            "wetlab/record_sample.html",
            {"reprocess_result": sample_recorded["reprocess_result"]},
        )

    # display the form to show the samples in pre-defined state that user requested to complete
    elif (
        request.method == "POST"
        and request.POST["action"] == "select_samples_pre_defined"
    ):
        if "samples_in_list" in request.POST:
            pre_defined_samples_id = request.POST.getlist("samples")
        sample_recorded = core.utils.samples.prepare_sample_project_input_table(pre_defined_samples_id)
        return render(
            request,
            "wetlab/record_sample.html",
            {"sample_recorded": sample_recorded},
        )

    # Add the additional information related to the project
    elif request.method == "POST" and request.POST["action"] == "sampleprojectdata":
        sample_recorded = core.utils.samples.analyze_input_sample_project_fields(request.POST)

        if request.POST["pending_pre_defined"] != "":
            sample_recorded.update(
                core.utils.samples.prepare_sample_project_input_table(
                    request.POST["pending_pre_defined"].split(",")
                )
            )
            return render(
                request,
                "wetlab/record_sample.html",
                {"sample_recorded": sample_recorded},
            )
        else:
            return render(
                request,
                "wetlab/record_sample.html",
                {"sample_recorded": sample_recorded},
            )
    # Load batch file
    elif request.method == "POST" and request.POST["action"] == "defineBatchSamples":
        sample_information = core.utils.samples.prepare_sample_input_table(__package__)
        if "samplesExcel" in request.FILES:
            samples_batch_df = core.utils.load_batch.read_batch_sample_file(request.FILES["samplesExcel"])
            if "ERROR" in samples_batch_df:
                return render(
                    request,
                    "wetlab/record_sample.html",
                    {
                        "sample_information": sample_information,
                        "error_message": samples_batch_df["ERROR"],
                    },
                )
            valid_file_result = core.utils.load_batch.valid_sample_batch_file(samples_batch_df, __package__)
            if valid_file_result != "OK":
                return render(
                    request,
                    "wetlab/record_sample.html",
                    {
                        "sample_information": sample_information,
                        "error_message": valid_file_result,
                    },
                )
            result_recorded = core.utils.load_batch.save_samples_in_batch_file(
                samples_batch_df, request.user.username, __package__
            )
            if result_recorded != "OK":
                return render(
                    request,
                    "wetlab/record_sample.html",
                    {
                        "sample_information": sample_information,
                        "error_message": result_recorded,
                    },
                )
            return render(
                request,
                "wetlab/record_sample.html",
                {
                    "sample_information": sample_information,
                    "successfuly_batch_load": "ok",
                },
            )
    # Form to get the new samples
    else:
        sample_information = core.utils.samples.prepare_sample_input_table(__package__)
        return render(
            request,
            "wetlab/record_sample.html",
            {"sample_information": sample_information},
        )


@login_required
def define_sample_projects_fields(request, sample_project_id):
    # Check user == WETLAB_MANAGER: if false,  redirect to 'login' page
    if request.user.is_authenticated:
        if not wetlab.utils.common.is_wetlab_manager(request):
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "You do not have enough privileges to see this page ",
                        "Contact with your administrator .",
                    ]
                },
            )
    else:
        # redirect to login webpage
        return redirect("/accounts/login")
    # get the list of defined sample Projects
    if (
        request.method == "POST"
        and request.POST["action"] == "defineSampleProjectFields"
    ):
        sample_project_field_data = core.utils.samples.set_sample_project_fields(request.POST)
        return render(
            request,
            "wetlab/defineSampleProjectFields.html",
            {"sample_project_field_data": sample_project_field_data},
        )

    elif request.method == "POST" and request.POST["action"] == "defineBatchFields":
        sample_project_data = core.utils.samples.define_table_for_sample_project_fields(
            request.POST["sample_project_id"]
        )
        # sample_information = prepare_sample_input_table(__package__)
        if "jsonSchema" in request.FILES:
            schema = core.utils.load_batch.read_json_schema(request.FILES["jsonSchema"])
            if "ERROR" in schema:
                return render(
                    request,
                    "wetlab/defineSampleProjectFields.html",
                    {
                        "sample_project_data": sample_project_data,
                        "error_message": schema["ERROR"],
                    },
                )
            result = core.utils.load_batch.store_schema(
                schema["schema"],
                request.POST["classification"],
                request.POST["subfilter"],
                request.POST["sample_project_id"],
            )
            if "ERROR" in result:
                return render(
                    request,
                    "wetlab/defineSampleProjectFields.html",
                    {
                        "error_message": result["ERROR"],
                        "sample_project_data": sample_project_data,
                    },
                )
            return render(
                request,
                "wetlab/defineSampleProjectFields.html",
                {"schema_result": result},
            )
        return render(
            request,
            "wetlab/defineSampleProjectFields.html",
            {"error_message": wetlab.config.ERROR_MESSAGE_UPLOAD_FILE_NOT_EXISTS},
        )
    else:
        if not core.utils.samples.check_if_sample_project_id_exists(sample_project_id):
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "The requested Protocol does not exist",
                        "Create the protocol name before assigning custom protocol parameters.",
                    ]
                },
            )

        sample_project_data = core.utils.samples.define_table_for_sample_project_fields(sample_project_id)
        return render(
            request,
            "wetlab/defineSampleProjectFields.html",
            {"sample_project_data": sample_project_data},
        )


@login_required
def modify_additional_kits(request, protocol_id):
    # Check user == WETLAB_MANAGER: if false,  redirect to 'login' page
    if request.user.is_authenticated:
        if not wetlab.utils.common.is_wetlab_manager(request):
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "You do not have enough privileges to see this page ",
                        "Contact with your administrator .",
                    ]
                },
            )
    else:
        # redirect to login webpage
        return redirect("/accounts/login")

    if request.method == "POST" and request.POST["action"] == "modifyAdditionalKits":
        additional_kits_data_saved = wetlab.utils.additional_kits.modify_fields_in_additional_kits(
            request.POST, request.user
        )
        return render(
            request,
            "wetlab/modifyAdditionalKits.html",
            {"additional_kits_data_saved": additional_kits_data_saved},
        )
    else:
        if not check_if_protocol_exists(protocol_id, __package__):
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "The requested additional kits do not exist",
                        "Create the addtional kits before.",
                    ]
                },
            )
        additional_kits_data = wetlab.utils.additional_kits.get_additional_kits_data_to_modify(protocol_id)
        return render(
            request,
            "wetlab/modifyAdditionalKits.html",
            {"additional_kits_data": additional_kits_data},
        )


@login_required
def modify_protocol_fields(request, protocol_id):
    # Check user == WETLAB_MANAGER: if false,  redirect to 'login' page
    if request.user.is_authenticated:
        if not wetlab.utils.common.is_wetlab_manager(request):
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "You do not have enough privileges to see this page ",
                        "Contact with your administrator .",
                    ]
                },
            )
    else:
        # redirect to login webpage
        return redirect("/accounts/login")

    if request.method == "POST" and request.POST["action"] == "modifyProtocolFields":
        protocol_field_saved = modify_fields_in_protocol(request.POST)
        return render(
            request,
            "wetlab/modifyProtocolFields.html",
            {"protocol_field_saved": protocol_field_saved},
        )
    else:
        if not check_if_protocol_exists(protocol_id, __package__):
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "The requested Protocol does not exist",
                        "Create the protocol name before assigning custom parameters.",
                    ]
                },
            )
        protocol_field = get_protocol_fields(protocol_id)
        return render(
            request,
            "wetlab/modifyProtocolFields.html",
            {"protocol_field": protocol_field},
        )


@login_required
def modify_sample_project_fields(request, sample_project_id):
    # Check user == WETLAB_MANAGER: if false,  redirect to 'login' page

    if request.user.is_authenticated:
        if not wetlab.utils.common.is_wetlab_manager(request):
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "You do not have enough privileges to see this page ",
                        "Contact with your administrator .",
                    ]
                },
            )
    else:
        # redirect to login webpage
        return redirect("/accounts/login")

    if (
        request.method == "POST"
        and request.POST["action"] == "modifySampleProjectFields"
    ):
        sample_project_field_saved = core.utils.samples.modify_fields_in_sample_project(request.POST)
        return render(
            request,
            "wetlab/modifySampleProjectFields.html",
            {"sample_project_field_saved": sample_project_field_saved},
        )

    else:
        if not core.utils.samples.check_if_sample_project_id_exists(sample_project_id):
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "The requested Sample project does not exist",
                        "Create the sample project name before assigning custom sample project parameters.",
                    ]
                },
            )
        sample_project_field = core.utils.samples.get_parameters_sample_project(sample_project_id)
        return render(
            request,
            "wetlab/modifySampleProjectFields.html",
            {"sample_project_field": sample_project_field},
        )


@login_required
def define_molecule_uses(request):
    """
    Functions:
        display_molecule_use
        record_molecule_use
    """
    molecule_use_data = core.utils.samples.display_molecule_use(__package__)
    if request.method == "POST" and request.POST["action"] == "record_molecule_use":
        molecule_use_data.update(core.utils.samples.record_molecule_use(request.POST, __package__))

    return render(
        request,
        "wetlab/defineMoleculeUses.html",
        {"molecule_use_data": molecule_use_data},
    )


@login_required
def define_type_of_samples(request):
    """
    Functions:
        display_sample_types
        save_type_of_sample
    """
    sample_types = core.utils.samples.display_sample_types(__package__)
    if request.method == "POST" and request.POST["action"] == "addNewSampleType":
        sample_types.update(core.utils.samples.save_type_of_sample(request.POST, __package__))

    return render(
        request,
        "wetlab/define_type_of_samples.html",
        {"sample_types": sample_types},
    )


@login_required
def display_sample(request, sample_id):
    """
    Functions:
        get_all_sample_information
        get_all_library_information
        get_additional_kits_used_in_sample
        get_sample_in_project_obj_from_sample_name
    """
    sample_information = core.utils.samples.get_all_sample_information(sample_id, True)
    if "Error" not in sample_information:
        sample_information.update(get_molecule_lot_kit_in_sample(sample_id))
        sample_information.update(get_all_library_information(sample_id))
        sample_information.update(wetlab.utils.additional_kits.get_additional_kits_used_in_sample(sample_id))
        sample_information.update(get_run_user_lot_kit_used_in_sample(sample_id))
    else:
        sample_information = {}
    sample_obj = core.utils.samples.get_sample_obj_from_id(sample_id)
    if sample_obj:
        sample_name = sample_obj.get_sample_name()
        run_sample_obj = get_sample_in_project_obj_from_sample_name(sample_name)
        if run_sample_obj:
            sample_information.update(get_info_sample_in_run(run_sample_obj))

    if len(sample_information) == 0:
        return render(
            request,
            "wetlab/error_page.html",
            {"content": ["No Sample was found"]},
        )
    else:
        return render(
            request,
            "wetlab/display_sample.html",
            {"sample_information": sample_information},
        )


@login_required
def display_sample_in_run(request, sample_run_id):
    """
    Functions:
        get_info_sample_in_run
    """
    sample_run_obj = get_sample_in_project_obj_from_id(sample_run_id)
    if not sample_run_obj:
        return render(
            request,
            "wetlab/error_page.html",
            {"content": ["No Sample was found"]},
        )
    sample_information = get_info_sample_in_run(sample_run_obj)
    return render(
        request,
        "wetlab/display_sample.html",
        {"sample_information": sample_information},
    )


@login_required
def display_type_of_sample(request, sample_type_id):
    """
    Functions:
        get_type_of_sample_information
    """
    type_of_sample_data = core.utils.samples.get_type_of_sample_information(sample_type_id)
    return render(
        request,
        "wetlab/display_type_of_sample.html",
        {"type_of_sample_data": type_of_sample_data},
    )


@login_required
def handling_library_preparations(request):
    """
    Functions:
        analyze_and_store_input_additional_kits
        get_samples_for_library_preparation
        create_library_preparation_instance
        extract_user_sample_sheet_data
        get_additional_kits_from_lib_prep
        get_type_of_sample_information
        get_library_preparation_heading_for_samples
        get_protocols_for_library_preparation
        get_samples_in_lib_prep_state
        validate_sample_sheet_data
    """
    # get the information for returning the uploaded file in case errors in the sample sheet
    samples_in_lib_prep = get_samples_for_library_preparation()

    if request.method == "POST" and request.POST["action"] == "assignProtocol":
        samples_in_lib_prep_protocol = extract_protocol_library_preparation_form(
            request.POST
        )
        if len(samples_in_lib_prep_protocol) == 0:
            return render(
                request,
                "wetlab/handlingLibraryPreparations.html",
                {"stored_lib_prep": samples_in_lib_prep},
            )
        library_preparation_objs = create_library_preparation_instance(
            samples_in_lib_prep_protocol, request.user
        )
        lib_prep_protocol_parameters = get_protocol_parameters_for_library_preparation(
            library_preparation_objs
        )
        return render(
            request,
            "wetlab/handlingLibraryPreparations.html",
            {"lib_prep_protocol_parameters": lib_prep_protocol_parameters},
        )

    # add protocol parameters for the user selected library preparation on defined state
    if request.method == "POST" and request.POST["action"] == "addProtocolParameter":
        lib_prep_ids = request.POST.getlist("libpreparation")
        library_preparation_objs = []
        for lib_prep_id in lib_prep_ids:
            library_preparation_objs.append(get_lib_prep_obj_from_id(lib_prep_id))
        lib_prep_protocol_parameters = get_protocol_parameters_for_library_preparation(
            library_preparation_objs
        )
        return render(
            request,
            "wetlab/handlingLibraryPreparations.html",
            {"lib_prep_protocol_parameters": lib_prep_protocol_parameters},
        )

    # store the parameter library preparation protocol
    if request.method == "POST" and request.POST["action"] == "recordProtocolParamters":
        stored_params = analyze_and_store_input_param_values(request.POST)
        if "ERROR" in stored_params:
            error_message = stored_params["ERROR"]
            lib_prep_ids = request.POST["lib_prep_ids"].split(",")
            library_preparation_objs = []
            for lib_prep_id in lib_prep_ids:
                library_preparation_objs.append(get_lib_prep_obj_from_id(lib_prep_id))
            lib_prep_protocol_parameters = (
                get_protocol_parameters_for_library_preparation(
                    library_preparation_objs
                )
            )
            # restore the user data
            lib_prep_protocol_parameters["data"] = json.loads(
                request.POST["protocol_data"]
            )
            return render(
                request,
                "wetlab/handlingLibraryPreparations.html",
                {
                    "ERROR": error_message,
                    "lib_prep_protocol_parameters": lib_prep_protocol_parameters,
                },
            )
        return render(
            request,
            "wetlab/handlingLibraryPreparations.html",
            {"stored_params": stored_params},
        )

    if request.method == "POST" and request.POST["action"] == "importsamplesheet":
        data = {}
        data["full_path_file"], data["file_name"] = store_user_input_file(
            request.FILES["uploadfile"]
        )
        file_read = read_user_iem_file(data["full_path_file"])
        if not valid_user_iem_file(file_read):
            # Error found when extracting data from sample sheet
            data["ERROR"] = wetlab.config.ERROR_INVALID_FILE_FORMAT
            if not delete_stored_file(data["full_path_file"]):
                data["ERROR"].append(wetlab.config.ERROR_UNABLE_TO_DELETE_USER_FILE)
            return render(
                request,
                "wetlab/handlingLibraryPreparations.html",
                {"ERROR": data["ERROR"], "samples_in_lib_prep": samples_in_lib_prep},
            )
        user_in_description = wetlab.utils.common.get_configuration_value(
            "DESCRIPTION_IN_SAMPLE_SHEET_MUST_HAVE_USERNAME"
        )
        if user_in_description == "TRUE":
            user_id_in_s_sheet = extract_userids_from_sample_sheet_data(file_read)
            if "ERROR" in user_id_in_s_sheet:
                if not delete_stored_file(data["full_path_file"]):
                    user_id_in_s_sheet["ERROR"].append(wetlab.config.ERROR_UNABLE_TO_DELETE_USER_FILE)
                return render(
                    request,
                    "wetlab/handlingLibraryPreparations.html",
                    {
                        "ERROR": user_id_in_s_sheet["ERROR"],
                        "samples_in_lib_prep": samples_in_lib_prep,
                    },
                )
        else:
            user_id_in_s_sheet = []
        sample_sheet_data = get_sample_sheet_data(file_read)

        valid_data = validate_sample_sheet_data(sample_sheet_data)
        if "ERROR" in valid_data:
            if not delete_stored_file(data["full_path_file"]):
                valid_data["ERROR"].append(wetlab.config.ERROR_UNABLE_TO_DELETE_USER_FILE)
            return render(
                request,
                "wetlab/handlingLibraryPreparations.html",
                {
                    "ERROR": valid_data["ERROR"],
                    "samples_in_lib_prep": samples_in_lib_prep,
                },
            )

        platform = request.POST["platform"]
        configuration = request.POST[request.POST["platform"]]
        sample_sheet_data["file_name"] = data["file_name"]
        sample_sheet_data["userid_names"] = user_id_in_s_sheet
        lib_prep_sample_sheet_obj = store_library_preparation_sample_sheet(
            sample_sheet_data, request.user, platform, configuration
        )
        display_sample_sheet = format_sample_sheet_to_display_in_form(sample_sheet_data)
        display_sample_sheet[
            "lib_prep_user_sample_sheet"
        ] = lib_prep_sample_sheet_obj.get_user_sample_sheet_id()
        display_sample_sheet["platform"] = platform
        display_sample_sheet["iem_version"] = sample_sheet_data["iem_version"]
        if user_in_description == "TRUE":
            display_sample_sheet["user_list"] = get_userid_list()
        return render(
            request,
            "wetlab/handlingLibraryPreparations.html",
            {"display_sample_sheet": display_sample_sheet},
        )

    if request.method == "POST" and request.POST["action"] == "storeIndexSample":
        store_data_result = store_confirmation_library_preparation_index(request.POST)
        if "ERROR" in store_data_result:
            return render(
                request,
                "wetlab/handlingLibraryPreparations.html",
                {
                    "ERROR": valid_data["ERROR"],
                    "samples_in_lib_prep": samples_in_lib_prep,
                },
            )
        stored_index = "True"
        return render(
            request,
            "wetlab/handlingLibraryPreparations.html",
            {"stored_index": stored_index},
        )

        # NOT WORKING. get_library_preparation_heading_for_samples does not exits.
        """
        if request.method == "POST" and request.POST["action"] == "libpreparationdefined":
        lib_prep_defined = request.POST.getlist("libpreparation")
        lib_protocols = get_protocol_from_library_id(lib_prep_defined[0])
        stored_lib_prep = get_library_preparation_heading_for_samples(
            lib_prep_defined, lib_protocols
        )
        return render(
            request,
            "wetlab/handlingLibraryPreparations.html",
            {"stored_lib_prep": stored_lib_prep},
        )
        """

    if request.method == "POST" and request.POST["action"] == "assignAdditionalKits":
        lib_prep_ids = request.POST.getlist("libpreparation")
        additional_kits = wetlab.utils.additional_kits.get_additional_kits_from_lib_prep(lib_prep_ids)
        return render(
            request,
            "wetlab/handlingLibraryPreparations.html",
            {"additional_kits": additional_kits},
        )

    if request.method == "POST" and request.POST["action"] == "storeAdditionalKits":
        stored_additional_kits = wetlab.utils.additional_kits.analyze_and_store_input_additional_kits(request.POST)
        if "ERROR" in stored_additional_kits:
            error_message = stored_additional_kits["ERROR"]
            lib_prep_ids = request.POST["lib_prep_ids"].split(",")
            additional_kits = wetlab.utils.additional_kits.get_additional_kits_from_lib_prep(lib_prep_ids)
            additional_kits["data"] = json.loads(request.POST["protocol_data"])
            return render(
                request,
                "wetlab/handlingLibraryPreparations.html",
                {"ERROR": error_message, "additional_kits": additional_kits},
            )

        return render(
            request,
            "wetlab/handlingLibraryPreparations.html",
            {"stored_additional_kits": stored_additional_kits},
        )
    else:
        return render(
            request,
            "wetlab/handlingLibraryPreparations.html",
            {"samples_in_lib_prep": samples_in_lib_prep},
        )


def handling_molecules(request):
    """
    Functions:
        get_samples_in_state
        create_table_to_select_molecules
        display_molecule_protocol_parameters
    """

    if request.method == "POST" and request.POST["action"] == "selectedMolecules":
        # If no samples are selected , call again this function to display again the sample list

        samples = core.utils.samples.get_selected_recorded_samples(request.POST["selected_samples"])
        if len(samples) == 0:
            return redirect("handling_molecules")
        molecule_protocol = core.utils.samples.get_table_record_molecule(samples, __package__)
        if "ERROR" in molecule_protocol:
            return render(
                request,
                "wetlab/handling_molecules.html",
                {"ERROR": "There was no valid sample selected "},
            )
        molecule_protocol["samples"] = ",".join(samples)

        return render(
            request,
            "wetlab/handling_molecules.html",
            {"molecule_protocol": molecule_protocol},
        )

    elif (
        request.method == "POST" and request.POST["action"] == "updateMoleculeProtocol"
    ):
        molecule_recorded = core.utils.samples.record_molecules(request.POST, request.user, __package__)

        if (
            "molecule_code_ids" in request.POST
            and request.POST["molecule_code_ids"] != ""
            and "molecule_code_ids" in molecule_recorded
        ):
            # Add the already recorded molecules to the new ones
            molecule_recorded["molecule_code_ids"] += (
                "," + request.POST["molecule_code_ids"]
            )
            molecule_recorded["molecule_ids"] += "," + request.POST["molecule_ids"]
        if "incomplete_sample_ids" in molecule_recorded:
            # collect the information to select in the option fields
            molecule_recorded.update(
                core.utils.samples.get_table_record_molecule(
                    molecule_recorded["incomplete_sample_ids"], __package__
                )
            )

            return render(
                request,
                "wetlab/handling_molecules.html",
                {"molecule_recorded": molecule_recorded},
            )

        show_molecule_parameters = core.utils.samples.display_molecule_protocol_parameters(
            molecule_recorded["molecule_ids"].split(","), request.user
        )
        return render(
            request,
            "wetlab/handling_molecules.html",
            {
                "molecule_recorded": molecule_recorded,
                "show_molecule_parameters": show_molecule_parameters,
            },
        )

    elif (
        request.method == "POST" and request.POST["action"] == "selectedOwnerMolecules"
    ):
        # If no samples are selected , call again this function to display again the sample list
        if "molecules" not in request.POST:
            return redirect("handling_molecules")
        molecules = request.POST.getlist("molecules")
        # Set to true to reuse the html Code
        molecule_recorded = True
        show_molecule_parameters = core.utils.samples.display_molecule_protocol_parameters(
            molecules, request.user
        )
        return render(
            request,
            "wetlab/handlingMolecules.html",
            {
                "molecule_recorded": molecule_recorded,
                "show_molecule_parameters": show_molecule_parameters,
            },
        )

    elif request.method == "POST" and request.POST["action"] == "addMoleculeParameters":
        molecule_parameters_updated = core.utils.samples.add_molecule_protocol_parameters(request.POST)
        if "pending" in request.POST:
            molecules = request.POST["pending"].split(",")
            show_molecule_parameters = core.utils.samples.display_molecule_protocol_parameters(
                molecules, request.user
            )
            return render(
                request,
                "wetlab/handling_molecules.html",
                {
                    "molecule_parameters_updated": molecule_parameters_updated,
                    "show_molecule_parameters": show_molecule_parameters,
                },
            )
        else:
            return render(
                request,
                "wetlab/handling_molecules.html",
                {"molecule_parameters_updated": molecule_parameters_updated},
            )

    elif request.method == "POST" and request.POST["action"] == "requestMoleculeUse":
        molecule_use = core.utils.samples.set_molecule_use(request.POST, __package__)
        return render(
            request,
            "wetlab/handling_molecules.html",
            {"molecule_use": molecule_use},
        )

    else:
        sample_availables, user_molecules, request_molecule_use = "", "", ""
        samples_list = core.utils.samples.get_samples_in_state("Defined")
        if samples_list:
            sample_availables = core.utils.samples.create_table_to_select_molecules(samples_list)

        user_owner_molecules = core.utils.samples.get_molecule_in_state("Defined", request.user)
        if len(user_owner_molecules) > 0:
            user_molecules = core.utils.samples.create_table_user_molecules(user_owner_molecules)

        samples_pending_use = core.utils.samples.get_samples_in_state("Pending for use")
        if samples_pending_use:
            request_molecule_use = core.utils.samples.create_table_pending_use(
                samples_pending_use, __package__
            )
        # check if there are defined the type
        molecule_use_defined = core.utils.samples.check_if_molecule_use_defined(__package__)

        return render(
            request,
            "wetlab/handling_molecules.html",
            {
                "sample_availables": sample_availables,
                "user_molecules": user_molecules,
                "molecule_use_defined": molecule_use_defined,
                "request_molecule_use": request_molecule_use,
            },
        )

    return


@login_required
def repeat_library_preparation(request):
    """
    Functions:
        analyze_reprocess_data  # located at utils/sample_functions.py
    """

    if (
        request.method == "POST"
        and request.POST["action"] == "repeat_library_preparation"
    ):
        molecule_code_id = request.POST["molecule_code_id"]
        sample_id = request.POST["sample_id"]
        result = analyze_reprocess_data(
            [molecule_code_id, "New Library Preparation"], sample_id, request.user
        )
        detail_description = {}
        if result == "Invalid options":
            detail_description["heading"] = wetlab.config.ERROR_UNABLE_SAVE_REQUEST
            detail_description["information"] = [
                wetlab.config.ERROR_INVALID_PARAMETERS_WHEN_REUSING_LIB_PREP
            ]
            return render(
                request,
                "wetlab/error_page.html",
                {"detail_description": detail_description},
            )
        detail_description["information"] = wetlab.config.SUCCESSFUL_REUSE_MOLECULE_EXTRACTION
        return render(
            request,
            "wetlab/successful_page.html",
            {"detail_description": detail_description},
        )
    # return to the main page because the page was not requested for the right page
    return redirect("")


@login_required
def repeat_molecule_extraction(request):
    """
    Functions:
    analyze_reprocess_data
    get_table_record_molecule
    """

    if request.method == "POST" and request.POST["action"] == "repeat_extraction":
        sample_id = request.POST["sample_id"]
        if analyze_reprocess_data(["New Extraction"], sample_id, request.user):
            molecule_protocol = core.utils.samples.get_table_record_molecule([sample_id], __package__)
            molecule_protocol["samples"] = sample_id
            # create a copy of the request, to allow to modify it
            # request.POST = request.POST.copy()
            # request.POST['action'] = 'selectedMolecules'
            # request.POST['samples'] = sample_id

            return render(
                request,
                "wetlab/handlingMolecules.html",
                {"molecule_protocol": molecule_protocol},
            )
    # return to the main page because the page was not requested for the right page
    return redirect("")


@login_required
def repeat_pool(request):
    """
    Functions:
    analyze_reprocess_data  : located at utils/sample_functions.py
    """
    if request.method == "POST" and request.POST["action"] == "repeat_pool":
        lib_prep_obj = get_lib_prep_obj_from_id(request.POST["lib_prep_id"])
        lib_prep_code_id = lib_prep_obj.get_lib_prep_code()
        molecule_code_id = lib_prep_obj.get_molecule_code_id()
        sample_id = lib_prep_obj.get_sample_id()

        result = analyze_reprocess_data(
            [molecule_code_id, lib_prep_code_id, "New Pool"], sample_id, request.user
        )
        detail_description = {}
        if result == "Invalid options":
            detail_description["heading"] = wetlab.config.ERROR_UNABLE_SAVE_REQUEST
            detail_description["information"] = [
                wetlab.config.ERROR_INVALID_PARAMETERS_WHEN_REUSING_LIB_PREP
            ]
            return render(
                request,
                "wetlab/error_page.html",
                {"detail_description": detail_description},
            )
        detail_description["information"] = wetlab.config.SUCCESSFUL_REUSE_LIB_PREP
        return render(
            request,
            "wetlab/successful_page.html",
            {"detail_description": detail_description},
        )
    # return to the main page because the page was not requested for the right page
    return redirect("")


@login_required
def search_sample(request):
    """
    Functions:
        get_sample_states
        check_valid_date_format
        search_samples
        search_run_samples
    """
    search_data = {}
    search_data["s_state"] = core.utils.samples.get_sample_states()

    if request.method == "POST" and request.POST["action"] == "searchsample":
        sample_name = request.POST["samplename"]
        start_date = request.POST["startdate"]
        end_date = request.POST["enddate"]
        user_name = request.POST["username"]
        sample_state = request.POST["sampleState"]

        # check that some values are in the request if not return the form
        if (
            user_name == ""
            and start_date == ""
            and end_date == ""
            and sample_name == ""
            and sample_state == ""
        ):
            return render(
                request,
                "wetlab/search_sample.html",
                {"search_data": search_data},
            )

        if user_name != "" and len(user_name) < 5:
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "The user name must contains at least 5 caracters ",
                        "ADVICE:",
                        "write the full user name to get a better match",
                    ]
                },
            )
        # check the right format of start and end date

        if start_date != "" and not wetlab.utils.common.check_valid_date_format(start_date):
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        'The format for the "Start Date Search" Field is incorrect ',
                        "ADVICE:",
                        "Use the format  (DD-MM-YYYY)",
                    ]
                },
            )
        if end_date != "" and not wetlab.utils.common.check_valid_date_format(end_date):
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        'The format for the "End Date Search" Field is incorrect ',
                        "ADVICE:",
                        "Use the format  (DD-MM-YYYY)",
                    ]
                },
            )
        # Get projects when sample name is not empty
        sample_list = core.utils.samples.search_samples(
            sample_name, user_name, sample_state, start_date, end_date
        )
        if sample_state == "":
            run_sample_list = search_run_samples(
                sample_name, user_name, start_date, end_date
            )
        else:
            run_sample_list = ""

        if len(sample_list) == 0 and len(run_sample_list) == 0:
            return render(
                request,
                "wetlab/search_sample.html",
                {
                    "no_samples": wetlab.config.ERROR_NO_SAMPLE_FOUND,
                    "sample_list": sample_list,
                    "search_data": search_data,
                },
            )
        elif len(sample_list) == 1 and len(run_sample_list) == 0:
            return redirect("display_sample", sample_id=sample_list[0])
        elif len(sample_list) == 0 and len(run_sample_list) == 1:
            # get_info_sample_in_run (run_sample_obj)
            return redirect("display_sample_in_run", sample_run_id=run_sample_list[0])
        else:
            # get the sample information to select it , because there are also matches on run_sample
            if len(sample_list) == 1:
                sample_obj = core.utils.samples.get_sample_obj_from_id(sample_list[0])
                sample_list = [sample_obj.get_info_for_searching()]
            if len(run_sample_list) == 1:
                run_sample_obj = get_sample_in_project_obj_from_id(run_sample_list[0])
                run_sample_list = [run_sample_obj.get_basic_info()]
            return render(
                request,
                "wetlab/search_sample.html",
                {"sample_list": sample_list, "run_sample_list": run_sample_list},
            )

    else:
        return render(
            request, "wetlab/search_sample.html", {"search_data": search_data}
        )


@login_required
def set_molecule_values(request):
    if request.method == "POST" and request.POST["action"] == "continueWithMolecule":
        if request.POST["samples"] == "":
            return render(
                request,
                "wetlab/error_page.html",
                {"content": ["There was no sample selected "]},
            )
        if "samples_in_list" in request.POST:
            samples = request.POST.getlist("samples")
        else:
            samples = request.POST["samples"].split(",")

        molecule_protocol = core.utils.samples.get_table_record_molecule(samples, __package__)
        if "ERROR" in molecule_protocol:
            return render(
                request,
                "wetlab/error_page.html",
                {"content": ["There was no valid sample selected "]},
            )

        molecule_protocol["samples"] = ",".join(samples)

        return render(
            request,
            "wetlab/setMoleculeValues.html",
            {"molecule_protocol": molecule_protocol},
        )

    elif (
        request.method == "POST" and request.POST["action"] == "updateMoleculeProtocol"
    ):
        molecule_recorded = core.utils.samples.record_molecules(request)

        if "heading" not in molecule_recorded:
            samples = request.POST["samples"].split(",")
            molecule_protocol = core.utils.samples.get_table_record_molecule(samples, __package__)
            molecule_protocol["data"] = molecule_recorded["incomplete_molecules"]
            molecule_protocol["samples"] = ",".join(samples)
            return render(
                request,
                "wetlab/setMoleculeValues.html",
                {"molecule_protocol": molecule_protocol},
            )
        else:
            if "incomplete_molecules" in molecule_recorded:
                samples = molecule_recorded["incomplete_molecules_ids"].split(",")
                molecule_recorded.update(core.utils.samples.get_table_record_molecule(samples))

            return render(
                request,
                "wetlab/setMoleculeValues.html",
                {"molecule_recorded": molecule_recorded},
            )

    elif (
        request.method == "POST"
        and request.POST["action"] == "displayMoleculeParameters"
    ):
        if "samples_in_list" in request.POST:
            molecules = request.POST.getlist("molecules")
        else:
            molecules = request.POST["molecules"].split(",")
        show_molecule_parameters = core.utils.samples.display_molecule_protocol_parameters(
            molecules, request.user
        )
        return render(
            request,
            "wetlab/setMoleculeValues.html",
            {"show_molecule_parameters": show_molecule_parameters},
        )

    elif request.method == "POST" and request.POST["action"] == "addMoleculeParameters":
        (
            added_molecule_protocol_parameters,
            sample_updated_list,
        ) = core.utils.samples.add_molecule_protocol_parameters(request)
        if "pending" in request.POST:
            molecules = request.POST["pending"].split(",")
            show_molecule_parameters = core.utils.samples.display_molecule_protocol_parameters(
                molecules, request.user
            )
            return render(
                request,
                "wetlab/setMoleculeValues.html",
                {
                    "added_molecule_protocol_parameters": added_molecule_protocol_parameters,
                    "show_molecule_parameters": show_molecule_parameters,
                },
            )
        else:
            return render(
                request,
                "wetlab/setMoleculeValues.html",
                {
                    "added_molecule_protocol_parameters": added_molecule_protocol_parameters
                },
            )

    else:
        register_user = request.user.username
        display_list = core.utils.samples.get_defined_samples(register_user)

        return render(
            request,
            "wetlab/setMoleculeValues.html",
            {"display_list": display_list},
        )
    return render(request, "wetlab/setMoleculeValues.html", {})


@login_required
def create_pool(request):
    # Check user == WETLAB_MANAGER: if false,  redirect to 'login' page
    if request.user.is_authenticated:
        if not wetlab.utils.common.is_wetlab_manager(request):
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "You do not have enough privileges to see this page ",
                        "Contact with your administrator .",
                    ]
                },
            )
    else:
        # redirect to login webpage
        return redirect("/accounts/login")
    # collect the information for collecting
    display_list = get_lib_prep_to_select_in_pool()

    if request.method == "POST" and request.POST["action"] == "createPool":
        new_pool = define_new_pool(request.POST, request.user)

        if not isinstance(new_pool, LibraryPool):
            display_list.update(new_pool)
            return render(
                request,
                "wetlab/createPool.html",
                {"display_list": display_list},
            )
        information_for_created_pool = get_info_to_display_created_pool(new_pool)
        return render(
            request,
            "wetlab/createPool.html",
            {"information_for_created_pool": information_for_created_pool},
        )

    else:
        return render(
            request, "wetlab/createPool.html", {"display_list": display_list}
        )


@login_required
def create_new_run(request):
    if request.user.is_authenticated:
        if not wetlab.utils.common.is_wetlab_manager(request):
            return render(
                request,
                "wetlab/error_page.html",
                {
                    "content": [
                        "You do not have enough privileges to see this page ",
                        "Contact with your administrator .",
                    ]
                },
            )
    else:
        # redirect to login webpage
        return redirect("/accounts/login")

    if request.method == "POST" and request.POST["action"] == "createNewRun":
        display_pools_for_run = display_available_pools()
        if "poolID" not in request.POST:
            error_message = wetlab.config.ERROR_NO_POOL_WAS_SELECTED_IN_FORM
            return render(
                request,
                "wetlab/CreateNewRun.html",
                {
                    "display_pools_for_run": display_pools_for_run,
                    "ERROR": error_message,
                },
            )
        compatibility = check_valid_data_for_creation_run(request.POST, request.user)
        if "ERROR" in compatibility:
            return render(
                request,
                "wetlab/CreateNewRun.html",
                {
                    "display_pools_for_run": display_pools_for_run,
                    "ERROR": compatibility["ERROR"],
                },
            )
        # compatible_index = check_index_compatible(lib_prep_ids)
        display_sample_information = (
            create_run_in_pre_recorded_and_get_data_for_confirmation(
                request.POST, request.user
            )
        )

        return render(
            request,
            "wetlab/CreateNewRun.html",
            {"display_sample_information": display_sample_information},
        )

    elif request.method == "POST" and request.POST["action"] == "continueWithRun":
        run_id = request.POST["run_ids"]
        experiment_name = wetlab.utils.common.get_experiment_name(run_id)
        pool_objs = LibraryPool.objects.filter(run_process_id__exact=run_id)
        pool_ids = []
        for pool in pool_objs:
            pool_ids.append(pool.get_id())
        lib_prep_ids = get_library_prep_in_pools(pool_ids)

        display_sample_information = get_library_preparation_data_in_run(
            lib_prep_ids, pool_ids
        )
        display_sample_information.update(get_stored_user_sample_sheet(lib_prep_ids))
        display_sample_information["experiment_name"] = experiment_name
        display_sample_information["run_process_id"] = run_id
        return render(
            request,
            "wetlab/CreateNewRun.html",
            {"display_sample_information": display_sample_information},
        )

    elif request.method == "POST" and request.POST["action"] == "storeDataNewRun":
        run_obj = get_run_obj_from_id(request.POST["run_process_id"])
        if run_obj.get_state() != "Pre-Recorded":
            exp_name = run_obj.get_run_name()
            error_message = wetlab.config.ERROR_RUN_NAME_CREATED_ALREADY.copy()
            error_message.insert(1, exp_name)
            display_pools_for_run = display_available_pools()
            return render(
                request,
                "wetlab/CreateNewRun.html",
                {
                    "display_pools_for_run": display_pools_for_run,
                    "ERROR": error_message,
                },
            )
        run_data = collect_data_and_update_library_preparation_samples_for_run(
            request.POST, request.user
        )

        projects_objs = create_new_projects_added_to_run(
            run_data["projects"], run_data["run_obj"], request.user
        )
        if "ERROR" in projects_objs:
            display_pools_for_run = display_available_pools()
            return render(
                request,
                "wetlab/CreateNewRun.html",
                {
                    "display_pools_for_run": display_pools_for_run,
                    "ERROR": projects_objs["ERROR"],
                },
            )

        run_obj.set_run_state("Recorded")

        sample_sheet_name = store_confirmation_sample_sheet(run_data)
        # update the sample state for each one in the run
        pools_obj = LibraryPool.objects.filter(run_process_id=run_obj)

        for pool_obj in pools_obj:
            pool_obj.set_pool_state("Used")
        update_batch_lib_prep_sample_state(run_data["lib_prep_ids"], "Sequencing")
        created_new_run = {}
        created_new_run["exp_name"] = run_data["exp_name"]
        created_new_run["run_process_id"] = request.POST["run_process_id"]
        created_new_run["sample_sheet"] = sample_sheet_name

        return render(
            request,
            "wetlab/CreateNewRun.html",
            {"created_new_run": created_new_run},
        )
    else:
        display_pools_for_run = display_available_pools()
        return render(
            request,
            "wetlab/CreateNewRun.html",
            {"display_pools_for_run": display_pools_for_run},
        )


@login_required
def pending_sample_preparations(request):
    pending = {}
    # get the samples in defined state
    pending["defined"] = core.utils.samples.get_samples_in_defined_state("")
    pending["extract_molecule"] = core.utils.samples.get_samples_in_extracted_molecule_state("")
    pending["create_library_preparation"] = get_samples_in_lib_prep_state()
    # pending['lib_prep_protocols'] = get_protocol_lib()
    # get the library preparation in defined state
    pending["add_lib_prep_parameters"] = get_lib_prep_to_add_parameters()
    pending["graphic_pending_samples"] = pending_samples_for_grafic(pending).render()
    return render(
        request, "wetlab/pendingSamplePreparations.html", {"pending": pending}
    )


@login_required
def compare_samples(request):
    user_is_wetlab_manager = wetlab.utils.common.is_wetlab_manager(request)

    import pdb; pdb.set_trace()
    samples_data = get_list_of_samples_in_projects(request.user, user_is_wetlab_manager)
    samples_data["user"] = request.user.username
    if request.method == "POST" and request.POST["action"] == "compareSamples":
        selected_sample_objs = analyze_compare_samples_form(request.POST["table_data"])
        if len(selected_sample_objs) == 0:
            error_message = wetlab.config.ERROR_NO_SAMPLES_SELECTED
            return render(
                request,
                "wetlab/compareSamples.html",
                {"ERROR": error_message, "samples_data": samples_data},
            )
        compared_data = get_comparation_sample_information(selected_sample_objs)

        return render(
            request,
            "wetlab/compareSamples.html",
            {"compared_data": compared_data},
        )
    else:
        return render(
            request,
            "wetlab/compareSamples.html",
            {"samples_data": samples_data},
        )


@login_required
def user_commercial_kit_inventory(request):
    expired_kit = get_expired_lot_user_kit(request.user)
    valid_kit = get_valid_lot_user_kit(request.user)
    if request.method == "POST" and request.POST["action"] == "runOutUserLotKit":
        selected_user_kits = request.POST.getlist("userKit")
        if len(selected_user_kits) == 0:
            return render(
                request,
                "wetlab/userCommercialKitInventory.html",
                {
                    "expired_kit": expired_kit,
                    "valid_kit": valid_kit,
                    "user_name": request.user.username,
                },
            )
        run_out_kits = set_user_lot_kit_to_run_out(selected_user_kits)
        return render(
            request,
            "wetlab/userCommercialKitInventory.html",
            {"run_out_kits": run_out_kits},
        )

    else:
        return render(
            request,
            "wetlab/userCommercialKitInventory.html",
            {
                "expired_kit": expired_kit,
                "valid_kit": valid_kit,
                "user_name": request.user.username,
            },
        )


@login_required
def search_user_lot_kit(request):
    protocol_list = core.utils.protocols.display_protocol_list()
    platform_list = core.utils.platforms.get_defined_platforms_and_ids("NGS")
    if request.method == "POST" and request.POST["action"] == "searchuserkit":
        if (
            request.POST["expired"] == ""
            and request.POST["lotNumber"] == ""
            and request.POST["commercial"] == ""
            and request.POST["protocol"] == ""
            and request.POST["platform"] == ""
        ) and "exclude_runout" not in request.POST:
            return render(
                request,
                "wetlab/searchUserLotKit.html",
                {
                    "protocol_list": protocol_list,
                    "platform_list": platform_list,
                },
            )

        if request.POST["expired"] != "" and not wetlab.utils.common.check_valid_date_format(
            request.POST["expired"]
        ):
            error_message = wetlab.config.ERROR_INVALID_FORMAT_FOR_DATES
            return render(
                request,
                "wetlab/searchUserLotKit.html",
                {
                    "protocol_list": protocol_list,
                    "platform_list": platform_list,
                    "ERROR": error_message,
                },
            )

        user_kits_objs = search_user_lot_kit_from_user_form(request.POST)
        if user_kits_objs == "No defined":
            error_message = wetlab.config.ERROR_NO_USER_LOT_KIT_DEFINED
            return render(
                request,
                "wetlab/searchUserLotKit.html",
                {
                    "protocol_list": protocol_list,
                    "platform_list": platform_list,
                    "ERROR": error_message,
                },
            )
        if len(user_kits_objs) > 1:
            display_user_kit_list = display_user_lot_kit_information_from_query_list(
                user_kits_objs
            )
            return render(
                request,
                "wetlab/searchUserLotKit.html",
                {"display_user_kit_list": display_user_kit_list},
            )
        elif len(user_kits_objs) == 0:
            error_message = wetlab.config.ERROR_NO_MATCHES_FOR_USER_LOT_KIT
            return render(
                request,
                "wetlab/searchUserLotKit.html",
                {
                    "protocol_list": protocol_list,
                    "platform_list": platform_list,
                    "ERROR": error_message,
                },
            )
        else:
            display_one_user_kit = user_kits_objs[0].get_user_lot_kit_id()
        return redirect("display_user_lot_kit", user_kit_id=display_one_user_kit)
    else:
        return render(
            request,
            "wetlab/searchUserLotKit.html",
            {"protocol_list": protocol_list, "platform_list": platform_list},
        )


@login_required
def display_user_lot_kit(request, user_kit_id):
    user_kit_obj = get_user_lot_commercial_kit_obj_from_id(user_kit_id)
    if user_kit_obj is None:
        return render(
            request,
            "wetlab/error_page.html",
            {"content": ["Invalid User Lot Commercial Kit"]},
        )
    user_lot_kit_data = get_user_lot_kit_data_to_display(user_kit_obj)
    return render(
        request,
        "wetlab/displayUserLotKit.html",
        {"user_lot_kit_data": user_lot_kit_data},
    )


@login_required
def sequencer_configuration(request):
    if not request.user.is_authenticated:
        # redirect to login webpage
        return redirect("/accounts/login")

    if not wetlab.utils.common.is_wetlab_manager(request):
        return render(
            request,
            "wetlab/error_page.html",
            {"content": wetlab.config.ERROR_USER_NOT_WETLAB_MANAGER},
        )
    sequencer_info = get_list_sequencer_configuration()
    sequencer_info["platforms"] = get_platform_data()
    sequencer_info["sequencer_names"] = get_defined_sequencers()

    if request.method == "POST" and request.POST["action"] == "addNewSequencer":
        new_sequencer = define_new_sequencer(request.POST)
        if "ERROR" in new_sequencer:
            return render(
                request,
                "wetlab/sequencerConfiguration.html",
                {"sequencer_info": sequencer_info, "ERROR": new_sequencer},
            )
        return render(
            request,
            "wetlab/sequencerConfiguration.html",
            {"sequencer_info": sequencer_info, "new_defined_sequencer": new_sequencer},
        )
    if request.method == "POST" and request.POST["action"] == "addNewConfiguration":
        new_defined_configuration = define_new_seq_configuration(request.POST)
        if "ERROR" in new_defined_configuration:
            return render(
                request,
                "wetlab/sequencerConfiguration.html",
                {"sequencer_info": sequencer_info, "ERROR": new_defined_configuration},
            )
        return render(
            request,
            "wetlab/sequencerConfiguration.html",
            {
                "sequencer_info": sequencer_info,
                "new_defined_configuration": new_defined_configuration,
            },
        )
    else:
        return render(
            request,
            "wetlab/sequencerConfiguration.html",
            {"sequencer_info": sequencer_info},
        )


@login_required
def sequencer_inventory(request):
    if not request.user.is_authenticated:
        # redirect to login webpage
        return redirect("/accounts/login")
    if request.method == "POST" and request.POST["action"] == "setRunOutDate":
        pass
    else:
        sequencer_data = get_sequencer_inventory_data()
        return render(
            request,
            "wetlab/sequencerInventory.html",
            {"sequencer_data": sequencer_data},
        )

    return
