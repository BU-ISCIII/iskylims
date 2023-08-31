# Generic imports
import json
import logging
import os
import re
import statistics
import time
from collections import OrderedDict  # noqa

import django.contrib.auth.models
import pandas as pd
import smb
from django.conf import settings
from django.contrib.auth.decorators import login_required
from django.core.files.storage import FileSystemStorage
from django.shortcuts import redirect, render

# Local imports
# import core.fusioncharts.fusioncharts
import core.utils.commercial_kits
import core.utils.common
import core.utils.load_batch
import core.utils.platforms
import core.utils.protocols
import core.utils.samples
import django_utils.models
import django_utils.views
import wetlab.config
import wetlab.models
import wetlab.utils.additional_kits
import wetlab.utils.collection_index
import wetlab.utils.common
import wetlab.utils.crontab_process
import wetlab.utils.fetch_info
import wetlab.utils.library
import wetlab.utils.pool
import wetlab.utils.reports
import wetlab.utils.run
import wetlab.utils.sample
import wetlab.utils.samplesheet
import wetlab.utils.sequencers
import wetlab.utils.statistics
import wetlab.utils.stats_graphs
import wetlab.utils.test_conf


def index(request):
    org_name = wetlab.utils.common.get_configuration_from_database("ORGANIZATION_NAME")
    return render(request, "wetlab/index.html", {"organization_name": org_name})


@login_required
def configuration_email(request):
    if request.user.username != "admin":
        return redirect("/wetlab")
    email_conf_data = core.utils.common.get_email_data()
    email_conf_data[
        "EMAIL_ISKYLIMS"
    ] = wetlab.utils.common.get_configuration_from_database("EMAIL_FOR_NOTIFICATIONS")
    if request.method == "POST" and (request.POST["action"] == "emailconfiguration"):
        result_email = core.utils.common.send_test_email(request.POST)
        if result_email != "OK":
            email_conf_data = core.utils.common.get_email_data()
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
        for field in wetlab.config.SAMBA_CONFIGURATION_FIELDS:
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
            conn = wetlab.utils.common.open_samba_connection()
        except OSError:
            error_message = wetlab.config.ERROR_WRONG_SAMBA_CONFIGURATION_SETTINGS
            return render(
                request,
                "wetlab/configuration_samba.html",
                {"samba_conf_data": samba_user_field, "error_message": error_message},
            )
        try:
            conn.listPath(samba_user_field["shared_folder_name"], "/")
            return render(
                request,
                "wetlab/configuration_samba.html",
                {"succesful_settings": True},
            )

        except smb.base.NotReadyError:
            error_message = wetlab.config.ERROR_WRONG_SAMBA_AUTHENTICATION_SETTINGS
        except Exception:
            error_message = wetlab.config.ERROR_WRONG_SAMBA_FOLDER_SETTINGS
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
        config_file = os.path.join(settings.BASE_DIR, "wetlab", "config.py")
        test_results[
            "iSkyLIMS_settings"
        ] = wetlab.utils.test_conf.get_iSkyLIMS_settings()
        test_results["config_file"] = wetlab.utils.test_conf.get_config_file(
            config_file
        )
        test_results["attr_files"] = wetlab.utils.test_conf.get_files_attribute(
            os.path.join(settings.MEDIA_ROOT, "wetlab")
        )
        test_results["database_access"] = wetlab.utils.test_conf.check_access_database()
        test_results[
            "samba_connection"
        ] = wetlab.utils.test_conf.check_samba_connection()

        test_results["basic_checks_ok"] = "OK"
        for result in test_results:
            if test_results[result] == "NOK":
                test_results["basic_checks_ok"] = "NOK"
                break
        available_run_test = []
        if wetlab.models.RunConfigurationTest.objects.all().exists():
            run_test_objs = wetlab.models.RunConfigurationTest.objects.all()
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
        if wetlab.models.RunConfigurationTest.objects.filter(
            pk__exact=request.POST["runTest"]
        ).exists():
            run_test_obj = wetlab.models.RunConfigurationTest.objects.filter(
                pk__exact=request.POST["runTest"]
            ).last()
            run_test_folder = run_test_obj.get_run_test_folder()
            run_test_name = run_test_obj.get_run_test_name()
        if not wetlab.utils.test_conf.folder_test_exists(run_test_folder):
            return render(
                request,
                "wetlab/configuration_test.html",
                {"error": wetlab.config.ERROR_NOT_FOLDER_RUN_TEST_WAS_FOUND},
            )
        run_test_result = wetlab.utils.test_conf.execute_test_for_testing_run(
            run_test_name
        )
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
        if wetlab.models.RunConfigurationTest.objects.filter(
            run_test_name__exact=request.POST["deleteRun"]
        ).exists():
            if wetlab.models.RunProcess.objects.filter(
                run_name__exact=request.POST["deleteRun"]
            ).exists():
                run_test_objs = wetlab.models.RunProcess.objects.filter(
                    run_name__exact=request.POST["deleteRun"]
                )
                for run_test_obj in run_test_objs:
                    wetlab.utils.test_conf.delete_test_run(run_test_obj)
                    delete_successful = {"run_name": request.POST["deleteRun"]}
                return render(
                    request,
                    "wetlab/configuration_test.html",
                    {"delete_successful": delete_successful},
                )
            return render(request, "wetlab/configuration_test.html")
    else:
        return render(request, "wetlab/configuration_test.html")


@login_required
def crontab_status(request):
    # check user privileges
    if request.user.is_authenticated:
        if request.user.username != "admin":
            return redirect("/wetlab")
    if request.method == "POST" and request.POST["action"] == "set_crontab_state":
        state = "activate" if "cron_status" in request.POST else "remove"
        c_status = wetlab.utils.crontab_process.set_crontab_status(state)
        c_status["updated"] = wetlab.config.SUCCESSFUL_CRONTAB_STATUS_CHANGED
        return render(request, "wetlab/crontab_status.html", {"c_status": c_status})
    else:
        c_status = wetlab.utils.crontab_process.get_crontab_status()
        return render(request, "wetlab/crontab_status.html", {"c_status": c_status})


@login_required
def initial_settings(request):
    if not wetlab.utils.common.is_wetlab_manager(request):
        return render(
            request,
            "wetlab/initial_settings.html",
            {"error_message": wetlab.config.ERROR_USER_NOT_WETLAB_MANAGER},
        )
    initial_data = core.utils.common.get_inital_sample_settings_values(__package__)
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
        new_inital_data = core.utils.common.save_inital_sample_setting_value(
            __package__, form_data
        )
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
            "wetlab/create_next_seq_run.html",
            {"error_message": wetlab.config.ERROR_USER_NOT_WETLAB_MANAGER},
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
                "wetlab/create_next_seq_run.html",
                {"error_message": "Uploaded file does not containt csv extension"},
            )
        ext_file = split_filename.group(2)

        # CHECK if file contains the csv extension.
        # Error page is shown if file does not contain the csv extension
        if ext_file != ".csv":
            return render(
                request,
                "wetlab/create_next_seq_run.html",
                {"error_message": "Sample Sheet must have a csv extension"},
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
        index_library_name = wetlab.utils.samplesheet.get_index_library_name(
            stored_file
        )
        run_name = wetlab.utils.common.get_experiment_name_from_file(stored_file)

        if run_name == "":
            # define an temporary unique value for the run name
            # until the real value is get from user FORM
            run_name = timestr

        # Check that runName is not already used in the database.
        # Error page is showed if runName is already  defined

        if (
            wetlab.models.RunProcess.objects.filter(run_name__iexact=run_name)
        ).exists():
            if wetlab.models.RunProcess.objects.filter(
                run_name__iexact=run_name, state__run_state_name__exact="Pre-Recorded"
            ).exists():
                # Delete the Sample Sheet file and the row in database
                delete_run_objs = wetlab.models.RunProcess.objects.filter(
                    run_name__iexact=run_name,
                    state__run_state_name__exact="Pre-Recorded",
                )
                for delete_run in delete_run_objs:
                    # sample_sheet_file = delete_run.get_sample_file()
                    # full_path_sample_sheet_file = os.path.join(settings.MEDIA_ROOT, sample_sheet_file)
                    # os.remove(full_path_sample_sheet_file)

                    if wetlab.models.Projects.objects.filter(
                        run_process=delete_run
                    ).exists():
                        project_objs = wetlab.models.Projects.objects.filter(
                            run_process=delete_run
                        )
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
                    "wetlab/create_next_seq_run.html",
                    {"error_message": "Run Name is already used. "},
                )

        # Fetch from the Sample Sheet file the projects included in
        # the run and the user. Error page is showed if not project/description
        # colunms are found

        project_list = wetlab.utils.samplesheet.get_projects_in_run(stored_file)

        if len(project_list) == 0:
            # delete sample sheet file
            fs.delete(file_name)
            return render(
                request,
                "wetlab/create_next_seq_run.html",
                {
                    "error_message": "Sample Sheet does not contain Sample project and/or Description fields"
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
            if wetlab.models.Projects.objects.filter(
                project_name__icontains=key
            ).exists():
                project_already_defined.append(key)

        if len(project_already_defined) > 0 and project_in_several_runs != "TRUE":
            # convert the list into string to display the user names on error page
            display_project = "  ".join(project_already_defined)
            # delete sample sheet file before showing the error page
            fs.delete(file_name)
            return render(
                request,
                "wetlab/create_next_seq_run.html",
                {"error_message": display_project + " already defined"},
            )

        # Once the information looks good. it will be stores in runProcess and projects table

        # store data in runProcess table, run is in pre-recorded state
        center_requested_id = django_utils.models.Profile.objects.get(
            profile_user_id=request.user
        ).profile_center.id
        center_requested_by = django_utils.models.Center.objects.get(
            pk=center_requested_id
        )
        new_run_obj = wetlab.models.RunProcess(
            run_name=run_name,
            sample_sheet=file_name,
            state=wetlab.models.RunStates.objects.get(
                run_state_name__exact="Pre-Recorded"
            ),
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
            if django.contrib.auth.models.User.objects.filter(
                username__exact=val
            ).exists():
                user_id = django.contrib.auth.models.User.objects.get(
                    username__exact=val
                )
            else:
                user_id = None

            if not wetlab.models.wetlab.models.Projects.objects.filter(
                project_name__iexact=key
            ).exists():
                data = {}
                data["user_id"] = user_id
                data["projectName"] = key
                project_obj = wetlab.models.Projects.objects.create_new_empty_project(
                    data
                )
            else:
                project_obj = wetlab.models.Projects.objects.filter(
                    project_name__iexact=key
                ).last()

            project_obj.add_run(new_run_obj)
            projects.append([key, val])

        run_info_values["projects_user"] = projects
        run_info_values["runname"] = run_name
        # Get the list of the library kit used (libraryKit)
        list_libraries = wetlab.models.LibraryKit.objects.order_by().values_list(
            "library_name", flat=True
        )
        run_info_values["used_libraryKit"] = list_libraries

        user_names = []
        all_users = django.contrib.auth.models.User.objects.all()
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
        if not wetlab.models.RunProcess.objects.filter(
            run_name__exact=run_name
        ).exists():
            return render(
                request,
                "wetlab/create_next_seq_run.html",
                {"error_message": "Reload again the page using the menu options"},
            )
        run_p = wetlab.models.RunProcess.objects.get(run_name__exact=run_name)
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
        wetlab.utils.samplesheet.create_unique_sample_id_values(in_file, index_file)
        # create the projects/users to update sample sheet
        user_names_in_projects = {}
        for p_index in range(len(projects)):
            user_names_in_projects[projects[p_index]] = user_name[p_index]

        wetlab.utils.samplesheet.set_user_names_in_sample_sheet(
            in_file, user_names_in_projects
        )
        # build the project list for each project_library kit
        for x in range(len(project_index_kit)):
            if project_index_kit[x] in library:
                library[project_index_kit[x]].append(projects[x])
            else:
                library[project_index_kit[x]] = [projects[x]]
        # save the project information on database
        for p in range(len(projects)):
            my_project = projects[p]
            library_kit_id = wetlab.models.LibraryKit.objects.filter(
                library_name__exact=library_kit[p]
            ).last()
            update_info_proj = wetlab.models.Projects.objects.get(
                project_name=my_project
            )
            update_info_proj.libraryKit = project_index_kit[p]
            # removed the link to base space file
            # update_info_proj.baseSpaceFile=bs_file[project_index_kit[p]]
            update_info_proj.baseSpaceFile = None
            update_info_proj.LibraryKit_id = library_kit_id
            update_info_proj.user_id = django.contrib.auth.models.User.objects.get(
                username__exact=user_name[p]
            )
            update_info_proj.save()
        results.append(["runname", experiment_name])
        run_p.set_run_state("Recorded")
        sample_sheet_lines = wetlab.utils.samplesheet.read_all_lines_in_sample_sheet(
            in_file
        )
        sample_names_and_data = wetlab.utils.samplesheet.get_samples_in_sample_sheet(
            sample_sheet_lines
        )
        wetlab.utils.run.increase_reuse_if_samples_exists(
            sample_names_and_data["samples"]
        )

        return render(
            request,
            "wetlab/create_next_seq_run.html",
            {"completed_form": results},
        )

    return render(request, "wetlab/create_next_seq_run.html")


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

    collection_indexes = wetlab.models.CollectionIndexKit.objects.all()
    if len(collection_indexes) > 0:
        for c_index in collection_indexes:
            collection_index_names.append(
                [c_index.get_id(), c_index.get_collection_index_name]
            )
    collection_index_information["collection_index"] = collection_index_names
    if request.method == "POST" and request.POST["action"] == "addCollectionIndexKit":
        # fetch the file from user form and  build the file name  including
        # the date and time on now to store in database

        file_name = request.FILES["newCollectionIndexFile"].name
        saved_file = wetlab.utils.collection_index.store_collection_kits_file(
            request.FILES["newCollectionIndexFile"]
        )

        # get the libary name to check if it is already defined
        if not wetlab.utils.collection_index.check_collection_index_file_format(
            os.path.join(settings.MEDIA_ROOT, saved_file)
        ):
            os.remove(saved_file)
            return render(
                request,
                "wetlab/add_collection_index_kit.html",
                {
                    "list_of_collection_index": collection_index_information,
                    "error_message": file_name + " does not have the right format",
                },
            )

        collection_name = wetlab.utils.collection_index.get_collection_index_name(
            os.path.join(settings.MEDIA_ROOT, saved_file)
        )
        if collection_name == "":
            # removing the uploaded file
            os.remove(saved_file)
            return render(
                request,
                "wetlab/add_collection_index_kit.html",
                {
                    "list_of_collection_index": collection_index_information,
                    "error_message": file_name + " does not contain kit  name",
                },
            )

        # check if library name is already defined on database
        if wetlab.utils.collection_index.check_collection_index_exists(collection_name):
            # removing the uploaded file
            os.remove(os.path.join(settings.MEDIA_ROOT, saved_file))
            return render(
                request,
                "wetlab/add_collection_index_kit.html",
                {
                    "list_of_collection_index": collection_index_information,
                    "error_message": file_name + " is already defined on iSkyLIMS",
                },
            )
        # Get the collection settings included in the file
        collection_settings = wetlab.utils.collection_index.get_collection_settings(
            os.path.join(settings.MEDIA_ROOT, saved_file)
        )
        new_collection_obj = wetlab.utils.collection_index.store_collection_settings(
            collection_settings, saved_file
        )
        # get the index name and index bases for the library
        collection_index = wetlab.utils.collection_index.get_index_values(
            os.path.join(settings.MEDIA_ROOT, saved_file)
        )
        wetlab.utils.collection_index.store_collection_indexes(
            collection_index, new_collection_obj
        )

        collection_index_information["collection_index_names"] = collection_settings[
            "name"
        ]

        return render(
            request,
            "wetlab/add_collection_index_kit.html",
            {"collection_index_information": collection_index_information},
        )
    else:
        return render(
            request,
            "wetlab/add_collection_index_kit.html",
            {"list_of_collection_index": collection_index_information},
        )


@login_required
def search_run(request):
    """The function is called from web, having 2 main parts:
            - User form with the information to search runs
            - Result information can be :
                - list of the matched runs
                - run information in case that only 1 match is found
    Parameters
    ----------
    request : _type_
        _description_

    Returns
    -------
    render: _type_
        Contains dictionnary with run data
    """
    # check user privileges
    if wetlab.utils.common.is_wetlab_manager(request):
        allowed_all_runs = True
    else:
        allowed_all_runs = False

    # Search for runs that fullfil the input values
    ########
    run_form_data = wetlab.utils.common.get_run_search_fields_form()
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
            return render(
                request, "wetlab/search_run.html", {"run_form_data": run_form_data}
            )

        # check the right format of start and end date

        if start_date != "":
            if not wetlab.utils.common.check_valid_date_format(start_date):
                error_message = wetlab.config.ERROR_INVALID_FORMAT_FOR_DATES
                return render(
                    request,
                    "wetlab/search_run.html",
                    {"run_form_data": run_form_data, "error_message": error_message},
                )
        if end_date != "":
            if not wetlab.utils.common.check_valid_date_format(end_date):
                error_message = wetlab.config.ERROR_INVALID_FORMAT_FOR_DATES
                return render(
                    request,
                    "wetlab/search_run.html",
                    {"run_form_data": run_form_data, "error_message": error_message},
                )

        # Get all the available runs to start the filtering
        if allowed_all_runs:
            runs_found = (
                wetlab.models.RunProcess.objects.all().order_by("run_date").reverse()
            )
        else:
            user_ids = wetlab.utils.common.get_allowed_user_for_sharing(request.user)
            user_projects = wetlab.models.wetlab.models.Projects.objects.filter(
                user_id__in=user_ids
            )
            run_list = []
            for user_project in user_projects:
                # run_list.append(user_project.runprocess_id.id)
                run_objs = user_project.run_process.all()
                for run_obj in run_objs:
                    run_list.append(run_obj.get_run_id())
            if wetlab.models.RunProcess.objects.filter(pk__in=run_list).exists():
                runs_found = wetlab.models.RunProcess.objects.filter(pk__in=run_list)
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
            if wetlab.models.RunProcess.objects.filter(
                run_name__iexact=run_name
            ).exists():
                run_name_found = wetlab.models.RunProcess.objects.filter(
                    run_name__iexact=run_name
                )
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
            sequencer_list = wetlab.utils.fetch_info.get_sequencer_names_from_platform(
                platform_name
            )
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
            s_state = wetlab.models.RunStates.objects.get(
                run_state_name__exact=run_state
            )
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
        elif len(runs_found) == 0:
            error_message = wetlab.config.ERROR_NO_MATCHES_FOR_RUN_SEARCH
            return render(
                request,
                "wetlab/search_run.html",
                {"run_form_data": run_form_data, "error_message": error_message},
            )

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
    project_form_data = wetlab.utils.common.get_project_search_fields_form()
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
            if not wetlab.utils.common.check_valid_date_format(start_date):
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
            if not wetlab.utils.common.check_valid_date_format(start_date):
                error_message = wetlab.config.ERROR_INVALID_FORMAT_FOR_DATES
                return render(
                    request,
                    "wetlab/search_project.html",
                    {
                        "project_form_data": project_form_data,
                        "error_message": error_message,
                    },
                )

        projects_found = wetlab.models.Projects.objects.all()

        if project_name != "":
            projects_found = projects_found.filter(project_name__icontains=project_name)
        if sequencer_name != "":
            run_objs = wetlab.models.RunProcess.objects.filter(
                used_sequencer__sequencer_name__exact=sequencer_name
            )
            projects_found = projects_found.filter(run_process__in=run_objs)
        if run_state != "":
            run_objs = wetlab.models.RunProcess.objects.filter(
                state__run_state_name__exact=run_state
            )
            projects_found = projects_found.filter(run_process__in=run_objs)
        if user_name != "":
            # check if user has a shared user
            if user_name == request.user.username:
                p_shared_list = wetlab.utils.common.get_allowed_user_for_sharing(
                    request.user
                )
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
    if not wetlab.utils.common.is_wetlab_manager(request):
        return render(
            request,
            "wetlab/display_run.html",
            {"error_message": "You do have the enough privileges to see this page "},
        )
    if request.method == "POST" and (request.POST["action"] == "retry_correct_error"):
        run_id = request.POST["run_id"]
        if wetlab.models.RunProcess.objects.filter(pk__exact=run_id).exists():
            run_name_found = wetlab.models.RunProcess.objects.get(pk__exact=run_id)
            previous_error_state = run_name_found.get_state_before_error()
            run_name_found.set_run_state(previous_error_state)
            detail_description = {}
            detail_description[
                "information"
            ] = wetlab.config.SUCCESSFUL_RUN_STATE_CHANGE_FOR_RETRY
            return render(
                request,
                "wetlab/display_run.html",
                {
                    "successful_correction": wetlab.config.SUCCESSFUL_RUN_STATE_CHANGE_FOR_RETRY,
                    "run_id": run_id,
                },
            )
        else:
            return render(
                request,
                "wetlab/display_run.html",
                {"error_message": "Run does not exist "},
            )
    else:
        return render(request, "wetlab/index.html")


@login_required
def skip_cancel_situation(request):
    # check user privileges
    if not wetlab.utils.common.is_wetlab_manager(request):
        return render(
            request,
            "wetlab/skip_cancel_situation.html",
            {"error_message": "You do have the enough privileges to see this page "},
        )
    if request.method == "POST" and (request.POST["action"] == "skip_cancel_situation"):
        run_id = request.POST["run_id"]
        if wetlab.models.RunProcess.objects.filter(pk__exact=run_id).exists():
            run_name_found = wetlab.models.RunProcess.objects.get(pk__exact=run_id)
            run_name_found.set_run_state("Sample Sent")
            run_name_found.set_forced_continue_on_error()
            detail_description = {}
            detail_description[
                "information"
            ] = wetlab.config.SUCCESSFUL_RUN_STATE_CHANGE_FOR_RETRY
            return render(
                request,
                "wetlab/skip_cancel_situation.html",
                {"detail_description": detail_description},
            )
        else:
            return render(
                request,
                "wetlab/skip_cancel_situation.html",
                {"error_message": ["Run does not exist"]},
            )
    else:
        # return redirect (request,'/')
        return render(request, "wetlab/index.html")


@login_required
def display_run(request, run_id):
    if not wetlab.models.wetlab.models.RunProcess.objects.filter(
        pk__exact=run_id
    ).exists():
        return render(
            request,
            "wetlab/display_run.html",
            {"error_message": wetlab.config.ERROR_RUN_DOES_NOT_EXIST},
        )
    # check user privileges
    if not wetlab.utils.common.is_wetlab_manager(request):
        # check if user is owner of the run or belongs to the shared user
        shared_user_ids = wetlab.utils.common.get_allowed_user_for_sharing(request.user)
        project_objs = wetlab.models.wetlab.models.Projects.objects.filter(
            run_process__exact=run_id
        )
        allowed = False
        for project_obj in project_objs:
            if int(project_obj.get_user_id()) in shared_user_ids:
                allowed = True
                break
        if not allowed:
            return render(
                request,
                "wetlab/display_run.html",
                {"error_message": wetlab.config.ERROR_USER_NOT_WETLAB_MANAGER},
            )
    run_name_found = wetlab.models.wetlab.models.RunProcess.objects.get(pk=run_id)
    r_data_display = wetlab.utils.fetch_info.get_information_run(run_name_found)
    return render(
        request,
        "wetlab/display_run.html",
        {"display_one_run": r_data_display},
    )


@login_required
def last_run_by_sequencer(request):
    # check user privileges
    if not wetlab.utils.common.is_wetlab_manager(request):
        return render(
            request,
            "wetlab/last_run_by_sequencer.html",
            {"error_message": "You do have the enough privileges to see this page "},
        )
    last_runs = wetlab.utils.fetch_info.get_last_runs_by_sequencer()
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
    if not wetlab.utils.common.is_wetlab_manager(request):
        return render(
            request,
            "wetlab/incompleted_runs.html",
            {"error_message": "You do have the enough privileges to see this page "},
        )
    if (
        wetlab.models.wetlab.models.RunProcess.objects.all()
        .exclude(state__run_state_name="Completed")
        .exists()
    ):
        display_incompleted_run = (
            wetlab.utils.fetch_info.get_information_for_incompleted_run()
        )
        return render(
            request,
            "wetlab/incompleted_runs.html",
            {"display_incompleted_run": display_incompleted_run},
        )
    else:
        return render(request, "wetlab/incompleted_runs.html")


@login_required
def display_project(request, project_id):
    if wetlab.models.Projects.objects.filter(pk=project_id).exists():
        project_obj = wetlab.models.Projects.objects.filter(pk=project_id).last()

        # check that user is allow to see the project
        groups = django.contrib.auth.models.Group.objects.get(
            name=wetlab.config.WETLAB_MANAGER
        )
        if groups not in request.user.groups.all():
            p_shared_list = wetlab.utils.common.get_allowed_user_for_sharing(
                request.user
            )
            if int(project_obj.get_user_id()) not in p_shared_list:
                return render(
                    request,
                    "wetlab/display_project.html",
                    {
                        "error_message": "You do have the enough privileges to see this page "
                    },
                )
        # Display the proyect information
        display_project_data = wetlab.utils.fetch_info.get_information_project(
            project_obj, request
        )
        return render(
            request,
            "wetlab/display_project.html",
            {"display_project_data": display_project_data},
        )
    else:
        return render(
            request,
            "wetlab/display_project.html",
            {"error_message": "No matches have been found for the project"},
        )


@login_required
def display_collection_index(request, collection_index_id):
    if wetlab.models.CollectionIndexKit.objects.filter(pk=collection_index_id).exists():
        collection_index_dict = (
            wetlab.utils.collection_index.get_collection_index_information(
                collection_index_id
            )
        )
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
                    "error_message": "There are recorded information for the collection index for your request"
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
                    {
                        "error_message": wetlab.config.ERROR_TOO_SHORT_INDEX_BASE_SEQUENCE
                    },
                )
            else:
                valid_seq_characters = ["a", "A", "c", "C", "g", "G", "t", "T"]
                for letter in index_sequence:
                    if letter not in valid_seq_characters:
                        return render(
                            request,
                            "wetlab/search_collection_index_library.html",
                            {
                                "error_message": wetlab.config.ERROR_INVALID_SEQUENCE_CHARACTERS
                            },
                        )

        collection_indexes = wetlab.models.CollectionIndexKit.objects.all()
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
            collection_values = wetlab.models.CollectionIndexValues.objects.all()
            if index_name != "":
                if collection_values.filter(
                    index_7__icontains=index_name,
                    collection_index_kit_id__in=collection_indexes,
                ).exists():
                    collection_values = collection_values.filter(
                        index_7__icontains=index_name,
                        collection_index_kit_id__in=collection_indexes,
                    )
                elif collection_values.filter(
                    index_5__icontains=index_name,
                    collection_index_kit_id__in=collection_indexes,
                ).exists():
                    collection_values = collection_values.filter(
                        index_5__icontains=index_name,
                        collection_index_kit_id__in=collection_indexes,
                    )
                else:
                    return render(
                        request,
                        "wetlab/search_collection_index_library.html",
                        {"not_found_matchs": "not_found_matchs"},
                    )

            if index_sequence != "":
                (
                    index_found,
                    sequence,
                ) = wetlab.utils.library.find_index_sequence_collection_values_kit(
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
    if wetlab.models.RunProcess.objects.filter(pk=run_id).exists():
        run_obj = wetlab.models.RunProcess.objects.get(pk=run_id)
        run_data = {}
        run_data["run_name"] = run_obj.get_run_name()
        run_data["run_id"] = run_id

        # check if user is allow to make the change
        groups = django.contrib.auth.models.Group.objects.filter(
            name="WetlabManager"
        ).last()
        # check if user belongs to WetlabManager . If true allow to see the page
        if groups not in request.user.groups.all():
            return render(
                request,
                "wetlab/change_run_name.html",
                {
                    "error_message": "You do have the enough privileges to see this page "
                },
            )
        if request.method == "POST" and request.POST["action"] == "change_run_name":
            new_run_name = request.POST["runName"]
            if wetlab.models.RunProcess.objects.filter(
                run_name__iexact=new_run_name
            ).exists():
                return render(
                    request,
                    "wetlab/change_run_name.html",
                    {
                        "error_message": "The given Run Name is already in use",
                        "run_data": run_data,
                    },
                )
            changed_name = {}
            run_obj.update_run_name(new_run_name)
            changed_name["new_name"] = new_run_name
            changed_name["run_id"] = run_id
            changed_name["old_name"] = run_data["run_name"]

            return render(
                request,
                "wetlab/change_run_name.html",
                {"changed_name": changed_name},
            )
        else:
            return render(
                request,
                "wetlab/change_run_name.html",
                {"run_data": run_data},
            )
    else:
        return render(
            request,
            "wetlab/change_run_name.html",
            {
                "error_message": "There is no Run for your query",
                "run_data": run_data,
            },
        )


@login_required
def stats_per_researcher(request):
    if request.method == "POST":
        r_name = request.POST["researchername"]
        start_date = request.POST["startdate"]
        end_date = request.POST["enddate"]

        researcher_statistics = wetlab.utils.statistics.get_researcher_statistics(
            r_name, start_date, end_date
        )
        if "ERROR" in researcher_statistics:
            error_message = researcher_statistics["ERROR"]
            return render(
                request,
                "wetlab/stats_per_researcher.html",
                {"error_message": error_message},
            )

        return render(
            request,
            "wetlab/stats_per_researcher.html",
            {"researcher_statistics": researcher_statistics},
        )
    else:
        return render(request, "wetlab/stats_per_researcher.html")


@login_required
def stats_per_sequencer(request):
    sequencer_names = wetlab.utils.fetch_info.get_sequencer_installed_names()
    if request.method == "POST":
        sequencer = request.POST["sequencer"]
        start_date = request.POST["startdate"]
        end_date = request.POST["enddate"]

        sequencer_statistics = wetlab.utils.statistics.get_sequencer_statistics(
            sequencer, start_date, end_date
        )
        if "ERROR" in sequencer_statistics:
            error_message = sequencer_statistics["ERROR"]
            return render(
                request,
                "wetlab/stats_per_sequencer.html",
                {
                    "error_message": error_message,
                    "sequencer_names": sequencer_names,
                },
            )
        return render(
            request,
            "wetlab/stats_per_sequencer.html",
            {"sequencer_statistics": sequencer_statistics},
        )
    else:
        return render(
            request,
            "wetlab/stats_per_sequencer.html",
            {"sequencer_names": sequencer_names},
        )


@login_required
def stats_per_time(request):
    if request.method == "POST":
        start_date = request.POST["startdate"]
        end_date = request.POST["enddate"]
        per_time_statistics = wetlab.utils.statistics.get_per_time_statistics(
            start_date, end_date
        )
        return render(
            request,
            "wetlab/stats_per_time.html",
            {"per_time_statistics": per_time_statistics},
        )
    return render(request, "wetlab/stats_per_time.html")


def get_list_of_libraries_values(
    library_found, q30_comparations, mean_comparations, n_bases_comparations
):
    for project_to_compare in library_found:
        library_to_compare_name = project_to_compare.get_index_library_name()
        # project_to_compare_id = project_to_compare.id
        q30_compare_lib, mean_compare_lib, yield_mb_compare_lib = [], [], []

        # This line must changed to handle project name is reused in several runs
        run_used_in_project = project_to_compare.runProcess.all().last()

        run_param_obj = wetlab.models.RunningParameters.objects.get(
            run_name_id=run_used_in_project
        )
        # get the number of lanes by quering the SequencerModel in the RunProcess
        # number_of_lanes = project_to_compare.runprocess_id.get_sequencing_lanes()
        number_of_lanes = int(run_param_obj.get_number_of_lanes())
        for lane_number in range(1, number_of_lanes + 1):
            try:
                lane_in_project = wetlab.models.StatsLaneSummary.objects.get(
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
def annual_report(request):
    # check user privileges
    if not wetlab.utils.common.is_wetlab_manager(request):
        return render(
            request,
            "wetlab/annual_report.html",
            {"error_message": "You do have the enough privileges to see this page "},
        )
    if request.method == "POST" and request.POST["action"] == "annualreport":
        annual_rep_data = wetlab.utils.reports.get_annual_report(
            request.POST["yearselected"]
        )
        if "ERROR" in annual_rep_data:
            return render(
                request,
                "wetlab/annual_report.html",
                {"error_message": annual_rep_data["ERROR"]},
            )
        return render(
            request,
            "wetlab/annual_report.html",
            {"annual_rep_data": annual_rep_data},
        )
    else:
        return render(request, "wetlab/annual_report.html")


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
    if not wetlab.utils.common.is_wetlab_manager(request):
        return render(
            request,
            "wetlab/create_protocol.html",
            {"error_message": "You do not have enough privileges to see this page "},
        )
    # get the list of defined protocols
    defined_protocols = core.utils.protocols.display_available_protocols(__package__)
    additional_kits = wetlab.utils.additional_kits.get_additional_kits_list(__package__)
    defined_protocol_types = core.utils.protocols.display_protocol_types(__package__)

    if request.method == "POST" and request.POST["action"] == "addNewProtocol":
        new_protocol = request.POST["newProtocolName"]
        protocol_type = request.POST["protocolType"]
        description = request.POST["description"]

        if core.utils.protocols.check_if_protocol_exists(new_protocol, __package__):
            return render(
                request,
                "wetlab/create_protocol.html",
                {"ERROR": "Protocol Name " + new_protocol + "Already exists."},
            )
        new_protocol_id = core.utils.protocols.create_new_protocol(
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
            },
        )

    return render(
        request,
        "wetlab/create_protocol.html",
        {
            "defined_protocols": defined_protocols,
            "defined_protocol_types": defined_protocol_types,
            "additional_kits": additional_kits,
        },
    )


@login_required
def define_sample_projects(request):
    if not wetlab.utils.common.is_wetlab_manager(request):
        return render(
            request,
            "wetlab/define_sample_projects.html",
            {"error_message": "You do not have enough privileges to see this page "},
        )
    # get the information of defined sample Projects
    defined_samples_projects = core.utils.samples.get_info_for_defined_sample_projects(
        __package__
    )

    if request.method == "POST" and request.POST["action"] == "addNewSampleProject":
        sample_project_name = request.POST["sampleProyectName"]
        # description = request.POST['description']

        if core.models.SampleProjects.objects.filter(
            sample_project_name__iexact=sample_project_name,
            apps_name__exact=__package__,
        ).exists():
            error_message = wetlab.config.ERROR_SAMPLE_PROJECT_ALREADY_EXISTS
            return render(
                request,
                "wetlab/define_sample_projects.html",
                {
                    "defined_samples_projects": defined_samples_projects,
                    "error_message": error_message,
                },
            )
        new_sample_project_id = core.utils.samples.create_new_sample_project(
            request.POST, __package__
        )
        new_defined_sample_project = sample_project_name
        return render(
            request,
            "wetlab/define_sample_projects.html",
            {
                "defined_samples_projects": defined_samples_projects,
                "new_sample_project_id": new_sample_project_id,
                "new_defined_sample_project": new_defined_sample_project,
            },
        )

    return render(
        request,
        "wetlab/define_sample_projects.html",
        {"defined_samples_projects": defined_samples_projects},
    )


@login_required
def define_additional_kits(request, protocol_id):
    if not wetlab.utils.common.is_wetlab_manager(request):
        return render(
            request,
            "wetlab/define_additional_kit.html",
            {"error_message": "You do not have enough privileges to see this page "},
        )
    additional_kits = wetlab.utils.additional_kits.define_table_for_additional_kits(
        protocol_id
    )
    if request.method == "POST" and request.POST["action"] == "defineAdditionalKits":
        recorded_additional_kits = wetlab.utils.additional_kits.set_additional_kits(
            request.POST, request.user
        )
        if len(recorded_additional_kits) == 0:
            return render(
                request,
                "wetlab/define_additional_kits.html",
                {"additional_kits": additional_kits},
            )

        return render(
            request,
            "wetlab/define_additional_kits.html",
            {"recorded_additional_kits": recorded_additional_kits},
        )

    else:
        if not core.utils.protocols.check_if_protocol_exists(protocol_id, __package__):
            return render(
                request,
                "wetlab/define_additional_kits.html",
                {"error_message": "The requested Protocol does not exist"},
            )

        return render(
            request,
            "wetlab/define_additional_kits.html",
            {"additional_kits": additional_kits},
        )


@login_required
def display_sample_project(request, sample_project_id):
    samples_project_data = core.utils.samples.get_info_to_display_sample_project(
        sample_project_id
    )
    if "ERROR" in samples_project_data:
        error_message = samples_project_data["ERROR"]
        return render(
            request,
            "wetlab/display_sample_project.html",
            {"error_message": error_message},
        )
    return render(
        request,
        "wetlab/display_sample_project.html",
        {"samples_project_data": samples_project_data},
    )


@login_required
def display_protocol(request, protocol_id):
    if not wetlab.utils.common.is_wetlab_manager(request):
        return render(
            request,
            "wetlab/display_protocol.html",
            {"error_message": "You do not have enough privileges to see this page "},
        )
    if not core.utils.protocols.check_if_protocol_exists(protocol_id, __package__):
        return render(
            request,
            "wetlab/display_protocol.html",
            {
                "error_message": "The protocol that you are trying to get DOES NOT exist."
            },
        )
    protocol_data = core.utils.protocols.get_all_protocol_info(protocol_id)
    kit_data = wetlab.utils.additional_kits.get_all_additional_kit_info(protocol_id)

    return render(
        request,
        "wetlab/display_protocol.html",
        {"protocol_data": protocol_data, "kit_data": kit_data},
    )


@login_required
def define_protocol_parameters(request, protocol_id):
    # Check user == WETLAB_MANAGER: if false,  redirect to 'login' page
    if not wetlab.utils.common.is_wetlab_manager(request):
        return render(
            request,
            "wetlab/define_protocol_parameters.html",
            {"error_message": "You do not have enough privileges to see this page "},
        )
    if (
        request.method == "POST"
        and request.POST["action"] == "define_protocol_parameters"
    ):
        recorded_prot_parameters = core.utils.protocols.set_protocol_parameters(request)

        return render(
            request,
            "wetlab/define_protocol_parameters.html",
            {"recorded_prot_parameters": recorded_prot_parameters},
        )

    else:
        if not core.utils.protocols.check_if_protocol_exists(protocol_id, __package__):
            return render(
                request,
                "wetlab/define_protocol_parameters.html",
                {"error_message": "The requested Protocol does not exist"},
            )

        prot_parameters = core.utils.protocols.define_table_for_prot_parameters(
            protocol_id
        )
        return render(
            request,
            "wetlab/define_protocol_parameters.html",
            {"prot_parameters": prot_parameters},
        )


@login_required
def add_commercial_kit(request):
    app_name = __package__.split(".")[0]
    defined_protocols = core.utils.protocols.get_defined_protocols(app_name, False)
    defined_platforms = core.utils.platforms.get_defined_platforms_and_ids("NGS")
    commercial_kits_data = core.utils.commercial_kits.get_data_for_commercial_kits(
        "NGS"
    )

    if request.method == "POST" and request.POST["action"] == "addCommercialKit":
        if core.utils.commercial_kits.get_commercial_kit_id(request.POST["kitName"]):
            return render(
                request,
                "wetlab/add_commercial_kit.html",
                {
                    "defined_protocols": defined_protocols,
                    "defined_platforms": defined_platforms,
                    "commercial_kits_data": commercial_kits_data,
                    "error_message": "Commercial kit "
                    + request.POST["kitName"]
                    + " is already defined",
                },
            )
        new_kit = core.utils.commercial_kits.store_commercial_kit(request.POST)
        new_kit_data = core.utils.commercial_kits.get_commercial_kit_basic_data(new_kit)
        return render(
            request,
            "wetlab/add_commercial_kit.html",
            {"new_kit_data": new_kit_data},
        )
    else:
        return render(
            request,
            "wetlab/add_commercial_kit.html",
            {
                "defined_protocols": defined_protocols,
                "defined_platforms": defined_platforms,
                "commercial_kits_data": commercial_kits_data,
            },
        )


@login_required
def add_user_lot_commercial_kit(request):
    defined_kits = core.utils.commercial_kits.get_defined_commercial_kits()
    if request.method == "POST" and request.POST["action"] == "addUserLotKit":
        if core.utils.commercial_kits.get_lot_user_commercial_kit_obj(
            request.POST["barCode"]
        ):
            return render(
                request,
                "wetlab/add_user_lot_commercial_kit.html",
                {
                    "defined_kits": defined_kits,
                    "error_message": "Lot barcode "
                    + request.POST["barCode"]
                    + " is already defined",
                },
            )
        new_lot_kit = core.utils.commercial_kits.store_lot_user_commercial_kit(
            request.POST, request.user
        )
        new_lot_kit_data = (
            core.utils.commercial_kits.get_lot_user_commercial_kit_basic_data(
                new_lot_kit
            )
        )
        return render(
            request,
            "wetlab/add_user_lot_commercial_kit.html",
            {"new_lot_kit_data": new_lot_kit_data},
        )
    else:
        return render(
            request,
            "wetlab/add_user_lot_commercial_kit.html",
            {"defined_kits": defined_kits},
        )


@login_required
def record_samples(request):
    """_summary_

    Parameters
    ----------
    request
        _description_

    Returns
    -------
        _description_
    """

    # Get Samples model fields with its verbose name, and dropdown lists for each of them
    fields_info = core.utils.samples.sample_table_fields(__package__)

    # Record new samples
    if request.method == "POST" and request.POST["action"] == "record_samples":
        # Get data from POST
        req_user = request.user.username
        excel_data = json.loads(request.POST["record_table_data"])
        header = json.loads(request.POST["record_table_header"])

        # Change excel header (verbose_name) to field names
        field_names = core.utils.common.sheet_header_to_field_name(
            header, fields_info["fields"]
        )
        # Convert excel list-list to dictionary with field_names
        excel_json_data = core.utils.common.jspreadsheet_to_dict(
            field_names, excel_data
        )

        # Test if json is empty and go back to table
        if len(excel_json_data) == 0:
            pre_def_samples = core.utils.samples.get_sample_objs_in_state("Pre-defined")
            return render(
                request,
                "wetlab/record_sample.html",
                {"fields_info": fields_info, "pre_def_samples": pre_def_samples},
            )
        else:
            # validate mandatory and redundant samples
            validation = core.utils.samples.validate_sample_data(
                excel_json_data, req_user, __package__
            )

            for val in validation:
                if not val["Validate"]:
                    return render(
                        request,
                        "wetlab/record_sample.html",
                        {
                            "fields_info": fields_info,
                            "validation": validation,
                            "excel_data": excel_data,
                        },
                    )

            # If all samples are validated
            try:
                # Record sample results.
                sample_record_result = core.utils.samples.save_recorded_samples(
                    excel_json_data, req_user, __package__
                )

                # If any samples couldn't be recorded for any error.
                for sample in sample_record_result:
                    if not sample["success"]:
                        # If some samples got an error. Some samples have been recorded and some of them don't.
                        return render(
                            request,
                            "wetlab/record_sample.html",
                            {
                                "fields_info": fields_info,
                                "error_message": "Some samples couldn't be recorded.",
                                "sample_record_result": sample_record_result,
                            },
                        )
            except Exception:
                # In case come uncatched error occurs
                error_message = (
                    "There was an unexpected error when recording the samples."
                )
                return render(
                    request,
                    "wetlab/record_sample.html",
                    {"error_message": error_message},
                )

        # If everything goes right, check if we need to add project data
        # Sort samples by project, and filter keeping only samples in "Pre-Defined" state.
        # Those only_recorded are already in complete state.
        sample_code_ids = [sample["sample_code_id"] for sample in sample_record_result]
        samples_query = core.models.Samples.objects.filter(
            sample_code_id__in=sample_code_ids
        )
        recorded_samples_info = wetlab.api.serializers.SampleSerializer(
            samples_query, many=True
        ).data
        project_ids = []
        for sample in recorded_samples_info:
            if sample["sample_project"] not in project_ids:
                project_ids.append(sample["sample_project"])

        # If no sample Pre-Defined just show result
        if not project_ids:
            return render(
                request,
                "wetlab/record_sample.html",
                {
                    "fields_info": fields_info,
                    "sample_record_result": sample_record_result,
                },
            )

        try:
            # if we have projects, get the fields for each projects associated with the recorded samples
            projects_fields = core.utils.samples.project_table_fields(project_ids)
            # Render the project data form
            return render(
                request,
                "wetlab/record_project_fields.html",
                {
                    "projects_fields": projects_fields,
                    "recorded_samples_info": recorded_samples_info,
                },
            )
        except Exception:
            # In case come uncatched error occurs
            error_message = (
                "There was an unexpected error when processing the project form."
            )
            return render(
                request,
                "wetlab/record_project_fields.html",
                {"error_message": error_message},
            )

    elif request.method == "POST" and request.POST["action"] == "select_samples_pre_defined":

        excel_data = json.loads(request.POST["predef_table_data"])
        header = json.loads(request.POST["predef_table_header"])

        # Change excel header (verbose_name) to field names
        field_names = core.utils.common.sheet_header_to_field_name(
            header, fields_info["fields"]
        )
        field_names.append("selected")
        # Convert excel list-list to dictionary with field_names
        excel_json_data = core.utils.common.jspreadsheet_to_dict(
            field_names, excel_data
        )

        selected_samples_info = [
            sample
            for sample in excel_json_data
            if sample["selected"]
        ]
        project_ids = []
        for sample in selected_samples_info:
            if sample["sample_project"] not in project_ids:
                project_ids.append(sample["sample_project"])

        try:
            # if we have projects, get the fields for each projects associated with the recorded samples
            projects_fields = core.utils.samples.project_table_fields(project_ids)
            # Render the project data form
            return render(
                request,
                "wetlab/record_project_fields.html",
                {
                    "projects_fields": projects_fields,
                    "recorded_samples_info": selected_samples_info,
                },
            )
        except Exception:
            # In case come uncatched error occurs
            error_message = (
                "There was an unexpected error when processing the project form."
            )
            return render(
                request,
                "wetlab/record_project_fields.html",
                {"error_message": error_message},
            )

    # Record project data
    elif request.method == "POST" and request.POST["action"] == "record_project_fields":
        not_validated_info = []
        not_saved_info = []
        projects_success = []
        json_data_all = []

        projects_fields = eval(request.POST["projects_fields"])

        for p_data in projects_fields:
            # Check if for any case there is no excel data sent for a project
            if not request.POST[p_data["sample_project_name"]]:
                # In case come uncatched error occurs
                error_message = f"Information for project {p_data['project'].sample_project_name} is missing."
                return render(
                    request,
                    "wetlab/record_project_fields.html",
                    {"error_message": error_message},
                )

            # Get field_names for the project
            field_names = [
                field["sample_project_field_name"]
                for field in p_data["sample_project_fields"]
            ]
            field_names.insert(0, "sample_name")
            field_names.insert(1, "sample_code_id")
            field_names.insert(2, "sample_project_name")
            field_names.insert(3, "only_recorded")

            # Get excel form data for this project from POST
            excel_data = json.loads(request.POST[f"{p_data['sample_project_name']}"])
            # Convert excel list-list to dictionary with field_names
            excel_json_data = core.utils.common.jspreadsheet_to_dict(
                field_names, excel_data
            )
            for json_data in excel_json_data:
                json_data_all.append(json_data)

            # validate types and option lists for projects
            validation = core.utils.samples.validate_project_data(
                excel_json_data, p_data["sample_project_name"]
            )
            # If some of the fields are not validated skip to next. DO NOT SAVE
            for val in validation:
                if not val["Validate"]:
                    not_validated_info.append(val)
                    next

            try:
                # save project data
                project_record_result = core.utils.samples.save_project_data(
                    excel_json_data, p_data
                )
                # Check if there was any error while saving
                if (
                    not project_record_result["success"]
                    and p_data["sample_project_name"] not in not_validated_info
                ):
                    not_saved_info.append(project_record_result)

            except Exception as e:
                # In case some uncatched error occurs
                error_message = (
                    "There was an unexpected error when recording the project data."
                    + str(e)
                )
                return render(
                    request,
                    "wetlab/record_sample.html",
                    {"fields_info": fields_info, "error_message": error_message},
                )

            if p_data["sample_project_name"] not in not_validated_info:
                projects_success.append(p_data["sample_project_name"])

        # Get recorded samples complete info
        sample_code_ids = [sample["sample_code_id"] for sample in json_data_all]
        samples_query = core.models.Samples.objects.filter(
            sample_code_id__in=sample_code_ids
        )
        recorded_samples_info = wetlab.api.serializers.SampleSerializer(
            samples_query, many=True, context={"output_label": "output_label"}
        ).data

        if not_validated_info:
            return render(
                request,
                "wetlab/record_project_fields.html",
                {
                    "projects_fields": projects_fields,
                    "projects_success": projects_success,
                    "not_validated_info": not_validated_info,
                    "json_data_all": json_data_all,
                    "recorded_samples_info": recorded_samples_info,
                },
            )
        elif not_saved_info:
            # If some project data got an error. Some project data have been recorded and some of them don't.
            return render(
                request,
                "wetlab/record_sample.html",
                {
                    "fields_info": fields_info,
                    "error_message": "Some project data couldn't be recorded.",
                },
            )
        # No validation nor saving errors. Show record summary.
        else:
            return render(
                request,
                "wetlab/record_project_fields.html",
                {
                    "recorded_samples_result": recorded_samples_info,
                },
            )
    # Record batch of samples
    elif request.method == "POST" and request.POST["action"] == "defineBatchSamples":
        req_user = request.user.username

        # Read excel file, remove empty rows and rename column names
        sample_batch_df = pd.read_excel(
            request.FILES["samplesExcel"], sheet_name=0, parse_dates=False
        )
        sample_batch_df = sample_batch_df.dropna(how="all")
        sample_batch_df = core.utils.load_batch.heading_refactor(sample_batch_df)

        # Test if column names are valid
        if core.utils.load_batch.validate_header(sample_batch_df):
            pre_def_samples = core.utils.samples.get_sample_objs_in_state("Pre-defined")
            return render(
                request,
                "wetlab/record_sample.html",
                {
                    "fields_info": fields_info,
                    "error_message": core.utils.load_batch.validate_header(
                        sample_batch_df
                    ),
                    "pre_def_samples": pre_def_samples,
                },
            )

        # Test if date columns have date format
        if core.utils.load_batch.check_format_date(sample_batch_df):
            pre_def_samples = core.utils.samples.get_sample_objs_in_state("Pre-defined")
            return render(
                request,
                "wetlab/record_sample.html",
                {
                    "fields_info": fields_info,
                    "error_message": core.utils.load_batch.check_format_date(
                        sample_batch_df
                    ),
                    "pre_def_samples": pre_def_samples,
                },
            )

        # Reformat date columns to correct format (mandatory befor converting to json)
        sample_batch_df = core.utils.load_batch.format_date(sample_batch_df)

        # Convert pandas dataframe to json list of dictionaries
        batch_json_data = core.utils.load_batch.read_batch_sample_file(sample_batch_df)

        # Test if json data is empty
        if len(batch_json_data) == 0:
            pre_def_samples = core.utils.samples.get_sample_objs_in_state("Pre-defined")
            return render(
                request,
                "wetlab/record_sample.html",
                {
                    "fields_info": fields_info,
                    "error_message": "Excell file is empty",
                    "pre_def_samples": pre_def_samples,
                },
            )
        else:
            # validate mandatory and redundant samples
            validation = core.utils.samples.validate_sample_data(
                batch_json_data, req_user, __package__
            )

            for val in validation:
                if not val["Validate"]:
                    pre_def_samples = core.utils.samples.get_sample_objs_in_state(
                        "Pre-defined"
                    )
                    return render(
                        request,
                        "wetlab/record_sample.html",
                        {
                            "fields_info": fields_info,
                            "validation": validation,
                            "pre_def_samples": pre_def_samples,
                        },
                    )

            # If all samples are validated
            try:
                # Record sample results.
                sample_record_result = core.utils.samples.save_recorded_samples(
                    batch_json_data, req_user, __package__
                )

                # If any samples couldn't be recorded for any error.
                for sample in sample_record_result:
                    if not sample["success"]:
                        # If some samples got an error. Some samples have been recorded and some of them don't.
                        return render(
                            request,
                            "wetlab/record_sample.html",
                            {
                                "fields_info": fields_info,
                                "error_message": "Some samples couldn't be recorded. Table summary:",
                                "sample_record_result": sample_record_result,
                            },
                        )
                    else:
                        return render(
                            request,
                            "wetlab/record_sample.html",
                            {
                                "fields_info": fields_info,
                                "sample_record_result": sample_record_result,
                            },
                        )

            except Exception:
                # In case come uncatched error occurs
                error_message = (
                    "There was an unexpected error when recording the samples."
                )
                return render(
                    request,
                    "wetlab/record_sample.html",
                    {"error_message": error_message},
                )

    # Form to get the new samples
    else:
        pre_def_samples = core.utils.samples.get_sample_objs_in_state("Pre-defined")
        return render(
            request,
            "wetlab/record_sample.html",
            {"fields_info": fields_info, "pre_def_samples": pre_def_samples},
        )


@login_required
def define_sample_projects_fields(request, sample_project_id):
    if not wetlab.utils.common.is_wetlab_manager(request):
        return render(
            request,
            "wetlab/define_sample_projects_fields.html",
            {"error_message": "You do not have enough privileges to see this page "},
        )
    # get the list of defined sample Projects
    if (
        request.method == "POST"
        and request.POST["action"] == "defineSampleProjectFields"
    ):
        sample_project_field_data = core.utils.samples.set_sample_project_fields(
            request.POST
        )
        return render(
            request,
            "wetlab/define_sample_project_fields.html",
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
                    "wetlab/define_sample_project_fields.html",
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
                    "wetlab/define_sample_project_fields.html",
                    {
                        "error_message": result["ERROR"],
                        "sample_project_data": sample_project_data,
                    },
                )
            return render(
                request,
                "wetlab/define_sample_project_fields.html",
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
                "wetlab/define_sample_project_fields.html",
                {
                    "error_message": "The requested Sample project does not exist",
                },
            )

        sample_project_data = core.utils.samples.define_table_for_sample_project_fields(
            sample_project_id
        )
        return render(
            request,
            "wetlab/define_sample_project_fields.html",
            {"sample_project_data": sample_project_data},
        )


@login_required
def modify_additional_kit(request, protocol_id):
    if not wetlab.utils.common.is_wetlab_manager(request):
        return render(
            request,
            "wetlab/modify_additional_kit.html",
            {"error_message": "You do not have enough privileges to see this page "},
        )

    if request.method == "POST" and request.POST["action"] == "modifyAdditionalKits":
        additional_kits_data_saved = (
            wetlab.utils.additional_kits.modify_fields_in_additional_kits(
                request.POST, request.user
            )
        )
        return render(
            request,
            "wetlab/modify_additional_kit.html",
            {"additional_kits_data_saved": additional_kits_data_saved},
        )
    else:
        if not core.utils.protocols.check_if_protocol_exists(protocol_id, __package__):
            return render(
                request,
                "wetlab/modify_additional_kit.html",
                {"error_message": "The requested additional kits do not exist"},
            )
        additional_kits_data = (
            wetlab.utils.additional_kits.get_additional_kits_data_to_modify(protocol_id)
        )
        return render(
            request,
            "wetlab/modify_additional_kit.html",
            {"additional_kits_data": additional_kits_data},
        )


@login_required
def modify_protocol_fields(request, protocol_id):
    if not wetlab.utils.common.is_wetlab_manager(request):
        return render(
            request,
            "wetlab/modify_protocol_fields.html",
            {
                "error_mes": [
                    "You do not have enough privileges to see this page ",
                    "Contact with your administrator .",
                ]
            },
        )
    if request.method == "POST" and request.POST["action"] == "modifyProtocolFields":
        protocol_field_saved = core.utils.protocols.modify_fields_in_protocol(
            request.POST
        )
        return render(
            request,
            "wetlab/modify_protocol_fields.html",
            {"protocol_field_saved": protocol_field_saved},
        )
    else:
        if not core.utils.protocols.check_if_protocol_exists(protocol_id, __package__):
            return render(
                request,
                "wetlab/modify_protocol_fields.html",
                {"error_message": "The requested Protocol does not exist"},
            )
        protocol_field = core.utils.protocols.get_protocol_fields(protocol_id)
        return render(
            request,
            "wetlab/modify_protocol_fields.html",
            {"protocol_field": protocol_field},
        )


@login_required
def modify_sample_project_fields(request, sample_project_id):
    if not wetlab.utils.common.is_wetlab_manager(request):
        return render(
            request,
            "wetlab/modify_sample_project_fields.html",
            {
                "content": [
                    "You do not have enough privileges to see this page ",
                    "Contact with your administrator .",
                ]
            },
        )
    if (
        request.method == "POST"
        and request.POST["action"] == "modifySampleProjectFields"
    ):
        sample_project_field_saved = core.utils.samples.modify_fields_in_sample_project(
            request.POST
        )
        return render(
            request,
            "wetlab/modify_sample_project_fields.html",
            {"sample_project_field_saved": sample_project_field_saved},
        )

    else:
        if not core.utils.samples.check_if_sample_project_id_exists(sample_project_id):
            return render(
                request,
                "wetlab/modify_sample_project_fields.html",
                {"error_message": "The requested Sample project does not exist"},
            )
        sample_project_field = core.utils.samples.get_parameters_sample_project(
            sample_project_id
        )
        return render(
            request,
            "wetlab/modify_sample_project_fields.html",
            {"sample_project_field": sample_project_field},
        )


@login_required
def define_molecule_uses(request):
    if not wetlab.utils.common.is_wetlab_manager(request):
        return render(
            request,
            "wetlab/define_molecule_uses.html",
            {
                "content": [
                    "You do not have enough privileges to see this page ",
                    "Contact with your administrator .",
                ]
            },
        )
    molecule_use_data = core.utils.samples.display_molecule_use(__package__)
    if request.method == "POST" and request.POST["action"] == "record_molecule_use":
        molecule_use_data.update(
            core.utils.samples.record_molecule_use(request.POST, __package__)
        )

    return render(
        request,
        "wetlab/define_molecule_uses.html",
        {"molecule_use_data": molecule_use_data},
    )


@login_required
def define_type_of_samples(request):
    if not wetlab.utils.common.is_wetlab_manager(request):
        return render(
            request,
            "wetlab/define_type_of_samples.html",
            {
                "content": [
                    "You do not have enough privileges to see this page ",
                    "Contact with your administrator .",
                ]
            },
        )
    sample_types = core.utils.samples.display_sample_types(__package__)
    if request.method == "POST" and request.POST["action"] == "addNewSampleType":
        sample_types.update(
            core.utils.samples.save_type_of_sample(request.POST, __package__)
        )

    return render(
        request,
        "wetlab/define_type_of_samples.html",
        {"sample_types": sample_types},
    )


@login_required
def display_sample(request, sample_id):
    sample_information = core.utils.samples.get_all_sample_information(sample_id, True)
    if "Error" not in sample_information:
        sample_information.update(
            core.utils.commercial_kits.get_molecule_lot_kit_in_sample(sample_id)
        )
        sample_information.update(
            wetlab.utils.library.get_all_library_information(sample_id)
        )
        sample_information.update(
            wetlab.utils.additional_kits.get_additional_kits_used_in_sample(sample_id)
        )
        sample_information.update(
            wetlab.utils.run.get_run_user_lot_kit_used_in_sample(sample_id)
        )
    else:
        sample_information = {}
    sample_obj = core.utils.samples.get_sample_obj_from_id(sample_id)
    if sample_obj:
        sample_name = sample_obj.get_sample_name()
        run_sample_obj = wetlab.utils.sample.get_sample_in_project_obj(sample_name)
        if run_sample_obj:
            sample_information.update(
                wetlab.utils.fetch_info.get_info_sample_in_run(run_sample_obj)
            )

    if len(sample_information) == 0:
        return render(
            request,
            "wetlab/display_sample.html",
            {"error_message": "Sample was not found"},
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
    sample_run_obj = wetlab.utils.sample.get_sample_in_project_obj(sample_run_id)
    if not sample_run_obj:
        return render(
            request,
            "wetlab/display_sample.html",
            {"error_message": "No Sample was found"},
        )
    # check user privileges
    if not wetlab.utils.common.is_wetlab_manager(request):
        # check if user is owner of the run or belongs to the shared user
        shared_user_ids = wetlab.utils.common.get_allowed_user_for_sharing(request.user)
        if sample_run_obj.get_user_id not in shared_user_ids:
            return render(
                request,
                "wetlab/display_sample.html",
                {"error_message": wetlab.config.ERROR_USER_NOT_WETLAB_MANAGER},
            )
    sample_information = wetlab.utils.fetch_info.get_info_sample_in_run(sample_run_obj)
    return render(
        request,
        "wetlab/display_sample.html",
        {"sample_information": sample_information},
    )


@login_required
def display_type_of_sample(request, sample_type_id):
    type_of_sample_data = core.utils.samples.get_type_of_sample_information(
        sample_type_id
    )
    if "ERROR" in type_of_sample_data:
        return render(
            request,
            "wetlab/display_type_of_sample.html",
            {"error_message": type_of_sample_data["ERROR"]},
        )
    return render(
        request,
        "wetlab/display_type_of_sample.html",
        {"type_of_sample_data": type_of_sample_data},
    )


@login_required
def handling_library_preparation(request):
    if wetlab.utils.common.is_wetlab_manager(request):
        samples_in_lib_prep = wetlab.utils.library.get_samples_for_library_preparation()
    else:
        samples_in_lib_prep = wetlab.utils.library.get_samples_for_library_preparation(
            request.user, True
        )

    if request.method == "POST" and request.POST["action"] == "assignProtocol":
        samples_in_lib_prep_protocol = (
            wetlab.utils.library.extract_protocol_library_preparation_form(request.POST)
        )
        if len(samples_in_lib_prep_protocol) == 0:
            return render(
                request,
                "wetlab/handling_library_preparation.html",
                {"stored_lib_prep": samples_in_lib_prep},
            )
        library_preparation_objs = (
            wetlab.utils.library.create_library_preparation_instance(
                samples_in_lib_prep_protocol, request.user
            )
        )
        lib_prep_protocol_parameters = (
            wetlab.utils.library.get_protocol_parameters_for_library_preparation(
                library_preparation_objs
            )
        )
        return render(
            request,
            "wetlab/handling_library_preparation.html",
            {"lib_prep_protocol_parameters": lib_prep_protocol_parameters},
        )

    # add protocol parameters for the user selected library preparation on defined state
    if request.method == "POST" and request.POST["action"] == "addProtocolParameter":
        lib_prep_ids = request.POST.getlist("libpreparation")
        if len(lib_prep_ids) == 0:
            error_message = "No selected library preparation was chosen"
            return render(
                request,
                "wetlab/handling_library_preparation.html",
                {
                    "samples_in_lib_prep": samples_in_lib_prep,
                    "error_message": error_message,
                },
            )
        # separate samples based on different protocols
        lib_prep_sep_protocol = wetlab.utils.library.separate_lib_prep_on_protocols(
            lib_prep_ids
        )

        lib_prep_protocol_parameters = (
            wetlab.utils.library.get_protocol_parameters_for_library_preparation(
                lib_prep_sep_protocol
            )
        )
        return render(
            request,
            "wetlab/handling_library_preparation.html",
            {"lib_prep_protocol_parameters": lib_prep_protocol_parameters},
        )

    # store the parameter library preparation protocol
    if request.method == "POST" and request.POST["action"] == "recordProtocolParamters":
        stored_params = wetlab.utils.library.analyze_and_store_prot_lib_param_values(
            request.POST
        )
        if "ERROR" in stored_params:
            error_message = stored_params["ERROR"]
            # TO DO in previoous form add lib_prep_ids
            lib_prep_ids = request.POST["lib_prep_ids"].split(",")

            library_preparation_objs = []
            for lib_prep_id in lib_prep_ids:
                library_preparation_objs.append(
                    wetlab.utils.library.get_lib_prep_obj_from_id(lib_prep_id)
                )
            lib_prep_protocol_parameters = (
                wetlab.utils.library.get_protocol_parameters_for_library_preparation(
                    library_preparation_objs
                )
            )
            # restore the user data
            lib_prep_protocol_parameters["data"] = json.loads(
                request.POST["protocol_data"]
            )
            return render(
                request,
                "wetlab/handling_library_preparation.html",
                {
                    "error_message": error_message,
                    "lib_prep_protocol_parameters": lib_prep_protocol_parameters,
                },
            )
        return render(
            request,
            "wetlab/handling_library_preparation.html",
            {"stored_params": stored_params},
        )

    if request.method == "POST" and request.POST["action"] == "importsamplesheet":
        data = {}
        (
            data["full_path_file"],
            data["file_name"],
        ) = wetlab.utils.samplesheet.store_user_input_file(request.FILES["uploadfile"])
        file_read = wetlab.utils.samplesheet.read_user_iem_file(data["full_path_file"])
        if not wetlab.utils.samplesheet.valid_user_iem_file(file_read):
            # Error found when extracting data from sample sheet
            data["ERROR"] = wetlab.config.ERROR_INVALID_FILE_FORMAT
            if not wetlab.utils.samplesheet.delete_stored_file(data["full_path_file"]):
                data["ERROR"].append(wetlab.config.ERROR_UNABLE_TO_DELETE_USER_FILE)
            return render(
                request,
                "wetlab/handling_library_preparation.html",
                {
                    "error_message": data["ERROR"],
                    "samples_in_lib_prep": samples_in_lib_prep,
                },
            )
        user_in_description = wetlab.utils.common.get_configuration_value(
            "DESCRIPTION_IN_SAMPLE_SHEET_MUST_HAVE_USERNAME"
        )
        if user_in_description == "TRUE":
            user_id_in_s_sheet = (
                wetlab.utils.library.extract_userids_from_sample_sheet_data(file_read)
            )
            if "ERROR" in user_id_in_s_sheet:
                if not wetlab.utils.samplesheet.delete_stored_file(
                    data["full_path_file"]
                ):
                    user_id_in_s_sheet["ERROR"].append(
                        wetlab.config.ERROR_UNABLE_TO_DELETE_USER_FILE
                    )
                return render(
                    request,
                    "wetlab/handling_library_preparation.html",
                    {
                        "error_message": user_id_in_s_sheet["ERROR"],
                        "samples_in_lib_prep": samples_in_lib_prep,
                    },
                )
        else:
            user_id_in_s_sheet = []
        sample_sheet_data = wetlab.utils.samplesheet.get_sample_sheet_data(file_read)

        valid_data = wetlab.utils.library.validate_sample_sheet_data(sample_sheet_data)
        if "ERROR" in valid_data:
            if not wetlab.utils.samplesheet.delete_stored_file(data["full_path_file"]):
                valid_data["ERROR"].append(
                    wetlab.config.ERROR_UNABLE_TO_DELETE_USER_FILE
                )
            return render(
                request,
                "wetlab/handling_library_preparation.html",
                {
                    "error_message": valid_data["ERROR"],
                    "samples_in_lib_prep": samples_in_lib_prep,
                },
            )

        platform = request.POST["platform"]
        configuration = request.POST[request.POST["platform"]]
        sample_sheet_data["file_name"] = data["file_name"]
        sample_sheet_data["userid_names"] = user_id_in_s_sheet
        lib_prep_sample_sheet_obj = (
            wetlab.utils.library.store_library_preparation_sample_sheet(
                sample_sheet_data, request.user, platform, configuration
            )
        )
        display_sample_sheet = (
            wetlab.utils.library.format_sample_sheet_to_display_in_form(
                sample_sheet_data
            )
        )
        display_sample_sheet[
            "lib_prep_user_sample_sheet"
        ] = lib_prep_sample_sheet_obj.get_user_sample_sheet_id()
        display_sample_sheet["platform"] = platform
        display_sample_sheet["iem_version"] = sample_sheet_data["iem_version"]
        if user_in_description == "TRUE":
            display_sample_sheet["user_list"] = wetlab.utils.common.get_userid_list()
        return render(
            request,
            "wetlab/handling_library_preparation.html",
            {"display_sample_sheet": display_sample_sheet},
        )

    if request.method == "POST" and request.POST["action"] == "storeIndexSample":
        store_data_result = (
            wetlab.utils.library.store_confirmation_library_preparation_index(
                request.POST
            )
        )
        if "ERROR" in store_data_result:
            return render(
                request,
                "wetlab/handling_library_preparation.html",
                {
                    "error_message": valid_data["ERROR"],
                    "samples_in_lib_prep": samples_in_lib_prep,
                },
            )
        stored_index = "True"
        return render(
            request,
            "wetlab/handling_library_preparation.html",
            {"stored_index": stored_index},
        )

    if request.method == "POST" and request.POST["action"] == "assignAdditionalKits":
        lib_prep_ids = request.POST.getlist("libpreparation")
        additional_kits = (
            wetlab.utils.additional_kits.get_additional_kits_from_lib_prep(lib_prep_ids)
        )
        return render(
            request,
            "wetlab/handling_library_preparation.html",
            {"additional_kits": additional_kits},
        )

    if request.method == "POST" and request.POST["action"] == "storeAdditionalKits":
        stored_additional_kits = (
            wetlab.utils.additional_kits.analyze_and_store_input_additional_kits(
                request.POST
            )
        )
        if "ERROR" in stored_additional_kits:
            error_message = stored_additional_kits["ERROR"]
            lib_prep_ids = request.POST["lib_prep_ids"].split(",")
            additional_kits = (
                wetlab.utils.additional_kits.get_additional_kits_from_lib_prep(
                    lib_prep_ids
                )
            )
            additional_kits["data"] = json.loads(request.POST["protocol_data"])
            return render(
                request,
                "wetlab/handling_library_preparation.html",
                {"error_message": error_message, "additional_kits": additional_kits},
            )

        return render(
            request,
            "wetlab/handling_library_preparation.html",
            {"stored_additional_kits": stored_additional_kits},
        )
    else:
        return render(
            request,
            "wetlab/handling_library_preparation.html",
            {"samples_in_lib_prep": samples_in_lib_prep},
        )


def handling_molecules(request):
    if request.method == "POST" and request.POST["action"] == "selectedMolecules":
        # If no samples are selected , call again this function to display again the sample list

        samples = core.utils.samples.get_selected_recorded_samples(
            request.POST["selected_samples"]
        )
        if len(samples) == 0:
            return redirect("handling_molecules")
        molecule_protocol = core.utils.samples.get_table_record_molecule(
            samples, __package__
        )
        if "ERROR" in molecule_protocol:
            return render(
                request,
                "wetlab/handling_molecules.html",
                {"error_message": "There was no valid sample selected "},
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
        molecule_recorded = core.utils.samples.record_molecules(
            request.POST, request.user, __package__
        )

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

        show_molecule_parameters = (
            core.utils.samples.display_molecule_protocol_parameters(
                molecule_recorded["molecule_ids"].split(","), request.user
            )
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
        show_molecule_parameters = (
            core.utils.samples.display_molecule_protocol_parameters(
                molecules, request.user
            )
        )
        return render(
            request,
            "wetlab/handling_molecules.html",
            {
                "molecule_recorded": molecule_recorded,
                "show_molecule_parameters": show_molecule_parameters,
            },
        )

    elif request.method == "POST" and request.POST["action"] == "addMoleculeParameters":
        molecule_parameters_updated = (
            core.utils.samples.add_molecule_protocol_parameters(request.POST)
        )
        if "pending" in request.POST:
            molecules = request.POST["pending"].split(",")
            show_molecule_parameters = (
                core.utils.samples.display_molecule_protocol_parameters(
                    molecules, request.user
                )
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
        sample_availables, molecules_availables, pending_to_use = "", "", ""
        if wetlab.utils.common.is_wetlab_manager(request):
            samples_list = core.utils.samples.get_sample_objs_in_state("Defined")
            samples_pending_use = core.utils.samples.get_sample_objs_in_state(
                "Pending for use"
            )
            molecule_list = core.utils.samples.get_molecule_objs_in_state("Defined")
        else:
            samples_list = core.utils.samples.get_sample_objs_in_state(
                "Defined", request.user, True
            )
            samples_pending_use = core.utils.samples.get_sample_objs_in_state(
                "Pending for use", request.user, True
            )
            molecule_list = core.utils.samples.get_molecule_objs_in_state(
                "Defined", request.user, True
            )
        if samples_list:
            sample_availables = core.utils.samples.create_table_to_select_molecules(
                samples_list
            )
        if molecule_list:
            molecules_availables = core.utils.samples.create_table_pending_molecules(
                molecule_list
            )
        if samples_pending_use:
            pending_to_use = core.utils.samples.create_table_molecule_pending_use(
                samples_pending_use, __package__
            )
        molecule_use_defined = core.utils.samples.check_if_molecule_use_defined(
            __package__
        )

        return render(
            request,
            "wetlab/handling_molecules.html",
            {
                "sample_availables": sample_availables,
                "molecules_availables": molecules_availables,
                "molecule_use_defined": molecule_use_defined,
                "pending_to_use": pending_to_use,
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
        result = wetlab.utils.sample.analyze_reprocess_data(
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
        detail_description[
            "information"
        ] = wetlab.config.SUCCESSFUL_REUSE_MOLECULE_EXTRACTION
        return render(
            request,
            "wetlab/successful_page.html",
            {"detail_description": detail_description},
        )
    # return to the main page because the page was not requested for the right page
    return redirect("")


@login_required
def repeat_molecule_extraction(request):
    if request.method == "POST" and request.POST["action"] == "repeat_extraction":
        sample_id = request.POST["sample_id"]
        if wetlab.utils.sample.analyze_reprocess_data(
            ["New Extraction"], sample_id, request.user
        ):
            molecule_protocol = core.utils.samples.get_table_record_molecule(
                [sample_id], __package__
            )
            molecule_protocol["samples"] = sample_id

            return render(
                request,
                "wetlab/handlingMolecules.html",
                {"molecule_protocol": molecule_protocol},
            )
    # return to the main page because the page was not requested for the right page
    return redirect("")


@login_required
def repeat_pool(request):
    if request.method == "POST" and request.POST["action"] == "repeat_pool":
        lib_prep_obj = wetlab.utils.library.get_lib_prep_obj_from_id(
            request.POST["lib_prep_id"]
        )
        lib_prep_code_id = lib_prep_obj.get_lib_prep_code()
        molecule_code_id = lib_prep_obj.get_molecule_code_id()
        sample_id = lib_prep_obj.get_sample_id()

        result = wetlab.utils.sample.analyze_reprocess_data(
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
    search_data = {}
    search_data["s_state"] = core.utils.samples.get_sample_states()

    if request.method == "POST" and request.POST["action"] == "searchsample":
        sample_name = request.POST["sample_name"]
        start_date = request.POST["startdate"]
        end_date = request.POST["enddate"]
        user_name = request.POST["username"]
        sample_state = request.POST["sampleState"]
        if request.POST["manager"] == "True":
            is_wetlab_manager = True
        else:
            is_wetlab_manager = False

        # check if all inputs are empty
        if (
            start_date == ""
            and end_date == ""
            and sample_name == ""
            and sample_state == ""
            and (
                (user_name != "" and not is_wetlab_manager)
                or (user_name == "" and is_wetlab_manager)
            )
        ):
            return render(
                request,
                "wetlab/search_sample.html",
                {"search_data": search_data},
            )
        if user_name != "" and len(user_name) < 5:
            return render(
                request,
                "wetlab/search_sample.html",
                {
                    "error_message": "The user name must contains at least 5 caracters ",
                    "search_data": search_data,
                },
            )
        # check the right format of start and end date
        if start_date != "" and not wetlab.utils.common.check_valid_date_format(
            start_date
        ):
            error_message = wetlab.config.ERROR_INVALID_FORMAT_FOR_DATES
            return render(
                request,
                "wetlab/search_sample.html",
                {"error_message": error_message, "search_data": search_data},
            )
        if end_date != "" and not wetlab.utils.common.check_valid_date_format(end_date):
            error_message = wetlab.config.ERROR_INVALID_FORMAT_FOR_DATES
            return render(
                request,
                "wetlab/search_sample.html",
                {
                    "error_message": error_message,
                    "search_data": search_data,
                },
            )
        # Get samples that are defined from collecting
        sample_list = core.utils.samples.search_samples(
            sample_name, user_name, sample_state, start_date, end_date
        )
        # Get samples that are defined in the run
        run_sample_list = wetlab.utils.sample.search_run_samples(
            sample_name, user_name, start_date, end_date, is_wetlab_manager
        )

        if len(sample_list) == 0 and len(run_sample_list) == 0:
            return render(
                request,
                "wetlab/search_sample.html",
                {
                    "error_message": wetlab.config.ERROR_NO_SAMPLE_FOUND,
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
                run_sample_obj = wetlab.utils.sample.get_sample_in_project_obj(
                    run_sample_list[0]
                )
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

        molecule_protocol = core.utils.samples.get_table_record_molecule(
            samples, __package__
        )
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
            molecule_protocol = core.utils.samples.get_table_record_molecule(
                samples, __package__
            )
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
                molecule_recorded.update(
                    core.utils.samples.get_table_record_molecule(samples)
                )

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
        show_molecule_parameters = (
            core.utils.samples.display_molecule_protocol_parameters(
                molecules, request.user
            )
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
            show_molecule_parameters = (
                core.utils.samples.display_molecule_protocol_parameters(
                    molecules, request.user
                )
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
    if not wetlab.utils.common.is_wetlab_manager(request):
        return render(
            request,
            "wetlab/create_pool.html",
            {"error_messge": "You do not have enough privileges to see this page "},
        )
    # collect the information for collecting
    display_list = wetlab.utils.pool.get_lib_prep_to_select_in_pool()

    if request.method == "POST" and request.POST["action"] == "createPool":
        new_pool = wetlab.utils.pool.define_new_pool(request.POST, request.user)

        if not isinstance(new_pool, wetlab.models.LibraryPool):
            display_list.update(new_pool)
            return render(
                request,
                "wetlab/create_pool.html",
                {"display_list": display_list},
            )
        information_for_created_pool = (
            wetlab.utils.pool.get_info_to_display_created_pool(new_pool)
        )
        return render(
            request,
            "wetlab/create_pool.html",
            {"information_for_created_pool": information_for_created_pool},
        )

    else:
        return render(
            request, "wetlab/create_pool.html", {"display_list": display_list}
        )


@login_required
def create_new_run(request):
    if not wetlab.utils.common.is_wetlab_manager(request):
        return render(
            request,
            "wetlab/create_new_run.html",
            {"error_message": "You do not have enough privileges to see this page "},
        )

    if request.method == "POST" and request.POST["action"] == "createNewRun":
        display_pools_for_run = wetlab.utils.run.display_available_pools()
        if "poolID" not in request.POST:
            error_message = wetlab.config.ERROR_NO_POOL_WAS_SELECTED_IN_FORM
            return render(
                request,
                "wetlab/create_new_run.html",
                {
                    "display_pools_for_run": display_pools_for_run,
                    "error_message": error_message,
                },
            )
        compatibility = wetlab.utils.run.check_valid_data_for_creation_run(
            request.POST, request.user
        )
        if "ERROR" in compatibility:
            return render(
                request,
                "wetlab/create_new_run.html",
                {
                    "display_pools_for_run": display_pools_for_run,
                    "error_message": compatibility["ERROR"],
                },
            )
        # compatible_index = check_index_compatible(lib_prep_ids)
        # check if run was created by crontab
        defined_run = wetlab.utils.run.check_run_already_defined_by_crontab(
            request.POST["experimentName"], request.POST["poolID"]
        )
        if "ERROR" in defined_run:
            return render(
                request,
                "wetlab/create_new_run.html",
                {
                    "display_pools_for_run": display_pools_for_run,
                    "error_message": defined_run["ERROR"],
                },
            )
        if not defined_run["defined"]:
            display_sample_information = wetlab.utils.run.create_run_in_pre_recorded_and_get_data_for_confirmation(
                request.POST, request.user
            )
            return render(
                request,
                "wetlab/create_new_run.html",
                {"display_sample_information": display_sample_information},
            )
        wetlab.utils.run.link_pool_with_existing_run(
            request.POST["experimentName"], request.POST["poolID"]
        )
        return render(
            request,
            "wetlab/create_new_run.html",
            {"info_message": wetlab.config.INFO_RUN_DEFINED_FROM_CRONTAB},
        )

    elif request.method == "POST" and request.POST["action"] == "continueWithRun":
        run_id = request.POST["run_ids"]
        experiment_name = wetlab.utils.run.get_experiment_name(run_id)
        pool_objs = wetlab.models.LibraryPool.objects.filter(
            run_process_id__exact=run_id
        )
        pool_ids = []
        for pool in pool_objs:
            pool_ids.append(pool.get_id())
        lib_prep_ids = wetlab.utils.run.get_library_prep_in_pools(pool_ids)

        display_sample_information = (
            wetlab.utils.run.get_library_preparation_data_in_run(lib_prep_ids, pool_ids)
        )
        display_sample_information.update(
            wetlab.utils.run.get_stored_user_sample_sheet(lib_prep_ids)
        )
        display_sample_information["experiment_name"] = experiment_name
        display_sample_information["run_process_id"] = run_id
        return render(
            request,
            "wetlab/create_new_run.html",
            {"display_sample_information": display_sample_information},
        )

    elif request.method == "POST" and request.POST["action"] == "storeDataNewRun":
        run_obj = wetlab.utils.run.get_run_obj_from_id(request.POST["run_process_id"])
        if run_obj.get_state() != "Pre-Recorded":
            exp_name = run_obj.get_run_name()
            error_message = str(exp_name + wetlab.config.ERROR_RUN_NAME_CREATED_ALREADY)
            display_pools_for_run = wetlab.utils.run.display_available_pools()
            return render(
                request,
                "wetlab/create_new_run.html",
                {
                    "display_pools_for_run": display_pools_for_run,
                    "error_message": error_message,
                },
            )
        run_data = wetlab.utils.run.collect_data_and_update_library_preparation_samples_for_run(
            request.POST, request.user
        )

        projects_objs = wetlab.utils.run.create_new_projects_added_to_run(
            run_data["projects"], run_data["run_obj"], request.user
        )
        if "ERROR" in projects_objs:
            display_pools_for_run = wetlab.utils.run.display_available_pools()
            return render(
                request,
                "wetlab/create_new_run.html",
                {
                    "display_pools_for_run": display_pools_for_run,
                    "error_message": projects_objs["ERROR"],
                },
            )

        run_obj.set_run_state("Recorded")

        wetlab.utils.run.store_confirmation_sample_sheet(run_data)
        # update the sample state for each one in the run
        pools_obj = wetlab.models.LibraryPool.objects.filter(run_process_id=run_obj)

        for pool_obj in pools_obj:
            pool_obj.set_pool_state("Used")
        wetlab.utils.library.update_batch_lib_prep_sample_state(
            run_data["lib_prep_ids"], "Sequencing"
        )
        created_new_run = {}
        created_new_run["exp_name"] = run_data["exp_name"]
        created_new_run["run_process_id"] = request.POST["run_process_id"]
        created_new_run["sample_sheet"] = run_obj.get_sample_file()
        return render(
            request,
            "wetlab/create_new_run.html",
            {"created_new_run": created_new_run},
        )
    else:
        display_pools_for_run = wetlab.utils.run.display_available_pools()
        return render(
            request,
            "wetlab/create_new_run.html",
            {"display_pools_for_run": display_pools_for_run},
        )


@login_required
def pending_sample_preparation(request):
    user_is_wetlab_manager = wetlab.utils.common.is_wetlab_manager(request)
    if user_is_wetlab_manager:
        req_user = None
        friend_list = False
    else:
        req_user = request.user
        friend_list = True
    pending_data = core.utils.samples.pending_sample_summary(req_user, friend_list)

    if len(pending_data["state"]) > 0:
        pending_data[
            "sample_heading"
        ] = wetlab.config.HEADING_FOR_PENDING_PROCESS_SAMPLES
        pending_data[
            "pending_sample_graphic"
        ] = wetlab.utils.statistics.get_pending_graphic_data(
            pending_data["state_number"],
            "Pending samples",
            "flint",
            "ex1",
            "500",
            "500",
            "chart-1",
        )
        # if wetlab manager create graphic for users on pending samples
        if user_is_wetlab_manager:
            pending_data[
                "pending_users_graphic"
            ] = wetlab.utils.statistics.get_pending_graphic_data(
                pending_data["users"],
                "Users with pending samples",
                "flint",
                "ex2",
                "500",
                "500",
                "chart-2",
            )

    return render(
        request,
        "wetlab/pending_sample_preparation.html",
        {"pending_data": pending_data},
    )


@login_required
def compare_samples(request):
    user_is_wetlab_manager = wetlab.utils.common.is_wetlab_manager(request)
    samples_data = wetlab.utils.sample.get_list_of_samples_in_projects(
        request.user, user_is_wetlab_manager
    )
    samples_data["user"] = request.user.username
    if request.method == "POST" and request.POST["action"] == "compareSamples":
        selected_sample_objs = wetlab.utils.sample.analyze_compare_samples_form(
            request.POST["table_data"]
        )
        if len(selected_sample_objs) == 0:
            error_message = wetlab.config.ERROR_NO_SAMPLES_SELECTED
            return render(
                request,
                "wetlab/compare_samples.html",
                {"error_message": error_message, "samples_data": samples_data},
            )
        compared_data = wetlab.utils.sample.get_comparation_sample_information(
            selected_sample_objs
        )

        return render(
            request,
            "wetlab/compare_samples.html",
            {"compared_data": compared_data},
        )
    else:
        return render(
            request,
            "wetlab/compare_samples.html",
            {"samples_data": samples_data},
        )


@login_required
def kit_inventory(request):
    if wetlab.utils.common.is_wetlab_manager(request):
        req_user = None
    else:
        req_user = request.user
    expired_kit = core.utils.commercial_kits.get_user_lot_kit_data(
        req_user, expired=True
    )
    valid_kit = core.utils.commercial_kits.get_user_lot_kit_data(request.user)
    if request.method == "POST" and request.POST["action"] == "runOutUserLotKit":
        selected_user_kits = request.POST.getlist("userKit")
        if len(selected_user_kits) == 0:
            return render(
                request,
                "wetlab/kit_inventory.html",
                {
                    "expired_kit": expired_kit,
                    "valid_kit": valid_kit,
                    "user_name": request.user.username,
                },
            )
        run_out_kits = core.utils.commercial_kits.set_user_lot_kit_to_run_out(
            selected_user_kits
        )
        return render(
            request,
            "wetlab/kit_inventory.html",
            {"run_out_kits": run_out_kits},
        )

    else:
        return render(
            request,
            "wetlab/kit_inventory.html",
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
                "wetlab/search_user_lot_kit.html",
                {
                    "protocol_list": protocol_list,
                    "platform_list": platform_list,
                },
            )

        if request.POST[
            "expired"
        ] != "" and not wetlab.utils.common.check_valid_date_format(
            request.POST["expired"]
        ):
            error_message = wetlab.config.ERROR_INVALID_FORMAT_FOR_DATES
            return render(
                request,
                "wetlab/search_user_lot_kit.html",
                {
                    "protocol_list": protocol_list,
                    "platform_list": platform_list,
                    "ERROR": error_message,
                },
            )

        user_kits_objs = core.utils.commercial_kits.search_user_lot_kit_from_user_form(
            request.POST
        )
        if user_kits_objs == "No defined":
            error_message = wetlab.config.ERROR_NO_USER_LOT_KIT_DEFINED
            return render(
                request,
                "wetlab/search_user_lot_kit.html",
                {
                    "protocol_list": protocol_list,
                    "platform_list": platform_list,
                    "ERROR": error_message,
                },
            )
        if len(user_kits_objs) > 1:
            display_user_kit_list = core.utils.commercial_kits.display_user_lot_kit_information_from_query_list(
                user_kits_objs
            )
            return render(
                request,
                "wetlab/search_user_lot_kit.html",
                {"display_user_kit_list": display_user_kit_list},
            )
        elif len(user_kits_objs) == 0:
            error_message = wetlab.config.ERROR_NO_MATCHES_FOR_USER_LOT_KIT
            return render(
                request,
                "wetlab/search_user_lot_kit.html",
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
            "wetlab/search_user_lot_kit.html",
            {"protocol_list": protocol_list, "platform_list": platform_list},
        )


@login_required
def display_user_lot_kit(request, user_kit_id):
    user_kit_obj = core.utils.commercial_kits.get_user_lot_commercial_kit_obj_from_id(
        user_kit_id
    )
    if user_kit_obj is None:
        return render(
            request,
            "wetlab/displayUserLotKit.html",
            {"error_message": ["Invalid User Lot Commercial Kit"]},
        )
    user_lot_kit_data = core.utils.commercial_kits.get_user_lot_kit_data_to_display(
        user_kit_obj
    )
    return render(
        request,
        "wetlab/displayUserLotKit.html",
        {"user_lot_kit_data": user_lot_kit_data},
    )


@login_required
def sequencer_configuration(request):
    if not wetlab.utils.common.is_wetlab_manager(request):
        return render(
            request,
            "wetlab/sequencer_configuration.html",
            {"error_message": wetlab.config.ERROR_USER_NOT_WETLAB_MANAGER},
        )
    sequencer_info = wetlab.utils.sequencers.get_list_sequencer_configuration()
    sequencer_info["platforms"] = wetlab.utils.sequencers.get_platform_data()
    sequencer_info["sequencer_names"] = wetlab.utils.sequencers.get_defined_sequencers()

    if request.method == "POST" and request.POST["action"] == "addNewSequencer":
        new_sequencer = wetlab.utils.sequencers.define_new_sequencer(request.POST)
        if "ERROR" in new_sequencer:
            return render(
                request,
                "wetlab/sequencer_configuration.html",
                {"sequencer_info": sequencer_info, "error_message": new_sequencer},
            )
        return render(
            request,
            "wetlab/sequencer_configuration.html",
            {"sequencer_info": sequencer_info, "new_defined_sequencer": new_sequencer},
        )
    if request.method == "POST" and request.POST["action"] == "addNewConfiguration":
        new_defined_configuration = (
            wetlab.utils.sequencers.define_new_seq_configuration(request.POST)
        )
        if "ERROR" in new_defined_configuration:
            return render(
                request,
                "wetlab/sequencer_configuration.html",
                {"sequencer_info": sequencer_info, "ERROR": new_defined_configuration},
            )
        return render(
            request,
            "wetlab/sequencer_configuration.html",
            {
                "sequencer_info": sequencer_info,
                "new_defined_configuration": new_defined_configuration,
            },
        )
    else:
        return render(
            request,
            "wetlab/sequencer_configuration.html",
            {"sequencer_info": sequencer_info},
        )


@login_required
def sequencer_details(request, seq_id):
    seq_obj = wetlab.utils.sequencers.get_sequencer_obj_from_id(seq_id)
    if seq_obj is None:
        return render(
            request,
            "wetlab/sequencer_inventory.html",
            {"error_message": "Sequencer not found"},
        )
    if request.method == "POST" and request.POST["action"] == "updateSequencer":
        seq_data = {}
        end_date = request.POST["seqEnd"]
        start_date = request.POST["seqStart"]
        if end_date == "" or end_date == "None":
            end_date = None
        else:
            if not wetlab.utils.common.check_valid_date_format(end_date):
                return render(
                    request,
                    "wetlab/sequencer_inventory.html",
                    {"error_message": "Format date is incorrect"},
                )
        if start_date == "" or start_date == "None":
            start_date = None
        else:
            if not wetlab.utils.common.check_valid_date_format(start_date):
                return render(
                    request,
                    "wetlab/sequencer_inventory.html",
                    {"error_message": "Format date is incorrect"},
                )

        seq_data["serial"] = request.POST["seqSerial"]
        seq_data["state"] = request.POST["seqState"]
        seq_data["lanes"] = request.POST["seqLanes"]
        seq_data["location"] = request.POST["seqLoc"]
        seq_data["description"] = request.POST["description"]
        seq_data["end_date"] = end_date
        seq_data["start_date"] = start_date

        seq_obj.update_sequencer_data(seq_data)

        return render(
            request,
            "wetlab/sequencer_inventory.html",
            {"updated_sequencer": request.POST["seq_name"]},
        )

    return render(
        request,
        "wetlab/sequencer_inventory.html",
        {"sequencer_detail_data": seq_obj.get_all_sequencer_data("%Y-%m-%d")},
    )


@login_required
def sequencer_inventory(request):
    sequencer_data = wetlab.utils.sequencers.get_sequencer_inventory_data("%d/%m/%Y")
    return render(
        request,
        "wetlab/sequencer_inventory.html",
        {"sequencer_data": sequencer_data},
    )
