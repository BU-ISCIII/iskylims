from iSkyLIMS_core.utils.common import get_friend_list
from iSkyLIMS_core.utils.handling_samples import *
from iSkyLIMS_wetlab.models import *
from iSkyLIMS_wetlab.wetlab_config import *

from ..fusioncharts.fusioncharts import FusionCharts
from .stats_graphics import *


def get_codeID_for_resequencing(sample_recorded):
    """
    Description:
        The function will get the already defined molecule and library preparation objects that
        are related to the sample to allow user to choose the difference alternatives for
        resequencing the sample
        Return a dictionary with all available possibilities
    Input:
        sample_recorded : sample id
    Functions:
        get_sample_obj_from_id : located at iSkyLIMS_core/utils/handling_samples.py
        get_molecule_objs_from_sample : located at iSkyLIMS_core/utils/handling_samples.py
        get_molecule_codeid_from_object : located at iSkyLIMS_core/utils/handling_samples.py
    Variables:
        lib_prep_available # list all possibilities for library preparation
        mol_lib_prep_available # list all possibilities for molecule
    Return:
        sample_recorded.
    """

    mol_lib_prep_available = {}
    lib_prep_available = ["New Library Preparation"]
    mol_lib_prep_available["New Extraction"] = [""]
    sample_obj = get_sample_obj_from_id(sample_recorded["sample_id_for_action"])
    molecule_objs = get_molecule_objs_from_sample(
        sample_recorded["sample_id_for_action"]
    )

    for molecule_obj in molecule_objs:
        molecule_id = get_molecule_codeid_from_object(molecule_obj)
        mol_lib_prep_available[molecule_id] = ["New Library Preparation"]
        if LibPrepare.objects.filter(
            molecule_id=molecule_obj, sample_id=sample_obj
        ).exists():
            libs_prep_obj = LibPrepare.objects.filter(
                molecule_id=molecule_obj, sample_id=sample_obj
            )
            for lib_prep_obj in libs_prep_obj:
                lib_prep_available.append(lib_prep_obj.get_lib_prep_code())
                mol_lib_prep_available[molecule_id].append(
                    lib_prep_obj.get_lib_prep_code()
                )

    sample_recorded["rep_filter_selection"] = []
    for key, value in mol_lib_prep_available.items():
        sample_recorded["rep_filter_selection"].append([key, value])
    sample_recorded["molecule_available"] = list(mol_lib_prep_available.keys())
    sample_recorded["lib_prep_available"] = lib_prep_available
    return sample_recorded


def analyze_reprocess_data(reprocess_data, reprocess_sample_id, reg_user):
    """
    Description:
        The function will get the option of reprocessing sample and it updates the sample state for reprocessing.
        In case that a new library preparation was required a new object is created.
    Input:
        reprocess_data           # data for creating the reuse
        reprocess_sample_id      # sample id to reprocess
        reg_user                 # register user
    Functions:
        update_sample_reused    # located at iSkyLIMS_core/utils/handling_samples.py
        update_molecule_reused  # located at iSkyLIMS_core/utils/handling_samples.
        get_sample_obj_from_id  # located at iSkyLIMS_core/utils/handling_samples.
    Return:
        True if user request were right, or Invalid options if user requests were wrong.
    """
    # options = json_data[-1]
    if "New Extraction" in reprocess_data:
        sample_obj = update_sample_reused(reprocess_sample_id)
        sample_obj.set_state("Defined")
        return True

    elif "New Library Preparation" in reprocess_data:
        molecule_code_id = reprocess_data[0]
        if molecule_code_id == "":
            return "Invalid options"
        else:
            sample_obj = get_sample_obj_from_id(reprocess_sample_id)
            if "Library preparation" != sample_obj.get_sample_state():
                molecule_obj = update_molecule_reused(
                    reprocess_sample_id, molecule_code_id
                )
                if not molecule_obj:
                    return "Invalid options"
                #
                sample_obj = update_sample_reused(reprocess_sample_id)
                sample_obj.set_state("Library preparation")

            return True
    elif "New Pool" in reprocess_data:
        molecule_code_id = reprocess_data[0]
        lib_prep_code_id = reprocess_data[1]

        if not LibPrepare.objects.filter(
            sample_id__pk__exact=reprocess_sample_id,
            lib_prep_code_id__exact=lib_prep_code_id,
        ).exists():
            return "Invalid options"
        lib_prep_obj = LibPrepare.objects.get(
            sample_id__pk__exact=reprocess_sample_id,
            lib_prep_code_id__exact=lib_prep_code_id,
        )
        sample_obj = update_sample_reused(reprocess_sample_id)
        molecule_obj = update_molecule_reused(reprocess_sample_id, molecule_code_id)
        lib_prep_obj.set_state("Reused pool")
        lib_prep_obj.set_increase_reuse()
        sample_obj.set_state("Pool preparation")

    else:
        return "Invalid options"


def analyze_compare_samples_form(form_data):
    """
    Description:
        The function get the selected samples and return the sample objs.
    Input:
        form_data               # data collected from the user
    Functions:
        get_sample_obj_from_id  # located at this file.
    Return:
        sample_objs
    """
    sample_objs = []
    cs_json_data = json.loads(form_data)
    for row_data in cs_json_data:
        if row_data[-1] is False:
            continue
        sample_objs.append(get_sample_in_project_obj_from_id(row_data[-2]))
    return sample_objs


def get_comparation_sample_information(sample_objs):
    """
    Description:
        The function get sample information to build the comparation data.
    Input:
        sample_objs       # samples objs
    Functions:
        get_sample_obj_from_id : located at iSkyLIMS_core/utils/handling_samples.
    Return:
        compared_data
    """
    compared_data = {}
    compared_data["table_data"] = []
    for sample_obj in sample_objs:
        run_obj = sample_obj.get_run_obj()
        stats_fl_obj = StatsFlSummary.objects.filter(runprocess_id=run_obj).last()
        data = sample_obj.get_sample_information()
        data.insert(2, sample_obj.get_run_name())
        data.insert(3, stats_fl_obj.get_sample_number())
        compared_data["table_data"].append(data)
    compared_data["table_heading"] = HEADING_COMPARATION_SAMPLE_INFORMATION
    return compared_data


def get_list_of_samples_in_projects(user, wetlab_manager):
    """
    Description:
        The function gets the list of the sampleInProject an returns sample_name, project, run and sample id
    Input:
        user                # user to filter the samples
        wetlab_manager      # boolean if user is wetlab_manager or not
    Constants:
        HEADING_COMPARATION_SAMPLE_LIST
    Functions:
        get_friend_list     # Located at iSkyLIMS_core.utils.common
    Return:
        samples_data.
    """
    samples_data = {}
    sample_objs = ""
    if wetlab_manager:
        if SamplesInProject.objects.all().exists():
            sample_objs = (
                SamplesInProject.objects.all().order_by("generated_at").reverse()
            )
    else:
        user_list_ids = get_friend_list(user)
        if SamplesInProject.objects.filter(user_id_id__in=user_list_ids).exists():
            sample_objs = (
                SamplesInProject.objects.filter(user_id_id__in=user_list_ids)
                .order_by("generated_at")
                .reverse()
            )
    if sample_objs != "":
        samples_data["data"] = []
        for sample_obj in sample_objs:
            run_obj = sample_obj.get_run_obj()
            if run_obj.get_state() != "Completed":
                continue
            data = []
            data.append(sample_obj.get_sample_name())
            data.append(sample_obj.get_project_name())
            data.append(sample_obj.get_run_name())
            data.append(run_obj.get_run_finish_date())
            data.append(sample_obj.get_sample_id())
            samples_data["data"].append(data)
        samples_data["heading"] = HEADING_COMPARATION_SAMPLE_LIST
    return samples_data


def get_sample_in_project_obj_from_id(sample_in_project_id):
    """
    Description:
        The function gets the sampleInProject id and return the object
        Return the if of sampleInProject
    Input:
        sample_name     # sample name to look at
    Return:
        sample_in_project_obj.
    """
    sample_in_project_obj = ""
    if SamplesInProject.objects.filter(pk__exact=sample_in_project_id).exists():
        sample_in_project_obj = SamplesInProject.objects.get(
            pk__exact=sample_in_project_id
        )

    return sample_in_project_obj


def get_sample_in_project_obj_from_sample_name(sample_name_in_project):
    """
    Description:
        The function gets the sampleInProject id and return the object
        Return the if of sampleInProject
    Input:
        sample_name     # sample name to look at
    Return:
        sample_in_project_obj.
    """
    sample_in_project_obj = ""
    if SamplesInProject.objects.filter(
        sample_name__exact=sample_name_in_project
    ).exists():
        sample_in_project_obj = SamplesInProject.objects.filter(
            sample_name__exact=sample_name_in_project
        ).last()

    return sample_in_project_obj


def get_run_sample_id(sample_name):
    """
    Description:
        The function gets the sample name and if found it returns the SamplesInProject objects.
        Return the if of sampleInProject
    Input:
        sample_name     # sample name to look at
    Return:
        run_sample_obj.
    """
    run_sample_obj = ""
    if SamplesInProject.objects.filter(sampleName__exact=sample_name).exists():
        sample_run_objs = SamplesInProject.objects.filter(
            sample_name__exact=sample_name
        )
        if len(sample_run_objs) > 1:
            pass
        else:
            run_sample_obj = SamplesInProject.objects.get(
                sample_name__exact=sample_name
            )
    return run_sample_obj


def search_run_samples(sample_name, user_name, start_date, end_date):
    """
    Description:
        The function search the run samples that matchs with the requested conditions.
        Return the if of sampleInProject
    Input:
        sample_name     # sample name to look at
        user_name       # user name
        start_date      # date from starting the search
        end_date        # date from ending the search
    Functions:
        get_friend_list # located at iSkyLIMS_core/utils/common.py
    Return:
        run_sample_obj.
    """
    run_sample_list = []

    if SamplesInProject.objects.all().exists():
        run_sample_founds = SamplesInProject.objects.all()
    else:
        return run_sample_list
    if user_name != "":
        if User.objects.filter(username__exact=user_name).exists():
            user_name_obj = User.objects.filter(username__exact=user_name).last()
            user_friend_list = get_friend_list(user_name_obj)
            if not run_sample_founds.filter(user_id__in=user_friend_list).exists():
                return run_sample_list
            else:
                run_sample_founds = run_sample_founds.filter(
                    user_id__in=user_friend_list
                )
        else:
            return run_sample_list
    if sample_name != "":
        if run_sample_founds.filter(sample_name__exact=sample_name).exists():
            run_sample_founds = run_sample_founds.filter(sample_name__exact=sample_name)
            if len(run_sample_founds) == 1:
                run_sample_list.append(run_sample_founds[0].pk)
                return run_sample_list

        elif run_sample_founds.filter(sample_name__icontains=sample_name).exists():
            run_sample_founds = run_sample_founds.filter(
                sample_name__icontains=sample_name
            )
        else:
            return run_sample_list
    if start_date != "" and end_date != "":
        run_sample_founds = run_sample_founds.filter(
            generated_at__range=(start_date, end_date)
        )

    if start_date != "" and end_date == "":
        run_sample_founds = run_sample_founds.filter(generated_at__gte=start_date)

    if start_date == "" and end_date != "":
        run_sample_founds = run_sample_founds.filter(generated_at__lte=end_date)

    if len(run_sample_founds) == 1:
        run_sample_list.append(run_sample_founds[0].pk)
        return run_sample_list

    for run_sample in run_sample_founds:
        run_sample_list.append(run_sample.get_basic_info())

    return run_sample_list


def pending_samples_for_grafic(pending):
    """
    Description:
        The function create the graphic for pending samples.
    Input:
        pending     # dictionary with the pending sample information
    Functions:
        graphic_3D_pie # located at iSkyLIMS_wetlab/utils/stats_graphics.py
    Return:
        graphic_pending_samples.
    """
    number_of_pending = {}
    number_of_pending["DEFINED"] = pending["defined"]["length"]
    number_of_pending["EXTRACTED MOLECULE"] = pending["extract_molecule"]["length"]
    if "create_library_preparation" in pending:
        number_of_pending["LIBRARY PREPARATION"] = pending[
            "create_library_preparation"
        ]["length"]

    data_source = graphic_3D_pie(
        "Number of Pending Samples", "", "", "", "fint", number_of_pending
    )
    graphic_pending_samples = FusionCharts(
        "pie3d", "ex1", "430", "450", "chart-1", "json", data_source
    )
    return graphic_pending_samples
