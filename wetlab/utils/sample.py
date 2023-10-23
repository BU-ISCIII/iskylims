# Generic imports
import json

# django imports
from django.db.models import F, Func, Value, CharField

# Local imports
import core.fusioncharts.fusioncharts
import core.utils.common
import core.utils.samples
import wetlab.config
import wetlab.models
import wetlab.utils.stats_graphs


def analyze_reprocess_data(reprocess_data, reprocess_sample_id):
    """
    Description:
        The function will get the option of reprocessing sample and it updates the sample state for reprocessing.
        In case that a new library preparation was required a new object is created.
    Input:
        reprocess_data           # data for creating the reuse
        reprocess_sample_id      # sample id to reprocess
        req_user                 # register user
    Functions:
        update_sample_reused
        update_molecule_reused
        get_sample_obj_from_id
    Return:
        True if user request were right, or Invalid options if user requests were wrong.
    """
    # options = json_data[-1]
    if "New Extraction" in reprocess_data:
        sample_obj = core.utils.samples.update_sample_reused(reprocess_sample_id)
        sample_obj.set_state("Defined")
        return True

    elif "New Library Preparation" in reprocess_data:
        molecule_code_id = reprocess_data[0]
        if molecule_code_id == "":
            return "Invalid options"
        else:
            sample_obj = core.utils.samples.get_sample_obj_from_id(reprocess_sample_id)
            if "Library preparation" != sample_obj.get_sample_state():
                molecule_obj = core.utils.samples.update_molecule_reused(
                    reprocess_sample_id, molecule_code_id
                )
                if not molecule_obj:
                    return "Invalid options"
                #
                sample_obj = core.utils.samples.update_sample_reused(
                    reprocess_sample_id
                )
                sample_obj.set_state("Library preparation")

            return True
    elif "New Pool" in reprocess_data:
        molecule_code_id = reprocess_data[0]
        lib_prep_code_id = reprocess_data[1]

        if not wetlab.models.LibPrepare.objects.filter(
            sample_id__pk__exact=reprocess_sample_id,
            lib_prep_code_id__exact=lib_prep_code_id,
        ).exists():
            return "Invalid options"
        lib_prep_obj = wetlab.models.LibPrepare.objects.get(
            sample_id__pk__exact=reprocess_sample_id,
            lib_prep_code_id__exact=lib_prep_code_id,
        )
        sample_obj = core.utils.samples.update_sample_reused(reprocess_sample_id)
        molecule_obj = core.utils.samples.update_molecule_reused(
            reprocess_sample_id, molecule_code_id
        )
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
        sample_objs.append(get_sample_in_project_obj(row_data[-2]))
    return sample_objs


def get_comparation_sample_information(sample_objs):
    """
    Description:
        The function get sample information to build the comparation data.
    Input:
        sample_objs       # samples objs
    Functions:
        get_sample_obj_from_id
    Return:
        compared_data
    """
    compared_data = {}
    compared_data["table_data"] = []
    for sample_obj in sample_objs:
        run_obj = sample_obj.get_run_obj()
        stats_fl_obj = wetlab.models.StatsFlSummary.objects.filter(
            runprocess_id=run_obj
        ).last()
        data = sample_obj.get_sample_information()
        data.insert(2, sample_obj.get_run_name())
        data.insert(3, stats_fl_obj.get_sample_number())
        compared_data["table_data"].append(data)
    compared_data[
        "table_heading"
    ] = wetlab.config.HEADING_COMPARATION_SAMPLE_INFORMATION
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
        get_friend_list     # Located at core.utils.common
    Return:
        samples_data.
    """
    samples_data = {}
    sample_objs = ""
    if wetlab_manager:
        if wetlab.models.SamplesInProject.objects.all().exists():
            sample_objs = (
                wetlab.models.SamplesInProject.objects.all()
                .order_by("generated_at")
                .reverse()
            )
    else:
        user_list_ids = core.utils.common.get_friend_list(user)
        if wetlab.models.SamplesInProject.objects.filter(
            user_id_id__in=user_list_ids
        ).exists():
            sample_objs = (
                wetlab.models.SamplesInProject.objects.filter(
                    user_id_id__in=user_list_ids
                )
                .order_by("generated_at")
                .reverse()
            )
    if sample_objs != "":
        samples_data["data"] = list(
            sample_objs.values_list(
                "run_process_id__run_name", "user_id__username", "sample_name"
            )
            .annotate(
                formated_date=Func(
                    F("run_process_id__run_finish_date"),
                    Value("%Y-%m-%d"),
                    function="DATE_FORMAT",
                    output_field=CharField(),
                )
            )
            .annotate(id=F("pk"))
        )
        samples_data["heading"] = wetlab.config.HEADING_COMPARATION_SAMPLE_LIST
    return samples_data


def get_sample_in_project_obj(sample_value, value_is_pk=False):
    """Return the sample instance for the request. Two input types are accepted
        sample name or the PK of the sample

    Args:
        sample_value (integer/string): If integer value is fetch then it is
            considered as the PK if variable contains a string it is considered
            as sample name.

    Returns:
        SamplesInProject: instance of the input value
    """
    sample_in_project_obj = ""
    if value_is_pk:
        try:
            _ = int(sample_value)
            if wetlab.models.SamplesInProject.objects.filter(
                pk__exact=sample_value
            ).exists():
                sample_in_project_obj = wetlab.models.SamplesInProject.objects.get(
                    pk__exact=sample_value
                )
        except ValueError:
            pass
    else:
        if wetlab.models.SamplesInProject.objects.filter(
            sample_name__exact=sample_value
        ).exists():
            sample_in_project_obj = wetlab.models.SamplesInProject.objects.filter(
                sample_name__exact=sample_value
            ).last()
    return sample_in_project_obj


def search_run_samples(sample_name, user_name, start_date, end_date, is_wetlab_manager):
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
        get_friend_list # located at core/utils/common.py
    Return:
        run_sample_obj.
    """
    run_sample_list = []

    if not wetlab.models.SamplesInProject.objects.all().exists():
        return run_sample_list
    if is_wetlab_manager:
        if user_name == "":
            run_sample_founds = wetlab.models.SamplesInProject.objects.all()
        else:
            run_sample_founds = wetlab.models.SamplesInProject.objects.filter(
                user_id__username__iexact=user_name
            )
    else:
        if wetlab.models.User.objects.filter(username__exact=user_name).exists():
            user_name_obj = wetlab.models.User.objects.filter(
                username__exact=user_name
            ).last()
            user_friend_list = core.utils.common.get_friend_list(user_name_obj)
            if not wetlab.models.SamplesInProject.objects.filter(
                user_id__in=user_friend_list
            ).exists():
                return run_sample_list
            else:
                run_sample_founds = wetlab.models.SamplesInProject.objects.filter(
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
    run_sample_list = list(
        run_sample_founds.values_list(
            "pk", "sample_name", "project_id__project_name", "run_process_id__run_name"
        ).annotate(
            formated_date=Func(
                F("run_process_id__run_finish_date"),
                Value("%Y-%m-%d"),
                function="DATE_FORMAT",
                output_field=CharField(),
            )
        )
    )
    return run_sample_list
