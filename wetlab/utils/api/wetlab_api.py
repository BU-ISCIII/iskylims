from django.db import connection

import wetlab.models


def get_runs_projects_samples_and_dates(user_list_ids):
    """
    Description:
        The function api return a dictionnary having as keys the run and the project and
        value a list of tupla sample_id, sample_name
    Input:
        user_list_ids     # user id list to get the list of samples
    Return:
        samples_data
    """
    samples_data = []
    if wetlab.models.SamplesInProject.objects.filter(user_id_id__in=user_list_ids).exists():
        cursor = connection.cursor()
        user_list_str = [str(id_str) for id_str in user_list_ids]
        user_list = "( " + ",".join(user_list_str) + " )"
        fetch_fields = " wetlab_run_process.run_name, wetlab_run_process.id ,wetlab_projects.project_name , wetlab_projects.id ,  wetlab_samples_in_project.sample_name, wetlab_samples_in_project.id, DATE_FORMAT(wetlab_run_process.run_completed_date,'%d/%m/%Y') AS niceDate , wetlab_running_parameters.run_id "
        q_tables = " FROM wetlab_samples_in_project inner join wetlab_projects ON wetlab_samples_in_project.project_id_id = wetlab_projects.id inner join wetlab_run_process ON wetlab_samples_in_project.run_process_id_id = wetlab_run_process.id   inner join wetlab_running_parameters ON wetlab_run_process.id =  wetlab_running_parameters.run_name_id_id "
        restrict_results = (
            "WHERE wetlab_samples_in_project.user_id_id IN " + user_list
        )
        query = (
            "SELECT "
            + fetch_fields
            + q_tables
            + restrict_results
            + "ORDER BY (wetlab_run_process.run_completed_date) DESC"
        )

        cursor.execute(query)
        samples_data = cursor.fetchall()

    return samples_data
