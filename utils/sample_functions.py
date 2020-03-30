from iSkyLIMS_wetlab.models import *
from iSkyLIMS_core.utils.handling_samples import *
from iSkyLIMS_core.utils.generic_functions import get_friend_list



def get_codeID_for_resequencing(sample_recorded):
    '''
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
    '''

    mol_lib_prep_available = {}
    lib_prep_available = ['New Library Preparation']
    mol_lib_prep_available['New Extraction'] =['']
    sample_obj = get_sample_obj_from_id (sample_recorded['sample_id_for_action'])
    molecule_objs = get_molecule_objs_from_sample(sample_recorded['sample_id_for_action'])

    for molecule_obj in molecule_objs:
        molecule_id = get_molecule_codeid_from_object (molecule_obj)
        mol_lib_prep_available[molecule_id] = ['New Library Preparation']
        if LibraryPreparation.objects.filter(molecule_id = molecule_obj, sample_id = sample_obj).exists():
            libs_prep_obj =  LibraryPreparation.objects.filter(molecule_id = molecule_obj, sample_id = sample_obj)
            for lib_prep_obj in libs_prep_obj :
                lib_prep_available.append(lib_prep_obj.get_lib_prep_code())
                mol_lib_prep_available[molecule_id].append(lib_prep_obj.get_lib_prep_code())

    sample_recorded['rep_filter_selection'] = []
    for key, value in mol_lib_prep_available.items():
        sample_recorded['rep_filter_selection'].append([key, value])
    sample_recorded['molecule_available'] = list(mol_lib_prep_available.keys())
    sample_recorded['lib_prep_available'] = lib_prep_available
    return sample_recorded

def analyze_reprocess_data(json_data, reprocess_id, reg_user):
    '''
    Description:
        The function will get the option of reprocessing sample and it updates the sample state for reprocessing.
        In case that a new library preparation was required a new object is created.
    Input:
        json_data       # data collected from the user
        reprocess_id    # sample id to reprocess
        reg_user        # register user
    Functions:
        update_sample_reused : located at iSkyLIMS_core/utils/handling_samples.py
        update_molecule_reused : located at iSkyLIMS_core/utils/handling_samples.
        get_sample_obj_from_id : located at iSkyLIMS_core/utils/handling_samples.
    Return:
        True if user request were right, or Invalid options if user requests were wrong.
    '''
    options = json_data[-3:]

    if "New Extraction" in options:
        sample_obj =  update_sample_reused(reprocess_id)
        sample_obj.set_state('Defined')
        return True

    elif 'New Library Preparation' in options:
        molecule_code_id = options[0]
        if molecule_code_id == '':
            return 'Invalid options'
        else:
            sample_obj = get_sample_obj_from_id (reprocess_id)
            if 'Library preparation' != sample_obj.get_sample_state():
                molecule_obj = update_molecule_reused(reprocess_id, molecule_code_id)
                if not molecule_obj:
                    return 'Invalid options'
                #
                sample_obj = update_sample_reused(reprocess_id)
                sample_obj.set_state('Library preparation')
                # create the new library preparation in "Created_for_reused" state
                new_library_preparation = LibraryPreparation.objects.create_reused_lib_preparation(reg_user, molecule_obj, sample_obj)
            return True
    elif 'New Pool' in options:
        molecule_code_id = options[0]
        lib_prep_code_id = options[1]

        if not LibraryPreparation.objects.filter(sample_id__pk__exact = reprocess_id, libPrepCodeID__exact = lib_prep_code_id).exists():
            return 'Invalid options'
        lib_prep_obj = LibraryPreparation.objects.get(sample_id__pk__exact = reprocess_id, libPrepCodeID__exact = lib_prep_code_id)
        sample_id = update_sample_reused(reprocess_id)
        molecule_obj = update_molecule_reused(reprocess_id, molecule_code_id)
        lib_prep_obj.set_state('Reused pool')
        reused = lib_prep_obj.set_increase_reuse()


    else:
        return 'Invalid options'


def get_sample_in_project_obj_from_id (sample_in_project_id):
    '''
    Description:
        The function gets the sampleInProject id and return the object
        Return the if of sampleInProject
    Input:
        sample_name     # sample name to look at
    Return:
        sample_in_project_obj.
    '''
    sample_in_project_obj = ''
    if SamplesInProject.objects.filter(pk__exact = sample_in_project_id).exists():
        sample_in_project_obj = SamplesInProject.objects.filter(sampleName__exact = sample_name)

    return sample_in_project_obj


def get_run_sample_id ( sample_name):
    '''
    Description:
        The function gets the sample name and if found it returns the SamplesInProject objects.
        Return the if of sampleInProject
    Input:
        sample_name     # sample name to look at
    Return:
        run_sample_obj.
    '''
    run_sample_obj = ''
    if SamplesInProject.objects.filter(sampleName__exact = sample_name).exists():
        sample_run_objs = SamplesInProject.objects.filter(sampleName__exact = sample_name)
        if len(sample_run_objs) > 1:
            pass
        else:
            run_sample_obj = SamplesInProject.objects.get(sampleName__exact = sample_name)
    return run_sample_obj

def search_run_samples(sample_name, user_name, start_date, end_date):
    '''
    Description:
        The function search the run samples that matchs with the requested conditions.
        Return the if of sampleInProject
    Input:
        sample_name     # sample name to look at
        user_name       # user name
        start_date      # date from starting the search
        end_date        # date from ending the search
    Functions:
        get_friend_list # located at iSkyLIMS_core/utils/generic_functions.py
    Return:
        run_sample_obj.
    '''
    run_sample_list = []

    if SamplesInProject.objects.all().exists():
        run_sample_founds = SamplesInProject.objects.all()
    else:
        return sample_list
    if user_name != '':
        user_name_obj = User.objects.get(username__exact = user_name)
        user_friend_list = get_friend_list(user_name_obj)
        if not run_sample_founds.filter(sampleUser__in = user_friend_list).exists():
            return run_sample_list
        else :
            run_sample_founds = run_sample_founds.filter(sampleUser__in = user_friend_list)
    if sample_name != '' :
        if run_sample_founds.filter(sampleName__exact = sample_name).exists():
            run_sample_founds = run_sample_founds.filter(sampleName__exact = sample_name)
            if len(run_sample_founds) == 1:
                run_sample_list.append(run_sample_founds[0].pk)
                return run_sample_list

        elif run_sample_founds.filter(sampleName__icontains = sample_name).exists():
            run_sample_founds = run_sample_founds.filter(sampleName__icontains = sample_name)
        else:
            return run_sample_list
    if start_date !='' and end_date != '':
        run_sample_founds = run_sample_founds.filter(generated_at___range=(start_date, end_date ))

    if start_date !='' and end_date  == '':
        run_sample_founds = run_sample_founds.filter(generated_at__gte = start_date)

    if start_date =='' and end_date  != '':
            run_sample_founds = run_sample_founds.filter(generated_at__lte = end_date )

    if len(run_sample_founds) == 1:
        sample_list.append(run_sample_founds[0].pk)
        return run_sample_list

    for run_sample in run_sample_founds :
        run_sample_list.append(run_sample.get_info_for_searching())

    return run_sample_list
