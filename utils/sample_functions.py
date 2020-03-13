from iSkyLIMS_wetlab.models import *
from iSkyLIMS_core.utils.handling_samples import *



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
    Return:
        True if user request were right, or Invalid options if user requests were wrong.
    '''
    options = json_data[-3:]

    if "New Extraction" in options:
        sample_obj =  update_sample_reused(reprocess_id)
        return True

    elif 'New Library Preparation' in options:
        molecule_code_id = options[0]
        if molecule_code_id == '':
            return 'Invalid options'
        else:
            sample_id = update_sample_reused(reprocess_id)
            molecule_obj = update_molecule_reused(reprocess_id, molecule_code_id)
            if not molecule_obj:
                return 'Invalid options'
            # create the new library preparation in "Created_for_reused" state
            new_library_preparation = LibraryPreparation.objects.create_reused_lib_preparation(reg_user, molecule_obj, sample_id)
            return True
    elif 'New Pool' in options:
        pass
    else:
        return 'Invalid options'
