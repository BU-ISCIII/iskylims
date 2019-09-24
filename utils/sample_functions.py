from iSkyLIMS_wetlab.models import *
from iSkyLIMS_core.utils.handling_samples import *



def get_available_codeID_for_resequencing(sample_recorded):
    mol_lib_prep_available = {}
    lib_prep_available = ['New Library Preparation']
    mol_lib_prep_available['New Extraction'] =['']
    sample_obj = get_sample_obj_from_id (sample_recorded['sample_id_for_action'])
    molecules_obj = get_molecule_obj_from_sample(sample_recorded['sample_id_for_action'])

    for molecule_obj in molecules_obj:
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
    options = json_data[-3:]

    if "New Extraction" in options:
        sample_id= update_sample_reused(reprocess_id)

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

def update_batch_lib_prep_sample_state(lib_prep_ids,  sample_state):
    for lib_id in lib_prep_ids:
        lib_obj = LibraryPreparation.objects.get(pk__exact = lib_id)
        lib_obj.get_sample_obj().set_state(sample_state)
    return
