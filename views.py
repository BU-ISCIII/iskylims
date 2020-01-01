from django.shortcuts import get_object_or_404, render, redirect
from django.http import HttpResponse
from django.contrib.auth.decorators import login_required
from django.conf import settings

#from django.core.urlresolvers import resolve

from iSkyLIMS_clinic.models import *
from iSkyLIMS_clinic.clinic_config import *
from iSkyLIMS_clinic.utils.managed_samples import *
from iSkyLIMS_clinic.utils.managed_results import *
from iSkyLIMS_clinic.utils.managed_patient import *

from iSkyLIMS_core.utils.handling_protocols import *
from iSkyLIMS_core.utils.handling_samples import *
from iSkyLIMS_core.utils.handling_commercial_kits import *

def index(request):
    #
    return render(request, 'iSkyLIMS_clinic/index.html')


@login_required
def add_commercial_kit (request):
    app_name = __package__.split('.')[0]
    # exclude the protocols that have no molecule involved
    defined_protocols = get_defined_protocols(app_name, True )
    commercial_kits_data = get_data_for_commercial_kits()

    if request.method == 'POST' and request.POST['action'] == 'addCommercialKit':
        if get_commercial_kit_id (request.POST['kitName']) :

            return render(request, 'iSkyLIMS_clinic/addCommercialKit.html',{'defined_protocols': defined_protocols, 'invalid_name': request.POST['kitName']})
        new_kit = store_commercial_kit(request.POST)
        new_kit_data = get_commercial_kit_basic_data(new_kit)
        return render(request, 'iSkyLIMS_clinic/addCommercialKit.html',{'new_kit_data': new_kit_data})
    else:
        return render(request, 'iSkyLIMS_clinic/addCommercialKit.html',{'defined_protocols': defined_protocols, 'commercial_kits_data': commercial_kits_data})

@login_required
def add_user_lot_commercial_kit (request):
    if request.method == 'POST' and request.POST['action'] == 'addUserLotKit':
        if get_lot_user_commercial_kit_id (request.POST['nickName']) :
            defined_kits = get_defined_commercial_kits()
            return render(request, 'iSkyLIMS_clinic/addUserLotCommercialKit.html',{'defined_kits': defined_kits, 'invalid_name': request.POST['nickName']})
        new_lot_kit = store_lot_user_commercial_kit(request.POST, request.user)
        new_lot_kit_data = get_lot_user_commercial_kit_basic_data(new_lot_kit)
        return render(request, 'iSkyLIMS_clinic/addUserLotCommercialKit.html',{'new_lot_kit_data':new_lot_kit_data})
    else:
        defined_kits = get_defined_commercial_kits()
        return render(request, 'iSkyLIMS_clinic/addUserLotCommercialKit.html',{'defined_kits':defined_kits})



@login_required
def add_result_data (request):

    if request.method == 'POST' and request.POST['action'] == 'updateClinicSampleProtocol':

        invalid, valid_c_samples_ids = record_result_protocol (request)
        if ('invalid_c_samples_ids' in invalid  and len(valid_c_samples_ids) == 0 ):
            c_samples_pending_protocol = get_clinic_sample_in_state('Pending protocol')
            result_protocol  = get_table_result_to_protocol(c_samples_pending_protocol)
            return render(request, 'iSkyLIMS_clinic/addResultData.html' ,{'result_protocol':result_protocol})

        show_result_parameters = define_table_for_result_parameters(valid_c_samples_ids)
        return render(request, 'iSkyLIMS_clinic/addResultData.html' ,{'show_result_parameters':show_result_parameters})
    elif request.method == 'POST' and request.POST['action'] == 'showTableResultParameters':
        if 'samples_in_list' in request.POST:
            valid_c_samples_ids = request.POST.getlist('c_samples')
        else:
            valid_c_samples_ids = request.POST['c_samples'].split(',')

        show_result_parameters = define_table_for_result_parameters(valid_c_samples_ids)
        return render(request, 'iSkyLIMS_clinic/addResultData.html' ,{'show_result_parameters':show_result_parameters})

    elif request.method == 'POST' and request.POST['action'] == 'addResultParameters':
        added_result_protocol_parameters = add_result_protocol_parameters(request)
        if 'pending' in request.POST :
            c_samples_pending = request.POST['pending'].split(',')
            show_result_parameters = define_table_for_result_parameters(c_samples_pending)
            return render(request, 'iSkyLIMS_clinic/addResultData.html' ,{'show_result_parameters':show_result_parameters})
        else:
            return render(request, 'iSkyLIMS_clinic/addResultData.html' ,{'added_result_protocol_parameters':added_result_protocol_parameters})



    else:
        if 'iSkyLIMS_wetlab' in settings.INSTALLED_APPS :
            for c_sample_state in ['Sequencing','Patient update']:
                check_if_need_update(c_sample_state)
        result_protocol = {}
        c_samples_pending_protocol = get_clinic_sample_in_state('Pending protocol')
        c_samples_pending_results = get_clinic_sample_in_state('Pending results')
        #pending['pending_results'] = get_clinic_samples_pending_results(request.user, 'Pending results')

        if  c_samples_pending_protocol :
            result_protocol['pending_protocol']  = get_table_result_to_protocol(c_samples_pending_protocol)
        if c_samples_pending_results:
            c_samples_ids = []
            for c_sample in c_samples_pending_results:
                c_samples_ids.append(c_sample.get_id())
            result_protocol['pending_results']  = define_table_for_result_parameters(c_samples_ids)

        if c_samples_pending_protocol or c_samples_pending_results :

            return render(request, 'iSkyLIMS_clinic/addResultData.html' ,{'result_protocol':result_protocol})
        else:
            return render(request, 'iSkyLIMS_clinic/addResultData.html' ,{'no_samples': True})

@login_required
def define_new_patient(request):
    if request.method == 'POST' and request.POST['action'] == 'defineNewPatient':
        import pdb; pdb.set_trace()
        defined_patient = create_new_patient(request.POST)
        if 'ERROR' in defined_patient:
            patient_definition_data = fields_for_new_patient()
            error_message = 'Patient Code already exists.'
            return render(request, 'iSkyLIMS_clinic/defineNewPatient.html' ,{'patient_definition_data': patient_definition_data,'error': error_message })
        return render(request, 'iSkyLIMS_clinic/defineNewPatient.html' ,{'defined_patient': defined_patient})
    else:
        patient_definition_data = fields_for_new_patient()


    return render(request, 'iSkyLIMS_clinic/defineNewPatient.html' ,{'patient_definition_data': patient_definition_data})

@login_required
def define_new_samples(request):
    if request.method == 'POST' and request.POST['action'] == 'recordsample':
        sample_recorded = analyze_input_samples (request)
        # if no samples are in any of the options, displays the inital page
        if (not 'valid_samples' in sample_recorded and not 'invalid_samples' in sample_recorded and not 'incomplete_samples'in sample_recorded) :
            sample_information = prepare_sample_input_table()
            return render(request, 'iSkyLIMS_clinic/defineNewSamples.html' ,{'sample_information':sample_information})

        if 'valid_samples' in sample_recorded :
            clinic_sample_list = []
            for sample_id in sample_recorded['valid_samples_ids']:
                new_clinic_sample = ClinicSampleRequest.objects.create(sampleCore = get_sample_obj_from_id(sample_id),
                            clinicSampleState = ClinicSampleState.objects.get(clinicState__exact = 'Defined'),
                            sampleRequestUser = request.user)
                new_clinic_sample.save()
                clinic_sample_list.append(new_clinic_sample.get_id())
            clinic_samples_ids = ','.join(clinic_sample_list)
            sample_recorded['clinic_samples_ids'] = clinic_samples_ids
        if 'incomplete_samples' in sample_recorded :
            sample_recorded.update(prepare_sample_input_table())
            sample_recorded['number_of_samples'] = len(sample_recorded['incomplete_samples'])

        return render(request, 'iSkyLIMS_clinic/defineNewSamples.html', {'sample_recorded':sample_recorded})
    else:
        sample_information = prepare_sample_input_table()
        return render(request, 'iSkyLIMS_clinic/defineNewSamples.html',{'sample_information':sample_information})

@login_required
def define_patient_information(request):
    if request.method == 'POST' and request.POST['action'] == 'continueWithPatient':
        #import pdb; pdb.set_trace()
        if 'samples_in_list' in request.POST:
            patient_information = prepare_patient_form(request.POST.getlist('c_samples'))
        else:
            patient_information = prepare_patient_form(request.POST['clinic_samples'].split(','))
        return render(request, 'iSkyLIMS_clinic/definePatientInformation.html',{'patient_information':patient_information})
    elif request.method == 'POST' and request.POST['action'] == 'storePatientInfo':
        updated_information = analyze_and_store_patient_data (request.POST, request.user)

        return render(request, 'iSkyLIMS_clinic/definePatientInformation.html',{'updated_information':updated_information})
    else:
        clinic_samples = get_clinic_samples_by_state('Defined')
        if not clinic_samples :
            return render(request, 'iSkyLIMS_clinic/definePatientInformation.html',{'no_samples':True})
        else:
            patient_information = prepare_patient_form(clinic_samples)

        return render(request, 'iSkyLIMS_clinic/definePatientInformation.html',{'patient_information':patient_information})


@login_required
def create_protocol (request):
    # get the list of defined protocols
    defined_protocols, other_protocol_list = display_available_protocols (__package__)
    defined_protocol_types = display_protocol_types (__package__)
    #import pdb; pdb.set_trace()

    if request.method == 'POST' and request.POST['action'] == 'addNewProtocol':
        new_protocol = request.POST['newProtocolName']
        protocol_type = request.POST['protocolType']
        description = request.POST['description']

        if check_if_protocol_exists (new_protocol, __package__):
            return render ( request,'iSkyLIMS_clinic/createProtocol.html',{'content':['Protocol Name ', new_protocol,
                            'Already exists.']})
        new_protocol_id = create_new_protocol(new_protocol, protocol_type, description, __package__)

        return render(request, 'iSkyLIMS_clinic/createProtocol.html',{'defined_protocols': defined_protocols,
                            'defined_protocol_types':defined_protocol_types, 'new_defined_protocol': new_protocol,
                            'new_protocol_id':new_protocol_id,  'other_protocol_list' :other_protocol_list})

    return render(request, 'iSkyLIMS_clinic/createProtocol.html',{'defined_protocols': defined_protocols,
                        'defined_protocol_types':defined_protocol_types, 'other_protocol_list' :other_protocol_list})


@login_required
def define_protocol_parameters (request, protocol_id):
    if request.method == 'POST' and request.POST['action'] == 'define_protocol_parameters':

        recorded_prot_parameters = set_protocol_parameters(request)

        return render(request, 'iSkyLIMS_clinic/defineProtocolParameters.html', {'recorded_prot_parameters':recorded_prot_parameters})

    else:
        if not check_if_protocol_exists(protocol_id, __package__):
            return render ( request,'iSkyLIMS_clinic/error_page.html',
                        {'content':['The requested Protocol does not exist',
                            'Create the protocol name before assigning custom protocol parameters.']})


        prot_parameters = define_table_for_prot_parameters(protocol_id)
        return render(request, 'iSkyLIMS_clinic/defineProtocolParameters.html', {'prot_parameters':prot_parameters})




@login_required
def define_result_protocol(request):
    defined_protocols, other_protocol_list = display_available_protocols (__package__)
    defined_protocol_types = display_protocol_types (__package__)
    if request.method == 'POST' and request.POST['action'] == 'addNewResultProtocol':
        new_result_protocol = request.POST['newResultProtocolName']
        protocol_type = request.POST['protocolType']
        description = request.POST['description']
        if check_if_protocol_exists (new_result_protocol, __package__):
            return render(request, 'iSkyLIMS_clinic/defineResultProtocol.html',{'other_protocol_list' :other_protocol_list,'defined_protocol_types':defined_protocol_types,'Error':new_result_protocol})

        new_result_protocol_id = create_new_protocol(new_result_protocol, protocol_type, description, __package__)

        return render(request, 'iSkyLIMS_clinic/defineResultProtocol.html',{'other_protocol_list' :other_protocol_list,
                            'defined_protocol_types':defined_protocol_types, 'new_defined_result_protocol': new_result_protocol,
                            'new_protocol_id':new_result_protocol_id})
        return render(request, 'iSkyLIMS_clinic/defineResultProcedure.html',{'recorded_result_parameters':recorded_result_parameters})
    else:
        return render(request, 'iSkyLIMS_clinic/defineResultProcedure.html',{'other_protocol_list' :other_protocol_list,'defined_protocol_types':defined_protocol_types})



@login_required
def define_result_protocol_parameters (request, result_protocol_id):
    if request.method == 'POST' and request.POST['action'] == 'defineResultProtocolParameters':
        recorded_result_prot_parameters = set_protocol_parameters(request)
        return render(request, 'iSkyLIMS_clinic/defineResultProtocolParameters.html',{'recorded_result_prot_parameters':recorded_result_prot_parameters})
    else:

        if not check_if_protocol_exists(result_protocol_id, __package__):
            return render ( request,'iSkyLIMS_clinic/error_page.html',
                        {'content':['The requested Protocol does not exist',
                            'Create the protocol name before assigning custom protocol parameters.']})
        result_parameters = define_table_for_prot_protocols(result_protocol_id)
        return render(request, 'iSkyLIMS_clinic/defineResultProtocolParameters.html',{'result_parameters':result_parameters})

@login_required
def display_patient_information (request, patient_id):

    display_patient_info = display_one_patient_info (patient_id)
    if 'ERROR' in display_patient_info:
        return render (equest, 'iSkyLIMS_clinic/displayPatientInformation.html', {'ERROR': 'ERROR'})
    else:
        return render(request, 'iSkyLIMS_clinic/displayPatientInformation.html', {'display_patient_info': display_patient_info })
    return

@login_required
def display_protocol (request, protocol_id):
    if not check_if_protocol_exists(protocol_id, __package__):
        return render (request,'iSkyLIMS_clinic/error_page.html',
            {'content':['The protocol that you are trying to get ',
                        'DOES NOT exists .']})
    protocol_data = get_all_protocol_info (protocol_id)
    #import pdb; pdb.set_trace()

    return render(request, 'iSkyLIMS_clinic/displayProtocol.html', {'protocol_data': protocol_data})



@login_required
def display_result_protocol (request, result_protocol_id):
    if not check_if_protocol_exists(result_protocol_id, __package__):
        return render (request,'iSkyLIMS_clinic/error_page.html',
            {'content':['The result protocol that you are trying to get ',
                        'DOES NOT exists .']})
    result_protocol_data = get_all_protocol_info (result_protocol_id)

    return render(request, 'iSkyLIMS_clinic/displayResultProtocol.html', {'result_protocol_data': result_protocol_data})


@login_required
def display_sample_info(request,sample_c_id):
    if check_if_sample_c_exists(sample_c_id):
        display_sample_info = display_one_sample_info (sample_c_id)
        #import pdb; pdb.set_trace()
        return render(request, 'iSkyLIMS_clinic/displaySampleInfo.html', {'display_sample_info': display_sample_info })
    else:
        search_sample_data = collect_data_for_search ()
        error_message = ['The clinic sample that you are trying to get', 'DOES NOT exists .']
        return render (request,'iSkyLIMS_clinic/searchSample.html', {'search_sample_data': search_sample_data,
                                'error_message' : error_message})


def pending_to_update(request):

    if 'iSkyLIMS_wetlab' in settings.INSTALLED_APPS :
        for c_sample_state in ['Sequencing','Patient update']:
            check_if_need_update(c_sample_state)
    pending = {}

    # get the samples in defined state
    pending['defined'] = get_clinic_samples_defined_state(request.user)

    pending['patient_update'] = get_clinic_samples_patient_sequencing_state(request.user, 'Patient update')

    pending['sequencing'] = get_clinic_samples_patient_sequencing_state(request.user, 'Sequencing')
    pending['pending_protocol'] = get_clinic_samples_pending_results(request.user, 'Pending protocol')
    pending['pending_results'] = get_clinic_samples_pending_results(request.user, 'Pending results')

    pending ['graphic_pending_samples'] = pending_clinic_samples_for_grafic(pending).render()

    return render(request, 'iSkyLIMS_clinic/pendingToUpdate.html', {'pending':pending})


@login_required
def search_sample(request):
    search_sample_data = collect_data_for_search ()
    if request.method == 'POST' and (request.POST['action'] == 'searchSample'):
        data_request = {}
        data_request['sample_name']  = request.POST['samplename']
        data_request['patient_name'] = request.POST['patientname']
        data_request['history_number'] = request.POST['historynumber']
        data_request['doctor_name'] = request.POST['doctor']
        data_request['requested_service_by'] = request.POST['requestedby']
        data_request['start_date'] = request.POST['startdate']
        data_request['end_date'] = request.POST['enddate']


        # check that some values are in the request if not return the form
        if len([v for v in data_request.values() if v !='']) == 0 :
            return render(request, 'iSkyLIMS_clinic/searchSample.html', {'search_sample_data': search_sample_data })
        # Check for valid date format
        if data_request['start_date'] != '':
            try:
                datetime.datetime.strptime(data_request['start_date'], '%Y-%m-%d')
            except:
                return render(request, 'iSkyLIMS_clinic/searchSample.html', {'search_sample_data': search_sample_data ,
                    'Error': ERROR_MESSAGE_FOR_INCORRECT_START_SEARCH_DATE })
        if data_request['end_date'] != '':
            try:
                datetime.datetime.strptime(end_date, '%Y-%m-%d')
            except:
                return render(request, 'iSkyLIMS_clinic/searchSample.html', {'search_sample_data': search_sample_data ,
                    'Error': ERROR_MESSAGE_FOR_INCORRECT_END_SEARCH_DATE })
        # Patient name length must be longer than 5 characters
        if data_request['patient_name'] !=''  and len(data_request['patient_name']) < 4 :
            return render(request, 'iSkyLIMS_clinic/searchSample.html', {'search_sample_data': search_sample_data ,
                'Error': ERROR_MESSAGE_FOR_SORT_PATIENT_NAME })
        sample_c_list = get_samples_clinic_in_search(data_request)

        if len(sample_c_list) == 0:
            return render(request, 'iSkyLIMS_clinic/searchSample.html', {'search_sample_data': search_sample_data ,
                'Error': ERROR_MESSAGE_FOR_NO_MATCH_IN_SEARCH })
        if len(sample_c_list) == 1:
            display_sample_info = display_one_sample_info (sample_c_list[0])

            return render(request, 'iSkyLIMS_clinic/displaySampleInfo.html', {'display_sample_info': display_sample_info })
        else:
            display_sample_list_info = display_sample_list(sample_c_list)
            return render(request, 'iSkyLIMS_clinic/displaySampleInfo.html', {'display_sample_list_info': display_sample_list_info })

    else:
        return render(request, 'iSkyLIMS_clinic/searchSample.html', {'search_sample_data': search_sample_data })


@login_required
def search_patient (request):
    if request.method == 'POST' and (request.POST['action'] == 'searchPatient'):
        data_request = {}
        data_request['p_name'] = request.POST['patientname']
        data_request['p_surname'] = request.POST['patientsurname']
        data_request['p_code'] = request.POST['patientcode']

        len_search_values = len(set(data_request.values()))

        if len_search_values == 1:
            return render(request, 'iSkyLIMS_clinic/searchPatient.html')

        if data_request['p_name'] !=''  and len(data_request['p_name']) < 4 :
            return render(request, 'iSkyLIMS_clinic/searchPatient.html', {'Error': ERROR_MESSAGE_FOR_SORT_PATIENT_NAME })
        patient_list = get_patients_in_search(data_request)

        if len(patient_list) == 0:
            return render(request, 'iSkyLIMS_clinic/searchPatient.html', {'Error': ERROR_MESSAGE_FOR_NO_MATCH_IN_SEARCH })
        if len(patient_list) == 1:
            return redirect ('display_patient_information', patient_id = patient_list[0] )
        else:
            display_patient_list_info = display_patient_list(patient_list)
            return render(request, 'iSkyLIMS_clinic/displayPatientInformation.html', {'display_patient_list_info': display_patient_list_info })

    else:
        return render(request, 'iSkyLIMS_clinic/searchPatient.html')

@login_required
def set_molecule_values(request):

    if request.method == 'POST' and request.POST['action'] == 'continueWithMolecule':
        if request.POST['c_samples'] == '':
            return render (request,'iSkyLIMS_clinic/error_page.html',
                {'content':['There was no sample selected ']})
        if  'samples_in_list' in request.POST:
            c_samples = request.POST.getlist('c_samples')
        else:
            c_samples = request.POST['c_samples'].split(',')

        molecule_protocol = get_table_record_molecule (c_samples, __package__)
        if 'ERROR' in molecule_protocol :
            return render (request, 'iSkyLIMS_clinic/error_page.html',
                {'content':['There was no valid sample selected ']})

        molecule_protocol['samples'] = ','.join(c_samples)
        #import pdb; pdb.set_trace()
        return render(request, 'iSkyLIMS_clinic/setMoleculeValues.html',{'molecule_protocol':molecule_protocol})

    elif request.method == 'POST' and request.POST['action'] == 'updateMoleculeProtocol':

        molecule_recorded = record_molecules (request)

        if not 'heading' in molecule_recorded:
            samples = request.POST['samples'].split(',')
            molecule_protocol = get_table_record_molecule (samples, __package__)
            molecule_protocol['data'] = molecule_recorded['incomplete_molecules']
            molecule_protocol['samples'] = ','.join(samples)
            return render(request, 'iSkyLIMS_clinic/setMoleculeValues.html',{'molecule_protocol':molecule_protocol})
        else:
            if 'incomplete_molecules' in molecule_recorded:
                samples = molecule_recorded['incomplete_molecules_ids'].split(',')
                molecule_recorded.update(get_table_record_molecule (samples, __package__))

            return render(request, 'iSkyLIMS_clinic/setMoleculeValues.html',{'molecule_recorded':molecule_recorded})

    elif request.method == 'POST' and request.POST['action'] == 'displayMoleculeParameters':
        if  'samples_in_list' in request.POST:
            molecules = request.POST.getlist('molecules')
        else:
            molecules = request.POST['molecules'].split(',')
        show_molecule_parameters = display_molecule_protocol_parameters(molecules, request.user)
        return render(request, 'iSkyLIMS_clinic/setMoleculeValues.html',{'show_molecule_parameters':show_molecule_parameters})

    elif request.method == 'POST' and request.POST['action'] == 'addMoleculeParameters':
        added_molecule_protocol_parameters , sample_updated_list = add_molecule_protocol_parameters(request)
        # Update the clinic sample request state
        for sample_updated in sample_updated_list:
            get_clinic_sample_obj_from_sample_id(sample_updated).set_state('Pending results')

        if 'pending' in request.POST :
            molecules = request.POST['pending'].split(',')
            show_molecule_parameters = display_molecule_protocol_parameters(molecules,request.user)
            return render(request, 'iSkyLIMS_clinic/setMoleculeValues.html',{'added_molecule_protocol_parameters':added_molecule_protocol_parameters, 'show_molecule_parameters':show_molecule_parameters})
        else:
            return render(request, 'iSkyLIMS_clinic/setMoleculeValues.html',{'added_molecule_protocol_parameters':added_molecule_protocol_parameters})

    else:
        register_user = request.user.username
        display_list = get_defined_samples (register_user)

        return render(request, 'iSkyLIMS_clinic/setMoleculeValues.html',{'display_list': display_list})
    return render(request, 'iSkyLIMS_clinic/setMoleculeValues.html',{})
