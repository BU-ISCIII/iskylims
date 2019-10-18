from django.shortcuts import get_object_or_404, render, redirect
from django.http import HttpResponse
from django.contrib.auth.decorators import login_required

#from django.core.urlresolvers import resolve

from iSkyLIMS_clinic.models import *
from iSkyLIMS_clinic.clinic_config import *
from iSkyLIMS_clinic.utils.managed_samples import *
from iSkyLIMS_clinic.utils.managed_results import *

from iSkyLIMS_core.utils.handling_protocols import *
from iSkyLIMS_core.utils.handling_samples import *

def index(request):
    #
    return render(request, 'iSkyLIMS_clinic/index.html')

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
def define_result_protocol(request):
    defined_protocols, other_protocol_list = display_available_protocols (__package__)
    defined_protocol_types = display_protocol_types (__package__)
    if request.method == 'POST' and request.POST['action'] == 'addNewResultProtocol':
        new_result_protocol = request.POST['newResultProtocolName']
        protocol_type = request.POST['protocolType']
        description = request.POST['description']
        if check_if_protocol_exists (new_result_protocol, __package__):
            return render(request, 'iSkyLIMS_clinic/defineResultProtocol.html',{'other_protocol_list' :other_protocol_list,'defined_protocol_types':defined_protocol_types,'Error':new_result_protocol})

        new_result_protocol_id = create_new_protocol(new_result_protocol, protocol_type, description)
        import pdb; pdb.set_trace()
        return render(request, 'iSkyLIMS_clini/defineResultProtocol.html',{'other_protocol_list' :other_protocol_list,
                            'defined_protocol_types':defined_protocol_types, 'new_defined_protocol': new_result_protocol,
                            'new_protocol_id':new_result_protocol_id})



        return render(request, 'iSkyLIMS_clinic/defineResultProtocol.html',{'recorded_result_parameters':recorded_result_parameters})
    else:


        return render(request, 'iSkyLIMS_clinic/defineResultProtocol.html',{'other_protocol_list' :other_protocol_list,'defined_protocol_types':defined_protocol_types})

@login_required
def display_result_protocol (request, result_protocol_id):
    if not check_if_protocol_exists(result_protocol_id, __package__):
        return render (request,'iSkyLIMS_clinic/error_page.html',
            {'content':['The result protocol that you are trying to get ',
                        'DOES NOT exists .']})
    result_protocol_data = get_all_protocol_info (result_protocol_id)
    import pdb; pdb.set_trace()
    return render(request, 'iSkyLIMS_clinic/displayResultProtocol.html', {'result_protocol_data': result_protocol_data})

def display_sample_info(request,sample_c_id):
    if check_if_sample_c_exists(sample_c_id):
        display_sample_info = display_one_sample_info (sample_c_id)

        return render(request, 'iSkyLIMS_clinic/displaySampleInfo.html', {'display_sample_info': display_sample_info })
    else:
        return render (request,'iSkyLIMS_clinic/error_page.html',
            {'content':['The clinic sample that you are trying to get ',
                        'DOES NOT exists .']})

@login_required
def define_result_protocol_parameters (request, result_protocol_id):
    if request.method == 'POST' and request.POST['action'] == 'defineResultProtocolParameters':
        recorded_result_prot_parameters = set_protocol_parameters(request)
        return render(request, 'iSkyLIMS_clinic/defineResultProtocolParameters.html',{'recorded_result_prot_parameters':recorded_result_prot_parameters})
    else:

        if not check_if_protocol_exists(result_protocol_id, __package__):
            return render ( request,'iSkyLIMS_wetlab/error_page.html',
                        {'content':['The requested Protocol does not exist',
                            'Create the protocol name before assigning custom protocol parameters.']})
        result_parameters = define_table_for_prot_parameters(result_protocol_id)
        return render(request, 'iSkyLIMS_clinic/defineResultProtocolParameters.html',{'result_parameters':result_parameters})

@login_required
def search_sample(request):
    search_sample_data = {}
    search_sample_data['doctors'] = get_available_doctor()
    search_sample_data['requested_service_by'] = get_service_units()
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
        if data_request['patient_name'] !=''  and len(data_request['patient_name']) <5 :
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
