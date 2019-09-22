from django.shortcuts import get_object_or_404, render, redirect
from django.http import HttpResponse
from django.contrib.auth.decorators import login_required
from iSkyLIMS_clinic.models import *
from iSkyLIMS_clinic.clinic_config import *
from iSkyLIMS_clinic.utils.managed_samples import *

from iSkyLIMS_core.utils.handling_samples import *

def index(request):
    #
    return render(request, 'iSkyLIMS_clinic/index.html')

@login_required
def define_new_samples(request):
    if request.method == 'POST' and request.POST['action'] == 'recordsample':
        sample_recorded = analyze_input_samples (request)
        # if no samples are in any of the options, displays the inital page
        if (not 'valid_samples' in sample_recorded and not 'invalid_samples' in sample_recorded) :
            sample_information = prepare_sample_input_table()
            return render(request, 'iSkyLIMS_clinic/defineNewSamples.html' ,{'sample_information':sample_information})

        if 'valid_samples' in sample_recorded :
            clinic_sample_list = []
            for sample_id in sample_recorded['valid_samples_ids']:
                new_clinic_sample = ClinicSampleRequest.objects.create(sampleCore = get_sample_obj_from_id(sample_id),
                            clinicSampleState = ClinicSampleState.objects.get(clinicState__exact = 'Defined'))
                new_clinic_sample.save()
                clinic_sample_list.append(new_clinic_sample.get_id())
            clinic_samples_ids = ','.join(clinic_sample_list)
            sample_recorded['clinic_samples_ids'] = clinic_samples_ids

        return render(request, 'iSkyLIMS_clinic/defineNewSamples.html', {'sample_recorded':sample_recorded})
    if request.method == 'POST' and request.POST['action'] == 'continueWithPatient':
        patient_information = define_patient_information(request.POST['clinic_samples'])
        return render(request, 'iSkyLIMS_clinic/defineNewSamples.html',{'patient_information':patient_information})
    else:
        sample_information = prepare_sample_input_table()
        import pdb; pdb.set_trace()
        return render(request, 'iSkyLIMS_clinic/defineNewSamples.html',{'sample_information':sample_information})
