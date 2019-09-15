from django.shortcuts import render
from django.shortcuts import get_object_or_404, render, redirect
from django.http import HttpResponse
from django.contrib.auth.decorators import login_required
from iSkyLIMS_clinic.models import *
from iSkyLIMS_clinic.clinic_config import *

from iSkyLIMS_core.utils.handling_samples import *

def index(request):
    #
    return render(request, 'iSkyLIMS_clinic/index.html')

@login_required
def record_samples(request):
    if request.method == 'POST' and request.POST['action'] == 'recordsample':
        sample_recorded = analize_input_samples (request)
        # if no samples are in any of the options, displays the inital page
        if (len(sample_recorded['valid_samples']) == 0 and
                len(sample_recorded['not_valid_samples_same_user']) == 0 and
                len(sample_recorded['not_valid_samples_other_user']) == 0):
            sample_information = prepare_sample_input_table()
            return render(request, 'iSkyLIMS_wetlab/recordSample.html',{'sample_information':sample_information})
        else :
            return render(request, 'iSkyLIMS_wetlab/recordSample.html',{'sample_recorded':sample_recorded})
    else:
        sample_information = prepare_sample_input_table()
        sample_information['heading'] += ADDITIONAL_HEADING_FOR_RECORDING_SAMPLES
        return render(request, 'iSkyLIMS_wetlab/recordSample.html',{'sample_information':sample_information})

