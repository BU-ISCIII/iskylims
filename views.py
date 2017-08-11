# -*- coding: utf-8 -*-
from django.shortcuts import get_object_or_404, render, redirect
from django.http import HttpResponse
from django.template import loader
import time
from django.conf import settings
from django.core.files.storage import FileSystemStorage
from .models import *
from .forms import *
from django.core.files.storage import FileSystemStorage
import re

####### Import libraries for static files
#from django.shortcuts import render_to_response
#from django.shortcuts import RequestContext
# import pdb; pdb.set_trace()


def index(request):
    #latest_question_list = Question.objects.order_by('-pub_date')[:5]
    #context = {'latest_question_list': latest_question_list}
    return render(request, 'drylab/index.html')

def service_request(request):
# 	if request.method == "POST":
# 		form = ServiceRequestForm(data.request.POST,
# 		files=request.FILES,
# 		)
# 		if form.is_valid():
# 			return redirect("service request done")
	form = ServiceRequestForm("Genomic Data Analysis")	
	return render(request, 'drylab/ServiceRequestForm.html' , { 'form' : form })
