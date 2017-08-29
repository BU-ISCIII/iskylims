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
import pdb 
#pdb.set_trace()


def index(request):
    return render(request, 'drylab/index.html')

def service_request(request):
	if request.method == "POST":
		#pdb.set_trace()  
		form = ServiceRequestForm(data=request.POST,files=request.FILES,serviceFilter="Genomic data analysis")
		if form.is_valid():
			#pdb.set_trace()   
			form.save()
			return render(request,'drylab/index.html')
	else:
		form = ServiceRequestForm(serviceFilter="Genomic data analysis")	
	
	return render(request, 'drylab/ServiceRequestForm.html' , { 'form' : form })
