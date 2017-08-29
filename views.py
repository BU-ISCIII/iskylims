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
		form = ServiceRequestForm(data=request.POST,files=request.FILES,serviceFilter="Genomic data analysis")
		if form.is_valid():
			new_service = form.save(commit=False)
			new_service.serviceStatus = "Recorded"
			new_service.save()
			form.save_m2m()
			return render(request,'drylab/info_page.html',{'content':['Your service request has been successfully recorded.','You will be contacted shortly.']})
	else:
		form = ServiceRequestForm(serviceFilter="Genomic data analysis")	
	
	return render(request, 'drylab/RequestForm.html' , { 'form' : form })

def counseling_request(request):                                                                                                                                    
 	if request.method == "POST":                                                                                                                                 
 		form = ServiceRequestForm(data=request.POST,files=request.FILES,serviceFilter="Bioinformatics consulting and training")                                                   
 		if form.is_valid():                                                                                                                                      
 			new_service = form.save(commit=False)                                                                                                                
 			new_service.serviceStatus = "Recorded"                                                                                                               
 			new_service.save()                                                                                                                                   
 			form.save_m2m()                                                                                                                                      
 			return render(request,'drylab/info_page.html',{'content':['Your service request has been successfully recorded.','You will be contacted shortly.']}) 
 	else:                                                                                                                                                        
 		form = ServiceRequestForm(serviceFilter="Bioinformatics consulting and training")	                                                                                     
 	                                                                                                                                                             
 	return render(request, 'drylab/RequestForm.html' , { 'form' : form })                                                                                        

def infrastructure_request(request):                                                                                                                                     
  	if request.method == "POST":                                                                                                                                  
  		form = ServiceRequestForm(data=request.POST,files=request.FILES,serviceFilter="User support")                                                    
  		if form.is_valid():                                                                                                                                       
  			new_service = form.save(commit=False)                                                                                                                 
  			new_service.serviceStatus = "Recorded"                                                                                                                
  			new_service.save()                                                                                                                                    
  			form.save_m2m()                                                                                                                                       
  			return render(request,'drylab/info_page.html',{'content':['Your service request has been successfully recorded.','You will be contacted shortly.']})  
  	else:                                                                                                                                                         
  		form = ServiceRequestForm(serviceFilter="User support")	                                                                                      
  	                                                                                                                                                              
  	return render(request, 'drylab/RequestForm.html' , { 'form' : form })                                                                                         
