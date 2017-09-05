# -*- coding: utf-8 -*-
from django.shortcuts import get_object_or_404, render, redirect
from .models import *
from .forms import *
from django.contrib.auth.decorators import login_required

####### Import libraries for static files
#from django.shortcuts import render_to_response
#from django.shortcuts import RequestContext
import pdb 
#pdb.set_trace()

def user_edit(request):                                                                                       
 	if request.method == "POST":                                                                                  
 		form1 = UserCreationForm(data=request.POST,instance=request.user)                                                               
 		#pdb.set_trace()
 		form2 = ProfileCreationForm(data=request.POST,instance=UserInfo.objects.get(userInfoUserID=request.user))                                                            
 		#pdb.set_trace()                                                                                           
 		if form1.is_valid() or form2.is_valid():                                                                  
 			user = form1.save()                                                                                   
 			profile = form2.save(commit=False)                                                                    
 			profile.userInfoUserID = user                                                                         
 			profile.save()                                                                                        
 			return render(request,'utils/info_page.html',{'content':["Your user has been successfully created"]}) 
 	else:                                                                                                         
 		form1 = UserCreationForm(instance=request.user)                                                                                
 		form2 = ProfileCreationForm(instance=UserInfo.objects.get(userInfoUserID=request.user))                                                                             
 		                                                                                                          
 		return render(request,'registration/user_creation.html',{'form1' : form1 ,'form2':form2})                 



def user_creation(request):
	if request.method == "POST":
		form1 = UserCreationForm(data=request.POST)
		form2 = ProfileCreationForm(data=request.POST)
		pdb.set_trace()
		if form1.is_valid() or form2.is_valid():
			user = form1.save()
			profile = form2.save(commit=False)
			profile.userInfoUserID = user
			profile.save()
			return render(request,'utils/info_page.html',{'content':["Your user has been successfully created"]})
	else:
		form1 = UserCreationForm()
		form2 = ProfileCreationForm() 
		
		return render(request,'registration/user_creation.html',{'form1' : form1 ,'form2':form2})

@login_required
def index(request):
    return render(request, 'drylab/index.html')

@login_required
def service_request(request):
	if request.method == "POST":
		form = ServiceRequestForm(data=request.POST,files=request.FILES)
		if form.is_valid():
			new_service = form.save(commit=False)
			new_service.serviceStatus = "recorded"
			new_service.save()
			form.save_m2m()
			return render(request,'utils/info_page.html',{'content':['Your service request has been successfully recorded.','You will be contacted shortly.']})
	else:
		form = ServiceRequestForm()
	
	form.fields['serviceAvailableService'].queryset = AvailableService.objects.filter(availServiceDescription__exact="Genomic data analysis").get_descendants(include_self=True) 
	return render(request, 'drylab/RequestForm.html' , { 'form' : form })

@login_required
def counseling_request(request):                                                                                                                                    
 	if request.method == "POST":                                                                                                                                 
 		form = ServiceRequestForm_extended(data=request.POST,files=request.FILES)                                                   
 		if form.is_valid():                                                                                                                                      
 			new_service = form.save(commit=False)                                                                                                                
 			new_service.serviceStatus = "recorded"                                                                                                               
 			new_service.save()                                                                                                                                   
 			form.save_m2m()                                                                                                                                      
 			return render(request,'utils/info_page.html',{'content':['Your service request has been successfully recorded.','You will be contacted shortly.']}) 
 	else:                                                                                                                                                        
 		form = ServiceRequestForm_extended()	                                                                                     
 	
 	form.fields['serviceAvailableService'].queryset = AvailableService.objects.filter(availServiceDescription__exact="Bioinformatics consulting and training").get_descendants(include_self=True)                                                                                                                                                              
 	return render(request, 'drylab/RequestForm.html' , { 'form' : form })                                                                                        

@login_required
def infrastructure_request(request):                                                                                                                                     
  	if request.method == "POST":                                                                                                                                  
  		form = ServiceRequestForm_extended(data=request.POST,files=request.FILES)                                                    
  		#pdb.set_trace()
  		if form.is_valid():                                                                                                                                       
  			new_service = form.save(commit=False)                                                                                                                 
  			new_service.serviceStatus = "recorded"                                                                                                                
  			#pdb.set_trace()
  			new_service.save()                                                                                                                                    
  			form.save_m2m()                                                                                                                                       
  			return render(request,'utils/info_page.html',{'content':['Your service request has been successfully recorded.','You will be contacted shortly.']})  
  	else:                                                                                                                                                         
  		form = ServiceRequestForm_extended()	                                                                                      
  	
  	form.fields['serviceAvailableService'].queryset = AvailableService.objects.filter(availServiceDescription__exact="User support").get_descendants(include_self=True)                                                                                                                                                                                                                                                                                                                            
  	#pdb.set_trace()
  	#form.helper[1].update_atrributes(hidden="true")
  	return render(request, 'drylab/RequestForm.html' , { 'form' : form })                                                                                         
