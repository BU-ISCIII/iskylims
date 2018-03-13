# -*- coding: utf-8 -*-
from django.shortcuts import get_object_or_404, render, redirect
from .models import *
from .forms import *
from django.contrib.auth.decorators import login_required

####### Import libraries for static files
#from django.shortcuts import render_to_response
#from django.shortcuts import RequestContext

@login_required
def index(request):
    return render(request, 'drylab/index.html')

@login_required
def service_request(request):
	if request.method == "POST":
		##Create a form instance with POST data
		# more information in: 
                # https://docs.djangoproject.com/en/2.0/topics/forms/modelforms/
		form = ServiceRequestForm(data=request.POST,files=request.FILES)

		if form.is_valid():
			# Create but dont save for following modif		
			new_service = form.save(commit=False)
			# Modification of 'serviceStatus'
			new_service.serviceStatus = "recorded"
			# Local ("in these premises :-)") sequencing
			new_service.serviceSeqCenter = "GENOMIC_SEQ_UNIT"
			# Previous 'serviceUsername'refactored to 'serviceUserId' which describes better its real nature
			new_service.serviceUserId = User.objects.get(id=request.user.id)
			# Save the new instance
			new_service.save()
			# Save the many-to-many data for the form
			# (required for cases like this one where form.save(commit=False)) 
			form.save_m2m()
			return render(request,'utils/info_page.html',{'content':['Your service request has been successfully recorded.','You will be contacted shortly.']})
	else:
		form = ServiceRequestForm()
	
	form.fields['serviceAvailableService'].queryset = AvailableService.objects.filter(availServiceDescription__exact="Genomic data analysis").get_descendants(include_self=True) 
	## Addition of serviceProjectName for
	# implementation of drop down menu to choose a project name of a list of projects
	# belonging to the logged-in user in the service request form
	form.fields['serviceProjectNames'].queryset = Projects.objects.filter(user_id__exact = request.user.id)
	
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
