from django.shortcuts import get_object_or_404, render, redirect 
from .models import *                                            
from .forms import *                                             

# Create your views here.

def user_edit(request):                                                                                                                                                       
  	if request.method == "POST":                                                                                                                                              
  		form1 = UserCreationForm(data=request.POST,instance=request.user)                                                                                                     
  		#pdb.set_trace()                                                                                                                                                      
  		form2 = ProfileCreationForm(data=request.POST,instance=request.user.profile)                                                             
  		#pdb.set_trace()                                                                                                                                                      
  		if form1.is_valid() and form2.is_valid():                                                                                                                              
  			user = form1.save()                                                                                                                                               
  			profile = form2.save(commit=False)                                                                                                                                
  			return render(request,'utils/info_page.html',{'content':["Your user has been successfully created"]})
  		else:
  			return render(request,'utils/error_page.html',{'content':[form1.errors,form2.errors]})
  		
  	else:                                                                                                                                                                     
  		form1 = UserCreationForm(instance=request.user)                                                                                                                       
  		form2 = ProfileCreationForm(instance=request.user.profile)                                                                               
  		                                                                                                                                                                      
  		return render(request,'registration/user_creation.html',{'form1' : form1 ,'form2':form2})                                                                             
                                                                                                                                                                               
                                                                                                                                                                               
def user_creation(request):                                                                                                                                                   
 	if request.method == "POST":                                                                                                                                              
 		form1 = UserCreationForm(data=request.POST)                                                                                                                           
 		form2 = ProfileCreationForm(data=request.POST)                                                                                                                        
 		if form1.is_valid() and form2.is_valid():                                                                                                                              
 			#form1.save()
 			import pdb ; pdb.set_trace()
 			#form2.save(commit=False)
 			#form2.save()
 			return render(request,'utils/info_page.html',{'content':["Your user has been successfully created"]})
 		else:
 			return render(request,'utils/error_page.html',{'content':[form1.errors,form2.errors]}) 
 	else:                                                                                                                                                                     
 		form1 = UserCreationForm()                                                                                                                                            
 		form2 = ProfileCreationForm()                                                                                                                                         
 		                                                                                                                                                                      
 		return render(request,'registration/user_creation.html',{'form1' : form1 ,'form2':form2})                                                                             
