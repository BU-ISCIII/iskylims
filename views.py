from django.shortcuts import get_object_or_404, render, redirect 
from .models import *  
from .forms import *
from django.db import transaction
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
            return render(request,'django_utils/info_page.html',{'content':["Your user has been successfully created"]})
        else:
            return render(request,'django_utils/error_page.html',{'content':[form1.errors,form2.errors]})

    else:
        form1 = UserCreationForm(instance=request.user)  
        form2 = ProfileCreationForm(instance=request.user.profile) 

    return render(request,'registration/user_creation.html',{'form1' : form1 ,'form2':form2})  

@transaction.atomic
def user_creation(request):
    if request.method == "POST": 
        form1 = UserCreationForm(request.POST)
        form2 = ProfileCreationForm(request.POST)
        if form1.is_valid() and form2.is_valid():
            form1.save()
            profile = form2.save(commit=False)
            user_name = form1.cleaned_data['username']
            profileUserID = User.objects.get(username__exact = user_name)
            profile.profileUserID = profileUserID
          
            profile.save()
            return render(request,'django_utils/info_page.html',{'content':["Your user has been successfully created"]})
        else:
            return render(request,'django_utils/error_page.html',{'content':[form1.errors,form2.errors]}) 
    else:
        form1 = UserCreationForm()  
        form2 = ProfileCreationForm()  

    return render(request,'registration/user_creation.html',{'form1' : form1 ,'form2':form2})  
