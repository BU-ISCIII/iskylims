from django.shortcuts import get_object_or_404, render, redirect 
from .models import *  
from .forms import *
from django.db import transaction
from django.contrib.auth.models import User , Group


def user_edit(request):
    if request.method == "POST":
        
        items = request.POST.keys()
        shared_list = []
        for item in items:
            if 'shared' in item :
                shared_list.append(request.POST[item])
        form1 = UserCreationForm(data=request.POST,instance=request.user)  
        #pdb.set_trace()
        form2 = ProfileCreationForm(data=request.POST,instance=request.user.profile) 
        #pdb.set_trace()
        if form1.is_valid() and form2.is_valid():
            user = form1.save()  
            profile = form2.save(commit=False)  
            # check the sharing list
            if len(shared_list) > 0 :
                #import pdb; pdb.set_trace()
                # create the group for the user to share their projects
                if Group.objects.filter(name__exact = request.user.username).exists():
                    new_group = Group.objects.get(name__exact = request.user.username)
                else:
                    new_group = Group (name = request.user.username)
                    new_group.save()
                groups = Group.objects.get(name = request.user.username)
                for item in shared_list :
                    if User.objects.filter(username__exact = item).exists():
                        user = User.objects.get(username__exact = item)
                        # check if user has already in the shared list
                        #import pdb; pdb.set_trace()
                        
                        if groups not in user.groups.all():
                            new_group.user_set.add(user)
                        
                import pdb; pdb.set_trace()
                
                
            return render(request,'django_utils/info_page.html',{'content':["Your user has been successfully updated"]})
        else:
            return render(request,'django_utils/error_page.html',{'content':[form1.errors,form2.errors]})

    else:
        sharing_list = []
        form1 = UserCreationForm(instance=request.user)  
        form2 = ProfileCreationForm(instance=request.user.profile) 
        #form3 = UserShared(instance=request.user)
        # get the list of users that are sharing 
        #import pdb; pdb.set_trace()
        if Group.objects.filter(name__exact = request.user.username).exists():
            import pdb; pdb.set_trace()
            group = Group.objects.get(name__exact = request.user.username )
            users = group.user_set.all().exclude(username = request.user.username)
            
            for user in users:
                sharing_list.append(user.username)
        
        #form3.fields['username'].queryset = User.objects.all()
        #import pdb; pdb.set_trace()
    return render(request,'registration/user_creation.html',{'form1' : form1 ,'form2':form2, 'sharing_list':sharing_list })  

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
