from django.conf import settings
from django.contrib.auth.models import Group, User
from django.db import transaction
from django.shortcuts import render

from .forms import *
from .models import *


def user_edit(request):
    if request.method == "POST":
        form1 = UserCreationForm(data=request.POST, instance=request.user)
        form2 = ProfileCreationForm(data=request.POST, instance=request.user.profile)
        if form1.is_valid() and form2.is_valid():
            email_address = form1.cleaned_data["email"]

            domain = email_address.split("@")[1]
            if len(settings.ALLOWED_EMAIL_DOMAINS) > 0:
                if domain not in settings.ALLOWED_EMAIL_DOMAINS:
                    error_description = str(
                        "Invalid email address . Your email domain "
                        + domain
                        + " is not allowed"
                    )
                    allowed_domains = " or ".join(settings.ALLOWED_EMAIL_DOMAINS)
                    return render(
                        request,
                        "django_utils/error_page.html",
                        {
                            "content": [
                                error_description,
                                "Only ",
                                allowed_domains,
                                " are allowed",
                            ]
                        },
                    )

            user = form1.save()
            form2.save(commit=False)
            # check the sharing list
            if Group.objects.filter(name__exact=request.user.username).exists():
                group = Group.objects.get(name__exact=request.user.username)
                if "shared_list" in request.POST:
                    existing_shared_list = []
                    shared_list = request.POST.getlist("shared_list")
                    existing_shared_users = group.user_set.all().exclude(
                        username=request.user.username
                    )
                    for existing_shared_user in existing_shared_users:
                        existing_shared_list.append(existing_shared_user.username)
                    delete_users = list(
                        set(existing_shared_list).difference(shared_list)
                    )
                    added_users = list(
                        set(shared_list).difference(existing_shared_list)
                    )
                    for delete_user in delete_users:
                        if User.objects.filter(username__exact=delete_user).exists():
                            user = User.objects.get(username__exact=delete_user)
                            group.user_set.remove(user)
                    for added_user in added_users:
                        if User.objects.filter(username__exact=added_user).exists():
                            user = User.objects.get(username__exact=added_user)
                            group.user_set.add(user)

                else:  # Remove all existing shared users
                    users = (
                        group.user_set.all()
                    )  # .exclude(username = request.user.username)
                    for user in users:
                        group.user_set.remove(user)
                    # delete group
                    group.delete()

            else:  # creating a new group and add user on it
                if "shared_list" in request.POST:
                    shared_list = request.POST.getlist("shared_list")
                    new_group = Group(name=request.user.username)
                    new_group.save()
                    # adding request user to its own group
                    new_group.user_set.add(
                        User.objects.get(username__exact=request.user.username)
                    )
                    for user in shared_list:
                        if User.objects.filter(username__exact=user).exists():
                            user = User.objects.get(username__exact=user)
                            new_group.user_set.add(user)

            return render(
                request,
                "django_utils/info_page.html",
                {"content": ["Your user has been successfully updated"]},
            )
        else:
            return render(
                request,
                "django_utils/error_page.html",
                {"content": [form1.errors, form2.errors]},
            )

    else:
        username_list = []
        form1 = UserCreationForm(instance=request.user)
        form2 = ProfileCreationForm(instance=request.user.profile)

        # get the list of users that are sharing
        # import pdb; pdb.set_trace()
        sharing_users = []
        sharing_list = []
        if Group.objects.filter(name__exact=request.user.username).exists():
            # import pdb; pdb.set_trace()
            group = Group.objects.get(name__exact=request.user.username)
            users = group.user_set.all().exclude(username=request.user.username)

            for user in users:
                sharing_users.append(user.username)
                sharing_list.append(
                    [user.username, user.first_name + " " + user.last_name]
                )
        # import pdb; pdb.set_trace()
        user_list = (
            User.objects.all()
            .exclude(username=request.user.username)
            .exclude(username__in=sharing_users)
        )
        for user in user_list:
            # do not include the user that are not well defined
            if user.first_name == "":
                continue
            username_list.append(
                [user.username, user.first_name + " " + user.last_name]
            )
        # import pdb; pdb.set_trace()
        if len(sharing_list) == 0:
            return render(
                request,
                "registration/user_creation.html",
                {"form1": form1, "form2": form2, "username_list": username_list},
            )
        else:
            return render(
                request,
                "registration/user_creation.html",
                {
                    "form1": form1,
                    "form2": form2,
                    "sharing_list": sharing_list,
                    "username_list": username_list,
                },
            )


@transaction.atomic
def user_creation(request):
    if request.method == "POST":
        form1 = UserCreationForm(request.POST)
        form2 = ProfileCreationForm(request.POST)
        if form1.is_valid() and form2.is_valid():
            email_address = form1.cleaned_data["email"]
            domain = email_address.split("@")[1]
            if len(settings.ALLOWED_EMAIL_DOMAINS) > 0:
                if domain not in settings.ALLOWED_EMAIL_DOMAINS:
                    error_description = str(
                        "Invalid email address . Your email domain "
                        + domain
                        + " is not allowed"
                    )
                    allowed_domains = " or ".join(settings.ALLOWED_EMAIL_DOMAINS)
                    return render(
                        request,
                        "django_utils/error_page.html",
                        {
                            "content": [
                                error_description,
                                "Only  ",
                                allowed_domains,
                                "are allowed",
                            ]
                        },
                    )

            form1.save()
            profile = form2.save(commit=False)
            user_name = form1.cleaned_data["username"]
            profileUserID = User.objects.get(username__exact=user_name)
            profile.profileUserID = profileUserID

            profile.save()
            return render(
                request,
                "django_utils/info_page.html",
                {"content": ["Your user has been successfully created"]},
            )
        else:
            return render(
                request,
                "django_utils/error_page.html",
                {"content": [form1.errors, form2.errors]},
            )
    else:
        form1 = UserCreationForm()
        form2 = ProfileCreationForm()

    return render(
        request, "registration/user_creation.html", {"form1": form1, "form2": form2}
    )


def check_user_group(request, group_name):
    """
    Description:
        The function is used to check if the loging user belongs to a group
    Input:
        request     # contains the request dictionary sent by django
        group_name  # contains the group name
    Return:
        True    # if logged user belongs to the group
        False   # if user does not belong to the group
    """
    group = Group.objects.get(name=group_name)
    if group not in request.user.groups.all():
        return False
    return True
