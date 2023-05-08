<<<<<<< HEAD
# from django.shortcuts import render
import json

from django.core.mail import BadHeaderError, send_mail
from django.http import HttpResponse
from django.shortcuts import redirect, render

from .forms import ContactForm
from .utils.common import *

# from django.conf import settings

# Create your views here.


=======
from django.core.mail import send_mail, BadHeaderError
from django.http import HttpResponse
from django.shortcuts import render, redirect
from .forms import ContactForm
from django.core.mail import send_mail

import core.utils.common 
import core.models
>>>>>>> Removing the iSkyLIMS_ references in html to statics and organize the statics to clean duplicated/triplicated files
def index(request):
    apps_in_iskylims = get_installed_apps()
    return render(request, "core/index.html", {"apps_in_iskylims": apps_in_iskylims})

<<<<<<< HEAD
=======
    apps_in_iskylims = core.utils.common.get_installed_apps ()
    return render(request, 'core/index.html',{'apps_in_iskylims': apps_in_iskylims})
>>>>>>> Removing the iSkyLIMS_ references in html to statics and organize the statics to clean duplicated/triplicated files

def add_new_contacts(request):
    """
    Description:
        The function will use to add new user detail contacts that are showed
        in the contact information.
        This function is only available for admin
    Input:
        request     # contains the request dictionary sent by django
    Variables:

    Return:

<<<<<<< HEAD
    """

    """
    library_kit_objects = BaseSpaceLibraryName.objects.all()
    if len(library_kit_objects) >0 :
        for l_kit in library_kit_objects :
            library_kits.append(l_kit.libraryName)
    """
    apps_installed = {}
    apps_installed["apps_names"] = get_installed_apps()

    if request.method == "POST" and request.POST["action"] == "addNewContacts":
=======
    '''
    apps_installed = {}
    apps_installed['apps_names'] = core.utils.common.get_installed_apps ()

    if request.method == 'POST' and request.POST['action'] == 'addNewContacts':
>>>>>>> Removing the iSkyLIMS_ references in html to statics and organize the statics to clean duplicated/triplicated files
        pass

    return render(
        request, "core/addNewContacts.html", {"apps_installed": apps_installed}
    )


def contact(request):
    contact_info = {}
    contacts = core.models.Contact.objects.all()
    for contact in contacts:
        contact_info[contact.get_contact_name()] = contact.get_contact_email()

    if request.method == "GET":
        form = ContactForm()
    else:
        form = ContactForm(request.POST)
        if form.is_valid():
            subject = form.cleaned_data["subject"]
            # from_email = form.cleaned_data['from_email']
            email_contact = form.cleaned_data["email_contact"]
            message = form.cleaned_data["message"]
            try:
                send_mail(subject, message, email_contact, ["bioinformatica@isciii.es"])
            except BadHeaderError:
                return HttpResponse("Invalid header found.")
            return redirect("thanks")
    return render(
        request, "core/contact_email.html", {"form": form, "contact_info": contact_info}
    )


def thanks(request):
    return render(request, "core/thanks.html")
    # return HttpResponse('Thank you for your message.')
