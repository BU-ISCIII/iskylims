# from django.shortcuts import render
from django.core.mail import BadHeaderError, send_mail
from django.http import HttpResponse
from django.shortcuts import redirect, render
from .forms import ContactForm
from core.models import Contact
from .utils.common import *

# from django.conf import settings

# Create your views here.


def index(request):
    apps_in_iskylims = get_installed_apps()
    return render(request, "core/index.html", {"apps_in_iskylims": apps_in_iskylims})


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

    """
    apps_installed = {}
    apps_installed["apps_names"] = get_installed_apps()

    if request.method == "POST" and request.POST["action"] == "addNewContacts":
        pass

    return render(
        request, "core/addNewContacts.html", {"apps_installed": apps_installed}
    )


def contact(request):
    contact_info = {}
    contacts = Contact.objects.all()
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
