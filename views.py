#from django.shortcuts import render
import json

from django.core.mail import send_mail, BadHeaderError
from django.http import HttpResponse, HttpResponseRedirect
from django.shortcuts import render, redirect
from .forms import ContactForm
from django.core.mail import send_mail

from .utils.generic_functions import *
#from django.conf import settings

# Create your views here.

def index(request):

    apps_in_iskylims = get_installed_apps ()

    #import pdb; pdb.set_trace()
    return render(request, 'iSkyLIMS_core/index.html',{'apps_in_iskylims': apps_in_iskylims})

def add_new_contacts (request):
    '''
    Description:
        The function will use to add new user detail contacts that are showed
        in the contact information.
        This function is only available for admin
    Input:
        request     # contains the request dictionary sent by django
    Variables:

    Return:

    '''

    '''
    library_kit_objects = BaseSpaceLibraryName.objects.all()
    if len(library_kit_objects) >0 :
        for l_kit in library_kit_objects :
            library_kits.append(l_kit.libraryName)
    '''
    apps_installed = {}
    apps_installed['apps_names'] = get_installed_apps ()

    if request.method == 'POST' and request.POST['action'] == 'addNewContacts':
        import pdb; pdb.set_trace()
        pass

    return render(request, "iSkyLIMS_core/addNewContacts.html",{'apps_installed':apps_installed})

def contact(request):
    if request.method == 'GET':
        form = ContactForm()
    else:
        form = ContactForm(request.POST)
        if form.is_valid():
            subject = form.cleaned_data['subject']
            #from_email = form.cleaned_data['from_email']
            email_contact = form.cleaned_data['email_contact']
            message = form.cleaned_data['message']
            try:
                #send_mail(subject, message, email_contact, ['bioinformatica@isciii.es'])
                send_mail(subject, message, email_contact, ['luis.chapado@amgitt.es',])
                #send_mail(subject, message, email_contact, ['chapado.l@gmail.com',])
            except BadHeaderError:
                return HttpResponse('Invalid header found.')
            return redirect('thanks')
    return render(request, "iSkyLIMS_core/contact_email.html", {'form': form})

def thanks(request):
    return render(request, 'iSkyLIMS_core/thanks.html')
    #return HttpResponse('Thank you for your message.')
