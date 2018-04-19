#from django.shortcuts import render

from django.core.mail import send_mail, BadHeaderError
from django.http import HttpResponse, HttpResponseRedirect
from django.shortcuts import render, redirect
from .forms import ContactForm
from django.core.mail import send_mail

from django.conf import settings

# Create your views here.

def index(request):
     return render(request, 'iSkyLIMS_home/index.html')



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
    return render(request, "iSkyLIMS_home/contact_email.html", {'form': form})

def thanks(request):
    return render(request, 'iSkyLIMS_home/thanks.html')
    #return HttpResponse('Thank you for your message.')
