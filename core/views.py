# Generic imports
import django.contrib.auth.models
from django.contrib.auth.decorators import login_required
from django.core.mail import BadHeaderError, send_mail
from django.http import HttpResponse
from django.shortcuts import redirect, render

# Local imports
import core.forms
import core.models
import core.utils.common


def index(request):
    apps_in_iskylims = core.utils.common.get_installed_apps()
    return render(request, "core/index.html", {"apps_in_iskylims": apps_in_iskylims})


def contact(request):
    contact_info = {}
    contacts = core.models.Contact.objects.all()
    for contact in contacts:
        contact_info[contact.get_contact_name()] = contact.get_contact_email()

    if request.method == "GET":
        form = core.forms.ContactForm()
    else:
        form = core.forms.ContactForm(request.POST)
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


@login_required
def user_login(request):
    if not request.user.is_authenticated:
        return redirect("/accounts/login")
    if request.user.username != "admin":
        return redirect("")

    user_data = []
    login_data = {}
    user_list = django.contrib.auth.models.User.objects.all().order_by("-last_login")
    for user in user_list:
        user_data.append(
            [
                user.username,
                user.first_name,
                user.last_name,
                user.email,
                user.last_login,
            ]
        )
    login_data["user_data"] = user_data

    return render(request, "core/user_login.html", {"login_data": login_data})
