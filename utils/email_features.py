from django.core.mail import send_mail

def request_send_mail( subject, body_message, from_user, to_user):
    send_mail (subject, body_message, from_user, to_user)
