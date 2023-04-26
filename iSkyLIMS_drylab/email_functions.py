#!/usr/bin/env python3

from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from subprocess import PIPE, Popen


# Create the message
def create_boby_msg(service, status):
    body_html = "<html><head><title>Skylims-Drylab</title></head>"
    body_html = body_html + "<body><div><h3>" + service + "</h3><br><br>"
    body_html = body_html + "<p>" + status + "</p>"
    body_html = body_html + "</div></body></html>"
    html = MIMEText(body_html, "html")
    return html


def create_mail(to_email, from_email, subject_email, service, status):
    """
    Description:
        create message  to send by email
    Input:
        msg_email # text  message to send by email
        to_email # who recieves the email
        from_email # who send the email
        subject_email #  mail reason
    Variables:

    Return:
        msg # list with the components of the mail
    """
    # html = MIMEText("<html><head><title>Skylims-Drylab</title></head><body>Servicio</body>", "html")

    html = create_boby_msg(service, status)
    msg = MIMEMultipart("alternative")
    msg["To"] = to_email
    msg["From"] = from_email
    msg["Subject"] = subject_email
    msg.attach(html)

    return msg


def send_mail_iskylims(msg):
    p = Popen(["/usr/sbin/sendmail", "-t", "-oi"], stdin=PIPE, universal_newlines=True)
    p.communicate(msg.as_string())


if __name__ == "__main__":
    to_email = "smzuniga@incliva.es"
    from_email = "iskylims@incliva.es"
    service = "Genomic data analysis- data analysis - Data Analysis-Tumor - Tumor-only samples"
    status = "In process"
    subject_email = "Prueba envío mail desde Python y servidor iSkyLIMS"
    msg = create_mail(to_email, from_email, subject_email, service, status)

    send_mail_iskylims(msg)
