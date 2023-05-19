from crispy_forms import bootstrap, layout
from crispy_forms.helper import FormHelper
from django import forms
from django.utils.translation import ugettext_lazy as _


class ContactForm(forms.Form):
    email_contact = forms.EmailField(label="Write your email ", required=True)
    subject = forms.CharField(
        label="Write the Subject of the email ", max_length=40, required=True
    )
    message = forms.CharField(
        label="Describe your request ", widget=forms.Textarea, required=True
    )

    def __init__(self, *args, **kwargs):
        super(ContactForm, self).__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.form_class = "container-md"
        self.helper.form_action = ""
        self.helper.form_method = "POST"

        self.helper.layout = layout.Layout(
            layout.Div(
                layout.Div(
                    layout.Field("email_contact"),
                    css_class="col-md-12",
                ),
                css_class="row",
            ),
            layout.Div(
                layout.Div(
                    layout.Field("subject"),
                    css_class="col-md-12",
                ),
                css_class="row",
            ),
            layout.Div(
                layout.Div(
                    layout.Field("message"),
                    css_class="col-md-12",
                ),
                css_class="row",
            ),
            layout.Div(
                layout.Div(
                    bootstrap.FormActions(
                        layout.Submit(
                            ("submit"), _("Submit")
                        ),
                    ),
                    css_class="col-md-auto",
                ),
                css_class="row justify-content-end"
            ),
        )
