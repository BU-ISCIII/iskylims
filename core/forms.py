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
        self.helper.form_class = "form-horizontal"
        self.helper.label_class = "col-lg-5"
        self.helper.field_class = "col-lg-7"
        self.helper.form_action = ""
        self.helper.form_method = "POST"

        self.helper.layout = layout.Layout(
            layout.Div(
                layout.HTML(
                    """<div class="card-header">
                        <h3 class="panel-title">Write your email and your request </h3>
                        </div>"""
                ),
                layout.Div(
                    layout.Div(
                        layout.Field("email_contact"),
                        css_class="col-md-10",
                    ),
                    layout.Div(
                        layout.Field("subject"),
                        css_class="col-md-10",
                    ),
                    layout.Div(
                        layout.Field("message"),
                        css_class="col-md-10",
                    ),
                    layout.Div(
                        bootstrap.FormActions(
                            layout.Reset(("Reset"), _("Reset")),
                            layout.Submit(
                                ("submit"), _("Submit"), style="margin-left: 80px"
                            ),
                        ),
                        css_class="col-md-10",
                    ),
                    css_class="row card-body",
                ),
                css_class="card ",
            ),
        )
