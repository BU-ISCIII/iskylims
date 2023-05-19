from crispy_forms import bootstrap, layout
from crispy_forms.helper import FormHelper
from django import forms
from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.models import User
from django.utils.translation import ugettext_lazy as _

from .models import *


class ProfileCreationForm(forms.ModelForm):
    class Meta:
        model = Profile
        fields = [
            "profile_position",
            "profile_center",
            "profile_area",
            "profile_extension",
        ]

    def __init__(self, *args, **kwargs):
        super(ProfileCreationForm, self).__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.form_action = ""
        self.helper.form_method = "POST"
        self.helper.form_tag = False
        self.helper.csrf = False
        self.fields["profile_position"].required = True
        self.fields["profile_center"].required = True
        self.fields["profile_area"].required = True

        self.helper.layout = layout.Layout(
            layout.Div(
                layout.HTML(
                    """<div class="card-header"><h3 class="panel-title">User data</h3></div>"""
                ),
                layout.Div(
                    layout.Div(
                        layout.Field("profile_position"),
                        layout.Field("profile_center"),
                        css_class="col-md-6",
                    ),
                    layout.Div(
                        layout.Field("profile_area"),
                        layout.Field("profile_extension"),
                        css_class="col-md-6",
                    ),
                    css_class="row card-body",
                ),
                css_class="card ",
            ),
        )


class UserCreationForm(UserCreationForm):
    class Meta(UserCreationForm.Meta):
        model = User
        fields = [
            "username",
            "email",
            "password1",
            "password2",
            "first_name",
            "last_name",
        ]
        help_texts = {
            "email": _(
                "Enter your email with domain @isciii.es or @externos.isciii.es"
            ),
        }

    def __init__(self, *args, **kwargs):
        super(UserCreationForm, self).__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.form_action = ""
        self.helper.form_method = "POST"
        self.helper.form_tag = False
        self.helper.csrf = False

        self.fields["email"].required = True
        self.fields["first_name"].required = True
        self.fields["last_name"].required = True

        self.helper.layout = layout.Layout(
            layout.Div(
                layout.HTML(
                    """<div class="card-header"><h3 class="panel-title">Researcher data</h3></div>"""
                ),
                layout.Div(
                    layout.Div(
                        layout.Field("username"),
                        layout.Field("first_name"),
                        layout.Field("password1"),
                        css_class="col-md-6",
                    ),
                    layout.Div(
                        bootstrap.PrependedText(
                            "email",
                            "@",
                            css_class="input-block-level",
                            placeholder="contact@example.com",
                        ),
                        layout.Field("last_name"),
                        layout.Field("password2"),
                        css_class="col-md-6",
                    ),
                    css_class="row card-body",
                ),
                css_class="card ",
            ),
        )
