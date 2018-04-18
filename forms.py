from django import forms
from django.utils.translation import ugettext_lazy as _, ugettext                                                                                  
from crispy_forms.helper import FormHelper                                                                                                         
from crispy_forms import layout, bootstrap 
from crispy_forms.layout import Field
from crispy_forms.layout import Layout, Div, Submit, HTML, Button, Row, Field

class ContactForm(forms.Form):

    #from_email = forms.EmailField(required=True)
    

    #subject = forms.CharField(required=True)
    #message = forms.CharField(widget=forms.Textarea,required=True)
    
    email_contact = forms.EmailField( label = "Write your email ",  required = True)
    subject = forms.CharField( label = "Write the Subject of the email ", max_length = 40, required = True)
    message = forms.CharField( label = "Describe your request ", widget=forms.Textarea, required = True)
    
    def __init__(self,*args, **kwargs):
        super(ContactForm, self).__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.form_class = 'form-horizontal'
        self.helper.label_class = 'col-lg-5'
        self.helper.field_class = 'col-lg-7'
        self.helper.form_action=""
        self.helper.form_method="POST"

        self.helper.layout = layout.Layout(
            layout.Div(
                layout.HTML(u"""<div class="panel-heading"><h3 class="panel-title">Write your email and your request </h3></div>"""),
                layout.Div(
                    layout.Div(
                        layout.Field('email_contact'),
                        css_class="col-md-10",
                    ),
                    layout.Div(
                        layout.Field('subject'),
                        css_class="col-md-10",
                    ),
                    layout.Div(
                        layout.Field('message'),
                        css_class="col-md-10",
                    ),
                    layout.Div(
                        bootstrap.FormActions( layout.Reset(('Reset'),_('Reset')),
                            layout.Submit(('submit'),_('Submit'),style='margin-left: 80px')),
                        css_class="col-md-10"
                    ),
                    css_class="row panel-body",
                ),
                css_class = "panel panel-default"
            ),
        )
