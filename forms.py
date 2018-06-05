from django import forms
from django.utils.translation import ugettext_lazy as _, ugettext
from crispy_forms.helper import FormHelper
from crispy_forms import layout, bootstrap
#from django_utils.fields import MultipleChoiceTreeField
from mptt.forms import TreeNodeMultipleChoiceField
from .models import *
#import pdb
from django.contrib.auth.forms import UserCreationForm


class ProfileCreationForm(forms.ModelForm):
 	class Meta:
 		model = Profile
 		fields = ['profilePosition','profileCenter','profileArea','profileExtension',]

 	def __init__(self,*args,**kwargs):
 		super(ProfileCreationForm,self).__init__(*args,**kwargs)
 		self.helper = FormHelper()
 		self.helper.form_action=""
 		self.helper.form_method="POST"
 		self.helper.form_tag=False
 		self.helper.csrf=False

 		self.helper.layout = layout.Layout(
          		layout.Div(
          			layout.HTML(u"""<div class="panel-heading"><h3 class="panel-title">User data</h3></div>"""),
          			layout.Div(
          				layout.Div(
          					layout.Field("profilePosition"),
          					layout.Field("profileCenter"),
          					css_class="col-md-6",
          				),
          				layout.Div(
          					layout.Field("profileArea"),
          					layout.Field("profileExtension"),
          					css_class="col-md-6",
          				),
          				css_class = 'row panel-body',
          			),
 		 			css_class = 'panel panel-default',
          			),
          		)


class UserCreationForm(UserCreationForm):

 	class Meta(UserCreationForm.Meta):
 		model = User
 		fields = ['username','email','password1','password2','first_name','last_name',]
 		help_texts ={
 				'email': _('Enter your email with domain @isciii.es or @externos.isciii.es'),
 				}

 	def __init__(self,*args, **kwargs):
 		super(UserCreationForm, self).__init__(*args, **kwargs)
 		self.helper = FormHelper()
 		self.helper.form_action=""
 		self.helper.form_method="POST"
 		self.helper.form_tag=False
 		self.helper.csrf=False

 		self.fields['email'].required=True
 		self.fields['first_name'].required=True
 		self.fields['last_name'].required=True

 		self.helper.layout = layout.Layout(
    				layout.Div(
    					layout.HTML(u"""<div class="panel-heading"><h3 class="panel-title">Researcher data</h3></div>"""),
    					layout.Div(
    						layout.Div(
    							layout.Field("username"),
    							layout.Field("first_name"),
    							layout.Field("password1"),
    							css_class="col-md-6",
    						),
    						layout.Div(
    							bootstrap.PrependedText("email","@",css_class="input-block-level",placeholder="contact@example.com"),
    							layout.Field("last_name"),
    							layout.Field("password2"),
    							css_class="col-md-6",
    						),
    						css_class = 'row panel-body',
    					),
    					css_class = 'panel panel-default'
    					),
      	)
