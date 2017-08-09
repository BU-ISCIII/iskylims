## drylab/forms.py

from django import forms
from django.utils.translation import ugettext_lazy as _, ugettext
from crispy_forms.helper import FormHelper
from crispy_forms import layout, bootstrap 
from .models import *


class ServiceRequestForm(forms.ModelForm):
 	class Meta:
 		model = Service
 		fields = ['serviceName','servicePosition','serviceCenter','serviceArea','serviceExtension','serviceEmail','serviceSeqCenter','servicePlatform','serviceRunSpecs','serviceFileExt','serviceFile','serviceNotes']
 		
 	def __init__(self, *args, **kwargs):
 		super(ServiceRequestForm, self).__init__(*args, **kwargs)
 		self.helper = FormHelper()
 		self.helper.form_action=""
 		self.helper.form_method="POST"
 		#self.fields['serviceCenter'] = forms.ModelChoiceField(queryset=Center.objects.all())
 		#self.fields['serviceFileExt'] = forms.ModelChoiceField(queryset=Center.objects.all())
 			
 		self.helper.layout = layout.Layout(
 				layout.Fieldset(
 					_("Researcher data"),
 					layout.Field("serviceName"),
 					layout.Field("servicePosition"),
 					layout.Field("serviceCenter"),
 					layout.Field("serviceArea"),
 					layout.Field("serviceExtension"),
 					bootstrap.PrependedText("serviceEmail","@",css_class="input-block-level",placeholder="contact@example.com"),
 					),
 				layout.Fieldset(
 					_("Sequencing data"),
 					layout.Field('serviceSeqCenter'),
 					layout.Field('servicePlatform'),
 					layout.Field('serviceRunSpecs'),
					layout.Field('serviceFileExt'),
					),
				layout.Fieldset(
					_("Service Description"),
					layout.Field('serviceFile'),
					layout.Field('serviceNotes'),						
					),
				bootstrap.FormActions(
					layout.Submit(('submit'),_('Save')),
                    )
				)

