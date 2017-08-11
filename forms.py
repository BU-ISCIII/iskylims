## drylab/forms.py

from django import forms
from django.utils.translation import ugettext_lazy as _, ugettext
from crispy_forms.helper import FormHelper
from crispy_forms import layout, bootstrap 
from utils.fields import MultipleChoiceTreeField
from .models import *


class ServiceRequestForm(forms.ModelForm):
 	class Meta:
 		model = Service
 		fields = ['serviceName','servicePosition','serviceCenter','serviceArea','serviceExtension','serviceEmail','serviceSeqCenter','servicePlatform','serviceRunSpecs','serviceFileExt','serviceFile','serviceNotes']
 		
 	def __init__(self, serviceFilter ,*args, **kwargs):
 		super(ServiceRequestForm, self).__init__(*args, **kwargs)
 		self.serviceFilter = serviceFilter
 		self.helper = FormHelper()
 		self.helper.form_action=""
 		self.helper.form_method="POST"
 		
 		#self.fields['serviceCenter'] = forms.ModelChoiceField(queryset=Center.objects.all())
 		#self.fields['serviceFileExt'] = forms.ModelChoiceField(queryset=Center.objects.all())
 		
 		servicesAvailable = MultipleChoiceTreeField(
 				label=_("Services Available"),
 				required = True,
 				queryset=AvailableService.objects.filter(availServiceDescription__exact=self.serviceFilter).get_descendants(include_self=True)
 				)                                                                                                                         		
 		
 		self.helper.layout = layout.Layout(
 				layout.Div(
 					layout.HTML(u"""<div class="panel-heading"><h3 class="panel-title">Researcher data</h3></div>"""), 
 					layout.Div(
 						layout.Div(
 							layout.Field("serviceName"),
 							layout.Field("servicePosition"),
 							layout.Field("serviceCenter"),
 							css_class="col-md-6",
 						),
 						layout.Div(
 							layout.Field("serviceArea"),
 							layout.Field("serviceExtension"),
 							bootstrap.PrependedText("serviceEmail","@",css_class="input-block-level",placeholder="contact@example.com"),
 							css_class="col-md-6",
 						),
 						css_class = 'row panel-body',
 					),
 					css_class = 'panel panel-default'
 					),
				layout.Div(               
 					layout.HTML(u"""<div class="panel-heading"><h3 class="panel-title">Sequencing Data</h3></div>"""),
 					layout.Div(
 						layout.Div(
 							layout.Field('serviceSeqCenter'),
 							layout.Field('servicePlatform'),
 							css_class="col-md-6",
 						),
 						layout.Div(
 							layout.Field('serviceRunSpecs'),
							layout.Field('serviceFileExt'),
							css_class="col-md-6",
						),
						css_class="row panel-body"
						),
					css_class = "panel panel-default"
					),
                layout.Div(                                                                                                
                	layout.HTML(u"""<div class="panel-heading"><h3 class="panel-title">Service selection</h3></div>"""), 
                	layout.Div(                                                                                            
                		layout.Div(                                                                                        
                			layout.Field('servicesAvailable',template="utils/checkbox_select_multiple_tree.html"),                                                                   
                			css_class="col-md-12"                                                                          
                			),                                                                                             
                    	css_class="row panel-body"                                                                         
                    	),                                                                                                 
                	css_class = "panel panel-default"                                                                      
                	),                                                                                                     
				layout.Div(
					layout.HTML(u"""<div class="panel-heading"><h3 class="panel-title">Service Description</h3></div>"""),
					layout.Div(
						layout.Div(
							layout.Field('serviceFile'),
							layout.Field('serviceNotes'),
							css_class="col-md-12"
							),
                    	css_class="row panel-body"
                    	),
					css_class = "panel panel-default"
					),
				bootstrap.FormActions(
					layout.Submit(('submit'),_('Save')),
                    )
				)

