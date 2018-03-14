## drylab/forms.py
from django import forms                                                                                                                           
from django.utils.translation import ugettext_lazy as _, ugettext                                                                                  
from crispy_forms.helper import FormHelper                                                                                                         
from crispy_forms import layout, bootstrap                                                                                                         
#from utils.fields import MultipleChoiceTreeField                                                                                                  
from mptt.forms import TreeNodeMultipleChoiceField                                                                                                 
from .models import *                                                                                                                              
from django.contrib.auth.forms import UserCreationForm                                                                                             

class ServiceRequestForm(forms.ModelForm):
 	class Meta:
 		model = Service
		##Addition of serviceProjectNames for
		# implementation of drop down menu to choose a project name of a list of projects
		# belonging to the logged-in user in the service request form
 		fields = ['serviceProjectNames','servicePlatform','serviceRunSpecs','serviceFileExt','serviceAvailableService','serviceFile','serviceNotes']
 		field_classes = {
				'serviceAvailableService': TreeNodeMultipleChoiceField,
				}

 	def __init__(self,*args, **kwargs):
 		super(ServiceRequestForm, self).__init__(*args, **kwargs)
 		self.helper = FormHelper()
 		self.helper.form_action=""
 		self.helper.form_method="POST"
 		
 		self.helper.layout = layout.Layout(
				layout.Div(              
 					layout.HTML(u"""<div class="panel-heading"><h3 class="panel-title">Sequencing Data</h3></div>"""),
 					layout.Div(
 						layout.Div(
 							layout.Field('serviceProjectNames'),
 							css_class="col-md-6",
 						),
 						layout.Div(
 							layout.Field('serviceRunSpecs'),
							layout.Field('serviceFileExt'),
							css_class="col-md-6",
						),
						css_class="row panel-body"
						)
,
 					layout.Div(
 						layout.Div(
 							layout.Field('servicePlatform'),
 							css_class="col-md-12",
 						),
						css_class="row panel-body"
						),

					css_class = "panel panel-default"
					),



                layout.Div(                                                                                                
                	layout.HTML(u"""<div class="panel-heading"><h3 class="panel-title">Service selection</h3></div>"""), 
                	layout.Div(                                                                                            
                		layout.Div(                                                                                        
                			layout.Field('serviceAvailableService',template="utils/checkbox_select_multiple_tree.html"),                                                                   
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
class ServiceRequestForm_extended(ServiceRequestForm):

	class Meta:
		model = Service
		exclude = [
				'serviceSeqCenter',
				'servicePlatform',
				'serviceRunSpecs',
				'serviceFileExt',
				'serviceStatus',
				]

	def __init__(self,*args, **kwargs):
		super(ServiceRequestForm_extended, self).__init__(*args, **kwargs)
		self.helper.layout.pop(0)
