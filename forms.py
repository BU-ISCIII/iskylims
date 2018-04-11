## drylab/forms.py
from django import forms                                                                                                                           
from django.utils.translation import ugettext_lazy as _, ugettext                                                                                  
from crispy_forms.helper import FormHelper                                                                                                         
from crispy_forms import layout, bootstrap 
from crispy_forms.layout import Field
from crispy_forms.layout import Layout, Div, Submit, HTML, Button, Row, Field
#from crispy_forms.bootstrap import InlineField                                                                                                        
#from utils.fields import MultipleChoiceTreeField                                                                                                  
from mptt.forms import TreeNodeMultipleChoiceField                                                                                                 
from .models import *                                                                                                                              
#from django.contrib.auth.forms import UserCreationForm                                                                                             


## Management of a 'Services request':
# case a) (internal) GENOMICS_UNIT_SEQUENCING

class ServiceRequestFormInternalSequencing(forms.ModelForm):
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
 		super(ServiceRequestFormInternalSequencing, self).__init__(*args, **kwargs)
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

# case b) (external) EXTERNAL_SEQUENCING
class ServiceRequestFormExternalSequencing(forms.ModelForm):
 	class Meta:
 		model = Service
 		fields = ['serviceSeqCenter','servicePlatform','serviceRunSpecs','serviceFileExt','serviceAvailableService','serviceFile','serviceNotes']
 		field_classes = {
				'serviceAvailableService': TreeNodeMultipleChoiceField,
				}

 	def __init__(self,*args, **kwargs):
 		super(ServiceRequestFormExternalSequencing, self).__init__(*args, **kwargs)
 		self.helper = FormHelper()
 		self.helper.form_action=""
 		self.helper.form_method="POST"
 		self.helper.layout = layout.Layout(
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


## Management of request forms for services OTHER THAN 'Services Request' :
# This form 1) extends the one used for 'Services request' and 2) excludes fields not filled-in through
# the corresponding service form (they are filled in 'by hand' in views.py)
class ServiceRequestForm_extended(ServiceRequestFormExternalSequencing):

	class Meta:
		model = Service
		exclude = [
				'serviceSeqCenter',
				'serviceProjectNames',
				'serviceUserId',
				'serviceRequestNumber',
				'servicePlatform',
				'serviceRunSpecs',
				'serviceFileExt',
				'serviceStatus',
				'serviceOnApprovedDate',
				'serviceOnRejectedDate',
				'serviceRequestID',
				'serviceOnDeliveredDate'
				]

	def __init__(self,*args, **kwargs):
		super(ServiceRequestForm_extended, self).__init__(*args, **kwargs)
		self.helper.layout.pop(0)


class DateInput(forms.DateInput):
    input_type = 'date'


class AddResolutionService(forms.ModelForm):

	class Meta:
		model = Resolution
		fields = ['resolutionEstimatedDate','resolutionFullNumber','resolutionAsignedUser','resolutionNotes']
		widgets = {'resolutionEstimatedDate': DateInput(),}
		exclude = ['resolutionNumber',
			    'resolutionDate',
			    'resolutionOnQueuedDate',
			    'resolutionOnInProgressDate',
			    'resolutionServiceID',
			    ]
	radio_buttons = forms.ChoiceField(
		label = 'Action to be done for the service',
		choices = (
			('accepted', "Service is Accepted"), 
			('rejected', "Service is Rejected")
		),
		widget = forms.RadioSelect,
		initial = 'accepted',
		required = True,
	)

	def __init__(self,*args, **kwargs):
		super(AddResolutionService, self).__init__(*args, **kwargs)
		self.fields['resolutionAsignedUser'].queryset = User.objects.filter( groups__name='Admin_iSkyLIMS')
		self.helper = FormHelper()
		self.helper.form_class = 'form-horizontal'
		self.helper.label_class = 'col-lg-4'
		self.helper.field_class = 'col-lg-7'
		self.helper.form_action=""
		self.helper.form_method="POST"
		
		self.helper.layout = layout.Layout(

				layout.Div(
					layout.HTML(u"""<div class="panel-heading"><h3 class="panel-title">Resolution Form for Service </h3></div>"""),
 					layout.Div(
 						layout.Div(
 							layout.Field('resolutionEstimatedDate'),
 							css_class="col-md-10",
 						),
						layout.Div(
 							layout.Field('resolutionFullNumber'),
 							css_class="col-md-10",
 						),
						layout.Div(
 							layout.Field('resolutionAsignedUser'),
 							css_class="col-md-10",
 						),
						layout.Div(
 							layout.Field('resolutionNotes'),
 							css_class="col-md-10",
 						),
						layout.Div(
							Field ('radio_buttons' , ), 
							css_class="col-md-10",
						),
						layout.Div(
						    bootstrap.FormActions( layout.Reset(('Reset'),_('Reset')),
									    layout.Submit(('submit'),_('Save'),style='margin-left: 80px')),
						    css_class="col-md-10"
						),
						css_class="row panel-body",
					),
					css_class = "panel panel-default"
					),
				)

class AddDeliveryService (forms.ModelForm) :
	
	class Meta:
		model = Delivery
		fields = ['deliveryNotes']
		exclude = ['deliveryResolutionID', 'deliveryDate'
			    ]


	def __init__(self,*args, **kwargs):
		super(AddDeliveryService, self).__init__(*args, **kwargs)
		self.helper = FormHelper()
		self.helper.form_class = 'form-horizontal'
		self.helper.label_class = 'col-lg-4'
		self.helper.field_class = 'col-lg-7'
		self.helper.form_action=""
		self.helper.form_method="POST"
		
		self.helper.layout = layout.Layout(
				layout.Div(
					layout.HTML(u"""<div class="panel-heading"><h3 class="panel-title">Delivery Form for Service </h3></div>"""),
 					layout.Div(
						layout.Div(
							layout.Field('deliveryNotes'),
							css_class="col-md-10",
						),
						layout.Div(
							bootstrap.FormActions( layout.Reset(('Reset'),_('Reset')),
										layout.Submit(('submit'),_('Save'),style='margin-left: 80px')),
										#layout.Submit(('submit'),_('Save'), css_class= "col-md-2 offset-md-3")),
							css_class="col-md-10"
						),
						css_class="row panel-body",
					),
					css_class = "panel panel-default"
					),
				)


class ByDateUserStats (forms.Form):

	user_name = forms.CharField(
		label = "User Name",
		max_length = 40,
		required = True,
		)
	start_date = forms.DateField(widget=forms.TextInput( attrs={'type': 'date'} ) ,
		label = 'Start Date',
		required = False,
		)
	end_date = forms.DateField(widget=forms.TextInput( attrs={'type': 'date'} ) ,
		label = 'End Date',
		required = False,
		)
	
	
	def __init__(self,*args, **kwargs):
		super(ByDateUserStats, self).__init__(*args, **kwargs)
		self.helper = FormHelper()
		self.helper.form_class = 'form-horizontal'
		self.helper.label_class = 'col-lg-4'
		self.helper.field_class = 'col-lg-7'
		self.helper.form_action=""
		self.helper.form_method="POST"
		
		self.helper.layout = layout.Layout(
				layout.Div(
					layout.HTML(u"""<div class="panel-heading"><h3 class="panel-title">Form for getting User statatistics </h3></div>"""),
 					layout.Div(
						layout.Div(
							layout.Field('user_name'),
							css_class="col-md-10",
						),
						layout.Div(
							layout.Field('start_date'),
							css_class="col-md-10",
						),
						layout.Div(
							layout.Field('end_date'),
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

class ByUserStats (forms.Form):


	start_date = forms.DateField(widget=forms.TextInput( attrs={'type': 'date'} ) ,
		label = 'Start Date',
		required = False,
		)
	end_date = forms.DateField(widget=forms.TextInput( attrs={'type': 'date'} ) ,
		label = 'End Date',
		required = False,
		)
	
	
	def __init__(self,*args, **kwargs):
		super(ByDateUserStats, self).__init__(*args, **kwargs)
		self.helper = FormHelper()
		self.helper.form_class = 'form-horizontal'
		self.helper.label_class = 'col-lg-4'
		self.helper.field_class = 'col-lg-7'
		self.helper.form_action=""
		self.helper.form_method="POST"
		
		self.helper.layout = layout.Layout(
				layout.Div(
					layout.HTML(u"""<div class="panel-heading"><h3 class="panel-title">Form for getting User statatistics </h3></div>"""),
 					layout.Div(
						layout.Div(
							layout.Field('user_name'),
							css_class="col-md-10",
						),
						layout.Div(
							layout.Field('start_date'),
							css_class="col-md-10",
						),
						layout.Div(
							layout.Field('end_date'),
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
	
