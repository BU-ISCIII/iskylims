## drylab/forms.py

from django import forms
from django.utils.translation import ugettext_lazy as _, ugettext
from crispy_forms.helper import FormHelper
from crispy_forms import layout, bootstrap 
#from utils.fields import MultipleChoiceTreeField
from mptt.forms import TreeNodeMultipleChoiceField
from .models import *
#import pdb
from django.contrib.auth.forms import UserCreationForm


class ProfileCreationForm(forms.ModelForm):
	class Meta:
		model = UserInfo
		fields = ['userInfoPosition','userInfoCenter','userInfoArea','userInfoExtension',]
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
         					layout.Field("userInfoPosition"),                                                                                    
         					layout.Field("userInfoCenter"),                                                                                
         					css_class="col-md-6",                                                                                           
         				),                                                                                                                  
         				layout.Div(                                                                                                         
         					layout.Field("userInfoArea"),                                                                                    
         					layout.Field("userInfoExtension"),                                                                               
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

class ServiceRequestForm(forms.ModelForm):
 	class Meta:
 		model = Service
 		fields = ['serviceSeqCenter','servicePlatform','serviceRunSpecs','serviceFileExt','serviceAvailableService','serviceFile','serviceNotes']
 		field_classes = {
				'serviceAvailableService': TreeNodeMultipleChoiceField,
				}

 	def __init__(self,*args, **kwargs):
 		super(ServiceRequestForm, self).__init__(*args, **kwargs)
 		self.helper = FormHelper()
 		self.helper.form_action=""
 		self.helper.form_method="POST"
 		#pdb.set_trace()
 		
 		self.helper.layout = layout.Layout(
#  				layout.Div(
#  					layout.HTML(u"""<div class="panel-heading"><h3 class="panel-title">Researcher data</h3></div>"""), 
#  					layout.Div(
#  						layout.Div(
#  							layout.Field("serviceName"),
#  							layout.Field("servicePosition"),
#  							layout.Field("serviceCenter"),
#  							css_class="col-md-6",
#  						),
#  						layout.Div(
#  							layout.Field("serviceArea"),
#  							layout.Field("serviceExtension"),
#  							bootstrap.PrependedText("serviceEmail","@",css_class="input-block-level",placeholder="contact@example.com"),
#  							css_class="col-md-6",
#  						),
#  						css_class = 'row panel-body',
#  					),
#  					css_class = 'panel panel-default'
#  					),
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
