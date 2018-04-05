#utils/fields.py
from django import forms

## For Multiple choice tree we need
# - Change default widget for checkbox select multiple 
# - Change label_from_instance to return object instead of choices value and choices text
class MultipleChoiceTreeField(forms.ModelMultipleChoiceField):
	widget = forms.CheckboxSelectMultiple
	def label_from_instance(self,obj):
		return obj
