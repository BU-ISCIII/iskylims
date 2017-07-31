from django import forms

from .models import Document

class DocumentForm(forms.ModelForm):
    class Meta:
        model = Document
        fields = ('run_name', 'project_name','description', 'name','user_id','email','csv_file','convert')
    
class Docutres(forms.ModelForm):
    class Meta:
        model=Document
        fields = ('run_name','user_id', 'email', 'csv_file')

