from django import forms

from .models import Document

class DocumentForm(forms.ModelForm):
    class Meta:
        model = Document
        fields = ('description', 'name','user_id','email','csv_file','convert')

