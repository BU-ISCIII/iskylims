# -*- coding: utf-8 -*-
from django.shortcuts import get_object_or_404, render, redirect
from django.http import HttpResponse
from django.template import loader
from .models import Question
from .models import Document
from .forms import DocumentForm
from .utils.sample_convertion import *


from django.core.files.storage import FileSystemStorage

####### Import libraries for static files
#from django.shortcuts import render_to_response
#from django.shortcuts import RequestContext

def index(request):
    latest_question_list = Question.objects.order_by('-pub_date')[:5]
    context = {'latest_question_list': latest_question_list}
    return render(request, 'polls/index.html', context)

def detail(request, question_id):
    return HttpResponse("You're looking at question %s." % question_id)

def results(request, question_id):
    response = "You're looking at the results of question %s."
    return HttpResponse(response % question_id)

def vote(request, question_id):
    return HttpResponse("You're voting on question %s." % question_id)
    
#############################
#### Form creation File System field
#############################
from django.shortcuts import render, redirect
from django.conf import settings
from django.core.files.storage import FileSystemStorage

from .models import Document
from .forms import DocumentForm

def home(request):
    form_sample_sheet = Document.objects.all()
    fs= FileSystemStorage()
    #filename= fs.
    #uploaded_file_url= 
    return render(request, 'polls/home.html', { 'documents': form_sample_sheet })

def error_page(request):
    
    return render(request, 'polls/error_page.html', {'error':error_text})

def downloadPdf(request):
    #from urllib.parse import urlparse
    #from os.path import splitext, basename
    #filename = object_name.file.name.split('/')[-1]
    #path_to_file = '/home/bioinfo/web_carlosIII/polls/documents/test.pdf'
    #f = open(path_to_file, 'r',encoding='utf-8')
    #myfile = File(f)
    file_tmp_name=str(request).split('/')[-1]
    temp_1=file_tmp_name.replace(">","")
    file_name=temp_1.replace("'","")
    #import pdb; pdb.set_trace()
    #disassembled = urlparse(str(request))
    #filename, file_ext = splitext(basename(disassembled.path))
    #file_name=str(filename + file_ext)
    #import pdb; pdb.set_trace()
    #to open the file as text but pdf is binary. Change 'r' to 'rb' in open
    with open(os.path.join(settings.MEDIA_ROOT, file_name), 'rb') as fh:
    #with open(os.path.join(settings.MEDIA_ROOT, 'test.pdf'), 'rb') as fh:
        response = HttpResponse(fh.read(), content_type="application/pdf")
        #### para ficheros en excel utilizar content_type='application/vnd.ms-excel'
        response['Content-Disposition'] = 'attachment; filename=donwload.pdf'
        return response


def result_form (request):
    doc= Document.object.all()[:5]
    return render(request , 'polls/result_form.html', {'form':doc})
    
def simple_upload(request):
    if request.method == 'POST' and request.FILES['myfile']:
        myfile = request.FILES['myfile']
        fs = FileSystemStorage()
        filename = fs.save(myfile.name,  myfile)
        uploaded_file_url = fs.url(filename)
        return render(request, 'polls/simple_upload.html', {
            'uploaded_file_url': uploaded_file_url
        })
    return render(request, 'polls/simple_upload.html')


def model_form_upload(request):
    if request.method == 'POST':
        form = DocumentForm(request.POST, request.FILES)
        if form.is_valid():
            doc_tmp=str(form.cleaned_data.get('csv_file'))
            #import pdb; pdb.set_trace() 
            tmp_name = re.search('(.*\.csv$)',doc_tmp)
            if (tmp_name):
                form.save()
                name_in_file=tmp_name.group(1)
                if (re.search('\.csv$',name_in_file)):
                    doc=str(name_in_file)
                    #doc=str('documents/'+ name_in_file)
                    #mapped_file=sample_sheet_map_basespace(doc)
                    mapped_file=True ## checking the file can upload and download
                    if (mapped_file):
                        return render( request, 'polls/result_form.html', {'text':['texto']})
                    else:
                        form.delete() ## using django-clenup app to delete the uploaded file
                        return redirect('polls/error_page.html', {'content':['Sample Sheet does not meet with the format']})
                else:
                    return redirect ('polls/error_page.html',  {'content':['invalid extension of Sample Sheet file' , 'Extension must be csv']})
            else:
                return render (request, 'polls/error_page.html', {'content':['invalid extension of Sample Sheet file', 'Extension must be csv']})
    else:
        form = DocumentForm()
    return render(request, 'polls/model_form_upload.html', {'form': form })

