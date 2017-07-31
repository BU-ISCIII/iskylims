from django.conf.urls import include , url
from django.conf import settings

from . import views
from django.conf.urls.static import static
from django.views.generic import ListView, DetailView
from polls.models import Document


urlpatterns = [
    # ex: /polls/
    url(r'^$', views.index, name='index'),
    # ex: /polls/5/
    #url(r'^(?P<question_id>[0-9]+)/$', views.detail, name='detail'),
    # ex: /polls/5/results/
    #url(r'^(?P<question_id>[0-9]+)/results/$', views.results, name='results'),
    # ex: /polls/5/vote/
    #url(r'^(?P<question_id>[0-9]+)/vote/$', views.vote, name='vote'),
    #################### formulario
    url(r'^home$', views.home, name='home'),
    url(r'^searchNextSeq', views.search_nextSeq, name='search_nextSeq'),
    url(r'^getSampleSheet', views.get_sample_file, name='get_sample_file'),

    #url(r'^uploads/simple/$', views.simple_upload, name='simple_upload'),
    #url(r'^uploadsForm/$', views.model_form_upload, name='model_form_upload'),
    #url(r'^error_page',views.error_page, name='error_page'),
    #url(r'^result_form', views.result_form, name='result_form'),
    #url(r'^list_forms', ListView.as_view(queryset=Document.objects.all().order_by("-uploaded_at")[:5],
    #                    template_name='polls/listForms.html')),
    url(r'^documents/', views.downloadFile, name='downloadFile'),
    #url(r'^results_run_folder', views.results_run_folder, name='results_run_folder'),
    #url(r'^get_run_data', views.get_run_data, name='get_run_data'),
    #url(r'^nextSeq/listStatistics',views.nextSeq_listStatistics, name='nextSeq_listStatistics'),
    #url(r'^test_stats',views.test_stats, name='test_stats'),

    url('^', include('django.contrib.auth.urls')),

] 
#+ static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)


