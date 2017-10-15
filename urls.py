from django.conf.urls import include , url
from django.conf import settings

from . import views
from django.conf.urls.static import static
from django.views.generic import ListView, DetailView



urlpatterns = [
    url(r'^$', views.index, name='index'),
    url(r'^searchNextSeq', views.search_nextSeq, name='search_nextSeq'),
    url(r'^searchNextProject', views.search_nextProject, name ='search_nextProject'),
    url(r'^searchNextSample', views.search_nextSample, name = 'search_nextSample'),
    url(r'^getSampleSheet', views.get_sample_file, name='get_sample_file'),
    #url(r'^documents/', views.downloadFile, name='downloadFile'),
    url(r'^search_run=(?P<run_id>[0-9]+)/$', views.search_run, name='search_run'),
    url(r'^latest_run/',views.latest_run, name='latest_run'),
    url(r'^search_project=(?P<project_id>[0-9]+)$', views.search_project, name='search_project'),
    url(r'^search_sample=(?P<sample_id>[0-9]+)$', views.search_sample, name= 'search_sample'),
    url(r'^NextSeqStatistics', views.next_seq_statistics, name ='next_seq_statistics'),
    url(r'^NextSeqStatsPerResearcher',views.nextSeqStats_per_researcher, name='nextSeqStats_per_researcher'),
    url('^', include('django.contrib.auth.urls')),
    url('^graphic', views.test_graphic, name='test_graphic'),
    url('^NextSeqStatsPerTime', views.nextSeqStats_per_time, name ='nextSeqStats_per_time'),
    #url(r'^documents/images_plot$', 'django.views.static.serve',{'document_root': settings.MEDIA_ROOT}),
    #url(r'^documents/(?P<path>.*)$', 'django.views.static.serve',{'document_root': settings.MEDIA_ROOT,}),

] 
#+ static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)

