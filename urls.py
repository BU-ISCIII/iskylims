#from django.conf.urls import include , url
from django.urls import path
from django.conf import settings

from . import views
from django.conf.urls.static import static
from django.views.generic import ListView, DetailView



urlpatterns = [
    path('',views.index, name = 'index'),
    path('getSampleSheet/', views.get_sample_file, name='get_sample_file'),
    path('latest_run/',views.latest_run, name='latest_run'),
    path('searchNextSeq/', views.search_nextSeq, name='search_nextSeq'),
    path('searchNextProject/', views.search_nextProject, name ='search_nextProject'),
    path('searchNextSample/', views.search_nextSample, name = 'search_nextSample'),
    path('search_run=<int:run_id>/', views.search_run, name='search_run'),
    path('search_project=<int:project_id>/', views.search_project, name='search_project'),
    path('search_sample=<int:sample_id>/', views.search_sample, name= 'search_sample'),
    path('NextSeqStatsExperiment/', views.next_seq_stats_experiment, name ='next_seq_stats_experiment'),
    path('NextSeqStatsPerResearcher/',views.nextSeqStats_per_researcher, name='nextSeqStats_per_researcher'),
    path('NextSeqStatsPerTime/', views.nextSeqStats_per_time, name ='nextSeqStats_per_time'),
    path('NextSeqStatsLibrary/', views.nextSeqStats_per_library , name ='nextSeqStats_per_library'),
    path('AnualReport',views.anual_report, name='anual_report'),
    path('MonthlyReport', views.monthly_report, name='montly_report'),
    path('QuarterReport', views.quarter_report, name='quarter_report'),
    path('updateTables',views.update_tables, name='update_tables'),
    path('register_wetlab',views.register_wetlab, name='register_wetlab')
    
    #url(r'^documents/images_plot$', 'django.views.static.serve',{'document_root': settings.MEDIA_ROOT}),
    #url(r'^documents/(?P<path>.*)$', 'django.views.static.serve',{'document_root': settings.MEDIA_ROOT,}),
    #url(r'^mail',views.email, name='email'),

] 
#+ static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)

