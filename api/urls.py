from django.urls import path
from iSkyLIMS_drylab.api.views import *
from . import views

    
urlpatterns = [
#     path('jobslist',views.jobs_list, name='jobs_list'),
#     path('<int:state>/jobsliststate',views.jobs_list_state, name='jobs_list_state'),
     path('drylab/api/',views.jobs_list, name='jobs_list'),
     path('drylab/api/list/<int:state>',views.jobs_list_state, name='jobs_list_state'),
     path('drylab/api/update/<str:service>',views.api_update_job, name='api_update_job'),
     path('drylab/api/updatestate/<str:service>',views.api_update_state, name='api_update_state'),
]
