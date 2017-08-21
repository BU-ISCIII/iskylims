from django.contrib import admin
from wetlab.models import *




class AppAdmin(admin.ModelAdmin):
    list_display = ('run_name','description','csv_file','name')
    #import pdb; pdb.set_trace()
    def file_link(self, obj):
        if obj.file:
            #import pdb; pdb.set_trace() 
            return "<a href='%s' download>Download</a>" % (obj.file.url,)
        else:
            return "No attachment"
    file_link.allow_tags = True
    file_link.short_description = 'File Download'

class RunningParametersAdmin(admin.ModelAdmin):
    list_display = ('runName_id','RunID','ExperimentName','RunStartDate')
    #import pdb; pdb.set_trace()

'''
class RawNextSeqStatsXml(admin.ModelAdmin):
    list_display =('document','project','rawYieldQ30','generated_at')
'''

class UserInfoAdmin (admin.ModelAdmin):
    list_display = ('userid','userFirstName','userLastName','userArea','userEmail')    
   
class RunProcessAdmin (admin.ModelAdmin):
    list_display = ('runName','sampleSheet','generatedBSFile', 'requestedCenter')


class ProjectsAdmin (admin.ModelAdmin):
    list_display= ('runprocess_id','user_id','projectName','procState','libraryKit','baseSpaceFile')
    #list_display= ('runprocess_id','projectName','procState','libraryKit','baseSpaceFile')
    
#admin.site.register(Document , AppAdmin)
admin.site.register(RunningParameters , RunningParametersAdmin)
admin.site.register(RunProcess , RunProcessAdmin)
#admin.site.register(RawNextSeqStatisticsXml , RawNextSeqStatsXml)
#admin.site.register(UserInfo,UserInfoAdmin)

admin.site.register(Projects, ProjectsAdmin)




