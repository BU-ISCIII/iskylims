from django.contrib import admin
from iSkyLIMS_wetlab.models import *




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

class RunProcessAdmin (admin.ModelAdmin):
    list_display = ('runName','sampleSheet','generatedat','run_date','runState','index_library','samples', 'useSpaceImgMb','useSpaceFastaMb','useSpaceOtherMb')

class ProjectsAdmin (admin.ModelAdmin):
    list_display= ('runprocess_id','user_id','LibraryKit_id','projectName','procState','libraryKit','baseSpaceFile','generatedat','project_run_date')
    #list_display= ('runprocess_id','projectName','procState','libraryKit','baseSpaceFile')

class LibraryKitAdmin(admin.ModelAdmin):
    list_display = ('libraryName','generatedat')


class SequencingPlatformAdmin(admin.ModelAdmin):
    list_display= ('platformName',)

class MachinesAdmin (admin.ModelAdmin):
    list_display=(
        'platformID','machineName','machineDescription','machineLocation','machineProvider',
        'machineSerialNumber','machineState','machineOperationStart','machineOperationEnd')



#admin.site.register(Document , AppAdmin)
admin.site.register(RunningParameters , RunningParametersAdmin)
admin.site.register(RunProcess , RunProcessAdmin)
admin.site.register(LibraryKit,LibraryKitAdmin)
admin.site.register(Projects, ProjectsAdmin)
admin.site.register(SequencingPlatform, SequencingPlatformAdmin)
admin.site.register(Machines, MachinesAdmin)
