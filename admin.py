from django.contrib import admin
from polls.models import *




class AppAdmin(admin.ModelAdmin):
    list_display = ('description','csv_file','name')
    #import pdb; pdb.set_trace()
    def file_link(self, obj):
        if obj.file:
            import pdb; pdb.set_trace() 
            #return "<a href='%s' download>Download</a>" % (obj.file.url,)
        else:
            return "No attachment"
    file_link.allow_tags = True
    file_link.short_description = 'File Download'

admin.site.register(Document , AppAdmin)

