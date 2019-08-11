from django.contrib import admin

class MoleculeTypeAdmin( admin.ModelAdmin):
    list_display = ('moleculeType',)

class ProtocolsAdmin(admin.ModelAdmin):
    list_display = ('type', 'name', 'description')

class ProtocolTypeAdmin( admin.ModelAdmin):
    list_display = ('molecule','protocol_type')


class ProtocolParametersAdmin (admin.ModelAdmin):
    list_display = ('protocol_id', 'parameterName', 'parameterOrder', 'parameterUsed', 'parameterMinValue', 'parameterMaxValue', 'parameterDescription')
