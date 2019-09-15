from django.db import models
from iSkyLIMS_core.models import Samples
import datetime
'''

class AdditionalParameters (models.Model):
    parameterName = models.CharField(max_length=255)
    parameterDescription = models.CharField(max_length= 400, null=True, blank=True)
    parameterOrder = models.IntegerField()
    parameterUsed = models.BooleanField()

    def __str__ (self) :
        return '%s' %(self.parameterName)



        
class SampleClinic (models.Manager):
    sampleCore = models.ForeignKey(
                Samples,
                on_delete = models.CASCADE, null = True)
    addPar = models.ForeignKey(
                AdditionalParameters,
                on_delete = models.CASCADE, null = True)
    generated_at = models.DateTimeField(auto_now_add=True)
    
'''


'''
    request
CÓDIGO
PRIORIDAD
NUMERO DE HISTORIA
SERVICIO
MÉDICO
SOSPECHA
COMENTARIOS
'''''
