from django.http import HttpResponse
import json

# Common base class for FC
class FusionCharts:

   #constructorOptions = {}
   constructorTemplate = """
     <script type="text/javascript">
         FusionCharts.ready(function () {
             new FusionCharts(__constructorOptions__);
         });
     </script>"""
   renderTemplate = """
     <script type="text/javascript">
         FusionCharts.ready(function () {
             FusionCharts("__chartId__").render();
         });
     </script>
   """
   # constructor
   def __init__(self, type, id, width, height, renderAt, dataFormat, dataSource): 
      self.constructorOptions = {}  
      self.constructorOptions['type'] = type
      self.constructorOptions['id'] = id
      self.constructorOptions['width'] = width
      self.constructorOptions['height'] = height
      self.constructorOptions['renderAt'] = renderAt
      self.constructorOptions['dataFormat'] = dataFormat
      #dataSource = unicode(dataSource, errors='replace')
      self.constructorOptions['dataSource'] = dataSource
   # render the chart created
   # It prints a script and calls the FusionCharts javascript render method of created chart   
   def render(self):
    self.readyJson = json.dumps(self.constructorOptions)
    self.readyJson = FusionCharts.constructorTemplate.replace('__constructorOptions__', self.readyJson)
    self.readyJson = self.readyJson + FusionCharts.renderTemplate.replace('__chartId__', self.constructorOptions['id'])
    self.readyJson = self.readyJson.replace('\\n', '')
    self.readyJson = self.readyJson.replace('\\t', '')

    if(self.constructorOptions['dataFormat'] == 'json'):
      self.readyJson = self.readyJson.replace('\\', '')
      self.readyJson = self.readyJson.replace('"{', "{")
      self.readyJson = self.readyJson.replace('}"', "}")
      
    return self.readyJson
 