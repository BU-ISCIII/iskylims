from django.shortcuts import render
from django.views import View
import core.models
import core.core_config
import core.utils.samples

class HandlingFragmentation(View):
    template_name = "wetlab/hamdling_fragmentation.html"
    package = __package__.split(".")[0]

    def get(self, request, *args, **kwargs):
        # collect the data to be shown in the form
        return render(request, self.template_name, {"fragmentation_data": self.get_fragmentation_data})

    def post(self, request, *args, **kwargs):


        pass


    def get_fragmentation_data(self):
        # get the data to be shown in the form
        pass