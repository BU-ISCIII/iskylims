from django.shortcuts import render
from django.views import View
import core.models # import Protocols, ProtocolParameters, MoleculeParameterValue, SampleProjects
import core.core_config
import core.utils.samples

class ExternalFilePreparationView(View):
    template_name = "wetlab/external_file_preparation.html"
    package = __package__.split(".")[0]

    def get(self, request, *args, **kwargs):
        # collect the data to be shown in the form
        return render(request, self.template_name, {"search_data": self.get_search_data})

    def post(self, request, *args, **kwargs):
        
        if request.POST.get("action") == "externalFilePreparation":
            # get the customer request to get the list of samples which where
            # showed in the excel form.
            selecting_samples = {}
            if "protocol" in request.POST:
                protocol_id = request.POST["protocol"]
                # find out the protocol type from the protocol id 
                if core.models.Protocols.objects.filter(pk__exact=protocol_id).exists():
                    p_obj = core.models.Protocols.objects.get(pk__exact=protocol_id)
                    # get the protocol type
                    p_type = p_obj.get_type()
                    selecting_samples["protocol"] = protocol_id
                    selecting_samples["protocol_type"] = p_type
                    if "DNA" in p_type or "RNA" in p_type:
                        # get the extraction object and then  samples
                        selecting_samples["samples"] = self.get_extraction_samples(p_obj)
                        if len(selecting_samples["samples"]) > 0:
                            selecting_samples["heading"] = core.core_config.HEADING_FOR_SELECTING_SAMPLE_EXTRACTION_ON_EXTRACTION
                            selecting_samples["heading_string"] = ",".join(selecting_samples["heading"])
                            return render(request, self.template_name, {"selecting_samples": selecting_samples})

        elif request.POST.get("action") == "selectSamples":
            heading = request.POST["heading"].split(",")
            heading.insert(-1, "extract_id")
            down_param_obj = list(
                core.models.ProtocolParameters.objects.filter(
                    protocol_id__pk=request.POST["protocolID"], parameter_download=True
                )
            )
            down_param_list = list(
                core.models.ProtocolParameters.objects.filter(
                    protocol_id__pk=request.POST["protocolID"], parameter_download=True
                ).values_list("parameter_name", flat=True)
            )
            project_list = []
            if "DNA" in request.POST["protocolType"] or "RNA" in request.POST["protocolType"]:
                map_sample_extraction = {}
                extractions, _ = core.utils.samples.get_selection_from_excel_data(
                    request.POST["selected_samples"], heading, "Select extraction", "extract_id"
                )
                if len(extractions) > 0:
                    data_for_download = {}
                    sample_ids = []
                    # split the extraction that have different sample project
                    for extraction in extractions:
                        sample_id = core.models.MoleculePreparation.objects.get(pk__exact=extraction).sample
                        map_sample_extraction.setdefault(sample_id, []).append(extraction)
                        sample_ids.append(sample_id)
                    # need to check if using the list 
                    # sample_ids = list(core.models.MoleculePreparation.objects.filter(pk_in=extractions).order_by('id').values_list("sample", flat=True)
                    sample_by_project = self.split_samples_by_project(sample_ids)

                    if len(sample_by_project) > 6:
                        error_message = core.core_config.ERROR_TOO_MANY_SAMPLE_PROJECTS
                        return render(request, self.template_name, {"error_message": error_message, "search_data": self.get_search_data()})
                    for project, samples in sample_by_project.items():
                        project_name = project.sample_project_name
                        project_list.append(project_name)
                        data_for_download[project_name] = {}
                        
                        proj_fields = list(core.models.SampleProjectsFields.objects.filter(sample_projects_id=project, sample_project_downloadable=True).values_list("sample_project_field_name", flat=True))
                        data_for_download[project_name]["heading"]  = ["sample name"] + proj_fields + ["extraction code"] +  down_param_list
                        for sample in samples:
                            s_name = sample.sample_name
                            s_proj_values = list(core.models.SampleProjectsFieldsValue.objects.filter(sample_id=sample, sample_project_field_id__sample_project_field_name__in=proj_fields).values_list("sample_project_field_value", flat=True))
                            for extract_id in map_sample_extraction[sample]:
                                extract_field_values = list(
                                    core.models.MoleculeParameterValue.objects.filter(
                                        molecule_id__pk=extract_id, molecule_parameter_id__in=down_param_obj
                                    ).values_list("parameter_value", flat=True)
                                )
                                extract_code_id = core.models.MoleculePreparation.objects.get(pk__exact=extract_id).molecule_code_id
                                data_for_download[project_name].setdefault("values", []).append([s_name] + s_proj_values + [extract_code_id] + extract_field_values)
                return render(request, self.template_name, {"data_for_download": data_for_download, "project_list": project_list})
        return render(request, self.template_name, {"search_data": self.get_search_data()})

    def get_extraction_samples(self, protocol_obj):
        # get the extraction object and then  samples
        if protocol_obj.extractions.filter(external_use=True).exists():
            extract_objs = protocol_obj.extractions.filter(external_use=True)
            return list(extract_objs.values_list("sample__sample_name", "molecule_code_id", "pk"))
        return []

    def get_search_data(self):
        search_data = {}
        if core.models.SampleProjects.objects.filter(apps_name__iexact=self.package).exists():
            search_data["s_projects"] = list(
                core.models.SampleProjects.objects.filter(
                    apps_name__iexact=self.package
                ).values_list("sample_project_name", "id")
            )

        if core.models.Protocols.objects.filter(type__apps_name__iexact=self.package).exists():
            search_data["protocols"] = list(
                core.models.Protocols.objects.filter(
                    type__apps_name__iexact=self.package
                ).values_list("name", "id")
            )

        return search_data

    def split_samples_by_project(self, sample_objs: list) -> dict:
        """Get the samples and split them by project

        Args:
            sample_ids (list): List of sample ids

        Returns:
            dict: Dictionary with the project obj as key and the sample  objs as value
        """
        sample_projects = {}
        for sample_obj in sample_objs:
            # sample_obj = core.models.Samples.objects.get(pk=sample_obj)
            s_project = sample_obj.sample_project
            sample_projects.setdefault(s_project, []).append(sample_obj)
        return sample_projects
