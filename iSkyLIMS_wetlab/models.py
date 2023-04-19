import os
import re

from django.contrib.auth.models import User
from django.db import models
from django_utils.models import Center
from iSkyLIMS_core.models import (
    CommercialKits,
    MoleculePreparation,
    ProtocolParameters,
    Protocols,
    Samples,
    SequencerInLab,
    SequencingConfiguration,
    SequencingPlatform,
    UserLotCommercialKits,
)
from . import wetlab_config


class RunErrors(models.Model):
    error_code = models.CharField(max_length=10)
    error_text = models.CharField(max_length=255)

    class Meta:
        db_table = "wetlab_run_errors"

    def __str__(self):
        return "%s" % (self.error_text)


class RunStates(models.Model):
    run_state_name = models.CharField(max_length=50)

    class Meta:
        db_table = "wetlab_run_states"

    def __str__(self):
        return "%s" % (self.run_state_name)

    def get_run_state_name(self):
        return "%s" % (self.run_state_name)


class RunProcessManager(models.Manager):
    def create_new_run_from_crontab(self, run_data):
        run_state = RunStates.objects.get(run_state_name__exact="Recorded")
        new_run = self.create(state=run_state, run_name=run_data["experiment_name"])
        return new_run


class RunProcess(models.Model):
    used_sequencer = models.ForeignKey(
        SequencerInLab, on_delete=models.CASCADE, blank=True, null=True
    )
    run_error = models.ForeignKey(
        RunErrors, on_delete=models.CASCADE, null=True, blank=True
    )
    state_before_error = models.ForeignKey(
        RunStates, on_delete=models.CASCADE, null=True, blank=True
    )
    state = models.ForeignKey(
        RunStates,
        on_delete=models.CASCADE,
        related_name="state_of_run",
        null=True,
        blank=True,
    )
    center_requested_by = models.ForeignKey(
        Center, on_delete=models.CASCADE, null=True, blank=True
    )
    reagent_kit = models.ManyToManyField(UserLotCommercialKits, blank=True)
    run_name = models.CharField(max_length=45)
    sample_sheet = models.FileField(
        upload_to=wetlab_config.RUN_SAMPLE_SHEET_DIRECTORY, null=True, blank=True
    )
    generated_at = models.DateTimeField(auto_now_add=True)
    run_date = models.DateField(auto_now=False, null=True, blank=True)
    run_finish_date = models.DateTimeField(auto_now=False, null=True, blank=True)
    bcl2fastq_finish_date = models.DateTimeField(auto_now=False, null=True, blank=True)
    run_completed_date = models.DateTimeField(auto_now=False, null=True, blank=True)
    run_forced_continue_on_error = models.BooleanField(
        null=True, blank=True, default=False
    )
    index_library = models.CharField(max_length=85, null=True, blank=True)
    samples = models.CharField(max_length=45, blank=True)
    use_space_img_mb = models.CharField(max_length=10, blank=True)
    use_space_fasta_mb = models.CharField(max_length=10, blank=True)
    use_space_other_mb = models.CharField(max_length=10, blank=True)

    class Meta:
        db_table = "wetlab_run_process"

    def __str__(self):
        return "%s" % (self.run_name)

    def get_run_id(self):
        return "%s" % (self.id)

    def get_recorded_date_no_format(self):
        return self.generated_at

    def get_run_date_no_format(self):
        return self.run_date

    def get_run_completion_date_no_format(self):
        return self.run_completed_date

    def get_run_finish_date_no_format(self):
        return self.run_finish_date

    def get_run_finish_date(self):
        if self.run_finish_date is None:
            run_finish_date = "Run NOT finished"
        else:
            run_finish_date = self.run_finish_date.strftime("%B %d, %Y")
        return run_finish_date

    def get_run_date(self):
        if self.run_date is None:
            rundate = "Run NOT started"
        else:
            rundate = self.run_date.strftime("%B %d, %Y")
        return rundate

    def get_run_year(self):
        return "%s" % (self.run_date.timetuple().tm_year)

    def get_run_generated_date_no_format(self):
        return self.generated_at

    def get_run_generated_date(self):
        return self.generated_at.strftime("%B %d, %Y")

    def get_error_text(self):
        return "%s" % (self.run_error)

    def get_forced_continue_on_error(self):
        return self.run_forced_continue_on_error

    def get_state(self):
        return "%s" % (self.state.get_run_state_name())

    def get_state_before_error(self):
        return "%s" % (self.state_before_error)

    def get_index_library(self):
        return "%s" % (self.index_library)

    def get_info_process(self):
        generated_date = self.generated_at.strftime("%I:%M%p on %B %d, %Y")

        requested_center = str(self.center_requested_by)
        if self.run_date is None:
            run_date = "Run NOT started"
        else:
            run_date = self.run_date.strftime("%B %d, %Y")

        if self.run_finish_date is None:
            finish_date = "Run multiplexation is not completed"
        else:
            finish_date = self.run_finish_date.strftime("%I:%M%p on %B %d, %Y")

        if self.bcl2fastq_finish_date is None:
            bcl2fastq_date = "bcl2fastq process is not completed"
        else:
            bcl2fastq_date = self.bcl2fastq_finish_date.strftime("%I:%M%p on %B %d, %Y")

        if self.run_completed_date is None:
            completed_date = "Run process is not completed"
        else:
            completed_date = self.run_completed_date.strftime("%I:%M%p on %B %d, %Y")

        return "%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s" % (
            self.run_name,
            self.state,
            requested_center,
            generated_date,
            run_date,
            completed_date,
            bcl2fastq_date,
            finish_date,
            self.use_space_img_mb,
            self.use_space_fasta_mb,
            self.use_space_other_mb,
        )

    def get_run_name(self):
        if self.run_name is None:
            return "Not defined"
        return "%s" % (self.run_name)

    def get_disk_space_utilization(self):
        image_size_str = "0" if self.use_space_fasta_mb == "" else self.use_space_fasta_mb
        data_size_str = (
            "0" if self.use_space_fasta_mb == "" else self.use_space_fasta_mb
        )
        other_size_str = (
            "0" if self.use_space_other_mb == "" else self.use_space_other_mb
        )
        image_size = int(image_size_str.replace(",", ""))
        data_size = int(data_size_str.replace(",", ""))
        other_size = int(other_size_str.replace(",", ""))
        total_size = image_size + data_size + other_size
        return "%s" % (total_size)

    def get_run_used_sequencer_name(self):
        if self.used_sequencer is not None:
            return "%s" % (self.used_sequencer.get_sequencer_name())
        else:
            return ""

    def get_run_used_sequencer(self):
        return "%s" % (self.used_sequencer)

    def get_run_platform(self):
        if self.used_sequencer is not None:
            return "%s" % (self.used_sequencer.get_sequencing_platform_name())
        else:
            return "None"

    def get_number_of_samples_in_run(self):
        return "%s" % (self.samples)

    def get_sample_file(self):
        return "%s" % (self.sample_sheet)

    def update_index_library(self, index_library_name):
        self.index_library = index_library_name
        self.save()
        return True

    def update_sample_sheet(self, full_path, file_name):
        self.sample_sheet.save(file_name, open(full_path, "r"), save=True)
        return self

    def set_used_space(self, disk_utilization):
        self.use_space_fasta_mb = disk_utilization["useSpaceFastaMb"]
        self.use_space_img_mb = disk_utilization["useSpaceImgMb"]
        self.use_space_other_mb = disk_utilization["useSpaceOtherMb"]
        self.save()
        return True

    def set_run_state(self, new_state):
        if RunStates.objects.filter(run_state_name__exact=new_state).exists():
            self.state = RunStates.objects.get(run_state_name__exact=new_state)
            self.save()
        else:
            return False
        return True

    def set_run_date(self, run_date):
        self.run_date = run_date
        self.save()
        return self

    def set_run_bcl2fastq_finished_date(self, bcl2fastq_finish_date):
        self.bcl2fastq_finish_date = bcl2fastq_finish_date
        self.save()
        return True

    def set_run_error_code(self, error_code):
        if RunErrors.objects.filter(error_code__exact=error_code).exists():
            self.run_error = RunErrors.objects.get(error_code__exact=error_code)
        else:
            self.run_error = RunErrors.objects.get(error_text__exact="Undefined")
        self.state_before_error = self.state
        self.state = RunStates.objects.get(run_state_name__exact="Error")
        self.save()
        return True

    def set_run_completion_date(self, completion_date):
        self.run_completed_date = completion_date
        self.save()
        return self

    def set_run_finish_date(self, finish_date):
        self.run_finish_date = finish_date
        self.save()
        return self

    def set_run_sample_sheet(self, sample_sheet):
        self.sample_sheet = sample_sheet
        self.save()
        return self

    def set_used_sequencer(self, sequencer):
        self.used_sequencer = sequencer
        self.save()
        return self

    def set_forced_continue_on_error(self):
        self.run_forced_continue_on_error = True
        self.save()
        return self

    objects = RunProcessManager()


class LibraryKit(models.Model):
    library_name = models.CharField(max_length=125)
    generated_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        db_table = "wetlab_library_kit"

    def get_bs_lib_name(self):
        return "%s" % (self.library_name)


class ProjectsManager(models.Manager):
    def create_new_empty_project(self, project_data):
        new_project = self.create(
            user_id=project_data["user_id"], project_name=project_data["projectName"]
        )
        return new_project

    def create_new_project(self, project_data):
        new_project = self.create(
            user_id=project_data["user_id"],
            project_name=project_data["projectName"],
            library_kit=project_data["libraryKit"],
            base_space_file=project_data["baseSpaceFile"],
            Base_space_library=project_data["BaseSpaceLibrary"],
        )
        return new_project


class Projects(models.Model):
    user_id = models.ForeignKey(User, on_delete=models.CASCADE, null=True, blank=True)
    library_kit_id = models.ForeignKey(
        LibraryKit, on_delete=models.CASCADE, null=True, blank=True
    )
    run_process = models.ManyToManyField(RunProcess)

    base_space_library = models.CharField(max_length=45, null=True, blank=True)
    project_name = models.CharField(max_length=45)
    library_kit = models.CharField(max_length=125, null=True, blank=True)
    base_space_file = models.CharField(max_length=255, null=True, blank=True)
    generated_at = models.DateTimeField(auto_now_add=True)
    project_run_date = models.DateField(auto_now=False, null=True, blank=True)

    class Meta:
        db_table = "wetlab_projects"

    def __str__(self):
        return "%s" % (self.project_name)

    def get_base_space_file(self):
        return "%s" % (self.base_space_file)

    def get_state(self):
        return "%s" % (self.runprocess_id.state.runStateName)

    def get_project_dates(self):
        if self.project_run_date is None:
            projectdate = "Run NOT started"
        else:
            projectdate = self.project_run_date.strftime("%B %d, %Y")
        p_info = []
        p_info.append(self.generated_at.strftime("%B %d, %Y"))
        p_info.append(projectdate)
        return p_info

    def get_p_info_change_library(self):
        run_name = self.runprocess_id.runName
        user_name = self.user_id.username
        if self.project_run_date is None:
            project_date = "Run NOT started"
        else:
            project_date = self.project_run_date.strftime("%B %d, %Y")

        return "%s;%s;%s;%s;%s" % (
            run_name,
            self.project_name,
            project_date,
            user_name,
            self.library_kit,
        )

    def get_run_id(self):
        return "%s" % (self.runprocess_id.get_run_id())

    def get_run_obj(self):
        return self.runprocess_id

    def get_user_name(self):
        user_name = self.user_id.username
        return "%s" % (user_name)

    def get_user_id(self):
        return "%s" % (self.user_id.pk)

    def get_project_name(self):
        return "%s" % (self.project_name)

    def get_index_library_name(self):
        return "%s" % (self.library_kit)

    def get_date(self):
        if self.project_run_date is None:
            project_date = "No Date"
        else:
            project_date = self.project_run_date.strftime("%B %d, %Y")
        return "%s" % project_date

    def get_project_id(self):
        return "%s" % (self.id)

    def add_run(self, run_obj):
        self.run_process.add(run_obj)
        return self

    def set_project_run_date(self, date):
        self.project_run_date = date
        self.save()
        return self

    objects = ProjectsManager()


class RunningParametersManager(models.Manager):
    def create_running_parameters(self, running_data, run_object):
        if running_data["PlannedRead1Cycles"] == "":
            running_data["PlannedRead1Cycles"] = 0
        if running_data["PlannedRead2Cycles"] == "":
            running_data["PlannedRead2Cycles"] = 0
        if running_data["PlannedIndex1ReadCycles"] == "":
            running_data["PlannedIndex1ReadCycles"] = 0
        if running_data["PlannedIndex2ReadCycles"] == "":
            running_data["PlannedIndex2ReadCycles"] = 0

        if running_data["LibraryID"] is None:
            running_data["LibraryID"] = " "
        if running_data["AnalysisWorkflowType"] is None:
            running_data["AnalysisWorkflowType"] = " "

        running_parameters = self.create(
            run_name_id=run_object,
            run_iD=running_data["RunID"],
            experiment_name=running_data["ExperimentName"],
            rta_version=running_data["RTAVersion"],
            system_suite_version=running_data["SystemSuiteVersion"],
            library_id=running_data["LibraryID"],
            chemistry=running_data["Chemistry"],
            run_start_date=running_data["RunStartDate"],
            analysis_workflow_type=running_data["AnalysisWorkflowType"],
            run_management_type=running_data["RunManagementType"],
            planned_read1_cycles=running_data["PlannedRead1Cycles"],
            planned_read2_cycles=running_data["PlannedRead2Cycles"],
            planned_index1_read_cycles=running_data["PlannedIndex1ReadCycles"],
            planned_index2_read_cycles=running_data["PlannedIndex2ReadCycles"],
            application_version=running_data["ApplicationVersion"],
            num_tiles_per_swath=running_data["NumTilesPerSwath"],
            image_channel=running_data["ImageChannel"],
            flowcell=running_data["Flowcell"],
            image_dimensions=running_data["ImageDimensions"],
            flowcell_layout=running_data["FlowcellLayout"],
        )

        return running_parameters


class RunningParameters(models.Model):
    run_name_id = models.OneToOneField(
        RunProcess, on_delete=models.CASCADE, primary_key=True
    )
    run_id = models.CharField(max_length=255)
    experiment_name = models.CharField(max_length=255)
    rta_version = models.CharField(max_length=255, null=True, blank=True)
    system_suite_version = models.CharField(max_length=255, null=True, blank=True)
    library_id = models.CharField(max_length=255, null=True, blank=True)
    chemistry = models.CharField(max_length=255, null=True, blank=True)
    run_start_date = models.CharField(max_length=255, null=True, blank=True)
    analysis_workflow_type = models.CharField(max_length=255, null=True, blank=True)
    run_management_type = models.CharField(max_length=255, null=True, blank=True)
    planned_read1_cycles = models.CharField(max_length=255, null=True, blank=True)
    planned_read2_cycles = models.CharField(max_length=255, null=True, blank=True)
    planned_index1_read_cycles = models.CharField(max_length=255, null=True, blank=True)
    planned_index2_read_cycles = models.CharField(max_length=255, null=True, blank=True)
    application_version = models.CharField(max_length=255, null=True, blank=True)
    num_tiles_per_swath = models.CharField(max_length=255, null=True, blank=True)
    image_channel = models.CharField(max_length=255, null=True, blank=True)
    flowcell = models.CharField(max_length=255, null=True, blank=True)
    image_dimensions = models.CharField(max_length=255, null=True, blank=True)
    flowcell_layout = models.CharField(max_length=255, null=True, blank=True)

    class Meta:
        db_table = "wetlab_running_parameters"

    def __str__(self):
        return "%s" % (self.run_id)
        # return '%s' %(self.runName_id)

    def get_run_chemistry(self):
        return "%s" % (self.chemistry)

    def get_run_parameters_info(self):
        # str_run_start_date=self.RunStartDate.strftime("%I:%M%p on %B %d, %Y")
        run_parameters_data = []
        if self.image_channel is None:
            img_channel = "None"
        else:
            img_channel = self.image_channel.strip("[").strip("]").replace("'", "")

        if self.image_dimensions is None:
            image_dimensions = ["None"]
        else:
            image_dimensions = (
                self.image_dimensions.strip("{").strip("}").replace("'", "").split(",")
            )

        flowcell_layout = (
            self.flowcell_layout.strip("{").strip("}").replace("'", "").split(",")
        )

        run_parameters_data.append(self.run_id)
        run_parameters_data.append(self.experiment_name)
        run_parameters_data.append(self.rta_version)
        run_parameters_data.append(self.system_suite_version)
        run_parameters_data.append(self.library_id)
        run_parameters_data.append(self.chemistry)
        run_parameters_data.append(self.run_start_date)
        run_parameters_data.append(self.analysis_workflow_type)
        run_parameters_data.append(self.run_management_type)
        run_parameters_data.append(self.planned_read1_cycles)
        run_parameters_data.append(self.planned_read2_cycles)
        run_parameters_data.append(self.planned_index1_read_cycles)
        run_parameters_data.append(self.planned_index2_read_cycles)
        run_parameters_data.append(self.application_version)
        run_parameters_data.append(self.num_tiles_per_swath)
        run_parameters_data.append(img_channel)
        run_parameters_data.append(self.flowcell)
        run_parameters_data.append(image_dimensions)
        run_parameters_data.append(flowcell_layout)

        return run_parameters_data

    def get_number_of_lanes(self):
        match_flowcell = re.match(
            r".*LaneCount.*'(\d+)'.*SurfaceCount.*", self.flowcell_layout
        )
        return match_flowcell.group(1)

    def get_number_of_reads(self):
        count = 0
        if self.planned_read1_cycles != "0":
            count += 1
        if self.planned_read2_cycles != "0":
            count += 1
        if self.planned_index1_read_cycles != "0":
            count += 1
        if self.planned_index2_read_cycles != "0":
            count += 1
        return count

    def get_number_of_cycles(self):
        cycles_number = (
            int(self.planned_read1_cycles)
            + int(self.planned_read2_cycles)
            + int(self.planned_index1_read_cycles)
            + int(self.planned_index2_read_cycles)
        )
        return cycles_number

    def get_run_folder(self):
        return "%s" % (self.run_id)

    objects = RunningParametersManager()


class StatsRunSummaryManager(models.Manager):
    def create_stats_run_summary(self, stats_run_summary, run_process_obj):
        s_run_summary = self.create(
            runprocess_id=run_process_obj,
            level=stats_run_summary["level"],
            yield_total=stats_run_summary["yieldTotal"],
            projected_total_yield=stats_run_summary["projectedTotalYield"],
            aligned=stats_run_summary["aligned"],
            error_rate=stats_run_summary["errorRate"],
            intensity_cycle=stats_run_summary["intensityCycle"],
            bigger_q30=stats_run_summary["biggerQ30"],
            stats_summary_run_date=run_process_obj.run_date,
        )
        return s_run_summary


class StatsRunSummary(models.Model):
    runprocess_id = models.ForeignKey(RunProcess, on_delete=models.CASCADE)
    level = models.CharField(max_length=20)
    yield_total = models.CharField(max_length=10)
    projected_total_yield = models.CharField(max_length=10)
    aligned = models.CharField(max_length=10)
    error_rate = models.CharField(max_length=10)
    intensity_cycle = models.CharField(max_length=10)
    bigger_q30 = models.CharField(max_length=10)
    generated_at = models.DateTimeField(auto_now_add=True)
    stats_summary_run_date = models.DateField(auto_now=False, null=True)

    class Meta:
        db_table = "wetlab_stats_run_summary"

    def __str__(self):
        return "%s" % (self.level)

    def get_bin_run_summary(self):
        return "%s;%s;%s;%s;%s;%s" % (
            self.yield_total,
            self.projected_total_yield,
            self.aligned,
            self.error_rate,
            self.intensity_cycle,
            self.bigger_q30,
        )

    objects = StatsRunSummaryManager()


class StatsRunReadManager(models.Manager):
    def create_stats_run_read(self, stats_run_read, experiment_name):
        run_process = RunProcess.objects.get(run_name__exact=experiment_name)
        s_run_read = self.create(
            runprocess_id=run_process,
            read=stats_run_read["read"],
            lane=stats_run_read["lane"],
            tiles=stats_run_read["tiles"],
            density=stats_run_read["density"],
            cluster_pf=stats_run_read["cluster_PF"],
            phas_prephas=stats_run_read["phas_prephas"],
            reads=stats_run_read["reads"],
            reads_pf=stats_run_read["reads_PF"],
            q30=stats_run_read["q30"],
            yields=stats_run_read["yields"],
            cycles_err_rated=stats_run_read["cyclesErrRated"],
            aligned=stats_run_read["aligned"],
            error_rate=stats_run_read["errorRate"],
            error_rate_35=stats_run_read["errorRate35"],
            error_rate_50=stats_run_read["errorRate50"],
            error_rate_75=stats_run_read["errorRate75"],
            error_rate_100=stats_run_read["errorRate100"],
            intensity_cycle=stats_run_read["intensityCycle"],
            stats_read_run_date=run_process.run_date,
        )

        return s_run_read


class StatsRunRead(models.Model):
    runprocess_id = models.ForeignKey(RunProcess, on_delete=models.CASCADE)
    read = models.CharField(max_length=10)
    lane = models.CharField(max_length=10)
    tiles = models.CharField(max_length=10)
    density = models.CharField(max_length=40)
    cluster_pf = models.CharField(max_length=40)
    phas_prephas = models.CharField(max_length=40)
    reads = models.CharField(max_length=40)
    reads_pf = models.CharField(max_length=40)
    q30 = models.CharField(max_length=40)
    yields = models.CharField(max_length=40)
    cycles_err_rated = models.CharField(max_length=40)
    aligned = models.CharField(max_length=40)
    error_rate = models.CharField(max_length=40)
    error_rate_35 = models.CharField(max_length=40)
    error_rate_50 = models.CharField(max_length=40)
    error_rate_75 = models.CharField(max_length=40)
    error_rate_100 = models.CharField(max_length=40)
    intensity_cycle = models.CharField(max_length=40)
    generated_at = models.DateTimeField(auto_now_add=True)
    stats_read_run_date = models.DateField(auto_now=False, null=True)

    class Meta:
        db_table = "wetlab_stats_run_read"

    def __str__(self):
        return "%s" % (self.read)

    def get_bin_run_read(self):
        return "%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s" % (
            self.lane,
            self.tiles,
            self.density,
            self.cluster_pf,
            self.phas_prephas,
            self.reads,
            self.reads_pf,
            self.q30,
            self.yields,
            self.cycles_err_rated,
            self.aligned,
            self.error_rate,
            self.error_rate_35,
            self.error_rate_50,
            self.error_rate_75,
            self.error_rate_100,
            self.intensity_cycle,
        )

    objects = StatsRunReadManager()


class RawDemuxStatsManager(models.Manager):
    def create_stats_run_read(self, raw_demux_stats, run_object_name):
        raw_stats = self.create(
            runprocess_id=run_object_name,
            project_id=raw_demux_stats["project_id"],
            default_all=raw_demux_stats["defaultAll"],
            raw_yield=raw_demux_stats["rawYield"],
            raw_yield_q30=raw_demux_stats["rawYieldQ30"],
            raw_quality=raw_demux_stats["rawQuality"],
            pf_yield=raw_demux_stats["PF_Yield"],
            pf_yield_q30=raw_demux_stats["PF_YieldQ30"],
            pf_quality_score=raw_demux_stats["PF_QualityScore"],
        )

        return raw_stats


class RawDemuxStats(models.Model):
    runprocess_id = models.ForeignKey(
        RunProcess,
        on_delete=models.CASCADE,
    )
    project_id = models.ForeignKey(Projects, on_delete=models.CASCADE, null=True)
    default_all = models.CharField(max_length=40, null=True)
    raw_yield = models.CharField(max_length=255)
    raw_yield_q30 = models.CharField(max_length=255)
    raw_quality = models.CharField(max_length=255)
    pf_yield = models.CharField(max_length=255)
    pf_yield_q30 = models.CharField(max_length=255)
    pf_quality_score = models.CharField(max_length=255)
    generated_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        db_table = "wetlab_raw_demux_stats"

    def __str__(self):
        return "%s" % (self.runprocess_id)

    objects = RawDemuxStatsManager()


class RawTopUnknowBarcodesManager(models.Manager):
    def create_unknow_barcode(self, unknow_barcode):
        unknow_barcode = self.create(
            runprocess_id=unknow_barcode["runprocess_id"],
            lane_number=unknow_barcode["lane_number"],
            top_number=unknow_barcode["top_number"],
            count=unknow_barcode["count"],
            sequence=unknow_barcode["sequence"],
        )

        return unknow_barcode


class RawTopUnknowBarcodes(models.Model):
    runprocess_id = models.ForeignKey(RunProcess, on_delete=models.CASCADE)
    lane_number = models.CharField(max_length=4)
    top_number = models.CharField(max_length=4)
    count = models.CharField(max_length=40)
    sequence = models.CharField(max_length=40)
    generated_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        db_table = "wetlab_raw_top_unknown_barcodes"

    def __str__(self):
        return "%s" % (self.lane_number)

    def get_unknow_barcodes(self):
        return "%s;%s" % (self.count, self.sequence)

    objects = RawTopUnknowBarcodesManager()


class StatsFlSummaryManager(models.Manager):
    def create_fl_summary(self, stats_fl_summary):
        fl_summary = self.create(
            runprocess_id=stats_fl_summary["runprocess_id"],
            project_id=stats_fl_summary["project_id"],
            default_all=stats_fl_summary["defaultAll"],
            flow_raw_cluster=stats_fl_summary["flowRawCluster"],
            flow_pf_cluster=stats_fl_summary["flowPfCluster"],
            flow_yield_mb=stats_fl_summary["flowYieldMb"],
            sample_number=stats_fl_summary["sampleNumber"],
        )

        return fl_summary


class StatsFlSummary(models.Model):
    runprocess_id = models.ForeignKey(RunProcess, on_delete=models.CASCADE)
    project_id = models.ForeignKey(
        Projects, on_delete=models.CASCADE, null=True, blank=True
    )
    default_all = models.CharField(max_length=40, null=True)
    flow_raw_cluster = models.CharField(max_length=40)
    flow_pf_cluster = models.CharField(max_length=40)
    flow_yield_mb = models.CharField(max_length=40)
    sample_number = models.CharField(max_length=40)
    generated_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        db_table = "wetlab_stats_fl_summary"

    def get_fl_summary(self):
        data = []
        data.append(self.flow_raw_cluster)
        data.append(self.flow_pf_cluster)
        data.append(self.flow_yield_mb)
        data.append(self.sample_number)
        return data

    def get_sample_number(self):
        return "%s" % (self.sample_number)

    objects = StatsFlSummaryManager()


class StatsLaneSummaryManager(models.Manager):
    def create_lane_summary(self, l_summary):
        lane_summary = self.create(
            runprocess_id=l_summary["runprocess_id"],
            project_id=l_summary["project_id"],
            default_all=l_summary["defaultAll"],
            lane=l_summary["lane"],
            pf_cluster=l_summary["pfCluster"],
            percent_lane=l_summary["percentLane"],
            perfect_barcode=l_summary["perfectBarcode"],
            one_mismatch=l_summary["oneMismatch"],
            yield_mb=l_summary["yieldMb"],
            bigger_q30=l_summary["biggerQ30"],
            mean_quality=l_summary["meanQuality"],
        )

        return lane_summary


class StatsLaneSummary(models.Model):
    runprocess_id = models.ForeignKey(RunProcess, on_delete=models.CASCADE)
    project_id = models.ForeignKey(
        Projects, on_delete=models.CASCADE, null=True, blank=True
    )
    default_all = models.CharField(max_length=40, null=True)
    lane = models.CharField(max_length=10)
    pf_cluster = models.CharField(max_length=64)
    percent_lane = models.CharField(max_length=64)
    perfect_barcode = models.CharField(max_length=64)
    one_mismatch = models.CharField(max_length=64)
    yield_mb = models.CharField(max_length=64)
    bigger_q30 = models.CharField(max_length=64)
    mean_quality = models.CharField(max_length=64)
    generated_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        db_table = "wetlab_stats_lane_summary"

    def __str__(self):
        return "%s" % (self.runprocess_id)

    def get_flow_cell_summary(self):
        return "%s" % (self.runprocess_id)

    def get_lane_summary(self):
        data = []
        data.append(self.lane)
        data.append(self.pf_cluster)
        data.append(self.percent_lane)
        data.append(self.perfect_barcode)
        data.append(self.one_mismatch)
        data.append(self.yield_mb)
        data.append(self.bigger_q30)
        data.append(self.mean_quality)
        return data

    def get_stats_info(self):
        stats_info = []
        stats_info.append(self.bigger_q30)
        stats_info.append(self.mean_quality)
        stats_info.append(self.yield_mb)
        stats_info.append(self.pf_cluster)
        return stats_info

    objects = StatsLaneSummaryManager()


class GraphicsStatsManager(models.Manager):
    def create_graphic_run_metrics(self, run_process_obj, run_folder):
        graphic_metric = self.create(
            runprocess_id=run_process_obj,
            folder_run_graphic=run_folder,
            cluster_count_graph="ClusterCount-by-lane.png",
            flowcell_graph="flowcell-Intensity.png",
            intensity_by_cycle_graph="Intensity-by-cycle.png",
            heatmap_graph="q-heat-map.png",
            histogram_graph="q-histogram.png",
            sample_qc_graph="sample-qc.png",
        )

        return graphic_metric


class GraphicsStats(models.Model):
    runprocess_id = models.ForeignKey(RunProcess, on_delete=models.CASCADE)
    folder_run_graphic = models.CharField(max_length=255)
    cluster_count_graph = models.CharField(max_length=255)
    flowcell_graph = models.CharField(max_length=255)
    intensity_by_cycle_graph = models.CharField(max_length=255)
    heatmap_graph = models.CharField(max_length=255)
    histogram_graph = models.CharField(max_length=255)
    sample_qc_graph = models.CharField(max_length=255)
    # Fix 24/07/2018: null=True added at definition of "generated_at" to allow json load
    generated_at = models.DateTimeField(auto_now_add=True, null=True)

    class Meta:
        db_table = "wetlab_graphics_stats"

    def __str__(self):
        return "%s" % (self.folder_run_graphic)

    def get_graphics(self):
        return "%s;%s;%s;%s;%s;%s" % (
            self.cluster_count_graph,
            self.flowcell_graph,
            self.intensity_by_cycle_graph,
            self.heatmap_graph,
            self.histogram_graph,
            self.sample_qc_graph,
        )

    def get_folder_graphic(self):
        return "%s" % (self.folder_run_graphic)

    objects = GraphicsStatsManager()


class SamplesInProjectManager(models.Manager):
    def create_sample_project(self, sample_p_data):
        try:
            user_obj = User.objects.get(username__exact=sample_p_data["user_id"])
        except Exception:
            user_obj = None
            sample_project = self.create(
                project_id=sample_p_data["project_id"],
                sample_name=sample_p_data["sampleName"],
                barcode_name=sample_p_data["barcodeName"],
                pf_clusters=sample_p_data["pfClusters"],
                percent_in_project=sample_p_data["percentInProject"],
                yield_mb=sample_p_data["yieldMb"],
                quality_q30=sample_p_data["qualityQ30"],
                mean_quality=sample_p_data["meanQuality"],
                run_process_id=sample_p_data["runProcess_id"],
                user_id=user_obj,
            )
        return sample_project


class SamplesInProject(models.Model):
    project_id = models.ForeignKey(Projects, on_delete=models.CASCADE)
    run_process_id = models.ForeignKey(
        RunProcess, on_delete=models.CASCADE, null=True, blank=True
    )
    user_id = models.ForeignKey(User, on_delete=models.CASCADE, null=True, blank=True)
    sample_in_core = models.ForeignKey(
        Samples, on_delete=models.CASCADE, null=True, blank=True
    )
    sample_name = models.CharField(max_length=255)
    barcode_name = models.CharField(max_length=255)
    pf_clusters = models.CharField(max_length=55)
    percent_in_project = models.CharField(max_length=25)
    yield_mb = models.CharField(max_length=55)
    quality_q30 = models.CharField(max_length=55)
    mean_quality = models.CharField(max_length=55)
    generated_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        db_table = "wetlab_samples_in_project"

    def __str__(self):
        return "%s" % (self.sample_name)

    def get_basic_info(self):
        sample_info = []
        sample_info.append(self.pk)
        sample_info.append(self.sample_name)
        if self.project_id is not None:
            sample_info.append(self.project_id.get_project_name())
        else:
            sample_info.append("None")
        if self.run_process_id is not None:
            sample_info.append(self.run_process_id.get_run_name())
        else:
            sample_info.append("None")
        sample_info.append(self.generated_at.strftime("%I:%M%p on %B %d, %Y"))
        return sample_info

    def get_sample_id(self):
        return "%s" % (self.id)

    def get_sample_information(self):
        data = []
        data.append(self.sample_name)
        data.append(self.barcode_name)
        data.append(self.pf_clusters)
        data.append(self.percent_in_project)
        data.append(self.yield_mb)
        data.append(self.quality_q30)
        data.append(self.mean_quality)
        return data

    def get_sample_information_with_project_run(self):
        data = []
        data.append(self.sample_name)
        if self.project_id is not None:
            data.append(self.project_id.get_project_name())
        else:
            data.append("None")
        if self.run_process_id is not None:
            data.append(self.run_process_id.get_run_name())
        else:
            data.append("None")
        data.append(self.generated_at.strftime("%B %d, %Y"))
        data.append(self.pf_clusters)
        data.append(self.yield_mb)
        data.append(self.quality_q30)
        data.append(self.mean_quality)
        return data

    def get_sample_name(self):
        return "%s" % (self.sample_name)

    def get_project_id(self):
        return "%s" % (self.project_id.pk)

    def get_project_name(self):
        if self.project_id is not None:
            project_name = self.project_id.get_project_name()
        else:
            project_name = "None"
        return "%s" % (project_name)

    def get_run_name(self):
        if self.run_process_id is not None:
            return "%s" % (self.run_process_id.get_run_name())
        else:
            return "None"

    def get_sequencer_used(self):
        return "%s" % (self.run_process_id.get_run_used_sequencer_name())

    def get_run_obj(self):
        return self.run_process_id

    def get_run_id(self):
        return self.run_process_id.pk

    def get_quality_sample(self):
        return "%s" % (self.quality_q30)

    def get_number_of_samples_in_run(self):
        return "%s" % (self.run_process_id.get_number_of_samples_in_run())

    def get_sample_user_name(self):
        if self.user_id is not None:
            return "%s" % (self.user_id.username)
        return "None"

    objects = SamplesInProjectManager()


class CollectionIndexKit(models.Model):
    collection_index_name = models.CharField(max_length=125)
    version = models.CharField(max_length=80, null=True)
    plate_extension = models.CharField(max_length=125, null=True)
    adapter_1 = models.CharField(max_length=125, null=True)
    adapter_2 = models.CharField(max_length=125, null=True)
    collection_index_file = models.FileField(
        upload_to=wetlab_config.COLLECTION_INDEX_KITS_DIRECTORY
    )
    generated_at = models.DateTimeField(auto_now_add=True, null=True)

    class Meta:
        db_table = "wetlab_collection_index_kit"

    def __str__(self):
        return "%s" % (self.collection_index_name)

    def get_collection_index_name(self):
        return "%s" % (self.collection_index_name)

    def get_id(self):
        return "%s" % (self.pk)

    def get_collection_index_information(self):
        if self.adapter_2 == "":
            adapter2 = "Not used on this collection"
        else:
            adapter2 = self.adapter_2
        collection_info = []
        collection_info.append(self.collection_index_name)
        collection_info.append(self.version)
        collection_info.append(self.plate_extension)
        collection_info.append(self.adapter_1)
        collection_info.append(adapter2)
        collection_info.append(self.collection_index_file)
        return collection_info


class CollectionIndexValues(models.Model):
    collection_index_kit_id = models.ForeignKey(
        CollectionIndexKit, on_delete=models.CASCADE
    )

    default_well = models.CharField(max_length=10, null=True)
    index_7 = models.CharField(max_length=25, null=True)
    i_7_seq = models.CharField(max_length=25, null=True)
    index_5 = models.CharField(max_length=25, null=True)
    i_5_seq = models.CharField(max_length=25, null=True)

    class Meta:
        db_table = "wetlab_collection_index_values"

    def get_index_value_information(self):
        return "%s;%s;%s;%s;%s" % (
            self.default_well,
            self.index_7,
            self.i_7_seq,
            self.index_5,
            self.i_5_seq,
        )

    def get_collection_index_id(self):
        return "%s" % (self.collection_index_kit_id.get_id())

    def get_collection_index_name(self):
        return "%s" % (self.collection_index_kit_id.get_collection_index_name())


class LibPreparationUserSampleSheetManager(models.Manager):
    def create_lib_prep_user_sample_sheet(self, user_sample_sheet_data):
        register_user_obj = User.objects.get(
            username__exact=user_sample_sheet_data["user"]
        )
        if user_sample_sheet_data["index_adapters"] == "":
            collection_index_kit_id = None
        else:
            try:
                collection_index_kit_id = CollectionIndexKit.objects.filter(
                    collection_index_name__exact=user_sample_sheet_data["index_adapters"]
                ).last()
            except Exception:
                collection_index_kit_id = None
        file_name = os.path.basename(user_sample_sheet_data["file_name"])
        configuration = SequencingConfiguration.objects.filter(
            platform_id__platformName__exact=user_sample_sheet_data["platform"],
            configuration_name__exact=user_sample_sheet_data["configuration"],
        ).last()
        new_lib_prep_user_sample_sheet = self.create(
            register_user=register_user_obj,
            collection_index_kit_id=collection_index_kit_id,
            reads=",".join(user_sample_sheet_data["reads"]),
            sample_sheet=file_name,
            application=user_sample_sheet_data["application"],
            instrument=user_sample_sheet_data["instrument type"],
            assay=user_sample_sheet_data["assay"],
            adapter1=user_sample_sheet_data["adapter1"],
            adapter2=user_sample_sheet_data["adapter2"],
            sequencing_configuration=configuration,
            iem_version=user_sample_sheet_data["iem_version"],
        )
        return new_lib_prep_user_sample_sheet


class LibUserSampleSheet(models.Model):
    register_user = models.ForeignKey(User, on_delete=models.CASCADE)

    collection_index_kit_id = models.ForeignKey(
        CollectionIndexKit, on_delete=models.CASCADE, null=True
    )

    sequencing_configuration = models.ForeignKey(
        SequencingConfiguration, on_delete=models.CASCADE, null=True
    )

    sample_sheet = models.FileField(
        upload_to=wetlab_config.LIBRARY_PREPARATION_SAMPLE_SHEET_DIRECTORY
    )
    generated_at = models.DateTimeField(auto_now_add=True, null=True)
    application = models.CharField(max_length=70, null=True, blank=True)
    instrument = models.CharField(max_length=70, null=True, blank=True)
    adapter_1 = models.CharField(max_length=70, null=True, blank=True)
    adapter_2 = models.CharField(max_length=70, null=True, blank=True)
    assay = models.CharField(max_length=70, null=True, blank=True)
    reads = models.CharField(max_length=10, null=True, blank=True)
    confirmed_used = models.BooleanField(default=False)
    iem_version = models.CharField(max_length=5, null=True, blank=True)

    class Meta:
        db_table = "wetlab_lib_user_samplesheet"

    def __str__(self):
        return "%s" % (self.sample_sheet)

    def get_assay(self):
        return "%s" % (self.assay)

    def get_adapters(self):
        adapters = []
        adapters.append(self.adapter1)
        adapters.append(self.adapter_2)
        return adapters

    def get_collection_index_kit(self):
        if self.collection_index_kit_id is not None:
            return "%s" % (self.collection_index_kit_id.get_collection_index_name())
        else:
            return "Not available"

    def get_iem_version(self):
        return "%s" % (self.iem_version)

    def get_sequencing_configuration_platform(self):
        return "%s" % (self.sequencing_configuration.get_platform_name())

    def get_sequencing_configuration_name(self):
        return "%s" % (self.sequencing_configuration.get_configuration_name())

    def get_user_sample_sheet_id(self):
        return "%s" % (self.pk)

    def get_all_data(self):
        if self.collection_index_kit_id is not None:
            collection_index = self.collection_index_kit_id.get_collection_index_name()
        else:
            collection_index = ""
        s_s_data = []
        s_s_data.append(collection_index)
        s_s_data.append(self.application)
        s_s_data.append(self.instrument)
        s_s_data.append(self.adapter1)
        s_s_data.append(self.adapter_2)
        s_s_data.append(self.assay)
        s_s_data.append(self.reads.split(",")[0])
        return s_s_data

    def update_confirm_used(self, confirmation):
        self.confirmed_used = confirmation
        self.save()
        return self

    objects = LibPreparationUserSampleSheetManager()


class PoolStates(models.Model):
    pool_state = models.CharField(max_length=50)

    class Meta:
        db_table = "wetlab_pool_states"

    def __str__(self):
        return "%s" % (self.pool_state)

    def get_pool_state(self):
        return "%s" % (self.pool_state)


class LibraryPoolManager(models.Manager):
    def create_lib_pool(self, pool_data):
        platform_obj = SequencingPlatform.objects.filter(
            platform_name__exact=pool_data["platform"]
        ).last()
        new_library_pool = self.create(
            register_user=pool_data["registerUser"],
            pool_state=PoolStates.objects.get(pool_state__exact="Defined"),
            pool_name=pool_data["poolName"],
            pool_code_id=pool_data["poolCodeID"],
            adapter=pool_data["adapter"],
            paired_end=pool_data["pairedEnd"],
            sample_number=pool_data["n_samples"],
            platform=platform_obj,
        )
        return new_library_pool


class LibraryPool(models.Model):
    register_user = models.ForeignKey(User, on_delete=models.CASCADE)
    pool_state = models.ForeignKey(PoolStates, on_delete=models.CASCADE)
    run_process_id = models.ForeignKey(
        RunProcess, on_delete=models.CASCADE, null=True, blank=True
    )

    platform = models.ForeignKey(
        SequencingPlatform, on_delete=models.CASCADE, null=True, blank=True
    )
    pool_name = models.CharField(max_length=50)
    sample_number = models.IntegerField(default=0)
    pool_code_id = models.CharField(max_length=50, blank=True)
    adapter = models.CharField(max_length=50, null=True, blank=True)
    paired_end = models.CharField(max_length=10, null=True, blank=True)
    generated_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        ordering = ("pool_name",)
        db_table = "wetlab_library_pool"

    def __str__(self):
        return "%s" % (self.pool_name)

    def get_adapter(self):
        return "%s" % (self.adapter)

    def get_id(self):
        return "%s" % (self.pk)

    def get_info(self):
        pool_info = []
        pool_info.append(self.pool_name)
        pool_info.append(self.pool_code_id)
        pool_info.append(self.sample_number)
        return pool_info

    def get_number_of_samples(self):
        return "%s" % (self.sample_number)

    def get_pool_name(self):
        return "%s" % (self.pool_name)

    def get_pool_code_id(self):
        return "%s" % (self.pool_code_id)

    def get_pool_single_paired(self):
        return "%s" % (self.paired_end)

    def get_platform_name(self):
        if self.platform is not None:
            return "%s" % (self.platform.get_platform_name())
        else:
            return "Not defined"

    def get_run_name(self):
        if self.run_process_id is not None:
            return "%s" % (self.run_process_id.get_run_name())
        else:
            return "Not defined yet"

    def get_run_id(self):
        if self.run_process_id is not None:
            return "%s" % (self.run_process_id.get_run_id())
        return None

    def get_run_obj(self):
        return self.run_process_id

    def set_pool_state(self, state):
        self.pool_state = PoolStates.objects.get(pool_state__exact=state)
        self.save()

    def update_number_samples(self, number_s_in_pool):
        self.sample_number = number_s_in_pool
        self.save()
        return self

    def update_run_name(self, run_name):
        self.run_process_id = run_name
        self.save()
        return self

    objects = LibraryPoolManager()


class LibPrepareStates(models.Model):
    lib_prep_state = models.CharField(max_length=50)

    class Meta:
        db_table = "wetlab_lib_prepare_states"

    def __str__(self):
        return "%s" % (self.lib_prep_state)

    def get_lib_state(self):
        return "%s" % (self.lib_prep_state)


class LibPrepareManager(models.Manager):
    def create_lib_preparation(self, lib_prep_data):
        register_user_obj = User.objects.get(
            username__exact=lib_prep_data["registerUser"]
        )
        sample_obj = Samples.objects.get(pk__exact=lib_prep_data["sample_id"])
        molecule_obj = MoleculePreparation.objects.get(
            pk__exact=lib_prep_data["molecule_id"]
        )
        lib_state_obj = LibPrepareStates.objects.get(lib_prep_state__exact="Defined")
        new_lib_prep = self.create(
            register_user=register_user_obj,
            molecule_id=molecule_obj,
            sample_id=sample_obj,
            protocol_id=lib_prep_data["protocol_obj"],
            lib_prep_state=lib_state_obj,
            lib_prep_code_id=lib_prep_data["lib_prep_code_id"],
            user_sample_id=lib_prep_data["user_sampleID"],
            unique_id=lib_prep_data["uniqueID"],
            prefix_protocol=lib_prep_data["prefixProtocol"],
            samplename_in_samplesheet=lib_prep_data["sampleNameInSampleSheet"],
        )
        return new_lib_prep


class LibPrepare(models.Model):
    register_user = models.ForeignKey(User, on_delete=models.CASCADE)
    molecule_id = models.ForeignKey(MoleculePreparation, on_delete=models.CASCADE)
    sample_id = models.ForeignKey(
        Samples, on_delete=models.CASCADE, null=True, blank=True
    )
    protocol_id = models.ForeignKey(
        Protocols, on_delete=models.CASCADE, null=True, blank=True
    )
    lib_prep_state = models.ForeignKey(
        LibPrepareStates, on_delete=models.CASCADE, null=True, blank=True
    )

    user_lot_kit_id = models.ForeignKey(
        UserLotCommercialKits, on_delete=models.CASCADE, null=True, blank=True
    )

    user_sample_sheet = models.ForeignKey(
        LibUserSampleSheet, on_delete=models.CASCADE, null=True, blank=True
    )

    pools = models.ManyToManyField(LibraryPool, blank=True)

    lib_prep_code_id = models.CharField(max_length=255, null=True, blank=True)
    user_sample_id = models.CharField(max_length=100, null=True, blank=True)
    project_in_samplesheet = models.CharField(max_length=80, null=True, blank=True)
    sample_plate = models.CharField(max_length=50, null=True, blank=True)
    sample_well = models.CharField(max_length=20, null=True, blank=True)
    index_plate_well = models.CharField(max_length=20, null=True, blank=True)
    i7_index_id = models.CharField(max_length=25, null=True, blank=True)
    i7_index = models.CharField(max_length=30, null=True, blank=True)
    i5_index_id = models.CharField(max_length=25, null=True, blank=True)
    i5_index = models.CharField(max_length=30, null=True, blank=True)
    genome_folder = models.CharField(max_length=180, null=True, blank=True)
    manifest = models.CharField(max_length=80, null=True, blank=True)
    reused_number = models.IntegerField(default=0)
    unique_id = models.CharField(max_length=16, null=True, blank=True)
    user_in_samplesheet = models.CharField(max_length=255, null=True, blank=True)
    samplename_in_samplesheet = models.CharField(max_length=255, null=True, blank=True)
    prefix_protocol = models.CharField(max_length=25, null=True, blank=True)

    class Meta:
        db_table = "wetlab_lib_prepare"
        ordering = ("lib_prep_code_id",)

    def __str__(self):
        return "%s" % (self.lib_prep_code_id)

    def get_adapters(self):
        return self.user_sample_sheet.get_adapters()

    def get_id(self):
        return "%s" % (self.pk)

    def get_info_for_display(self):
        lib_info = []
        lib_info.append(self.lib_prep_code_id)
        lib_info.append(self.molecule_id.get_molecule_code_id())
        lib_info.append(self.lib_prep_state.lib_prep_state)
        lib_info.append(self.protocol_id.get_name())
        lib_info.append(self.project_in_samplesheet)
        lib_info.append(self.i7_index_id)
        lib_info.append(self.i5_index_id)
        lib_info.append(self.reused_number)
        return lib_info

    def get_info_for_selection_in_pool(self):
        lib_info = []
        lib_info.append(self.lib_prep_code_id)
        lib_info.append(self.sample_id.get_sample_name())
        lib_info.append(self.register_user.username)
        lib_info.append(self.user_sample_sheet.get_collection_index_kit())
        lib_info.append(self.i7_index_id)
        lib_info.append(self.i5_index_id)
        lib_info.append(self.pk)
        return lib_info

    def get_info_for_run_paired_end(self):
        lib_info = []
        lib_info.append(self.unique_id)
        lib_info.append(self.samplename_in_samplesheet)
        lib_info.append(self.sample_plate)
        lib_info.append(self.sample_well)
        lib_info.append(self.i7_index_id)
        lib_info.append(self.i7_index)
        lib_info.append(self.i5_index_id)
        lib_info.append(self.i5_index)
        lib_info.append(self.project_in_samplesheet)
        lib_info.append(self.user_in_samplesheet)
        return lib_info

    def get_info_for_run_single_read(self):
        lib_info = []
        lib_info.append(self.unique_id)
        lib_info.append(self.sample_id.get_sample_name())
        lib_info.append(self.sample_plate)
        lib_info.append(self.sample_well)
        lib_info.append(self.i7_index_id)
        lib_info.append(self.i7_index)
        lib_info.append(self.project_in_samplesheet)
        lib_info.append(self.user_in_samplesheet)
        return lib_info

    def get_basic_data(self):
        lib_info = []
        lib_info.append(self.sample_id.get_sample_name())
        lib_info.append(self.unique_id)
        lib_info.append(self.protocol_id.get_name())
        lib_info.append(self.pk)
        return lib_info

    def get_info_for_display_pool(self):
        if self.register_user.username is None:
            user_name = ""
        else:
            user_name = self.register_user.username

        lib_info = []
        lib_info.append(self.lib_prep_code_id)
        lib_info.append(self.get_sample_name())
        lib_info.append(user_name)
        lib_info.append(self.get_sample_id())
        return lib_info

    def get_i7_index(self):
        if self.i7_index_id is not None:
            return "%s" % (self.i7_index_id)
        else:
            return ""

    def get_i5_index(self):
        if self.i5_index_id is not None:
            return "%s" % (self.i5_index_id)
        else:
            return ""

    def get_iem_version(self):
        if self.user_sample_sheet is not None:
            return "%s" % (self.user_sample_sheet.get_iem_version())
        else:
            return "None"

    def get_indexes(self):
        index = {}
        if self.i5_index_id is None:
            index["i5_indexID"] = ""
            index["i5_seq"] = ""
        else:
            index["i5_indexID"] = self.i5_index_id
            index["i5_seq"] = self.i5_index
        index["i7_indexID"] = self.i7_index_id
        index["i7_seq"] = self.i7_index
        return index

    def get_index_plate_well(self):
        if self.index_plate_well is not None:
            return "%s" % (self.index_plate_well)
        else:
            return ""

    def get_lib_prep_code(self):
        return "%s" % (self.lib_prep_code_id)

    def get_lib_prep_id(self):
        return "%s" % (self.pk)

    def get_manifest(self):
        if self.manifest is not None:
            return "%s" % (self.manifest)
        else:
            return ""

    def get_genome_folder(self):
        if self.genome_folder is not None:
            return "%s" % (self.genome_folder)
        else:
            return ""

    def get_molecule_code_id(self):
        return "%s" % (self.molecule_id.get_molecule_code_id())

    def get_molecule_obj(self):
        return self.molecule_id

    def get_protocol_used(self):
        return "%s" % (self.protocol_id.get_name())

    def get_protocol_obj(self):
        return self.protocol_id

    def get_protocol_id(self):
        return "%s" % (self.protocol_id.pk)

    def get_reused_value(self):
        return "%s" % (self.reused_number)

    def get_sample_name(self):
        return "%s" % (self.sample_id.get_sample_name())

    def get_sample_id(self):
        return self.sample_id.pk

    def get_sample_obj(self):
        return self.sample_id

    def get_sample_well(self):
        return "%s" % (self.sample_well)

    def get_state(self):
        return "%s" % (self.lib_prep_state.lib_prep_state)

    def get_unique_id(self):
        return "%s" % (self.unique_id)

    def get_user_obj(self):
        return self.register_user

    def get_user_sample_sheet_obj(self):
        return self.user_sample_sheet

    def set_increase_reuse(self):
        self.reused_number += 1
        self.save()
        return

    def set_pool(self, pool_obj):
        self.pools.add(pool_obj)
        self.save()
        return

    def set_state(self, state_value):
        self.lib_prep_state = LibPrepareStates.objects.get(
            lib_prep_state__exact=state_value
        )
        self.save()
        return

    def update_i7_index(self, i7_seq_value):
        self.i7_index = i7_seq_value
        self.save()
        return

    def update_i5_index(self, i5_seq_value):
        self.i5_index = i5_seq_value
        self.save()
        return

    def update_library_preparation_with_indexes(self, lib_prep_data):
        self.user_sample_sheet = lib_prep_data["user_sample_sheet"]
        self.project_in_samplesheet = lib_prep_data["projectInSampleSheet"]
        self.sample_plate = lib_prep_data["samplePlate"]
        self.sample_well = lib_prep_data["sampleWell"]
        self.i7_index_id = lib_prep_data["i7IndexID"]
        self.i7_index = lib_prep_data["i7Index"]
        if "i5IndexID" in lib_prep_data:
            self.i5_index_id = lib_prep_data["i5IndexID"]
        if "i5Index" in lib_prep_data:
            self.i5_index = lib_prep_data["i5Index"]
        if "indexPlateWell" in lib_prep_data:
            self.index_plate_well = lib_prep_data["indexPlateWell"]
        if "genomeFolder" in lib_prep_data:
            self.genome_folder = lib_prep_data["genomeFolder"]
        if "manifest" in lib_prep_data:
            self.manifest = lib_prep_data["manifest"]
        self.user_in_samplesheet = lib_prep_data["userInSampleSheet"]
        self.user_sample_id = lib_prep_data["userSampleID"]
        self.save()
        return self

    def update_lib_preparation_info_in_reuse_state(self, lib_prep_data):
        # Update the instance with the sample sheet information
        self.register_user = User.objects.get(
            username__exact=lib_prep_data["registerUser"]
        )
        self.protocol_id = lib_prep_data["protocol_obj"]
        self.user_sample_sheet = lib_prep_data["user_sample_sheet"]
        self.lib_prep_code_id = lib_prep_data["lib_prep_code_id"]
        self.project_in_samplesheet = lib_prep_data["projectInSampleSheet"]
        self.user_sample_id = lib_prep_data["userSampleID"]
        self.sample_plate = lib_prep_data["samplePlate"]
        self.sample_well = lib_prep_data["sampleWell"]
        self.index_plate_well = lib_prep_data["indexPlateWell"]
        self.i7_index_id = lib_prep_data["i7IndexID"]
        self.i7_index = lib_prep_data["i7Index"]
        self.i5_index_id = lib_prep_data["i5_index_id"]
        self.i5_index = lib_prep_data["i5Index"]
        self.unique_id = lib_prep_data["uniqueID"]
        self.save()
        return self

    objects = LibPrepareManager()


class LibParameterValueManager(models.Manager):
    def create_library_parameter_value(self, parameter_value):
        new_parameter_data = self.create(
            parameter_id=parameter_value["parameter_id"],
            library_id=parameter_value["library_id"],
            parameter_value=parameter_value["parameterValue"],
        )
        return new_parameter_data


class LibParameterValue(models.Model):
    parameter_id = models.ForeignKey(ProtocolParameters, on_delete=models.CASCADE)
    library_id = models.ForeignKey(LibPrepare, on_delete=models.CASCADE)
    parameter_value = models.CharField(max_length=255)
    generated_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        db_table = "wetlab_lib_parameter_value"

    def __str__(self):
        return "%s" % (self.parameter_value)

    def get_parameter_information(self):
        return "%s" % (self.parameter_value)

    objects = LibParameterValueManager()


class AdditionaKitsLibPrepareManager(models.Manager):
    def create_additional_kit(self, kit_data):
        new_additional_kit = self.create(
            register_user=kit_data["user"],
            protocol_id=kit_data["protocol_id"],
            commercial_kit_id=kit_data["commercialKit_id"],
            kit_name=kit_data["kitName"],
            description=kit_data["description"],
            kit_order=kit_data["kitOrder"],
            kit_used=kit_data["kitUsed"],
        )
        return new_additional_kit


class AdditionaKitsLibPrepare(models.Model):
    register_user = models.ForeignKey(User, on_delete=models.CASCADE)
    protocol_id = models.ForeignKey(Protocols, on_delete=models.CASCADE)
    commercial_kit_id = models.ForeignKey(CommercialKits, on_delete=models.CASCADE)
    kit_name = models.CharField(max_length=255)
    description = models.CharField(max_length=400, null=True, blank=True)
    kit_order = models.IntegerField()
    kit_used = models.BooleanField()
    generated_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        db_table = "wetlab_lib_additional_kits_lib_prepare"

    def __str__(self):
        return "%s" % (self.kit_name)

    def get_kit_name(self):
        return "%s" % (self.kit_name)

    def get_add_kit_id(self):
        return "%s" % (self.pk)

    def get_all_kit_info(self):
        data = []
        data.append(self.kit_name)
        data.append(self.kit_order)
        data.append(self.kit_used)
        data.append(self.commercial_kit_id.get_name())
        data.append(self.description)
        return data

    def get_add_kit_data_for_javascript(self):
        if self.kit_used:
            used = "true"
        else:
            used = "false"
        add_kit_data = []
        add_kit_data.append(self.kit_name)
        add_kit_data.append(self.kit_order)
        add_kit_data.append(used)
        add_kit_data.append(self.commercial_kit_id.get_name())
        add_kit_data.append(self.description)
        return add_kit_data

    def get_commercial_kit_obj(self):
        return self.commercial_kit_id

    def get_commercial_kit_name(self):
        return "%s" % (self.commercial_kit_id.get_name())

    def update_add_kit_fields(self, kit_data):
        self.register_user = kit_data["user"]
        self.commercial_kit_id = kit_data["commercialKit_id"]
        self.kit_name = kit_data["kitName"]
        self.description = kit_data["description"]
        self.kit_order = kit_data["kitOrder"]
        self.kit_used = kit_data["kitUsed"]
        self.save()
        return self

    objects = AdditionaKitsLibPrepareManager()


class AdditionalUserLotKitManager(models.Manager):
    def create_additional_user_lot_kit(self, user_additional_kit):
        new_user_additional_kit = self.create(
            lib_prep_id=user_additional_kit["lib_prep_id"],
            additional_lot_kits=user_additional_kit["additionalLotKits"],
            user_lot_kit_id=user_additional_kit["userLotKit_id"],
            value=0,
        )
        return new_user_additional_kit


class AdditionalUserLotKit(models.Model):
    lib_prep_id = models.ForeignKey(LibPrepare, on_delete=models.CASCADE)
    additional_lot_kits = models.ForeignKey(
        AdditionaKitsLibPrepare, on_delete=models.CASCADE
    )
    user_lot_kit_id = models.ForeignKey(
        UserLotCommercialKits, on_delete=models.CASCADE, null=True, blank=True
    )
    value = models.CharField(max_length=255)
    generated_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        db_table = "wetlab_lib_additional_user_lot_kit"

    def __str__(self):
        return "%s" % (self.userLotKit_id)

    def get_additional_kit_info(self):
        data = []
        data.append(self.additional_lot_kits.get_kit_name())
        data.append(self.additional_lot_kits.get_commercial_kit_name())
        if self.user_lot_kit_id is None:
            data.append("Not set")
        else:
            data.append(self.user_lot_kit_id.get_lot_number())
        data.append(self.generated_at.strftime("%d %B %Y"))
        return data

    objects = AdditionalUserLotKitManager()


class SambaConnectionData(models.Model):
    samba_folder_name = models.CharField(
        max_length=80, null=True, blank=True
    )
    domain = models.CharField(max_length=80, null=True, blank=True)
    host_name = models.CharField(max_length=80, null=True, blank=True)
    ip_server = models.CharField(max_length=20, null=True, blank=True)
    port_server = models.CharField(max_length=10, null=True, blank=True)
    remote_server_name = models.CharField(max_length=80, null=True, blank=True)
    shared_folder_name = models.CharField(max_length=80, null=True, blank=True)
    user_id = models.CharField(max_length=80, null=True, blank=True)
    user_password = models.CharField(max_length=20, null=True, blank=True)
    is_direct_tcp = models.BooleanField(default=True)
    ntlm_used = models.BooleanField(default=True)

    class Meta:
        db_table = "wetlab_samba_connection_data"

    def __str__(self):
        return "%s" % (self.remote_server_name)

    def get_samba_data(self):
        samba_data = {}
        for field in wetlab_config.SAMBA_CONFIGURATION_FIELDS:
            samba_data[field] = getattr(self, field)
        return samba_data

    def get_samba_application_folder_name(self):
        if self.samba_folder_name is not None:
            return "%s" % (self.samba_folder_name)
        else:
            return ""

    def get_samba_folder_name(self):
        if self.shared_folder_name is not None:
            return "%s" % (self.shared_folder_name)
        else:
            return ""

    def update_data(self, data):
        for field in wetlab_config.SAMBA_CONFIGURATION_FIELDS:
            setattr(self, field, data[field])
        self.save()
        return self


class ConfigSettingManager(models.Manager):
    def create_config_setting(self, configuration_name, configuration_value):
        new_config_settings = self.create(
            configuration_name=configuration_name, configurationValue=configuration_value
        )
        return new_config_settings


class ConfigSetting(models.Model):
    configuration_name = models.CharField(max_length=80)
    configuration_value = models.CharField(max_length=255, null=True, blank=True)
    generated_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        db_table = "wetlab_config_setting"

    def __str__(self):
        return "%s" % (self.configuration_name)

    def get_configuration_value(self):
        return "%s" % (self.configuration_value)

    def set_configuration_value(self, new_value):
        self.configuration_value = new_value
        self.save()
        return self

    objects = ConfigSettingManager()


class RunConfigurationTest(models.Model):
    run_test_name = models.CharField(max_length=80)
    run_test_folder = models.CharField(max_length=200)

    class Meta:
        db_table = "wetlab_lib_run_configuration_test"

    def __str__(self):
        return "%s" % (self.run_test_name)

    def get_run_test_id(self):
        return "%s" % (self.pk)

    def get_run_test_name(self):
        return "%s" % (self.run_test_name)

    def get_run_test_folder(self):
        return "%s" % (self.run_test_folder)
