import datetime, os

from django.db import models
from django import forms
from django.utils import timezone
#from django.utils.encoding import python_2_unicode_compatible
from django.contrib.auth.models import User
from django_utils.models import Center
from django.utils.translation import ugettext_lazy as _

from .  import wetlab_config
from iSkyLIMS_core.models import MoleculePreparation , Samples , ProtocolType, Protocols, ProtocolParameters, UserLotCommercialKits, CommercialKits, SequencerInLab, SequencingConfiguration, SequencingPlatform

class FlexibleConfSettings(models.Model):
    confParameterName = models.CharField(max_length=255)
    confParameterValue = models.CharField(max_length=255)

    def __str__ (self):
        return '%s' %(self.confParameterName)

    def get_parameter_value(self):
        return '%s' %(self.confParamValue)

class RunErrors (models.Model):
    errorCode = models.CharField(max_length=10)
    errorText = models.CharField(max_length=255)

    def __str__ (self):
        return '%s' %(self.errorText)

class RunStates (models.Model):
    runStateName = models.CharField(max_length=50)

    def __str__ (self):
        return '%s' %(self.runStateName)

    def get_run_state_name(self):
        return '%s' %(self.runStateName)

class RunProcess(models.Model):
    usedSequencer = models.ForeignKey(
                        SequencerInLab,
                        on_delete=models.CASCADE, blank=True, null = True)
    runError = models.ForeignKey(
                        RunErrors,
                        on_delete=models.CASCADE, null = True, blank = True)
    stateBeforeError = models.ForeignKey (
                        RunStates,
                        on_delete = models.CASCADE, null = True, blank = True)
    state = models.ForeignKey (
                        RunStates,
                        on_delete = models.CASCADE, related_name = 'state_of_run', null = True, blank = True)
    centerRequestedBy = models.ForeignKey (
                        Center,
                        on_delete=models.CASCADE, null = True, blank = True )
    reagent_kit = models.ManyToManyField(UserLotCommercialKits, blank = True)


    runName = models.CharField(max_length=45)
    sampleSheet = models.FileField(upload_to = wetlab_config.RUN_SAMPLE_SHEET_DIRECTORY, null = True, blank = True)
    generatedat = models.DateTimeField(auto_now_add=True)
    run_date = models.DateField(auto_now = False, null=True, blank = True)
    run_finish_date = models.DateTimeField(auto_now = False, null=True, blank=True)
    bcl2fastq_finish_date = models.DateTimeField(auto_now = False, null=True, blank=True)
    run_completed_date = models.DateTimeField(auto_now = False, null=True, blank=True)
    #runState = models.CharField(max_length=25)

    index_library = models.CharField(max_length=85, null = True, blank = True)
    samples= models.CharField(max_length=45,blank=True)
    useSpaceImgMb=models.CharField(max_length=10, blank=True)
    useSpaceFastaMb=models.CharField(max_length=10, blank=True)
    useSpaceOtherMb=models.CharField(max_length=10, blank=True)

    # sequencerModel = models.ForeignKey ('iSkyLIMS_drylab.Machines', on_delete=models.CASCADE, null=True, blank=True)


    def __str__(self):
        return '%s' %(self.runName)

    def get_run_id (self):
        return '%s' %(self.id)

    def get_recorded_date_no_format(self):
        return self.generatedat

    def get_run_date_no_format(self):
        return self.run_date

    def get_run_finish_date_no_format(self):
        return self.run_finish_date

    def get_run_finish_date(self):
        return self.run_finish_date.strftime("%B %d, %Y")

    def get_run_date (self):
        if self.run_date is None :
            rundate = 'Run NOT started'
        else :
            rundate=self.run_date.strftime("%B %d, %Y")
        return rundate

    def get_run_year(self):
        return '%s' %(self.run_date.timetuple().tm_year)

    def get_error_text (self):
        return '%s' %(self.runError)

    def get_state(self):
        return '%s' %(self.state.get_run_state_name())

    def get_state_before_error(self):
        return '%s' %(self.stateBeforeError)

    def get_info_process (self):
        generated_date=self.generatedat.strftime("%I:%M%p on %B %d, %Y")

        requested_center = str(self.centerRequestedBy)
        if self.run_date is None :
            rundate = 'Run NOT started'
        else :
            rundate=self.run_date.strftime("%B %d, %Y")

        if self.run_finish_date is None:
            finish_date = 'Run multiplexation is not completed'
        else:
            finish_date = self.run_finish_date.strftime("%I:%M%p on %B %d, %Y")

        if self.bcl2fastq_finish_date is None:
            bcl2fastq_date = 'bcl2fastq process is not completed'
        else:
            bcl2fastq_date = self.bcl2fastq_finish_date.strftime("%I:%M%p on %B %d, %Y")

        if self.run_completed_date is None:
            completed_date = 'Run process is not completed'
        else:
            completed_date = self.run_completed_date.strftime("%I:%M%p on %B %d, %Y")


        return '%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s'  %(self.runName, self.state,
                requested_center, generated_date, rundate,completed_date, bcl2fastq_date,
                finish_date, self.useSpaceImgMb, self.useSpaceFastaMb,
                self.useSpaceOtherMb )


    def get_run_name (self):
        return '%s' %(self.runName)

    def get_disk_space_utilization (self):
        image_size_str = '0' if self.useSpaceImgMb == '' else self.useSpaceImgMb
        data_size_str = '0' if self.useSpaceFastaMb == '' else self.useSpaceFastaMb
        other_size_str = '0' if self.useSpaceOtherMb == '' else self.useSpaceOtherMb
        image_size = int(image_size_str.replace(',',''))
        data_size = int(data_size_str.replace(',',''))
        other_size = int(other_size_str.replace(',',''))
        total_size = image_size + data_size + other_size
        return '%s'%(total_size)

    def get_run_used_sequencer_name (self):
        return '%s' %(self.usedSequencer.get_sequencer_name())

    def get_run_used_sequencer (self):
        return '%s' %(self.usedSequencer)

    def get_run_platform (self):
        return '%s' %(self.usedSequencer.get_sequencing_platform_name())
    '''
    def get_machine_lanes(self):
        number_of_lanes = self.sequencerModel.get_number_of_lanes()
        return int(number_of_lanes)
    '''
    def get_sequencing_lanes(self):
        number_of_lanes = self.usedSequencer.get_number_of_lanes()
        return int(number_of_lanes)


    def get_sample_file (self):
        return '%s' %(self.sampleSheet)

    def update_index_library (self, index_library_name):
        self.index_library = index_library_name
        self.save()
        return True

    def update_sample_sheet(self, sample_file):
        self.sampleSheet = sample_file
        self.save()
        return True

    def set_used_space (self, disk_utilization):
        self.useSpaceFastaMb  = disk_utilization ['useSpaceFastaMb']
        self.useSpaceImgMb  = disk_utilization ['useSpaceImgMb']
        self.useSpaceOtherMb  = disk_utilization ['useSpaceOtherMb']
        self.save()
        return True

    def set_run_state (self, new_state):
        if RunStates.objects.filter(runStateName__exact = new_state).exists():
            self.state = RunStates.objects.get(runStateName__exact = new_state)
            self.save()
        else:
            return False
        return True

    def set_run_date (self, run_date):
        self.run_date = run_date
        self.save()
        return True

    def set_run_bcl2fastq_finished_date (self, bcl2fastq_finish_date):
        self.bcl2fastq_finish_date = bcl2fastq_finish_date
        self.save()
        return True

    def set_run_error_code (self, error_code):
        if RunErrors.objects.filter(errorCode__exact = error_code).exists():
            self.runError = RunErrors.objects.get(errorCode__exact = error_code)
        else:
            self.runError = RunErrors.objects.get(errorText__exact = 'Undefined')
        self.stateBeforeError = self.state
        self.state = RunStates.objects.get(runStateName__exact = 'Error')
        self.save()
        return True

    def set_run_completion_date (self, completion_date):
        self.run_completed_date = completion_date
        self.save()
        return True

    def set_run_finish_date (self, finish_date):
        self.run_finish_date = finish_date
        self.save()
        return True

    def set_run_sample_sheet(self, sample_sheet):
        self.sampleSheet = sample_sheet
        self.save()
        return True

    def set_used_sequencer (self, sequencer):
        self.usedSequencer = sequencer
        self.save()
        return True

class LibraryKit (models.Model):
    libraryName = models.CharField(max_length=125)
    #sampleNumber = models.CharField(max_length = 25)
    #indexNumber = models.CharField(max_length = 25)
    generatedat = models.DateTimeField(auto_now_add=True)

    def get_bs_lib_name(self):
        return '%s' %(self.libraryName)


class ProjectsManager(models.Manager):
    def create_new_empty_project(self,project_data):
        new_project = self.create(user_id = project_data['user_id'], projectName = project_data['projectName'])
        return new_project

    def create_new_project(self,project_data):
        new_project = self.create(user_id = project_data['user_id'], projectName = project_data['projectName'],
                    libraryKit = project_data['libraryKit'], baseSpaceFile = project_data['baseSpaceFile'],
                    BaseSpaceLibrary = project_data['BaseSpaceLibrary'])
        return new_project


class Projects(models.Model):
    #runprocess_id = models.ForeignKey(
    #        RunProcess,
    #        on_delete=models.CASCADE, null = True) # added null for new lab process functionality
    user_id= models.ForeignKey(
                    User,
                    on_delete = models.CASCADE, null = True)
    LibraryKit_id = models.ForeignKey(
                    LibraryKit,
                    on_delete=models.CASCADE , null=True, blank = True)
    runProcess = models.ManyToManyField(RunProcess)

    BaseSpaceLibrary = models.CharField(max_length=45, null=True, blank=True)
    projectName= models.CharField(max_length=45)
    libraryKit=models.CharField(max_length=125)
    baseSpaceFile = models.CharField(max_length=255)
    generatedat = models.DateTimeField(auto_now_add=True)
    project_run_date = models.DateField(auto_now = False, null=True)
    #sampleSheet = models.FileField(upload_to = wetlab_config.SAMPLE_SHEET_CREATED_ON_LAB, null = True) # added null for new lab process functionality
    #pairedEnd  = models.BooleanField(null = True, blank = True)
    #read_length = models.CharField(max_length = 10)

    def __str__(self):
        return '%s' %(self.projectName)

    def get_base_space_file (self):
        return '%s' %(self.baseSpaceFile)

    def get_state(self):
        return '%s' %(self.runprocess_id.state.runStateName)

    def get_project_dates (self):
        if self.project_run_date is None:
            projectdate = 'Run NOT started'
        else :
            projectdate=self.project_run_date.strftime("%B %d, %Y")
        p_info = []
        p_info.append(self.generatedat.strftime("%B %d, %Y"))
        p_info.append(projectdate)
        return p_info

    def get_p_info_change_library (self):
        run_name = self.runprocess_id.runName
        user_name = self.user_id.username
        if self.project_run_date is None:
            projectdate = 'Run NOT started'
        else :
            projectdate=self.project_run_date.strftime("%B %d, %Y")

        return '%s;%s;%s;%s;%s'%(run_name, self.projectName, projectdate, user_name, self.libraryKit)

    #def get_run_name(self):
    #   return '%s' %(self.runprocess_id.get_run_name())

    def get_run_id(self):
        return '%s' %(self.runprocess_id.get_run_id())

    def get_run_obj(self):
        return self.runprocess_id

    def get_user_name (self):
        user_name = self.user_id.username
        return '%s' %(user_name)

    def get_project_name (self):
        return '%s' %(self.projectName)

    def get_index_library_name (self):
        return '%s' %(self.libraryKit)

    def get_date (self):
        if self.project_run_date is None:
            projectdate = 'No Date'
        else :
            projectdate=self.project_run_date.strftime("%B %d, %Y")
        return '%s' %(projectdate)

    def get_project_id(self):
        return '%s' %(self.id)

    def add_run(run_obj):
        self.runProcess.add(run_obj)
        return self

    def set_project_run_date(self, date):
        self.project_run_date = date
        self.save()
        return self

    objects = ProjectsManager()

class RunningParametersManager (models.Manager) :

    def create_running_parameters (self, running_data, run_object) :

        if running_data['PlannedRead1Cycles'] == "":
            running_data['PlannedRead1Cycles'] = 0
        if running_data['PlannedRead2Cycles'] == "":
            running_data['PlannedRead2Cycles'] = 0
        if running_data['PlannedIndex1ReadCycles'] == "":
            running_data['PlannedIndex1ReadCycles'] = 0
        if running_data['PlannedIndex2ReadCycles'] == "":
            running_data['PlannedIndex2ReadCycles'] = 0

        if running_data['LibraryID'] == None:
            running_data['LibraryID'] = " "
        if running_data['AnalysisWorkflowType'] == None :
            running_data['AnalysisWorkflowType'] = " "

        running_parameters = self.create (runName_id = run_object,
                         RunID=running_data['RunID'], ExperimentName=running_data['ExperimentName'],
                         RTAVersion=running_data['RTAVersion'], SystemSuiteVersion= running_data['SystemSuiteVersion'],
                         LibraryID= running_data['LibraryID'], Chemistry= running_data['Chemistry'],
                         RunStartDate= running_data['RunStartDate'], AnalysisWorkflowType= running_data['AnalysisWorkflowType'],
                         RunManagementType= running_data['RunManagementType'], PlannedRead1Cycles= running_data['PlannedRead1Cycles'],
                         PlannedRead2Cycles= running_data['PlannedRead2Cycles'], PlannedIndex1ReadCycles= running_data['PlannedIndex1ReadCycles'],
                         PlannedIndex2ReadCycles= running_data['PlannedIndex2ReadCycles'], ApplicationVersion= running_data['ApplicationVersion'],
                         NumTilesPerSwath= running_data['NumTilesPerSwath'], ImageChannel= running_data['ImageChannel'],
                         Flowcell= running_data['Flowcell'], ImageDimensions= running_data['ImageDimensions'],
                         FlowcellLayout= running_data['FlowcellLayout'])

        return running_parameters

class RunningParameters (models.Model):
    runName_id = models.OneToOneField(
            RunProcess,
            on_delete=models.CASCADE,
            primary_key=True,
            )
    RunID= models.CharField(max_length=255)
    ExperimentName= models.CharField(max_length=255)
    RTAVersion= models.CharField(max_length=255)
    SystemSuiteVersion= models.CharField(max_length=255, null=True)
    LibraryID= models.CharField(max_length=255, null=True)
    Chemistry= models.CharField(max_length=255)
    RunStartDate= models.CharField(max_length=255)
    AnalysisWorkflowType= models.CharField(max_length=255)
    RunManagementType= models.CharField(max_length=255)
    PlannedRead1Cycles= models.CharField(max_length=255)
    PlannedRead2Cycles= models.CharField(max_length=255)
    PlannedIndex1ReadCycles= models.CharField(max_length=255)
    PlannedIndex2ReadCycles= models.CharField(max_length=255)
    ApplicationVersion= models.CharField(max_length=255)
    NumTilesPerSwath= models.CharField(max_length=255)
    ImageChannel= models.CharField(max_length=255,null=True)
    Flowcell= models.CharField(max_length=255)
    ImageDimensions= models.CharField(max_length=255,null=True)
    FlowcellLayout= models.CharField(max_length=255)

    def __str__(self):
        return '%s' %(self.RunID)
        #return '%s' %(self.runName_id)

    def get_run_chemistry(self):
        return '%s' %(self.Chemistry)

    def get_run_parameters_info (self):
        #str_run_start_date=self.RunStartDate.strftime("%I:%M%p on %B %d, %Y")
        run_parameters_data = []
        if self.ImageChannel == None:
            img_channel = 'None'
        else:
            img_channel=self.ImageChannel.strip('[').strip(']').replace("'","")

        if self.ImageDimensions == None:
            image_dimensions = ['None']
        else:
            image_dimensions = self.ImageDimensions.strip('{').strip('}').replace("'","").split(',')

        flowcell_layout = self.FlowcellLayout.strip('{').strip('}').replace("'","").split(',')

        run_parameters_data.append(self.RunID); run_parameters_data.append(self.ExperimentName)
        run_parameters_data.append(self.RTAVersion); run_parameters_data.append(self.SystemSuiteVersion)
        run_parameters_data.append(self.LibraryID); run_parameters_data.append(self.Chemistry)
        run_parameters_data.append(self.RunStartDate); run_parameters_data.append(self.AnalysisWorkflowType)
        run_parameters_data.append(self.RunManagementType); run_parameters_data.append(self.PlannedRead1Cycles)
        run_parameters_data.append(self.PlannedRead2Cycles); run_parameters_data.append(self.PlannedIndex1ReadCycles)
        run_parameters_data.append(self.PlannedIndex2ReadCycles); run_parameters_data.append(self.ApplicationVersion)
        run_parameters_data.append(self.NumTilesPerSwath); run_parameters_data.append(img_channel)
        run_parameters_data.append(self.Flowcell); run_parameters_data.append(image_dimensions)
        run_parameters_data.append(flowcell_layout)

        return run_parameters_data


    def get_number_of_reads (self):
        count = 0
        if self.PlannedRead1Cycles != "0" :
            count +=1
        if self.PlannedRead2Cycles != "0" :
            count +=1
        if self.PlannedIndex1ReadCycles != "0" :
            count +=1
        if self.PlannedIndex2ReadCycles != "0" :
            count +=1
        return count

    def get_number_of_cycles (self):

        number_of_cycles = int(self.PlannedRead1Cycles) + int(self.PlannedRead2Cycles) + int(self.PlannedIndex1ReadCycles) + int(self.PlannedIndex2ReadCycles)
        return number_of_cycles

    def get_run_folder (self):
        return '%s' %(self.RunID)

    objects = RunningParametersManager ()


class StatsRunSummaryManager (models.Manager):

    def create_stats_run_summary (self, stats_run_summary , experiment_name) :
        run_process = RunProcess.objects.get(runName__exact = experiment_name)
        s_run_summary = self.create(runprocess_id = run_process,
                                level=stats_run_summary['level'], yieldTotal = stats_run_summary['yieldTotal'],
                                projectedTotalYield= stats_run_summary['projectedTotalYield'],
                                aligned= stats_run_summary['aligned'], errorRate= stats_run_summary['errorRate'],
                                intensityCycle= stats_run_summary['intensityCycle'], biggerQ30= stats_run_summary['biggerQ30'],
                                stats_summary_run_date = run_process.run_date)
        return s_run_summary

class StatsRunSummary (models.Model):
    runprocess_id = models.ForeignKey(
            RunProcess,
            on_delete=models.CASCADE)
    level= models.CharField(max_length=20)
    yieldTotal= models.CharField(max_length=10)
    projectedTotalYield = models.CharField(max_length=10)
    aligned= models.CharField(max_length=10)
    errorRate=models.CharField(max_length=10)
    intensityCycle= models.CharField(max_length=10)
    biggerQ30= models.CharField(max_length=10)
    generatedat = models.DateTimeField(auto_now_add=True)
    stats_summary_run_date = models.DateField(auto_now = False, null=True)

    def __str__(self):
        return '%s' %(self.level)

    def get_bin_run_summary(self):
        return '%s;%s;%s;%s;%s;%s' %( self.yieldTotal,
                self.projectedTotalYield, self.aligned, self.errorRate,
                self.intensityCycle, self.biggerQ30)

    objects = StatsRunSummaryManager ()

class StatsRunReadManager (models.Manager):

    def create_stats_run_read (self, stats_run_read, experiment_name):
        run_process = RunProcess.objects.get(runName__exact = experiment_name)
        s_run_read = self.create (runprocess_id = run_process,
                                    read= stats_run_read['read'], lane = stats_run_read['lane'] ,
                                    tiles= stats_run_read['tiles'], density= stats_run_read['density'] ,
                                    cluster_PF= stats_run_read['cluster_PF'], phas_prephas= stats_run_read['phas_prephas'],
                                    reads= stats_run_read['reads'], reads_PF= stats_run_read['reads_PF'],
                                    q30= stats_run_read['q30'], yields= stats_run_read['yields'],
                                    cyclesErrRated= stats_run_read['cyclesErrRated'], aligned= stats_run_read['aligned'],
                                    errorRate= stats_run_read['errorRate'], errorRate35= stats_run_read['errorRate35'],
                                    errorRate50= stats_run_read['errorRate50'] , errorRate75= stats_run_read['errorRate75'] ,
                                    errorRate100= stats_run_read['errorRate100'] , intensityCycle= stats_run_read['intensityCycle'] ,
                                    stats_read_run_date = run_process.run_date)

        return s_run_read

class StatsRunRead (models.Model):
    runprocess_id = models.ForeignKey(
            RunProcess,
            on_delete=models.CASCADE)
    read= models.CharField(max_length=10)
    lane = models.CharField(max_length=10)
    tiles = models.CharField(max_length=10)
    density = models.CharField(max_length=40)
    cluster_PF = models.CharField(max_length=40)
    phas_prephas = models.CharField(max_length=40)
    reads = models.CharField(max_length=40)
    reads_PF = models.CharField(max_length=40)
    q30 = models.CharField(max_length=40)
    yields = models.CharField(max_length=40)
    cyclesErrRated = models.CharField(max_length=40)
    aligned = models.CharField(max_length=40)
    errorRate = models.CharField(max_length=40)
    errorRate35 = models.CharField(max_length=40)
    errorRate50 = models.CharField(max_length=40)
    errorRate75 = models.CharField(max_length=40)
    errorRate100 = models.CharField(max_length=40)
    intensityCycle = models.CharField(max_length=40)
    generatedat = models.DateTimeField(auto_now_add=True)
    stats_read_run_date = models.DateField(auto_now = False, null=True)

    def __str__(self):
        return '%s' %(self.read)

    def get_bin_run_read(self):
        return '%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s' %( self.lane,self.tiles, self.density,
                self.cluster_PF, self.phas_prephas, self.reads, self.reads_PF,
                self.q30, self.yields, self.cyclesErrRated, self.aligned,
                self.errorRate, self.errorRate35, self.errorRate50, self.errorRate75,
                self.errorRate100, self.intensityCycle)

    objects = StatsRunReadManager ()


class RawDemuxStatsManager(models.Manager) :
    def create_stats_run_read (self, raw_demux_stats, run_object_name):

        raw_stats = self.create(runprocess_id = run_object_name,
                                project_id = raw_demux_stats ['project_id'], defaultAll = raw_demux_stats ['defaultAll'],
                                rawYield= raw_demux_stats ['rawYield'], rawYieldQ30= raw_demux_stats ['rawYieldQ30'],
                                rawQuality= raw_demux_stats ['rawQuality'], PF_Yield= raw_demux_stats ['PF_Yield'],
                                PF_YieldQ30= raw_demux_stats ['PF_YieldQ30'], PF_QualityScore =raw_demux_stats ['PF_QualityScore'] )

        return raw_stats

class RawDemuxStats (models.Model):
    runprocess_id = models.ForeignKey(
            RunProcess,
            on_delete=models.CASCADE,)
    project_id = models.ForeignKey(Projects, on_delete=models.CASCADE, null=True)
    defaultAll =models.CharField(max_length=40, null=True)
    rawYield = models.CharField(max_length=255)
    rawYieldQ30= models.CharField(max_length=255)
    rawQuality= models.CharField(max_length=255)
    PF_Yield= models.CharField(max_length=255)
    PF_YieldQ30= models.CharField(max_length=255)
    PF_QualityScore= models.CharField(max_length=255)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return '%s' %(self.runprocess_id)

    def get_raw_xml_stats(self):
        return '%s;%s;%s;%s;%s;%s;%s;%s;%s;%s' %(self.project, self.barcodeCount,
                self.perfectBarcodeCount, self.rawYield, self.PF_Yield,
                self.rawYieldQ30, self.PF_YieldQ30, self.rawQuality,
                self.PF_Quality, self.sampleNumber)


    objects = RawDemuxStatsManager ()

class RawTopUnknowBarcodesManager(models.Manager) :
    def create_unknow_barcode (self, unknow_barcode):
        unknow_barcode = self.create (runprocess_id = unknow_barcode['runprocess_id'] , lane_number = unknow_barcode['lane_number'] ,
                                top_number =  unknow_barcode['top_number'], count = unknow_barcode['count'],
                                sequence = unknow_barcode['sequence']   )

        return unknow_barcode


class RawTopUnknowBarcodes (models.Model):
    runprocess_id = models.ForeignKey(
            RunProcess,
            on_delete=models.CASCADE)
    lane_number= models.CharField(max_length=4)
    top_number = models.CharField(max_length=4)
    count= models.CharField(max_length=40)
    sequence = models.CharField(max_length=40)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__ (self):
        return '%s' %(self.lane_number)

    def get_unknow_barcodes (self):
        return '%s;%s' %(self.count, self.sequence)

    objects = RawTopUnknowBarcodesManager ()

class StatsFlSummaryManager (models.Manager) :
    def create_fl_summary (self, stats_fl_summary):
        fl_summary = self.create (runprocess_id = stats_fl_summary ['runprocess_id'], project_id = stats_fl_summary ['project_id'],
                            defaultAll = stats_fl_summary ['defaultAll'], flowRawCluster = stats_fl_summary ['flowRawCluster'],
                            flowPfCluster = stats_fl_summary ['flowPfCluster'], flowYieldMb = stats_fl_summary ['flowYieldMb'],
                            sampleNumber = stats_fl_summary ['sampleNumber'] )

        return fl_summary

class StatsFlSummary(models.Model):
    runprocess_id = models.ForeignKey(
            RunProcess,
            on_delete=models.CASCADE)
    project_id = models.ForeignKey(Projects, on_delete=models.CASCADE, null=True, blank=True)
    defaultAll =models.CharField(max_length=40, null=True)
    flowRawCluster = models.CharField(max_length=40)
    flowPfCluster= models.CharField(max_length=40)
    flowYieldMb= models.CharField(max_length=40)
    sampleNumber= models.CharField(max_length=40)
    generated_at = models.DateTimeField(auto_now_add=True)

    def get_fl_summary(self):
        data = []
        data.append(self.flowRawCluster)
        data.append(self.flowPfCluster)
        data.append(self.flowYieldMb)
        data.append(self.sampleNumber)
        return data

    objects = StatsFlSummaryManager ()

class StatsLaneSummaryManager (models.Manager) :
    def create_lane_summary (self, l_summary) :
        lane_summary = self.create (runprocess_id = l_summary['runprocess_id'], project_id = l_summary['project_id'],
                            defaultAll = l_summary['defaultAll'], lane = l_summary['lane'],
                            pfCluster = l_summary['pfCluster'], percentLane = l_summary['percentLane'],
                            perfectBarcode = l_summary['perfectBarcode'], oneMismatch = l_summary['oneMismatch'],
                            yieldMb = l_summary['yieldMb'], biggerQ30 = l_summary['biggerQ30'],
                            meanQuality = l_summary['meanQuality'] )

        return lane_summary

class StatsLaneSummary (models.Model):
    runprocess_id = models.ForeignKey(
            RunProcess,
            on_delete=models.CASCADE)
    project_id = models.ForeignKey(Projects, on_delete=models.CASCADE,  null=True, blank=True)
    defaultAll =models.CharField(max_length=40, null=True)
    lane= models.CharField(max_length=10)
    pfCluster = models.CharField(max_length=64)
    percentLane = models.CharField(max_length=64)
    perfectBarcode = models.CharField(max_length=64)
    oneMismatch = models.CharField(max_length=64)
    yieldMb = models.CharField(max_length=64)
    biggerQ30 = models.CharField(max_length=64)
    meanQuality = models.CharField(max_length=64)
    generated_at = models.DateTimeField(auto_now_add=True)


    def __str__(self):
        return '%s' %(self.runprocess_id)

    def get_flow_cell_summary (self):
        return '%s' %(self.runprocess_id)

    def get_lane_summary(self):
        data = []
        data.append(self.lane)
        data.append(self.pfCluster)
        data.append(self.percentLane)
        data.append(self.perfectBarcode)
        data.append(self.oneMismatch)
        data.append(self.yieldMb)
        data.append(self.biggerQ30)
        data.append(self.meanQuality)
        return data

    def get_stats_info (self):
        stats_info = []
        stats_info.append(self.biggerQ30)
        stats_info.append(self.meanQuality)
        stats_info.append(self.yieldMb)
        stats_info.append(self.pfCluster)
        return stats_info

    objects = StatsLaneSummaryManager ()


class GraphicsStats (models.Model):
    runprocess_id = models.ForeignKey(
            RunProcess,
            on_delete=models.CASCADE)
    folderRunGraphic= models.CharField(max_length=255)
    cluserCountGraph= models.CharField(max_length=255)
    flowCellGraph= models.CharField(max_length=255)
    intensityByCycleGraph= models.CharField(max_length=255)
    heatMapGraph= models.CharField(max_length=255)
    histogramGraph= models.CharField(max_length=255)
    sampleQcGraph= models.CharField(max_length=255)
    ## Fix 24/07/2018: null=True added at definition of "generated_at" to allow json load
    generated_at = models.DateTimeField(auto_now_add=True, null=True)

    def __str__(self):
        return '%s' %(self.folderRunGraphic)

    def get_graphics(self):
        return '%s;%s;%s;%s;%s;%s'%(self.cluserCountGraph,
                self.flowCellGraph, self.intensityByCycleGraph,
                self.heatMapGraph, self.histogramGraph,
                self.sampleQcGraph)

    def get_folder_graphic(self):
        return '%s' %(self.folderRunGraphic)

class SamplesInProjectManager (models.Manager):
    def create_sample_project (self, s_project):
        sample_project = self.create( project_id = s_project['project_id'] , sampleName =  s_project['sampleName'],
                            barcodeName =  s_project['barcodeName'], pfClusters =  s_project['pfClusters'],
                            percentInProject =  s_project['percentInProject'], yieldMb =  s_project['yieldMb'],
                            qualityQ30 =  s_project['qualityQ30'], meanQuality =  s_project['meanQuality'] )

class SamplesInProject (models.Model):
    project_id = models.ForeignKey(
                Projects,
                on_delete= models.CASCADE)
    runProcess_id = models.ForeignKey(
                RunProcess,
                on_delete = models.CASCADE, null = True, blank = True)
    user_id = models.ForeignKey(
                User,
                on_delete=models.CASCADE, null = True, blank = True)
    sampleName = models.CharField(max_length=255)
    barcodeName = models.CharField(max_length=255)
    pfClusters = models.CharField(max_length=55)
    percentInProject = models.CharField(max_length=25)
    yieldMb = models.CharField(max_length=55)
    qualityQ30 = models.CharField(max_length=55)
    meanQuality = models.CharField(max_length=55)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__ (self):
        return '%s' %(self.sampleName)

    def get_info_for_searching(self):
        sample_info = []
        sample_info.append(self.pk)
        sample_info.append(self.sampleName)
        sample_info.append(self.project_id.get_project_name())
        sample_info.append(self.runProcess_id.get_run_name())
        sample_info.append(self.generated_at.strftime("%I:%M%p on %B %d, %Y"))
        return sample_info

    def get_sample_id(self):
        return '%s' %(self.id)

    def get_sample_information(self):
        data = []
        data.append(self.sampleName)
        data.append(self.barcodeName)
        data.append(self.pfClusters)
        data.append(self.percentInProject)
        data.append(self.yieldMb)
        data.append(self.qualityQ30)
        data.append(self.meanQuality)
        return data

    def get_sample_information_with_project_run(self):
        data = []
        data.append(self.sampleName)
        data.append(self.project_id.get_project_name())
        data.append(self.runProcess_id.get_run_name())
        data.append(self.barcodeName)
        data.append(self.pfClusters)
        data.append(self.percentInProject)
        data.append(self.yieldMb)
        data.append(self.qualityQ30)
        data.append(self.meanQuality)
        return data

    def get_sample_name(self):
        return '%s' %(self.sampleName)

    def get_project_id (self) :
        return '%s' %(self.project_id.pk)

    def get_project_name (self) :
        #p_id = self.project_id
        #project_name =Projects.objects.get(projectName=p_id).get_project_name()
        project_name = self.project_id.get_project_name()
        #Projects.objects.prefetch_related('user_id').filter(user_id = user_id)
        return '%s' %(project_name)

    def get_run_name(self):
        return '%s' %(self.runProcess_id.get_run_name())

    def get_sequencer_used(self):
        return '%s' %(self.runProcess_id.get_run_used_sequencer_name())

    def get_run_obj(self):
        return self.runProcess_id

    def get_run_id(self):
        return self.runProcess_id.pk

    def get_quality_sample (self):
        return '%s' %(self.qualityQ30)

    objects = SamplesInProjectManager()


################################################
################################################
# New objets for version 2.0.0
################################################


class CollectionIndexKit (models.Model):
    collectionIndexName = models.CharField(max_length=125)
    version = models.CharField(max_length=80,null=True)
    plateExtension = models.CharField(max_length=125, null=True)
    adapter1 = models.CharField(max_length=125,null=True)
    adapter2 = models.CharField(max_length=125, null=True)
    collectionIndexFile =  models.FileField(upload_to = wetlab_config.COLLECTION_INDEX_KITS_DIRECTORY )
    generatedat = models.DateTimeField(auto_now_add=True, null=True)

    def __str__ (self):
        return '%s' %(self.collectionIndexName)

    def get_collection_index_name(self):
        return '%s' %(self.collectionIndexName)

    def get_id(self):
        return '%s' %(self.pk)

    def get_collection_index_information (self):
        if self.adapter2 =='':
            adapter2 = 'Not used on this collection'
        else:
            adapter2 = self.adapter2
        collection_info = []
        collection_info.append(self.collectionIndexName)
        collection_info.append(self.version)
        collection_info.append(self.plateExtension)
        collection_info.append(self.adapter1)
        collection_info.append(self.adapter2)
        collection_info.append(self.collectionIndexFile)
        return collection_info




class CollectionIndexValues (models.Model):
    collectionIndexKit_id = models.ForeignKey(
        CollectionIndexKit,
        on_delete=models.CASCADE)

    defaultWell = models.CharField(max_length=10,null=True)
    index_7 = models.CharField(max_length=25,null=True)
    i_7_seq = models.CharField(max_length=25,null=True)
    index_5 = models.CharField(max_length=25,null=True)
    i_5_seq = models.CharField(max_length=25,null=True)


    def get_index_value_information (self):
        return '%s;%s;%s;%s;%s' %(self.defaultWell, self.index_7, self.i_7_seq,
                        self.index_5, self.i_5_seq)

    def get_collection_index_id(self):
        return '%s' %(self.collectionIndexKit_id.get_id())

    def get_collection_index_name(self):
        return '%s' %(self.collectionIndexKit_id.get_collection_index_name())

class libPreparationUserSampleSheetManager (models.Manager):

    def create_lib_prep_user_sample_sheet (self, user_sample_sheet_data):
        register_user_obj = User.objects.get(username__exact = user_sample_sheet_data['user'])
        if user_sample_sheet_data['index_adapter'] == '':
            collection_index_kit_id = None
        else:
            collection_index_kit_id = CollectionIndexKit.objects.get(collectionIndexName__exact = user_sample_sheet_data['index_adapter'])
        file_name =  os.path.basename(user_sample_sheet_data['file_name'])
        configuration = SequencingConfiguration.objects.filter( platformID__platformName__exact = user_sample_sheet_data['platform'], configurationName__exact = user_sample_sheet_data['configuration']).last()
        new_lib_prep_user_sample_sheet = self.create(registerUser = register_user_obj,
                    collectionIndexKit_id  = collection_index_kit_id, reads = ','.join(user_sample_sheet_data['reads']),
                    sampleSheet = file_name, application = user_sample_sheet_data['application'],
                    instrument = user_sample_sheet_data ['instrument'], assay = user_sample_sheet_data['assay'],
                    adapter1 = user_sample_sheet_data['adapter1'], adapter2 = user_sample_sheet_data['adapter2'],
                    sequencingConfiguration = configuration)
        return new_lib_prep_user_sample_sheet

class libPreparationUserSampleSheet (models.Model):
    registerUser = models.ForeignKey(
            User,
            on_delete=models.CASCADE)

    collectionIndexKit_id = models.ForeignKey(
                CollectionIndexKit,
                on_delete= models.CASCADE, null = True)

    sequencingConfiguration = models.ForeignKey(
                SequencingConfiguration,
                on_delete= models.CASCADE, null = True)

    sampleSheet = models.FileField(upload_to = wetlab_config.LIBRARY_PREPARATION_SAMPLE_SHEET_DIRECTORY)
    generatedat = models.DateTimeField(auto_now_add=True, null=True)
    application = models.CharField(max_length=70, null = True, blank = True)
    instrument = models.CharField(max_length=70, null = True, blank = True)
    adapter1 = models.CharField(max_length=70, null = True, blank = True)
    adapter2 = models.CharField(max_length=70, null = True, blank = True)
    assay = models.CharField(max_length=70, null = True, blank = True)
    reads = models.CharField(max_length=10, null = True, blank = True)
    confirmedUsed = models.BooleanField(default = False)

    def __str__ (self):
        return '%s' %(self.sampleSheet)

    def get_assay(self):
        return '%s' %(self.assay)

    def get_adapters(self):
        adapters = []
        adapters.append(self.adapter1)
        adapters.append(self.adapter2)
        return adapters



    def get_collection_index_kit (self):
        return '%s' %(self.collectionIndexKit_id.get_collection_index_name())

    def get_sequencing_configuration_platform (self):
        return '%s'  %(self.sequencingConfiguration.get_platform_name())
    def get_user_sample_sheet_id (self):
        return '%s' %(self.pk)

    def get_all_data(self):
        s_s_data = []
        s_s_data.append(self.collectionIndexKit_id.get_collection_index_name())
        s_s_data.append(self.application)
        s_s_data.append(self.instrument)
        s_s_data.append(self.adapter1)
        s_s_data.append(self.adapter2)
        s_s_data.append(self.assay)
        s_s_data.append(self.reads.split(',')[0])
        return s_s_data

    def update_confirm_used (self, confirmation):
        self.confirmedUsed = confirmation
        self.save()
        return self

    objects = libPreparationUserSampleSheetManager()

class StatesForPool (models.Model):
    poolState = models.CharField(max_length=50)

    def __str__ (self):
        return '%s' %(self.poolState)

    def get_pool_state(self):
        return '%s' %(self.poolState)

class LibraryPoolManager (models.Manager):
    def create_lib_pool (self, pool_data):
        platform_obj = SequencingPlatform.objects.filter(platformName__exact = pool_data['platform']).last()
        new_library_pool = self.create(registerUser = pool_data['registerUser']  ,
                    poolState = StatesForPool.objects.get(poolState__exact = 'Defined'),
                    poolName = pool_data['poolName'], poolCodeID = pool_data['poolCodeID'],
                    adapter = pool_data['adapter'],  pairedEnd = pool_data['pairedEnd'],
                    numberOfSamples = pool_data['n_samples'], platform = platform_obj)
        return new_library_pool


class LibraryPool (models.Model):
    registerUser = models.ForeignKey(
            User,
            on_delete=models.CASCADE)
    poolState = models.ForeignKey(
                StatesForPool,
                on_delete= models.CASCADE)
    runProcess_id = models.ForeignKey(
            RunProcess,
            on_delete = models.CASCADE, null = True, blank = True)

    platform = models.ForeignKey(
            SequencingPlatform,
            on_delete = models.CASCADE, null = True, blank = True)
    poolName = models.CharField(max_length=50)
    #sampleSheet = models.FileField(upload_to = wetlab_config.RUN_SAMPLE_SHEET_DIRECTORY)
    #experiment_name = models.CharField(max_length=50)
    #plateName = models.CharField(max_length=50)
    #containerID = models.CharField(max_length=50, null =True, blank = True)
    #libUsedInBaseSpace = models.CharField(max_length=50)
    numberOfSamples = models.IntegerField(default=0)
    poolCodeID = models.CharField(max_length=50, blank = True)
    adapter = models.CharField(max_length=50, null = True, blank = True)
    pairedEnd = models.CharField(max_length=10, null = True, blank = True)
    #assay = models.CharField(max_length=50, null = True, blank = True)
    #collectionIndex = models.CharField(max_length=50, null = True, blank = True)
    generated_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        ordering = ('poolName',)

    def __str__ (self):
        return '%s' %(self.poolName)

    def get_adapter(self):
        return '%s' %(self.adapter)

    def get_id(self):
        return '%s' %(self.pk)

    def get_info(self):
        pool_info =[]
        pool_info.append(self.poolName)
        pool_info.append(self.poolCodeID)
        pool_info.append(self.numberOfSamples)
        return pool_info

    def get_number_of_samples (self):
        return '%s' %(self.numberOfSamples)

    def get_pool_name(self):
        return '%s' %(self.poolName)

    def get_pool_code_id(self):
        return '%s' %(self.poolCodeID)

    def get_pool_single_paired(self):
        return '%s' %(self.pairedEnd)

    def get_platform_name(self):
        return '%s' %(self.platform.get_platform_name())
    def get_run_name(self):
        if self.runProcess_id != None :
            return '%s' %(self.runProcess_id.get_run_name())
        else:
            return 'Not defined yet'

    def get_run_id(self):
        return '%s' %(self.runProcess_id.get_run_id())


    def set_pool_state(self, state):
        self.poolState = StatesForPool.objects.get(poolState__exact = state)
        self.save()


    def update_number_samples(self, number_s_in_pool):
        self.numberOfSamples = number_s_in_pool
        self.save()
        return self

    def update_run_name (self, run_name):
        self.runProcess_id = run_name
        self.save()
        return self

    objects = LibraryPoolManager()



class StatesForLibraryPreparation (models.Model):
    libPrepState = models.CharField(max_length=50)

    def __str__ (self):
        return '%s' %(self.libPrepState)

    def get_lib_state(self):
        return '%s' %(self.libPrepState)



class libraryPreparationManager(models.Manager):
    def create_lib_preparation (self, lib_prep_data):
        #import pdb; pdb.set_trace()
        registerUser_obj = User.objects.get(username__exact  = lib_prep_data['registerUser'])
        sample_obj = Samples.objects.get(pk__exact = lib_prep_data['sample_id'])
        molecule_obj = MoleculePreparation.objects.get(pk__exact = lib_prep_data['molecule_id'])
        lib_state_obj = StatesForLibraryPreparation.objects.get(libPrepState__exact =  'Defined')
        new_lib_prep = self.create(registerUser = registerUser_obj, molecule_id = molecule_obj, sample_id = sample_obj,
            protocol_id =   lib_prep_data['protocol_obj'], libPrepState = lib_state_obj,
            libPrepCodeID = lib_prep_data['lib_prep_code_id'], userSampleID = lib_prep_data['user_sampleID'],
            uniqueID = lib_prep_data['uniqueID'])
        return new_lib_prep

    def create_reused_lib_preparation (self, reg_user, molecule_obj, sample_id):
        lib_state = StatesForLibraryPreparation.objects.get(libPrepState =  'Created for Reuse')
        new_library_preparation = self.create(registerUser = reg_user, molecule_id = molecule_obj, sample_id = sample_id,libPrepState = lib_state)
        return new_library_preparation

class LibraryPreparation (models.Model):
    registerUser = models.ForeignKey(
            User,
            on_delete=models.CASCADE)
    molecule_id = models.ForeignKey(
                MoleculePreparation,
                on_delete= models.CASCADE)
    sample_id = models.ForeignKey(
                Samples,
                on_delete= models.CASCADE, null = True, blank = True)
    protocol_id = models.ForeignKey(
                Protocols,
                on_delete= models.CASCADE, null = True, blank = True)
    libPrepState = models.ForeignKey(
                StatesForLibraryPreparation,
                on_delete= models.CASCADE, null = True, blank = True)

    userLotKit_id = models.ForeignKey(
                UserLotCommercialKits,
                on_delete= models.CASCADE, null = True, blank = True)

    user_sample_sheet = models.ForeignKey(
                libPreparationUserSampleSheet,
                on_delete= models.CASCADE, null = True, blank = True)

    pools = models.ManyToManyField(LibraryPool, blank = True)

    libPrepCodeID = models.CharField(max_length=255, null = True, blank = True)
    userSampleID = models.CharField(max_length =20 , null = True, blank = True)
    projectInSampleSheet = models.CharField(max_length =50, null = True, blank = True)
    samplePlate = models.CharField(max_length =50, null = True, blank = True)
    sampleWell = models.CharField(max_length =20, null = True, blank = True)
    indexPlateWell = models.CharField(max_length =20, null = True, blank = True)
    i7IndexID = models.CharField(max_length =16, null = True, blank = True)
    i7Index = models.CharField(max_length =16, null = True, blank = True)
    i5IndexID = models.CharField(max_length =16, null = True, blank = True)
    i5Index = models.CharField(max_length =16, null = True, blank = True)
    ## Miseq fields
    genomeFolder = models.CharField(max_length =80, null = True, blank = True)
    manifest = models.CharField(max_length =80, null = True, blank = True)
    ##### End Miseq fields
    #singlePairedEnd  = models.CharField(max_length =20, null = True, blank = True)
    #lengthRead = models.CharField(max_length =5, null = True, blank = True)
    numberOfReused = models.IntegerField(default=0)
    uniqueID = models.CharField(max_length =16, null = True, blank = True)
    userInSampleSheet = models.CharField(max_length=255, null = True, blank = True)

    class Meta:
        ordering = ('libPrepCodeID',)


    def __str__ (self):
        return '%s' %(self.libPrepCodeID)

    def get_adapters(self):
        return self.user_sample_sheet.get_adapters()

    '''
    def get_collection_index_name (self):
        return '%s' %(self.collectionIndex_id.get_collection_index_name())
    '''

    def get_id (self):
        return '%s' %(self.pk)

    def get_info_for_display(self):
        lib_info = []
        lib_info.append(self.libPrepCodeID)
        lib_info.append(self.molecule_id.get_molecule_code_id())
        lib_info.append(self.libPrepState.libPrepState)
        lib_info.append(self.protocol_id.get_name())
        lib_info.append(self.projectInSampleSheet)
        lib_info.append(self.i7IndexID)
        lib_info.append(self.i5IndexID)
        lib_info.append(self.numberOfReused)
        return lib_info

    def get_info_for_selection_in_pool(self):
        lib_info = []
        lib_info.append(self.libPrepCodeID)
        lib_info.append(self.sample_id.get_sample_name())
        lib_info.append(self.registerUser.username)
        lib_info.append(self.user_sample_sheet.get_collection_index_kit())
        lib_info.append(self.i7IndexID)
        lib_info.append(self.i5IndexID)
        lib_info.append(self.pk)
        return lib_info



    def get_info_for_run_paired_end(self):
        lib_info = []
        lib_info.append(self.uniqueID)
        lib_info.append(self.sample_id.get_sample_name())
        lib_info.append(self.samplePlate)
        lib_info.append(self.sampleWell)
        lib_info.append(self.indexPlateWell)
        lib_info.append(self.i7IndexID)
        lib_info.append(self.i7Index)
        lib_info.append(self.i5IndexID)
        lib_info.append(self.i5Index)
        lib_info.append(self.projectInSampleSheet)
        lib_info.append(self.userInSampleSheet)
        #lib_info.append(self.collectionIndex_id.get_collection_index_name())
        return lib_info


    def get_info_for_run_single_read(self):
        lib_info = []
        lib_info.append(self.uniqueID)
        lib_info.append(self.sample_id.get_sample_name())
        lib_info.append(self.samplePlate)
        lib_info.append(self.sampleWell)
        lib_info.append(self.indexPlateWell)
        lib_info.append(self.i7IndexID)
        lib_info.append(self.i7Index)
        lib_info.append(self.projectInSampleSheet)
        lib_info.append(self.userInSampleSheet)
        #lib_info.append(self.collectionIndex_id.get_collection_index_name())
        return lib_info

    def get_basic_data(self):
        lib_info = []
        lib_info.append(self.sample_id.get_sample_name())
        lib_info.append(self.uniqueID)
        lib_info.append(self.protocol_id.get_name())
        lib_info.append(self.pk)
        return lib_info

    def get_info_for_display_pool (self):
        if self.registerUser.username == None:
            user_name = ''
        else :
            user_name = self.registerUser.username

        lib_info = []
        lib_info.append(self.libPrepCodeID)
        lib_info.append(self.get_sample_name())
        lib_info.append(user_name)
        lib_info.append(self.get_sample_id())
        return lib_info

    def get_collection_index_kit (self):
        return '%s' %(self.collectionIndex_id.get_collection_index_name())

    def get_i7_index(self):
        return '%s' %(self.i7IndexID)

    def get_i5_index(self):
        return '%s' %(self.i5IndexID)

    def get_indexes (self):
        index = {}
        if self.i5IndexID == None:
            index['i5_indexID'] = ''
            index['i5_seq'] = ''
        else:
            index['i5_indexID']= self.i5IndexID
            index['i5_seq'] = self.i5Index
        index['i7_indexID'] = self.i7IndexID
        index['i7_seq'] = self.i7Index
        return index

    def get_lib_prep_code (self):
        return '%s' %(self.libPrepCodeID)

    def get_lib_prep_id (self):
        return '%s' %(self.pk)

    def get_molecule_code_id(self):
        return '%s' %(self.molecule_id.get_molecule_code_id())

    def get_molecule_obj(self):
        return self.molecule_id

    def get_protocol_used (self):
        return '%s'  %(self.protocol_id.get_name())

    def get_protocol_obj(self):
        return self.protocol_id

    def get_protocol_id (self):
        return '%s' %(self.protocol_id.pk)

    def get_reagents_kit_used(self):
        return '%s' %(self.reagent_id.get_nick_name())

    def get_reused_value(self):
        return '%s' %(self.numberOfReused)

    def get_sample_name(self):
        return '%s' %(self.sample_id.get_sample_name())

    def get_sample_id(self):
        return self.sample_id.pk

    def get_sample_obj(self):
        return self.sample_id


    def get_state(self):
        return '%s' %(self.libPrepState.libPrepState)

    def get_unique_id(self):
        return '%s' %(self.uniqueID)

    def get_user_obj(self):
        return self.registerUser

    def get_user_sample_sheet_obj (self):
        return self.user_sample_sheet

    def set_increase_reuse(self):
        self.numberOfReused += 1
        self.save()
        return

    def set_pool (self , pool_obj):
        self.pools.add(pool_obj)
        self.save()
        return

    def set_state(self , state_value):
        self.libPrepState = StatesForLibraryPreparation.objects.get(libPrepState__exact = state_value)
        self.save()
        return
    '''
    def set_reagent_user_kit(self, kit_value):
        self.user_reagentKit_id = UserLotCommercialKits.objects.get(nickName__exact = kit_value)
        self.save()
        return
    '''
    def update_i7_index(self, i7_seq_value):
        self.i7Index = i7_seq_value
        self.save()
        return

    def update_i5_index(self, i5_seq_value):
        self.i5Index = i5_seq_value
        self.save()
        return

    def update_library_preparation_with_indexes (self,lib_prep_data ):
        self.user_sample_sheet = lib_prep_data['user_sample_sheet']
        self.projectInSampleSheet = lib_prep_data['projectInSampleSheet']
        self.samplePlate = lib_prep_data['samplePlate']
        self.sampleWell = lib_prep_data['sampleWell']
        self.i7IndexID = lib_prep_data['i7IndexID']
        self.i7Index = lib_prep_data['i7Index']
        self.i5IndexID = lib_prep_data['i5IndexID']
        self.i5Index = lib_prep_data['i5Index']
        self.indexPlateWell = lib_prep_data['indexPlateWell']
        self.genomeFolder = lib_prep_data['genomeFolder']
        self.manifest = lib_prep_data['manifest']
        self.userInSampleSheet = lib_prep_data['userInSampleSheet']
        self.userSampleID = lib_prep_data['userSampleID']
        self.save()
        return self


    def update_lib_preparation_info_in_reuse_state (self, lib_prep_data):
        # Update the instance with the sample sheet information
        self.registerUser = User.objects.get(username__exact  = lib_prep_data['registerUser'])
        self.protocol_id =   lib_prep_data['protocol_obj']
        lib_state = StatesForLibraryPreparation.objects.get(libPrepState__exact =  'Defined')
        self.user_sample_sheet = lib_prep_data['user_sample_sheet']
        self.libPrepCodeID = lib_prep_data['lib_prep_code_id']
        self.projectInSampleSheet = lib_prep_data['projectInSampleSheet']
        self.userSampleID = lib_prep_data['userSampleID']
        self.samplePlate = lib_prep_data['samplePlate']
        self.sampleWell = lib_prep_data['sampleWell']
        self.indexPlateWell =  lib_prep_data['indexPlateWell']
        self.i7IndexID = lib_prep_data['i7IndexID']
        self.i7Index = lib_prep_data['i7Index']
        self.i5IndexID = lib_prep_data['i5IndexID']
        self.i5Index = lib_prep_data['i5Index']
        self.uniqueID = lib_prep_data['uniqueID']

        self.save()
        return self

    objects = libraryPreparationManager()


class LibParameterValueManager (models.Manager):
    def create_library_parameter_value (self, parameter_value):
        new_parameter_data = self.create(parameter_id = parameter_value['parameter_id'],
                library_id = parameter_value['library_id'],
                parameterValue = parameter_value['parameterValue'])
        return new_parameter_data

class LibParameterValue (models.Model):
    parameter_id = models.ForeignKey(
                    ProtocolParameters,
                    on_delete= models.CASCADE)
    library_id = models.ForeignKey(
                LibraryPreparation,
                on_delete= models.CASCADE)
    parameterValue = models.CharField(max_length=255)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__ (self):
        return '%s' %(self.parameterValue)

    def get_parameter_information (self):
        return '%s' %(self.parameterValue)

    objects = LibParameterValueManager()


class AdditionaKitsLibraryPreparationManager(models.Manager):

    def create_additional_kit (self, kit_data):
        new_additional_kit = self.create(registerUser = kit_data['user'],
                protocol_id = kit_data['protocol_id'], commercialKit_id = kit_data['commercialKit_id'],
                kitName = kit_data['kitName'], description = kit_data['description'],
                kitOrder = kit_data['kitOrder'], kitUsed = kit_data['kitUsed'])
        return new_additional_kit

class AdditionaKitsLibraryPreparation (models.Model):
    registerUser = models.ForeignKey(
                    User,
                    on_delete=models.CASCADE)
    protocol_id =  models.ForeignKey(
                    Protocols,
                    on_delete= models.CASCADE)
    commercialKit_id =  models.ForeignKey(
                    CommercialKits,
                    on_delete= models.CASCADE)
    kitName = models.CharField(max_length=255)
    description = models.CharField(max_length= 400, null=True, blank=True)
    kitOrder = models.IntegerField()
    kitUsed = models.BooleanField()
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__ (self):
        return '%s' %(self.kitName)

    def get_all_kit_info(self):
        data = []
        data.append(self.kitName)
        data.append(self.kitOrder)
        data.append(self.kitUsed)
        data.append(self.commercialKit_id.get_name())
        data.append(self.description)
        return data

    def get_kit_name (self):
        return '%s' %(self.kitName)

    def get_commercial_kit_obj (self):
        return self.commercialKit_id

    def get_commercial_kit_name(self):
        return '%s' %(self.commercialKit_id.get_name())



    objects = AdditionaKitsLibraryPreparationManager()

class AdditionalUserLotKitManager(models.Manager):
    def create_additional_user_lot_kit(self, user_additional_kit):
        new_user_additional_kit = self.create(lib_prep_id  = user_additional_kit['lib_prep_id'],
                    additionalLotKits = user_additional_kit['additionalLotKits'],
                    userLotKit_id = user_additional_kit['userLotKit_id'], value = 0)



class AdditionalUserLotKit (models.Model):
    lib_prep_id =  models.ForeignKey(
                    LibraryPreparation,
                    on_delete= models.CASCADE)
    additionalLotKits = models.ForeignKey(
                    AdditionaKitsLibraryPreparation,
                    on_delete= models.CASCADE)
    userLotKit_id = models.ForeignKey(
                    UserLotCommercialKits,
                    on_delete= models.CASCADE, null=True , blank=True)
    value = models.CharField(max_length=255)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__ (self):
        return '%s' %(self.userLotKit_id)

    def get_additional_kit_info(self):
        data = []
        data.append(self.additionalLotKits.get_kit_name())
        data.append(self.additionalLotKits.get_commercial_kit_name())
        if self.userLotKit_id == None:
            data.append('Not set')
        else:
            data.append(self.userLotKit_id.get_lot_number())
        data.append(self.generated_at.strftime("%d %B %Y"))
        return data


    objects = AdditionalUserLotKitManager()



class SambaConnectionData (models.Model):
    SAMBA_APPLICATION_FOLDER_NAME = models.CharField(max_length=80, null=True , blank=True)
    SAMBA_DOMAIN = models.CharField(max_length=80, null=True , blank=True)
    SAMBA_HOST_NAME = models.CharField(max_length=80, null=True , blank=True)
    SAMBA_IP_SERVER = models.CharField(max_length=20, null=True , blank=True)
    SAMBA_PORT_SERVER = models.CharField(max_length=10,  null=True , blank=True)
    SAMBA_REMOTE_SERVER_NAME = models.CharField(max_length=80, null=True , blank=True)
    SAMBA_SHARED_FOLDER_NAME =models.CharField(max_length=80, null=True , blank=True)
    SAMBA_USER_ID = models.CharField(max_length=20, null=True , blank=True)
    SAMBA_USER_PASSWORD = models.CharField(max_length=20, null=True , blank=True)
    IS_DIRECT_TCP = models.BooleanField(default = True)
    SAMBA_NTLM_USED = models.BooleanField(default = True)

    def __str__ (self):
        return '%s' %(self.SAMBA_REMOTE_SERVER_NAME)

    def get_samba_data(self):
        samba_data = {}
        for field in wetlab_config.SAMBA_CONFIGURATION_FIELDS:
            samba_data[field] = getattr(self, field)
        return samba_data

    def update_data(self, data):
        for field in wetlab_config.SAMBA_CONFIGURATION_FIELDS:
            setattr(self, field , data[field])
        self.save()
        return self

class EmailDataMamager(models.Manager):
    def create_email_data(self, data):
        new_email_data = self.create( hostName = data['EMAIL_HOST'], emailPort = data['EMAIL_PORT'],
            emailOnError= data['SENT_EMAIL_ON_ERROR'], emailUseTLS = data['USE_TLS'], userEmail = data['USER_EMAIL'],
            userName = data['USER_NAME'],  userPassword = data['USER_PASSWORD'])

        return new_email_data

class EmailData (models.Model):
    hostName =  models.CharField(max_length=80, null=True , blank=True)
    emailPort =  models.CharField(max_length=8, null=True , blank=True)
    emailOnError =  models.BooleanField(default = False)
    emailUseTLS =  models.BooleanField(default = False)
    userEmail =  models.CharField(max_length=80, null=True , blank=True)
    userName =  models.CharField(max_length=80, null=True , blank=True)
    userPassword = models.CharField(max_length=20, null=True , blank=True)

    def __str__ (self):
        return '%s' %(self.userEmail)

    def get_email_data(self):
        email_data ={}
        email_data['EMAIL_HOST'] = self.hostName
        email_data['EMAIL_PORT'] = self.emailPort
        email_data['SENT_EMAIL_ON_ERROR'] = self.emailOnError
        email_data['USE_TLS'] = self.emailUseTLS
        email_data['USER_EMAIL'] = self.userEmail
        email_data['USER_NAME'] = self.userName
        email_data['USER_PASSWORD'] = self.userPassword
        return email_data

    def update_data(self, data):
        self.hostName = data['EMAIL_HOST']
        self.emailPort = data['EMAIL_PORT']
        self.emailOnError= data['SENT_EMAIL_ON_ERROR']
        self.emailUseTLS = data['USE_TLS']
        self.userEmail = data['USER_EMAIL']
        self.userName = data['USER_NAME']
        self.userPassword = data['USER_PASSWORD']
        self.save()
        return self

    objects = EmailDataMamager()
