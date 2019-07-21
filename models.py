import datetime

from django.db import models
from django import forms
from django.utils import timezone
from django.utils.encoding import python_2_unicode_compatible
from django.contrib.auth.models import User
from django_utils.models import Center
from django.utils.translation import ugettext_lazy as _

from .  import wetlab_config

class RunErrors (models.Model):
    errorCode = models.CharField(max_length=10)
    errorText = models.CharField(max_length=255)

    def __str__ (self):
        return '%s' %(self.errorText)

class RunStates (models.Model):
    runStateName = models.CharField(max_length=50)

    def __str__ (self):
        return '%s' %(self.runStateName)

class RunProcess(models.Model):
    runName = models.CharField(max_length=45)
    sampleSheet = models.FileField(upload_to = wetlab_config.RUN_SAMPLE_SHEET_DIRECTORY)
    generatedat = models.DateTimeField(auto_now_add=True)
    run_date = models.DateField(auto_now = False, null=True)
    run_finish_date = models.DateTimeField(auto_now = False, null=True, blank=True)
    bcl2fastq_finish_date = models.DateTimeField(auto_now = False, null=True, blank=True)
    run_completed_date = models.DateTimeField(auto_now = False, null=True, blank=True)
    runState = models.CharField(max_length=25)
    state = models.ForeignKey ( RunStates, on_delete = models.CASCADE, related_name = 'state_of_run', null = True, blank = True)
    index_library = models.CharField(max_length=85)
    samples= models.CharField(max_length=45,blank=True)
    useSpaceImgMb=models.CharField(max_length=10, blank=True)
    useSpaceFastaMb=models.CharField(max_length=10, blank=True)
    useSpaceOtherMb=models.CharField(max_length=10, blank=True)
    centerRequestedBy = models.ForeignKey (Center, on_delete=models.CASCADE)
    sequencerModel = models.ForeignKey ('iSkyLIMS_drylab.Machines', on_delete=models.CASCADE, null=True, blank=True)
    runError = models.ForeignKey( RunErrors, on_delete=models.CASCADE, null = True, blank = True)
    stateBeforeError = models.ForeignKey ( RunStates, on_delete = models.CASCADE, null = True, blank = True)

    def __str__(self):
        return '%s' %(self.runName)

    def get_run_id (self):
        return '%s' %(self.id)

    def get_run_date (self):
        if self.run_date is None :
            rundate = 'Run NOT started'
        else :
            rundate=self.run_date.strftime("%B %d, %Y")
        return rundate

    def get_error_text (self):
        return '%s' %(self.runError)

    def get_state(self):
        return '%s' %(self.state)

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

    def get_run_sequencerModel (self):
        return '%s' %(self.sequencerModel)

    def get_run_platform (self):
        return '%s' %self.sequencerModel.platformID

    def get_machine_lanes(self):
        number_of_lanes = self.sequencerModel.get_number_of_lanes()
        return int(number_of_lanes)

    def get_sample_file (self):
        return '%s' %(self.sampleSheet)

    def update_library (self, library_name):
        self.index_library = library_name
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

    def set_sequencer (self, sequencer):
        self.sequencerModel = sequencer
        self.save()
        return True

class BaseSpaceLibraryName (models.Model):
    libraryName = models.CharField(max_length=125)
    #sampleNumber = models.CharField(max_length = 25)
    #indexNumber = models.CharField(max_length = 25)
    generatedat = models.DateTimeField(auto_now_add=True)



class IndexLibraryKit (models.Model):
    indexLibraryName = models.CharField(max_length=125)
    version = models.CharField(max_length=80,null=True)
    plateExtension = models.CharField(max_length=125, null=True)
    adapter1 = models.CharField(max_length=125,null=True)
    adapter2 = models.CharField(max_length=125, null=True)
    indexLibraryFile =  models.FileField(upload_to=wetlab_config.LIBRARY_KITS_DIRECTORY )
    generatedat = models.DateTimeField(auto_now_add=True, null=True)

    def __srt__ (self):
        return '%s'(self.indexLibraryName)

    def get_index_library_information (self):
        if self.adapter2 =='':
            adapter2 = 'Not used on this library'
        else:
            adapter2 = self.adapter2
        return '%s;%s;%s;%s;%s;%s' %(self.indexLibraryName, self.version, self.plateExtension ,
                            self.adapter1  , adapter2, self.indexLibraryFile)



class IndexLibraryValues (models.Model):
    indexLibraryKit_id = models.ForeignKey(
        IndexLibraryKit,
        on_delete=models.CASCADE)
    defaultWell = models.CharField(max_length=10,null=True)
    index_7 = models.CharField(max_length=25,null=True)
    i_7_seq = models.CharField(max_length=25,null=True)
    index_5 = models.CharField(max_length=25,null=True)
    i_5_seq = models.CharField(max_length=25,null=True)


    def get_index_value_information (self):
        return '%s;%s;%s;%s;%s' %(self.defaultWell, self.index_7, self.i_7_seq,
                        self.index_5, self.i_5_seq)


class Projects(models.Model):
    runprocess_id = models.ForeignKey(
            RunProcess,
            on_delete=models.CASCADE, null = True) # added null for new lab process functionality
    user_id= models.ForeignKey(User,on_delete=models.CASCADE, null = True)
    LibraryKit_id = models.ForeignKey(
            BaseSpaceLibraryName,
            on_delete=models.CASCADE , null=True)
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

    def get_state(self):
        return '%s' %(self.runprocess_id.state.runStateName)

    def get_project_info (self):
        generated_date=self.generatedat.strftime("%I:%M%p on %B %d, %Y")
        if self.project_run_date is None:
            projectdate = 'Run NOT started'
        else :
            projectdate=self.project_run_date.strftime("%B %d, %Y")
        return '%s;%s;%s;%s;%s' %(self.projectName, self.libraryKit,
                        self.baseSpaceFile, generated_date, projectdate )

    def get_p_info_change_library (self):
        run_name = self.runprocess_id.runName
        user_name = self.user_id.username
        if self.project_run_date is None:
            projectdate = 'Run NOT started'
        else :
            projectdate=self.project_run_date.strftime("%B %d, %Y")

        return '%s;%s;%s;%s;%s'%(run_name, self.projectName, projectdate, user_name, self.libraryKit)

    def get_user_name (self):
        user_name = self.user_id.username
        return '%s' %(user_name)

    def get_project_name (self):
        return '%s' %(self.projectName)

    def get_library_name (self):
        return '%s' %(self.libraryKit)

    def get_date (self):
        if self.project_run_date is None:
            projectdate = 'No Date'
        else :
            projectdate=self.project_run_date.strftime("%B %d, %Y")
        return '%s' %(projectdate)

    def get_project_id(self):
        return '%s' %(self.id)


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
    project_id = models.ForeignKey(Projects, on_delete=models.CASCADE, null=True)
    defaultAll =models.CharField(max_length=40, null=True)
    flowRawCluster = models.CharField(max_length=40)
    flowPfCluster= models.CharField(max_length=40)
    flowYieldMb= models.CharField(max_length=40)
    sampleNumber= models.CharField(max_length=40)
    generated_at = models.DateTimeField(auto_now_add=True)

    def get_fl_summary(self):
        return '%s;%s;%s;%s' %(self.flowRawCluster, self.flowPfCluster,
            self.flowYieldMb, self.sampleNumber)

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
    project_id = models.ForeignKey(Projects, on_delete=models.CASCADE,  null=True)
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
        return '%s;%s;%s;%s;%s;%s;%s;%s' %(self.lane, self.pfCluster,
                self.percentLane, self.perfectBarcode, self.oneMismatch,
                self.yieldMb, self.biggerQ30, self.meanQuality)

    def get_stats_info (self):

        return'%s;%s;%s;%s' %(self.biggerQ30, self.meanQuality, self.yieldMb, self.pfCluster)

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



# New objets for handling process requests form lab user

REAGENT_USED_CHOICE = (
    ('Not used', _("NOT USED")),
    ('Partial used', _("PARTIAL USED")),
    ('Completed used', _("COMPLETED USED"))
    )

class ReagentsLibraryKit (models.Model):
#    batch_id = models.CharField(max_length= 255)
    provider = models.CharField(max_length =60)
    reagentLibraryName = models.CharField(max_length = 100)
    cuantityUsed = models.CharField(choices = REAGENT_USED_CHOICE, max_length = 50, default = 'Not used')
    usedDate = models.DateTimeField(null = True)
    chipLot = models.CharField(max_length = 100)
#    sampleNumber = models.CharField(max_length = 25)
#    indexNumber = models.CharField(max_length = 25)
    expirationDate = models.DateField(auto_now_add=False)
    generatedat = models.DateTimeField(auto_now_add=True, null=True)


class ProtocolInLab (models.Model):
    protocolName = models.CharField(max_length=255)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__ (self):
        return '%s' %(self.protocolName)

    def get_name (self):
        return '%s' %(self.protocolName)

# Q-Fluor ng/ul	Tamaño	nM Quantifluor

class ProtocolParameters (models.Model) :
    protocol_id = models.ForeignKey(
                    ProtocolInLab,
                    on_delete= models.CASCADE)
    parameterName = models.CharField(max_length=255)
    parameterDescription = models.CharField(max_length= 400, null=True, blank=True)
    parameterOrder = models.IntegerField()
    parameterUsed = models.BooleanField()
    parameterMaxValue = models.CharField(max_length = 50, null = True, blank = True)
    parameterMinValue = models.CharField(max_length = 50, null = True, blank = True)


    def __str__ (self):
        return '%s' %(self.parameterName)




class ReferenceGenome (models.Model):
    refGenomeName = models.CharField(max_length=255)
    refGenomeSize = models.CharField(max_length=100)
    refGenomeID = models.CharField(max_length=255)


    def __str__ (self):
        return '%s' %(self.refGenomeName)


class StatesForSample (models.Model):
    sampleStateName = models.CharField(max_length=50)

    def __str__ (self):
        return '%s' %(self.sampleStateName)

class SamplesInProjectManager (models.Manager):
    def create_sample_project (self, s_project):
        sample_project = self.create( project_id = s_project['project_id'] , sampleState__StatesForSample = 'Completed',
                            sampleProtocol = '', sampleType = '', registerUser = '', sampleName =  s_project['sampleName'],
                            barcodeName =  s_project['barcodeName'], pfClusters =  s_project['pfClusters'],
                            percentInProject =  s_project['percentInProject'], yieldMb =  s_project['yieldMb'],
                            qualityQ30 =  s_project['qualityQ30'], meanQuality =  s_project['meanQuality'] )
    def update_sample_project (self, s_project):
        sample_project = self.update( project_id = s_project['project_id'] , sampleState__StatesForSample = 'Completed',
                            barcodeName =  s_project['barcodeName'], pfClusters =  s_project['pfClusters'],
                            percentInProject =  s_project['percentInProject'], yieldMb =  s_project['yieldMb'],
                            qualityQ30 =  s_project['qualityQ30'], meanQuality =  s_project['meanQuality'] )

    def create_sample_from_investigator (self, sample_data):
        sample_from_investigator = self.create( project_id = '', sampleState__StatesForSample = 'Defined',
                            sampleProtocol__ProtocolInLab = sample_data['protocol'],
                            sampleType = sample_data['type'] , registerUser__username__exact = sample_data['user'],
                            sampleName =  sample_data['sample_name'],
                            barcodeName =  '', pfClusters = '', percentInProject = '',
                            yieldMb = '', qualityQ30 = '', meanQuality = '')


class SamplesInProject (models.Model):
    project_id = models.ForeignKey(
                Projects,
                on_delete= models.CASCADE, null = True)
    sampleState = models.ForeignKey(
                StatesForSample,
                on_delete = models.CASCADE, null = True)
    sampleProtocol = models.ForeignKey(
                ProtocolInLab,
                on_delete = models.CASCADE, null = True)
    sampleType = models.CharField(max_length=50, null = True)
    registerUser = models.ForeignKey(
                User,
                on_delete=models.CASCADE, null = True)
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

    def get_sample_information(self):
        return '%s;%s;%s;%s;%s;%s;%s' %(self.sampleName , self.barcodeName,
                self.pfClusters ,self.percentInProject, self.yieldMb ,
                self.qualityQ30 , self.meanQuality )
    def get_sample_name(self):
        return '%s' %(self.sampleName)

    def get_register_user(self):
        return '%s' %(self.registerUser)

    def get_project_name (self) :
        #p_id = self.project_id
        #project_name =Projects.objects.get(projectName=p_id).get_project_name()
        project_name = self.project_id.get_project_name()
        #Projects.objects.prefetch_related('user_id').filter(user_id = user_id)
        return '%s' %(project_name)

    def get_quality_sample (self):
        return '%s' %(self.qualityQ30)

    objects = SamplesInProjectManager()

class SampleProtocolParameterData (models.Model):
    parameter_id = models.ForeignKey(
                    ProtocolParameters,
                    on_delete= models.CASCADE)
    sample_id = models.ForeignKey(
                SamplesInProject,
                on_delete= models.CASCADE)
    parameterValue = models.CharField(max_length=255)
    generated_at = models.DateTimeField(auto_now_add=True)

    def __str__ (self):
        return '%s' %(self.parameterValue)
