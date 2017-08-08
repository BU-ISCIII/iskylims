import datetime

from django.db import models
from django import forms
from django.utils import timezone
from django.utils.encoding import python_2_unicode_compatible

class UserInfo(models.Model):
    userid=models.CharField(max_length=25)
    userFirstName=models.CharField(max_length=45)
    userLastName=models.CharField(max_length=45)
    userArea=models.CharField(max_length=25)
    userEmail=models.EmailField(max_length=45)
    
    def __str__ (self):
        return '%s' %(self.userid)


class BioInfo(models.Model):
    info1=models.CharField(max_length=45)
    info2=models.CharField(max_length=45, default='')
    serviceRegistration= models.CharField(max_length=45, default='not required')
    researcher =models.CharField(max_length=45, default='')
    department= models.CharField(max_length=45, default='not assigned')
    

    def __str__(self):
        return '%s' %(self.info1)

class RunProcess(models.Model):
    runName = models.CharField(max_length=45)
    sampleSheet = models.FileField(upload_to='documents/')
    generatedat = models.DateTimeField(auto_now_add=True)
    runState = models.CharField(max_length=25)
    generatedBSFile = models.BooleanField(default=False)
    samples= models.CharField(max_length=45)
    useSpaceImgMb=models.CharField(max_length=10)
    useSpaceFastaMb=models.CharField(max_length=10)
    useSpaceOtherMb=models.CharField(max_length=10)
    bioinfo_id=models.ForeignKey(
            BioInfo)
    requestedCenter= models.CharField(max_length=45)
    
    def __str__(self):
        return '%s' %(self.runName)
        
    def get_sample_file (self):
        return '%s' %(self.sampleSheet)
        
    def get_state(self):
        return '%s' %(self.runState)
        
    def get_info_process (self):
        str_date=self.generatedat.strftime("%I:%M%p on %B %d, %Y")
        if (self.runState == 'Recorded' or self.runState == 'SampleSent'):
            return '%s;%s;%s;%s;%s'  %(self.runName, self.runState, self.requestedCenter, self.sampleSheet, str_date )
        else:
            return '%s;%s;%s;%s;%s;%s;%s'  %(self.runName, self.requestedCenter, self.useSpaceImgMb, self.useSpaceFastaMb, self.useSpaceOtherMb, str_date)

class Projects(models.Model):
    runprocess_id = models.ForeignKey(
            RunProcess,
            on_delete=models.CASCADE)
    user_id= models.ForeignKey(
            UserInfo)
    projectName= models.CharField(max_length=45)
    procState=models.CharField(max_length=25, default='Not Started')
    libraryKit=models.CharField(max_length=45)
    baseSpaceFile = models.CharField(max_length=255)
    
    def __str__(self):
        return '%s' %(self.projectName)


class RunningParameters (models.Model):
    runName_id = models.OneToOneField(
            RunProcess,
            on_delete=models.CASCADE,
            primary_key=True,
            )            
    RunID= models.CharField(max_length=255)
    ExperimentName= models.CharField(max_length=255)
    RTAVersion= models.CharField(max_length=255)
    SystemSuiteVersion= models.CharField(max_length=255)
    LibraryID= models.CharField(max_length=255)
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
    ImageChannel= models.CharField(max_length=255)
    Flowcell= models.CharField(max_length=255)
    ImageDimensions= models.CharField(max_length=255)
    FlowcellLayout= models.CharField(max_length=255)
    
    def __str__(self):
        return '%s' %(self.RunID)
        
    def get_run_parameters_info (self):
        return '%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s' %(self.RunID, self.ExperimentName, self.RTAVersion,
            self.SystemSuiteVersion, self.LibraryID, self.Chemistry, self.RunStartDate,
            self.AnalysisWorkflowType, self.RunManagementType, self.PlannedRead1Cycles, self.PlannedRead2Cycles,
            self.PlannedIndex1ReadCycles, self.PlannedIndex2ReadCycles, self.ApplicationVersion, self.NumTilesPerSwath,
            self.ImageChannel, self.Flowcell, self.ImageDimensions, self.FlowcellLayout)
    


'''
class Question(models.Model):
    question_text = models.CharField(max_length=200)
    pub_date = models.DateTimeField('date published')
    def __str__(self):
        return self.question_text
    
    def was_published_recently(self):
        return self.pub_date >= timezone.now() - datetime.timedelta(days=1)


class Choice(models.Model):
    question = models.ForeignKey(Question, on_delete=models.CASCADE)
    choice_text = models.CharField(max_length=200)
    votes = models.IntegerField(default=0)
    def __str__(self):
        return self.choice_text

class luis(models.Model):
	votes = models.IntegerField(default=0)
	choice_text = models.CharField(max_length=200)
'''
########## class for upload field
class Document(models.Model):
    run_name = models.CharField(max_length=255, blank=True)
    project_name = models.CharField(max_length=255, blank=True) 
    description = models.CharField(max_length=255, blank=True)
    csv_file = models.FileField(upload_to='documents/')
    uploaded_at = models.DateTimeField(auto_now_add=True)
    name = models.CharField(max_length=255, blank=True)
    user_id = models.CharField(max_length=50,blank=True)
    email = models.EmailField(max_length=255,blank=True)
    convert = models.BooleanField(default=False)
    
    
    def __str__(self):
        return '%s' %(self.run_name)
    
    def get_run_info (self):
        return '%s;%s;%s;%s;%s;%s;%s' %(self.run_name, self.project_name, 
                    self.user_id, self.description,  self.name, 
                    self.csv_file, self.uploaded_at)
        
class BaseSpaceFile (models.Model):
    document = models.OneToOneField(
            Document,
            on_delete=models.CASCADE,
            primary_key=True,
            )
    baseSpace_file = models.CharField(max_length=255)
    generated_at = models.DateTimeField(auto_now_add=True)
    
    def __str__(self):
        return '%s' %(self.baseSpace_file)
'''   
class NextSeqStatisticsBin (models.Model):
    document = models.OneToOneField(
            Document,
            on_delete=models.CASCADE,
            primary_key=True,
            )
    lane = models.CharField(max_length=10)
    singleRead = models.CharField(max_length=10)
    tiles = models.CharField(max_length=10)
    density = models.CharField(max_length=10)
    cluster_PF = models.CharField(max_length=10)
    phas_prephas = models.CharField(max_length=10)
    reads = models.CharField(max_length=10)
    reads_PF = models.CharField(max_length=10)
    q30 = models.CharField(max_length=10)
    yields = models.CharField(max_length=10)
    cyclesErrRated = models.CharField(max_length=10)
    aligned = models.CharField(max_length=10)
    errorRate = models.CharField(max_length=10)
    errorRate35 = models.CharField(max_length=10)
    errorRate75 = models.CharField(max_length=10)
    errorRate100 = models.CharField(max_length=10)
    intensityCycle = models.CharField(max_length=10)
    comments = models.CharField(max_length=10)
    status = models.CharField(max_length=10)
    
    def __str__(self):
        return '%s' %(self.document)
'''
''' 
class RawNextSeqStatisticsXml (models.Model):
    document = models.ForeignKey(
            Document,
            on_delete=models.CASCADE)
    rawYield = models.CharField(max_length=64)
    rawYieldQ30= models.CharField(max_length=64)
    rawQuality= models.CharField(max_length=64)
    PF_Yield= models.CharField(max_length=64)
    PF_YieldQ30= models.CharField(max_length=64)
    PF_Quality= models.CharField(max_length=64)
    barcodeCount= models.CharField(max_length=64)
    perfectBarcodeCount= models.CharField(max_length=64)
    project=models.CharField(max_length=30,default='-')
    generated_at = models.DateTimeField(auto_now_add=True)
    
    def __str__(self):
        return '%s' %(self.document)
        
    def get_nextSeq_stats(self):
        return '%s;%s;%s;%s;%s;%s;%s;%s;%s' %(self.project, self.barcodeCount,
                self.perfectBarcodeCount, self.rawYield, self.PF_Yield,
                self.rawYieldQ30, self.PF_YieldQ30, self.rawQuality,
                self.PF_Quality)
'''
''' 
class NextSeqStatisticsXml (models.Model):
    document = models.ForeignKey(
            Document,
            on_delete=models.CASCADE)
    flowSummRaw = models.CharField(max_length=20)
    flowSummPF= models.CharField(max_length=20)
    rawQuality= models.CharField(max_length=64)
    PF_Yield= models.CharField(max_length=64)
    PF_YieldQ30= models.CharField(max_length=64)
    PF_Quality= models.CharField(max_length=64)
    barcodeCount= models.CharField(max_length=64)
    perfectBarcodeCount= models.CharField(max_length=64)
    project=models.CharField(max_length=30,default='-')
    generated_at = models.DateTimeField(auto_now_add=True)

'''
'''
class NextSeqGraphicsStats (models.Model):
    document = models.OneToOneField(
            Document,
            on_delete=models.CASCADE,
            primary_key=True,
            )
    folderGraphic= models.CharField(max_length=255)
    cluserCountGraph= models.CharField(max_length=255)
    flowCellGraph= models.CharField(max_length=255)
    intensityByCycleGraph= models.CharField(max_length=255)
    heatMapGraph= models.CharField(max_length=255)
    histogramGraph= models.CharField(max_length=255)
    sampleQcGraph= models.CharField(max_length=255)
    
    def __str__(self):
        return '%s' %(self.document)
    
class MiSeqStatisticsBin (models.Model):
    document = models.OneToOneField(
            Document,
            on_delete=models.CASCADE,
            primary_key=True,
            )
    lane = models.CharField(max_length=10)
    singleRead = models.CharField(max_length=10)
    tiles = models.CharField(max_length=10)
    density = models.CharField(max_length=10)
    cluster_PF = models.CharField(max_length=10)
    phas_prephas = models.CharField(max_length=10)
    reads = models.CharField(max_length=10)
    reads_PF = models.CharField(max_length=10)
    q30 = models.CharField(max_length=10)
    yields = models.CharField(max_length=10)
    cyclesErrRated = models.CharField(max_length=10)
    aligned = models.CharField(max_length=10)
    errorRate = models.CharField(max_length=10)
    errorRate35 = models.CharField(max_length=10)
    errorRate75 = models.CharField(max_length=10)
    errorRate100 = models.CharField(max_length=10)
    intensityCycle = models.CharField(max_length=10)
    comments = models.CharField(max_length=10)
    status = models.CharField(max_length=10)
    
    def __str__(self):
        return '%s' %(self.document)

'''
'''    
class MiSeqStatisticsXml (models.Model):
    document = models.ForeignKey(
            Document,
            on_delete=models.CASCADE)
    rawYield = models.CharField(max_length=20)
    rawYieldQ30= models.CharField(max_length=20)
    rawQuality= models.CharField(max_length=20)
    PF_Yield= models.CharField(max_length=20)
    PF_YieldQ30= models.CharField(max_length=20)
    PF_Quality= models.CharField(max_length=20)
    barcodeCount= models.CharField(max_length=20)
    perfectBarcodeCount= models.CharField(max_length=20)
    project=models.CharField(max_length=30,default=False)
    generated_at = models.DateTimeField(auto_now_add=True)
    
    def __str__(self):
        return '%s' %(self.document)
'''
'''
class MiSeqGraphicsStats (models.Model):
    document = models.ForeignKey(
            Document,
            on_delete=models.CASCADE)
    folderGraphic= models.CharField(max_length=255)
    cluserCountGraph= models.CharField(max_length=255)
    flowCellGraph= models.CharField(max_length=255)
    intensityByCycleGraph= models.CharField(max_length=255)
    heatMapGraph= models.CharField(max_length=255)
    histogramGraph= models.CharField(max_length=255)
    sampleQcGraph= models.CharField(max_length=255)
    
    
    def __str__(self):
        return '%s' %(self.document)
'''  

    
    



	
