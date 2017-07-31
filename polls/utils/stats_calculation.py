import sys
import xml.etree.ElementTree as ET
from polls.models import *


  


def get_statistics_xml(demux_file, conversion_file):
 
    stats_result={}

    #demux_file='example.xml'
    demux_stat=ET.parse(demux_file)
    root=demux_stat.getroot()
    projects=[]
    for child in root.iter('Project'):
        projects.append(child.attrib['name'])
    
    for i in range(len(projects)):
        p_temp=root[0][i]
        samples=p_temp.findall('Sample')
        sample_all_index=len(samples)-1
        barcodeCount=[]
        perfectBarcodeCount=[]
        b_count=[]
        p_b_count=[]
        dict_stats={}
        for c in p_temp[sample_all_index].iter('BarcodeCount'):
            #b_count.append(c.text)
            barcodeCount.append(c.text)
        for c in p_temp[sample_all_index].iter('PerfectBarcodeCount'):
            p_b_count.append(c.text)
        perfectBarcodeCount.append(p_b_count)
        dict_stats['BarcodeCount']=barcodeCount
        dict_stats['PerfectBarcodeCount']=perfectBarcodeCount
        stats_result[projects[i]]=dict_stats
    
    
    #demux_file='example.xml'
    conversion_stat=ET.parse(conversion_file)
    root_conv=conversion_stat.getroot()
    projects=[]
    
    
    for child in root_conv.iter('Project'):
        projects.append(child.attrib['name'])
    for i in range(len(projects)):
        p_temp=root_conv[0][i]
        samples=p_temp.findall('Sample')
        sample_all_index=len(samples)-1
        tiles=p_temp[sample_all_index][0][0].findall('Tile')
        tiles_index=len(tiles)-1
        list_raw_yield=[]
        list_raw_yield_q30=[]
        list_raw_quality=[]
        list_pf_yield=[]
        list_pf_yield_q30=[]
        list_pf_quality=[]
    
        for l_index in range(4):
            raw_yield_value = 0
            raw_yield_q30_value = 0
            raw_quality_value = 0   
            pf_yield_value = 0
            pf_yield_q30_value = 0
            pf_quality_value = 0
            for t_index in range(tiles_index):
                
                     # get the yield value for RAW and for read 1 and 2
                for c in p_temp[sample_all_index][0][l_index][t_index][0].iter('Yield'):
                    raw_yield_value +=int(c.text)
                    # get the yield Q30 value for RAW  and for read 1 and 2
                for c in p_temp[sample_all_index][0][l_index][t_index][0].iter('YieldQ30'):
                    raw_yield_q30_value +=int(c.text)
                for c in p_temp[sample_all_index][0][l_index][t_index][0].iter('QualityScoreSum'):
                    raw_quality_value +=int(c.text)
                 # get the yield value for PF and for read 1 and 2
                for c in p_temp[sample_all_index][0][l_index][t_index][1].iter('Yield'):
                    pf_yield_value +=int(c.text)
                # get the yield Q30 value for PF and for read 1 and 2
                for c in p_temp[sample_all_index][0][l_index][t_index][1].iter('YieldQ30'):
                    pf_yield_q30_value +=int(c.text)
                for c in p_temp[sample_all_index][0][l_index][t_index][1].iter('QualityScoreSum'):
                    pf_quality_value +=int(c.text)
            list_raw_yield.append(str(raw_yield_value))
            list_raw_yield_q30.append(str(raw_yield_q30_value))
            list_raw_quality.append(str(raw_quality_value))
            list_pf_yield.append(str(pf_yield_value))
            list_pf_yield_q30.append(str(pf_yield_q30_value))
            list_pf_quality.append(str(pf_quality_value))
                
        stats_result[projects[i]]['RAW_Yield']=list_raw_yield
        stats_result[projects[i]]['RAW_YieldQ30']=list_raw_yield_q30
        stats_result[projects[i]]['RAW_Quality']=list_raw_quality
        stats_result[projects[i]]['PF_Yield']=list_pf_yield
        stats_result[projects[i]]['PF_YieldQ30']=list_pf_yield_q30
        stats_result[projects[i]]['PF_Quality']=list_pf_quality
    return stats_result

def get_running_data(run_file, run_parameter):
    running_data={}
    image_channel=[]
    
    run_data=ET.parse(run_file)
    run_root=run_data.getroot()
    p_run=run_root[0]
    for i in run_root.iter('Name'):
        image_channel.append(i.text)
    running_data['ImageChannel']=image_channel
    running_data['Flowcell']=p_run.find('Flowcell').text
    running_data['ImageDimensions']=p_run.find('ImageDimensions').attrib
    running_data['FlowcellLayout']=p_run.find('FlowcellLayout').attrib
    
    parameter_data=ET.parse(run_parameter)
    parameter_data_root=parameter_data.getroot()
    p_parameter=parameter_data_root[1]
    running_data['RunID']=parameter_data_root.find('RunID').text
    running_data['ExperimentName']=parameter_data_root.find('ExperimentName').text
    running_data['RTAVersion']=parameter_data_root.find('RTAVersion').text
    running_data['SystemSuiteVersion']=parameter_data_root.find('SystemSuiteVersion').text
    running_data['LibraryID']=parameter_data_root.find('LibraryID').text
    running_data['Chemistry']=parameter_data_root.find('Chemistry').text
    running_data['RunStartDate']=parameter_data_root.find('RunStartDate').text
    running_data['AnalysisWorkflowType']=parameter_data_root.find('AnalysisWorkflowType').text
    running_data['RunManagementType']=parameter_data_root.find('RunManagementType').text
    running_data['PlannedRead1Cycles']=parameter_data_root.find('PlannedRead1Cycles').text
    running_data['PlannedRead2Cycles']=parameter_data_root.find('PlannedRead2Cycles').text
    running_data['PlannedIndex1ReadCycles']=parameter_data_root.find('PlannedIndex1ReadCycles').text
    running_data['PlannedIndex2ReadCycles']=parameter_data_root.find('PlannedIndex2ReadCycles').text
    running_data['ApplicationVersion']=p_parameter.find('ApplicationVersion').text
    running_data['NumTilesPerSwath']=p_parameter.find('NumTilesPerSwath').text
   
    return running_data



def store_in_db(statistics, used_table,run_name_value):
    
    if (used_table == 'running_table'):
        run_name_value = statistics['RunID']
        run_db= RunDataParameters(document=Document.objects.get(run_name__icontains = run_name_value),
            RunID= statistics['RunID'], ExperimentName= statistics['ExperimentName'] ,
            ImageChannel = statistics['ImageChannel'],  Flowcell = statistics['Flowcell'],
            ImageDimensions = statistics['ImageDimensions'],  FlowcellLayout = statistics['FlowcellLayout'],
            RTAVersion = statistics['RTAVersion'], SystemSuiteVersion = statistics['SystemSuiteVersion'],
            LibraryID = statistics['LibraryID'], Chemistry = statistics['Chemistry'],
            RunStartDate = statistics['RunStartDate'],  AnalysisWorkflowType= statistics['AnalysisWorkflowType'],
            RunManagementType = statistics['RunManagementType'], PlannedRead1Cycles = statistics['PlannedRead1Cycles'],
            PlannedRead2Cycles = statistics['PlannedRead2Cycles'],  PlannedIndex1ReadCycles= statistics['PlannedIndex1ReadCycles'],
            ApplicationVersion = statistics['ApplicationVersion'],  NumTilesPerSwath= statistics['NumTilesPerSwath'])
                        
        #run_db.save()
        return 
    if (used_table == 'nextSeqXml_table'):
        for project in statistics:
            stats_db = NextSeqStatisticsXml(document=Document.objects.get(run_name__icontains = run_name_value),
                RAW_Yield = xml_statistics[project]['RAW_Yield'], RAW_YieldQ30 = xml_statistics[project]['RAW_YieldQ30'],
                RAW_Quality= xml_statistics[project]['RAW_Quality'],  PF_Yield= xml_statistics[project]['PF_Yield'],
                PF_YieldQ30= xml_statistics[project]['PF_YieldQ30'],  PF_Quality= xml_statistics[project]['PF_Quality'],
                BarcodeCount= xml_statistics[project]['BarcodeCount'],  
                PerfectBarcodeCount= xml_statistics[project]['PerfectBarcodeCount'],
                project=project)
            #stats_db.save()
            #import pdb; pdb.set_trace()
        return
'''
def get_nextSeq_stats_fromDB(run_name_value):
    stats_index=NextSeqStatisticsXml.objects.filter(document__icontains =run_name_value)
    for index in stats_index:
        stats_values=index.
'''            
