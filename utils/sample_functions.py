from iSkyLIMS_wetlab.models import *


def get_data_from_excel_form (spreadsheet, n_columns):
    '''
    Description:
        The function split the spreadsheet data. Every row in the spreadsheet is stored
        in a list. Return a array with evey item list.
    Input:
        spreadsheet     # contains all data
        n_columns       # number of columns used to
    Variables:
        spread_data # array containing row data list
    Return:
        d_list # data splited in rows
    '''
    split_data = spreadsheet.split(',')
    d_list = []
    for i in range(0,len(split_data),n_columns):
        d_list.append(split_data[i:i+n_columns])
    return d_list

def get_protocols_name ():
    '''
    Description:
        The function will return the name of the protols defined in database.
    Input:
        none
    Variables:
        prot_names # list containing all protocol names
    Return:
        prot_names.
    '''
    prot_names = []
    if ProtocolInLab.objects.filter().exists():
        protocols = ProtocolInLab.objects.all()

        for protocol in protocols:
            prot_names.append(protocol.get_name())
    return prot_names


def sample_already_defined(sample_name):
    if SamplesInProject.objects.filter(sampleName__exact = sample_name, sampleState__StatesForSample = 'Defined').exists():
        return True
    else:
        return False
