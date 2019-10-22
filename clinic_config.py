################# SAMPLE SETTINGS ##############################
ADDITIONAL_HEADING_FOR_RECORDING_SAMPLES = ['Sample','Order', 'Confirmation Code', 'Priority', 'Requested Service Date',
                'History Number', 'Requested Service by', 'Doctor', 'Suspicious', 'Comments']

MAP_ADDITIONAL_HEADING_TO_DATABASE = [('Order','orderInEntry'), ('Confirmation Code', 'confirmationCode'), ('Priority', 'priority'),
                ('Requested Service Date','serviceDate' ), ('Comments', 'comments')]

HEADING_FOR_STORED_PATIENT_DATA = ['Sample', 'History Number', 'Requested Service by']

HEADING_FOR_DISPLAY_PATIENT_INFORMATION = ['Patient Name', 'History Number']
HEADING_FOR_DISPLAY_REQUESTED_BY_INFORMATION = ['Requested Service by', 'Doctor Name']

HEADING_FOR_DISPLAY_SAMPLE_CORE_INFORMATION =['Order', 'Confirmation Code', 'Priority','Laboratory', 'Type of Sample', 'Species', 'Date for entry in Lab']

HEADING_FOR_NOT_MATCH = ['Sample', 'Order', 'Confirmation' , 'Priority' ,'Requested Service Date', 'Wrong History Number']

ERROR_MESSAGE_FOR_INCORRECT_START_SEARCH_DATE = ['The format for the "Start Date Search" Field is incorrect ', 'Use the format  (DD-MM-YYYY)']

ERROR_MESSAGE_FOR_INCORRECT_END_SEARCH_DATE = ['The format for the "End Date Search" Field is incorrect ', 'Use the format  (DD-MM-YYYY)']

ERROR_MESSAGE_FOR_SORT_PATIENT_NAME = ['The Patient name must contains more than 5 characters']

ERROR_MESSAGE_FOR_NO_MATCH_IN_SEARCH = 'Your query did not return any results'

################ SHOWING LIST WHEN SEARCH MATCH MORE ONE RESULT ##############################
HEADING_SEARCH_LIST_SAMPLES_CLINIC = ['Sample Clinic name', 'History Number ' ,'Doctor name' , 'Priority', 'Requested Service Date', 'State']


################ RESULT PARAMETER SETTINGS ##############################
### Headings used when defining the custom result parameters
HEADING_FOR_RESULT_TO_PROTOCOL = ['Clinic Sample', 'History Number', 'Priority', 'Protocol Type', 'Protocol to be used']

################ MAPPING SAMPLE STATES ON CORE WITH CLINIC ######################
MAPPING_SAMPLES_CORE_VS_CLINIC = {'Defined': 'Patient update', 'Extract molecule': 'Sequencing',  'Library preparation': 'Sequencing', 'Pool Preparation' : 'Sequencing',
                                    'Sequencing': 'Sequencing', 'Completed':'Pending protocol', 'Error':'Invalid'}
