################# SAMPLE SETTINGS ##############################
'''
ADDITIONAL_HEADING_FOR_RECORDING_SAMPLES = ['Sample','Order', 'Confirmation Code', 'Priority', 'Requested Service Date',
                'History Number', 'Requested Service by', 'Doctor', 'Suspicion', 'Comments']

MAP_ADDITIONAL_HEADING_TO_DATABASE = [('Order','orderInEntry'), ('Confirmation Code', 'confirmationCode'), ('Priority', 'priority'),
                ('Requested Service Date','serviceDate' ), ('Comments', 'comments')]

HEADING_FOR_STORED_PATIENT_DATA = ['Sample', 'History Number', 'Requested Service by']
'''

HEADING_FOR_DISPLAY_PATIENT_BASIC_INFORMATION = ['Patient Name', 'Patient Surame', 'Patient Code' ,'Sex']
HEADING_FOR_DISPLAY_REQUESTED_BY_INFORMATION = ['Requested Service by', 'Doctor Name']
HEADING_FOR_DISPLAY_PATIENT_ADDITIONAL_INFORMATION = ['Address', 'Phone number', 'email', 'Birthday', 'Smoker', 'Notification Preferences', 'Comments']

HEADING_FOR_DISPLAY_SAMPLE_DATA_IN_PATIENT_INFO = ['Sample Name', 'Sample Origin', 'Entry Date', 'Sample type', 'Sample state']



HEADING_FOR_DISPLAY_SAMPLE_MAIN_INFORMATION = ['Order', 'Confirmation Code', 'Priority','Sample Origin', 'Type of Sample', 'Species', 'Entry Date', 'User name']

HEADING_FOR_DISPLAY_SAMPLE_DEFINED_STATE = ['Sample Name','Sample Origin', 'Type of Sample', 'Species', 'Date for entry in Lab']


#HEADING_FOR_DISPLAY_SAMPLE_PATIENT_SEQUENCING_STATE = ['Sample Name','Order', 'Confirmation Code', 'Priority', 'History Number', 'Doctor Name' , 'Service Requested by' ]

#HEADING_FOR_DISPLAY_SAMPLE_PENDING_RESULT_STATE = ['Sample Name','Order', 'Confirmation Code', 'Priority', 'History Number', 'Doctor Name' , 'Service Requested by', 'To be included' ]

#HEADING_FOR_NOT_MATCH = ['Sample', 'Order', 'Confirmation' , 'Priority' ,'Requested Service Date', 'Wrong History Number']

##############################################################
################## ERROR MESSAGES  ###########################
##############################################################
ERROR_MESSAGE_FOR_PATIENT_CODE_EXISTS = 'Patient Code already exists.'
ERROR_MESSAGE_FOR_PROJECT_NAME_EXISTS = 'Project Name already exists.'

ERROR_MESSAGE_FOR_INCORRECT_START_SEARCH_DATE = ['The format for the "Start Date Search" Field is incorrect ', 'Use the format  (DD-MM-YYYY)']

ERROR_MESSAGE_FOR_INCORRECT_END_SEARCH_DATE = ['The format for the "End Date Search" Field is incorrect ', 'Use the format  (DD-MM-YYYY)']

ERROR_MESSAGE_FOR_SORT_PATIENT_NAME = ['The Patient name must contains more than 4 characters']

ERROR_MESSAGE_FOR_NO_MATCH_IN_SEARCH = ['Your query did not return any results']

###################### FIELDS NAME USED IN USER FORM  ########################################
FORM_MAIN_DATA_PATIENT_DEFINITION = ['patientName', 'patientSurname', 'patientCode', 'patientSex']
FORM_OPT_DATA_PATIENT_DEFINITION = ['address', 'phone', 'email', 'birthday', 'smoker', 'notificationPreference' , 'comments']

################ SHOWING LIST WHEN SEARCH MATCH MORE ONE RESULT ##############################
HEADING_SEARCH_LIST_SAMPLES_CLINIC = ['Sample Clinic name', 'History Number ' ,'Doctor name' , 'Priority', 'Requested Service Date', 'State']


################ RESULT PARAMETER SETTINGS ##############################
### Headings used when defining the custom result parameters
HEADING_FOR_RESULT_TO_PROTOCOL = ['Clinic Sample', 'History Number', 'Priority', 'Protocol Type', 'Protocol to be used']

################ MAPPING SAMPLE STATES ON CORE WITH CLINIC ######################
MAPPING_SAMPLES_CORE_VS_CLINIC = {'Defined': 'Patient update', 'Extract molecule': 'Sequencing',  'Library preparation': 'Sequencing', 'Pool Preparation' : 'Sequencing',
                                    'Sequencing': 'Sequencing', 'Completed':'Pending protocol', 'Error':'Invalid'}
