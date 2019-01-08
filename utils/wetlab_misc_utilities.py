##General purpose utility functions used in different iSkyLIMS_wetlab modules

import datetime

import os, errno, re



from django.core.mail import send_mail

def send_error_email_to_user ( subject, body_message, from_user, to_user):
    '''
    Description:
        Send an email to users  defined in the user list "to_user"
    Input:
    '''
    send_mail (subject, body_message, from_user, to_user)


    








def normalized_data (set_data, all_data) :
    '''
    Description:
        The function is used to normalized data from diferent range of values
    Input:
        set_data    # contains a gruop of data
        all_data    # contains all data to be used for the normalization
    Variables:
        normalized_set_data # to keep the normalized set data
        normalized_all_data # to keep the normalized value of all data
    Return:
        normalized_set_data
        normalized_all_data.
    '''
    normalized_set_data, normalized_all_data = [] , []
    min_value = min(min(set_data),min(all_data))
    max_value = max(max(set_data), max(all_data))
    for value in set_data :
        normalized_set_data.append(format((value - min_value)/max_value,'.2f'))
    for value in all_data :
        normalized_all_data.append(format((value - min_value)/max_value,'.2f'))

    return normalized_set_data, normalized_all_data







