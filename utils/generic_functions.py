

def check_empty_fields (row_data):
    for data in row_data:
        if data == '':
            return True
    return False
