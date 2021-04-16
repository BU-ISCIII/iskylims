


def graphic_3D_pie (heading, sub_title, axis_x_description, axis_y_description, theme, source_data) :
    data_source = {}
    data_source['chart'] = {
                "caption": heading,
                "subCaption": sub_title,
                "xAxisName": axis_x_description,
                "yAxisName": axis_y_description,
                "theme": theme,
                "numberPrefix": "",
                "showlegend": "1",
                "placevaluesInside": "1",
                "showpercentvalues": "1",
                "rotatevalues": "1",
                #Showing canvas bg to apply background color
                "showCanvasBg": "1",
                #Shwoing canvas base to apply base color
                "showCanvasBase": "1",
                #Changing canvas base depth
                "canvasBaseDepth": "14",
                #Changing canvas background depth
                "canvasBgDepth": "5",
                #Changing canvas base color
                "canvasBaseColor": "#aaaaaa",
                #Changing canvas background color
                "canvasBgColor": "#eeeeee",
                "exportEnabled": "1"
            }

    data =[]

    for key , values in source_data.items() :
        data_dict = {}
        data_dict['label'] = key
        data_dict['value'] = values
        data.append(data_dict)
    data_source['data'] = data
    return data_source


def column_graphic_dict (heading, sub_caption, x_axis_name, y_axis_name, theme, source_data) :
    data_source = {}
    data_source['chart'] = {
                "caption": heading,
                "subCaption": sub_caption,
                "xAxisName": x_axis_name,
                "yAxisName": y_axis_name,
                "theme": theme,
                "numberPrefix": "",
                "placevaluesInside": "1",
                "rotatevalues": "1",
                #Showing canvas bg to apply background color
                "showCanvasBg": "1",
                #Shwoing canvas base to apply base color
                "showCanvasBase": "1",
                #Changing canvas base depth
                "canvasBaseDepth": "14",
                #Changing canvas background depth
                "canvasBgDepth": "5",
                #Changing canvas base color
                "canvasBaseColor": "#aaaaaa",
                #Changing canvas background color
                "canvasBgColor": "#eeeeee",
                "exportEnabled": "1"
            }

    data =[]

    for key , values in source_data.items() :
        data_dict = {}
        data_dict['label'] = key
        data_dict['value'] = values
        data.append(data_dict)
    data_source['data'] = data

    return data_source

def column_graphic_tupla (heading, sub_caption, x_axis_name, y_axis_name, theme, source_data) :
    data_source = {}
    data_source['chart'] = {
                "caption": heading,
                "subCaption": sub_caption,
                "xAxisName": x_axis_name,
                "yAxisName": y_axis_name,
                "theme": theme,
                "numberPrefix": "",
                "placevaluesInside": "1",
                "rotatevalues": "1",
                #Showing canvas bg to apply background color
                "showCanvasBg": "1",
                #Shwoing canvas base to apply base color
                "showCanvasBase": "1",
                #Changing canvas base depth
                "canvasBaseDepth": "14",
                #Changing canvas background depth
                "canvasBgDepth": "5",
                #Changing canvas base color
                "canvasBaseColor": "#aaaaaa",
                #Changing canvas background color
                "canvasBgColor": "#eeeeee",
                "exportEnabled": "1"
            }

    data =[]

    for key , values in source_data :
        data_dict = {}
        data_dict['label'] = key
        data_dict['value'] = int(values)
        data.append(data_dict)
    data_source['data'] = data

    return data_source

def column_graphic_per_time (heading, sub_caption, x_axis_name, y_axis_name, time_values , service_values):
    data_source = {}
    # Chart data is passed to the `dataSource` parameter, as hashes, in the form of
    # key-value pairs.
    data_source['chart'] = {
        "caption": heading,
        "subCaption": sub_caption,
        "xAxisName": x_axis_name,
        "yAxisName": y_axis_name,
        #"theme": theme,
        "palette": "2",
        "numberprefix": "",
        "showvalues": "0",
        "legendshadow": "0",
        "legendborderalpha": "0",
        #"yaxismaxvalue": "90",
        "legendbgcolor": "FFFFFF",
        "exportEnabled": "1"
    }

    category_list = []
    for key  in time_values :
        category_list.append({ "label": key, "stepSkipped": 'false',"appliedSmartLabel": 'true'})
    data_source["categories"] = [{"category": category_list}]
    #index_color = ["005476", "0054dc","a1fddc", "a1c74a", "9966ff", "ffcc66","cc8800","ccff66", "0086b3" , "5c5c8a", "cc6699", "006699"]
    index_color = ["005476", "0054dc","a1fddc", "a1c74a", "9966ff", "ffcc66","cc8800","ccff66", "0086b3"  , "5c5c8a"]
    data_set_list =[]
    counter = 0
    for key ,values in service_values.items():

        series_name_list =[]
        for date in time_values :
            series_name_list.append({"value" : service_values[key][date]})
        if counter > 10:
            import pdb; pdb.set_trace()
        data_set_list.append({"seriesname": key, "color": index_color[counter],'data' : series_name_list})
        counter +=1
        #import pdb; pdb.set_trace()
        if counter > len(index_color):
            counter = 0

    data_source ["dataset"] =data_set_list

    return data_source
