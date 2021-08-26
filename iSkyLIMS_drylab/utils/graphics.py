


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

def graphic_multi_level_pie (heading, plottooltext, tooltext, theme, colors, source_data):
    data_source = {}
    data_source['chart'] = {
            "caption": heading,
            "subCaption": '',
            "showplotborder": "1",
            "plotfillalpha": "60",
            "hoverfillcolor": "#CCCCCC",
            "numberprefix": "$",
            "plottooltext":  plottooltext ,
            "theme": theme,
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

    total_value = 0
    color_index = 0
    category_root = []
    root_data ={}
    for key in source_data.keys() :
        data = []
        value_in_cat = 0
        category_parent = {}
        for sub_key in source_data[key].keys():
            sub_data = {'label': sub_key}
            sub_data['color'] = colors[color_index]
            value_in_cat += source_data[key][sub_key]
            sub_data['value'] = str(source_data[key][sub_key])
            data.append(sub_data)

        total_value += value_in_cat
        category_parent['category'] = data
        category_parent['value'] = str(value_in_cat)
        category_parent['color'] = colors[color_index]
        category_parent['label'] = key
        category_root.append(category_parent)
        color_index += 1
    root_data['label'] = "Pending Services"
    root_data['tooltext'] = tooltext
    root_data['color'] = "#ffffff"
    root_data['value'] = total_value
    root_data['category'] = category_root
    data_source['category'] = [root_data]

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
    index_color = ["005476", "0054dc","a1fddc", "a1c74a", "9966ff", "ffcc66","cc8800","ccff66", "0086b3" , "5c5c8a", "cc6699", "006699"]
    data_set_list =[]
    counter = 0
    for key ,values in service_values.items():

        series_name_list =[]
        for date in time_values :
            series_name_list.append({"value" : service_values[key][date]})
        data_set_list.append({"seriesname": key, "color": index_color[counter],'data' : series_name_list})
        counter +=1

        if counter >= len(index_color):
            counter = 0

    data_source ["dataset"] =data_set_list

    return data_source
