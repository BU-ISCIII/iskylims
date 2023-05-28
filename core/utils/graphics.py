def preparation_3D_pie (heading, sub_title, theme, source_data) :
    """_summary_

    Parameters
    ----------
    heading : str
        text to insert on top of graphic
    sub_title : str
        additional text to include in graphic
    axis_x_description : str
        description of x axis
    axis_y_description : str
        description
    theme : str
        name of the available themes to be used for the graphic
    source_data : dict
        Dictionary containing value for the graphic

    Returns
    -------
    _type_
        _description_
    """    
    data_source = {}
    data_source['chart'] = {
                "caption": heading,
                "subCaption": sub_title,
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

def preparation_graphic_data (heading, sub_caption, x_axis_name, y_axis_name, theme, input_data, label_key=None, label_value=None) :
    data_source = {}
    data_source['chart'] = {
                "caption": heading,
                "subCaption": sub_caption,
                "xAxisName": x_axis_name,
                "yAxisName": y_axis_name,
                "labelDisplay": "rotate",
                "useEllipsesWhenOverflow": "0",
                "theme": theme,
                "numberPrefix": "",
                "placevaluesInside": "1",
                "rotatevalues": "1",
                "showlegend": "1",
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
    if isinstance(input_data, dict):
        for key , values in input_data.items() :
            data_dict = {}
            data_dict['label'] = key
            data_dict['value'] = values
            data.append(data_dict)
    # converting data when is a list of dictionnary
    else:
        if label_key is None:
            # process dictionnary when each item in list contains only
            # key / value  
            for dict_item in input_data:
                for key, values in dict_item.items():
                    data_dict = {}
                    if isinstance(values, float):
                        values = round(values,2)
                    data_dict['label'] = key
                    data_dict['value'] = values
                    data.append(data_dict)
        else:
            # process dictionary when item list contains 2 dictionnaries
            # one for the key and the other for the value
            # label_key and label_value are used to map correctly
            for dict_item in input_data:
                data_dict = {}
                values = dict_item[label_value]
                if isinstance(values, float):
                    values = round(values,2)
                data_dict["label"] = dict_item[label_key]
                data_dict["value"] = values
                data.append(data_dict)
    data_source['data'] = data
    return data_source