def preparation_3D_pie(
    heading: str, sub_title: str, theme: str, source_data: dict
) -> dict:
    """Join the parameters to create a dictionary with two keys: chart and data.
    The chart key contains the heading, sub_title, theme, numberPrefix, showlegend, placevaluesInside, showpercentvalues, rotatevalues, showCanvasBg, showCanvasBase, canvasBaseDepth, canvasBgDepth, canvasBaseColor, canvasBgColor, exportEnabled.
    The data key contains a list of dictionaries with the keys label and value.

    Args:
        heading (str): title of the chart
        sub_title (str): additional information
        theme (str): theme of the chart
        source_data (dict): dictionary with the data to be displayed

    Returns:
        dict: dictionary with the keys chart and data
    """
    data_source = {}
    data_source["chart"] = {
        "caption": heading,
        "subCaption": sub_title,
        "theme": theme,
        "numberPrefix": "",
        "showlegend": "1",
        "placevaluesInside": "1",
        "showpercentvalues": "1",
        "rotatevalues": "1",
        # Showing canvas bg to apply background color
        "showCanvasBg": "1",
        # Shwoing canvas base to apply base color
        "showCanvasBase": "1",
        # Changing canvas base depth
        "canvasBaseDepth": "14",
        # Changing canvas background depth
        "canvasBgDepth": "5",
        # Changing canvas base color
        "canvasBaseColor": "#aaaaaa",
        # Changing canvas background color
        "canvasBgColor": "#eeeeee",
        "exportEnabled": "1",
    }

    data = []

    for key, values in source_data.items():
        data_dict = {}
        data_dict["label"] = key
        data_dict["value"] = values
        data.append(data_dict)
    data_source["data"] = data
    return data_source


def preparation_graphic_data(
    heading: str,
    sub_caption: str,
    x_axis_name: str,
    y_axis_name: str,
    theme: str,
    input_data: dict,
    label_key: str = None,
    label_value: str = None,
) -> dict:
    """Join the parameters to create a dictionary with two keys: chart and data.

    Args:
        heading (str): heading of the chart
        sub_caption (str): sub_caption of the chart
        x_axis_name (str): name of the x axis
        y_axis_name (str): name of the y axis
        theme (str): theme of the chart
        input_data (dict): data to be displayed
        label_key (str, optional): key name from the input dict. Defaults to None.
        label_value (str, optional): value field in the innput dict. Defaults to None.

    Returns:
        dict: dictionary with the keys chart and data
    """
    data_source = {}
    data_source["chart"] = {
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
        # Showing canvas bg to apply background color
        "showCanvasBg": "1",
        # Shwoing canvas base to apply base color
        "showCanvasBase": "1",
        # Changing canvas base depth
        "canvasBaseDepth": "14",
        # Changing canvas background depth
        "canvasBgDepth": "5",
        # Changing canvas base color
        "canvasBaseColor": "#aaaaaa",
        # Changing canvas background color
        "canvasBgColor": "#eeeeee",
        "exportEnabled": "1",
    }

    data = []
    if isinstance(input_data, dict):
        for key, values in input_data.items():
            data_dict = {}
            data_dict["label"] = key
            data_dict["value"] = values
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
                        values = round(values, 2)
                    data_dict["label"] = key
                    data_dict["value"] = values
                    data.append(data_dict)
        else:
            # process dictionary when item list contains 2 dictionnaries
            # one for the key and the other for the value
            # label_key and label_value are used to map correctly
            for dict_item in input_data:
                data_dict = {}
                values = dict_item[label_value]
                if isinstance(values, float):
                    values = round(values, 2)
                data_dict["label"] = dict_item[label_key]
                data_dict["value"] = values
                data.append(data_dict)
    data_source["data"] = data
    return data_source
