def graphic_3D_pie(
    heading, sub_title, axis_x_description, axis_y_description, theme, source_data
):
    data_source = {}
    data_source["chart"] = {
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
