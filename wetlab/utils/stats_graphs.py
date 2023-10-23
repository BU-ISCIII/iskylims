def json_unknow_barcode_graphic(heading, barcode_data):
    data_source = {}

    # Chart data is passed to the `dataSource` parameter, as hashes, in the form of
    # key-value pairs.
    data_source["chart"] = {
        "caption": heading,
        "subcaption": "Top Sequence found in the Run",
        "startingangle": "120",
        "showlabels": "0",
        "showlegend": "1",
        "enablemultislicing": "0",
        "slicingdistance": "15",
        "showpercentvalues": "1",
        "showpercentintooltip": "0",
        "plottooltext": "Sequence : $label Total count : $datavalue",
        "theme": "ocean",
        "exportEnabled": "1",
    }
    data = []
    for key, values in barcode_data.items():
        data_dict = {}
        data_dict["label"] = key
        data_dict["value"] = values
        data.append(data_dict)
    data_source["data"] = data

    return data_source


def graphic_for_unbarcodes(heading, theme, lane_unbarcode):
    data_source = {}
    data_source["chart"] = {
        "caption": heading,
        "subCaption": "",
        "xAxisName": "Sequence",
        "yAxisName": "Number undetermined barcode",
        "theme": theme,
        "numberPrefix": "",
        "placevaluesInside": "1",
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
        "canvasBgColor": "#8c8c8c",
        "exportEnabled": "1",
    }

    data = []

    for key, values in lane_unbarcode.items():
        data_dict = {}
        data_dict["label"] = key
        data_dict["value"] = values
        data.append(data_dict)
    data_source["data"] = data

    return data_source


def graphic_for_library_kit(
    heading, sub_caption, x_axis_name, y_axis_name, theme, lane_quality
):
    data_source = {}
    data_source["chart"] = {
        "caption": heading,
        "subCaption": sub_caption,
        "xAxisName": x_axis_name,
        "yAxisName": y_axis_name,
        # "theme": "fint",
        "theme": theme,
        "numberPrefix": "",
        "placevaluesInside": "1",
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

    for key, values in lane_quality.items():
        data_dict = {}
        data_dict["label"] = key
        data_dict["value"] = values
        data.append(data_dict)
    data_source["data"] = data

    return data_source


def pie_graphic(heading, theme, top_count_sequence):
    data_source = {}

    data_source["chart"] = {
        "caption": heading,
        "subcaption": "Top Sequence found in the Run",
        "startingangle": "120",
        "showlabels": "0",
        "showlegend": "1",
        "enablemultislicing": "0",
        "slicingdistance": "15",
        "showpercentvalues": "1",
        "showpercentintooltip": "0",
        "plottooltext": "Sequence : $label Total count : $datavalue",
        "theme": theme,
        "exportEnabled": "1",
    }
    data = []
    for key, values in top_count_sequence.items():
        data_dict = {}
        data_dict["label"] = key
        data_dict["value"] = values
        data.append(data_dict)
    data_source["data"] = data

    return data_source


########
# JSON for the quality of sample/project/Run
########
def graphic_for_quality_angular(heading, value):
    data_source = {}
    data_source["chart"] = {
        "caption": heading,
        "lowerlimit": "0",
        "upperlimit": "100",
        "lowerlimitdisplay": "Bad",
        "upperlimitdisplay": "Good",
        "palette": "1",
        "numbersuffix": "%",
        "tickvaluedistance": "10",
        "showvalue": "0",
        # inner white radious
        "gaugeinnerradius": "45",
        "bgcolor": "FFFFFF",
        "pivotfillcolor": "333333",
        "pivotradius": "8",
        "pivotfillmix": "333333, 333333",
        "pivotfilltype": "radial",
        "pivotfillratio": "0,100",
        "showtickvalues": "1",
        "showborder": "0",
    }
    data_source["colorrange"] = {
        "color": [
            {"minvalue": "0", "maxvalue": "45", "code": "e44a00"},
            {"minvalue": "45", "maxvalue": "75", "code": "f8bd19"},
            {"minvalue": "75", "maxvalue": "100", "code": "6baa01"},
        ]
    }
    data_source["dials"] = {
        "dial": [
            {
                "value": value,
                "rearextension": "25",
                # length fo the arrow
                "radius": "80",
                "bgcolor": "333333",
                "bordercolor": "333333",
                "basewidth": "8",
            }
        ]
    }
    return data_source


def pie_graphic_standard(heading, subcaption, theme, input_values):
    data_source = {}
    data_source["chart"] = {
        "caption": heading,
        "subcaption": subcaption,
        "startingangle": "120",
        "showlabels": "0",
        "showlegend": "1",
        "enablemultislicing": "0",
        "slicingdistance": "15",
        "showpercentvalues": "1",
        "showpercentintooltip": "0",
        "plottooltext": "Sequence : $label Total count : $datavalue",
        "theme": theme,
        "exportEnabled": "1",
    }
    data = []
    for key, values in input_values.items():
        data_dict = {}
        data_dict["label"] = key
        data_dict["value"] = values
        data.append(data_dict)
    data_source["data"] = data

    return data_source


def column_graphic_for_year_report(
    heading, sub_caption, x_axis_name, y_axis_name, theme, year_report_data
):
    data_source = {}
    data_source["chart"] = {
        "caption": heading,
        "subCaption": sub_caption,
        "xAxisName": x_axis_name,
        "yAxisName": y_axis_name,
        # "theme": "fint",
        "theme": theme,
        "numberPrefix": "",
        "placevaluesInside": "1",
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

    for key, values in year_report_data.items():
        data_dict = {}
        data_dict["label"] = key
        data_dict["value"] = values
        data.append(data_dict)
    data_source["data"] = data

    return data_source


def researcher_project_column_graphic(
    heading, sub_caption, x_axis_name, y_axis_name, theme, lane_report_data
):
    data_source = {}
    data_source["chart"] = {
        "caption": heading,
        "subCaption": sub_caption,
        "xAxisName": x_axis_name,
        "yAxisName": y_axis_name,
        # "theme": "fint",
        "theme": theme,
        "numberPrefix": "",
        "placevaluesInside": "1",
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

    for key, values in lane_report_data.items():
        data_dict = {}
        data_dict["label"] = key
        data_dict["value"] = values
        data.append(data_dict)
    data_source["data"] = data

    return data_source


def column_graphic_one_column_highligthed(
    heading,
    sub_caption,
    x_axis_name,
    y_axis_name,
    theme,
    percentage_in_project,
    sample_name,
):
    data_source = {}
    data_source["chart"] = {
        "caption": heading,
        "subCaption": sub_caption,
        "xAxisName": x_axis_name,
        "yAxisName": y_axis_name,
        "theme": theme,
        "numberPrefix": "",
        "placevaluesInside": "1",
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

    for key, values in percentage_in_project.items():
        data_dict = {}
        data_dict["label"] = key
        data_dict["value"] = values
        if key == sample_name:
            data_dict["dashed"] = "1"
            data_dict["color"] = "ff0000"
        data.append(data_dict)
    data_source["data"] = data

    return data_source


def bloxplot_graphic(
    heading, sub_caption, x_axis_name, y_axis_name, theme, categories, series, data
):
    """
    Description:
        The function will prepare the input data into the json format
        to be used in the graphic.
    Input:
        heading      # contains the title name to be present in the graphic
        sub_caption  # contains the sub title name to be present in the graphic
        x_axis_name  # contains the title to display in X axis
        y_axis_name  # contains the title to display in Y axis
        theme        # contains the theme used in the graphic
        categories   # contains the categories name values
        series       # contains the series name
        data         # contains the data for the series

    Variables:
        category        # dictionary to contain the category list
        category_list   # list of the categories to display
        dataset         # list with all the values for each category
        label_category  # title name for each category
    return:
        data_source
    """
    data_source = {}
    data_source["chart"] = (
        {
            "theme": theme,
            "caption": heading,
            "subcaption": sub_caption,
            "xAxisName": x_axis_name,
            "YAxisName": y_axis_name,
            # "numberPrefix": "%",
            "legendBorderAlpha": "1",
            "legendShadow": "0",
            "legendPosition": "botom",
            "showValues": "0",
            "toolTipColor": "#ffffff",
            "toolTipBorderThickness": "0",
            "toolTipBgColor": "#000000",
            "toolTipBgAlpha": "80",
            "toolTipBorderRadius": "2",
            "toolTipPadding": "5",
            "yAxisMaxValue": "0.4",
            "showMean": "1",
            "exportEnabled": "1",
        },
    )

    category = {}
    category_list = []
    for item in categories:
        label_category = {}
        label_category["label"] = item
        category_list.append(label_category)
    category["category"] = category_list

    data_source["categories"] = [category]

    dataset = []

    for serie in range(len(series)):
        dataset_dict = {}
        dataset_dict["seriesname"] = series[serie][0]
        dataset_dict["lowerBoxColor"] = series[serie][1]
        dataset_dict["upperBoxColor"] = series[serie][2]
        serie_data = []
        for item in data[serie]:
            value = {}
            value["value"] = item
            serie_data.append(value)
        dataset_dict["data"] = serie_data
        dataset.append(dataset_dict)
    data_source["dataset"] = dataset

    return data_source


def column_graphic_simple(
    heading, sub_caption, x_axis_name, y_axis_name, theme, input_data
):
    data_source = {}
    data_source["chart"] = {
        "caption": heading,
        "subCaption": sub_caption,
        "xAxisName": x_axis_name,
        "yAxisName": y_axis_name,
        # "theme": "fint",
        "theme": theme,
        "numberPrefix": "",
        "placevaluesInside": "1",
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

    for key, values in input_data.items():
        data_dict = {}
        data_dict["label"] = key
        data_dict["value"] = values
        data.append(data_dict)
    data_source["data"] = data
    return data_source


def column_graphic_tupla(
    heading, sub_caption, x_axis_name, y_axis_name, theme, source_data, highlight_value
):
    data_source = {}
    data_source["chart"] = {
        "caption": heading,
        "subCaption": sub_caption,
        "xAxisName": x_axis_name,
        "yAxisName": y_axis_name,
        "theme": theme,
        "numberPrefix": "",
        "placevaluesInside": "1",
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

    for key, values in source_data:
        data_dict = {}
        data_dict["label"] = key
        data_dict["value"] = float(values)
        # data_dict['value'] = int(values)
        if key == highlight_value:
            data_dict["dashed"] = "1"
            data_dict["color"] = "ff0000"
        data.append(data_dict)
    data_source["data"] = data

    return data_source


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
