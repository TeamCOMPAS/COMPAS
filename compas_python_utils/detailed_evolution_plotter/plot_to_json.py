
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import _get_dash_pattern, _scale_dashes
from matplotlib.colors import to_hex, to_rgba

class NumpyEncoder(json.JSONEncoder):
    """ Special json encoder for numpy types """
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

def _artist_is_line(artist):
    return isinstance(artist, plt.Line2D) and not artist.get_label().startswith("_child")

def _artist_is_ref_line(artist):
    return isinstance(artist, plt.Line2D) and artist.get_label().startswith("_child")

def _artist_is_text(artist):
    return isinstance(artist, plt.Text) and artist.get_text()

def get_artist_colour(artist):
    """Get the colour of the artist in hex format

    Parameters
    ----------
    artist : pyplot.Line2D or pyplot.Text
        Pyplot artist for which to obtain the colour

    Returns
    -------
    str
        String with the colour in hex format
    """    
    return to_hex(to_rgba(artist.get_color(), artist.get_alpha()), keep_alpha=True)

def get_ref_line_data(ref_line, label):
    """Obtain the data and metadata needed to render reference lines on the GWLandscape service

    Parameters
    ----------
    ref_line : pyplot.Line2D
        Reference line object to be rendered
    label : str
        Name of the reference line

    Returns
    -------
    dict
        Dictionary with the data and metadata needed to render the reference line
    """    
    ref_line_meta = get_line_meta(ref_line, "y", label)

    xs, ys = ref_line.get_xdata(), ref_line.get_ydata()

    if xs[0] == xs[-1] and len(set(xs)) == 1:
        ref_line_meta["type"] = "vline"
    elif ys[0] == ys[-1] and len(set(ys)) == 1:
        ref_line_meta["type"] = "hline"
    else:
        ref_line_meta["type"] = "ref"

    ref_line_meta["points"] = [
        {"x": xs[0], "y": ys[0]},
        {"x": xs[-1], "y": ys[-1]},
    ]
    return ref_line_meta

def get_text_data(text):
    """Obtain the data needed to render text elements on the GWLandscape service

    Parameters
    ----------
    text : pyplot.Text
        Text object to be rendered

    Returns
    -------
    dict
        Dictionary containing the data and metadata necessary to render the text elements
    """    
    return {
        "label": text.get_text(),
        "x": text.get_position()[0],
        "y": text.get_position()[1],
        "colour": get_artist_colour(text)
    }

def get_line_dashes(line):
    """Obtain the dash pattern of a line. This uses a private attribute of the Line2D class or
    a private function from pyplot

    Parameters
    ----------
    line : _type_pyplot.Line2D
        Line object for which to obtain the dash pattern

    Returns
    -------
    str or None
        String of numbers separated by spaces representing the lengths of dashes and spaces in the pattern
    """    
    if hasattr(line, "_dash_pattern"):
        _, dashes = line._dash_pattern
    else:
        _, dashes = _scale_dashes(*_get_dash_pattern(line.get_linestyle()), line.get_linewidth())
    
    return " ".join(map(str, dashes)) if dashes else None

def get_line_meta(line, y_key, label=None):
    """Get the metadata for a line

    Parameters
    ----------
    line : pyplot.Line2D
        A line from a line plot
    y_key : str
        The key for the data, which will be used to identify it on the frontend
    label : str, optional
        The name of the line, by default None. If None, the label will be obtained from the Line2D object

    Returns
    -------
    dict
        Dictionary containing the metadata necessary to render the line properly on the frontend
    """    
    return {
        "colour": get_artist_colour(line),
        "dashes": get_line_dashes(line),
        "width": line.get_linewidth(),
        "xKey": "x",
        "yKey": y_key,
        "label": line.get_label() if label is None else label,
        "type": "data"
    }

def get_line_groups(lines):
    """Takes a list of Line2D objects and organises them into groups based on whether or not they have the same
    x data. Each group contains the x data for the group, a list of y axis data for each line in the group, and
    the metadata for each line.


    Parameters
    ----------
    lines : list
        List of pyplot Line2D objects

    Returns
    -------
    list
        List of groups containing the x data of the group and the y data and metadata for all lines in the group
    """    
    groups = []
    for i, line in enumerate(lines):
        x_data = line.get_xdata()
        y_data = line.get_ydata()
        meta = get_line_meta(line, f"y{i}")
        for group in groups:
            if np.array_equal(group["x_data"], x_data):
                group["y_data"].append(y_data)
                group["meta"].append(meta)
                break
        else:
            groups.append({"x_data": x_data, "y_data": [y_data], "meta": [meta]})
    return groups

def get_plot_data(fig):
    """Takes a pyplot Figure instance and outputs JSON data to render the plots on the GWLandscape service

    Parameters
    ----------
    fig : pyplot.Figure
        Pyplot Figure to replicate

    Returns
    -------
    dict
        Dictionary containing necessary data to render the plots on a webpage
    """    
    json_data = {
        "plots": {}
    }

    for ax in fig.get_axes():
        if not hasattr(ax, "tag"):
            continue

        if ax.xaxis_inverted():
            ax.invert_xaxis()

        artists = ax.get_children()
        ref_lines = [
            get_ref_line_data(ref_line, f"refLine{i}")
            for i, ref_line in enumerate(filter(_artist_is_ref_line, artists))
        ]
        texts = [get_text_data(text) for text in filter(_artist_is_text, artists)]

        line_groups = get_line_groups(filter(_artist_is_line, artists))

        groups = []
        for line_group in line_groups:
            group = {"data": [], "meta": line_group["meta"]}
            for i, x in enumerate(line_group["x_data"]):
                row = {line["yKey"]: line_group["y_data"][j][i] for j, line in enumerate(group["meta"])}
                row[group["meta"][0]["xKey"]] = x
                group["data"].append(row)
            groups.append(group)

        json_data["plots"][ax.tag] = {
            "meta": {
                "xAxis": {
                    "label": ax.get_xlabel(),
                    "ticks": ax.get_xticks(),
                    "limits": ax.get_xlim(),
                    "scale": ax.get_xscale()
                },
                "yAxis": {
                    "label": ax.get_ylabel(),
                    "ticks": ax.get_yticks(),
                    "limits": ax.get_ylim(),
                    "scale": ax.get_yscale()
                },
            },
            "groups": groups,
            "refLines": ref_lines,
            "texts": texts
        }
    
    return json_data

def get_events_data(events):
    """Uses a list of Events to generate a JSON structure for rendering Van Den Heuvel plots on the GWLandscape service

    Parameters
    ----------
    events : list
        List of Events

    Returns
    -------
    dict
        Dictionary containing data necessary to render VDH plot on a webpage
    """    
    return {
        "events": [
            {
                "eventChar": chr(ord('@') + 1 + i),
                "time": event.time,
                "a": [event.aprev, event.a],
                "m1": [event.m1prev, event.m1],
                "m2": [event.m2prev, event.m2],
                "eventString": event.eventString,
                "imageNum": event.image_num,
                "flipImage": event.rotate_image,
            } for i, event in enumerate(events)
        ]
    }