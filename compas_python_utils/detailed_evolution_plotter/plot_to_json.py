import json
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import _get_dash_pattern, _scale_dashes
from matplotlib.colors import to_hex, to_rgba
from compas_python_utils.detailed_evolution_plotter.plot_detailed_evolution import run_main_plotter


class NumpyEncoder(json.JSONEncoder):
    """JSON Encoder to help with numpy types"""

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return round(float(obj), 6)
        elif isinstance(obj, np.ndarray):
            if np.issubdtype(obj.dtype, np.floating):
                obj = np.round(obj, decimals=6)
            return obj.tolist()
        elif isinstance(obj, np.bool_):
            return bool(obj)
        return json.JSONEncoder.default(self, obj)


def _artist_is_line(artist):
    return isinstance(artist, plt.Line2D) and not artist.get_label().startswith(
        "_child"
    )


def _artist_is_ref_line(artist):
    return isinstance(artist, plt.Line2D) and artist.get_label().startswith("_child")


def _artist_is_text(artist):
    return isinstance(artist, plt.Text) and artist.get_text()


def replace_latex_commands(latex_expr):
    commands = [(r"\\odot", "⊙"), (r"\\;", " "), (r"\s+", " ")]
    for pattern, replacement in commands:
        latex_expr = re.sub(pattern, replacement, latex_expr)
    return latex_expr


def parse_latex_math(latex_expr):
    """Use regex to parse latex maths into tokens of subscript, superscript or neither

    Parameters
    ----------
    latex_expr : str
        The input latex maths expression

    Returns
    -------
    list
        List of tuples, where the first item of the tuple is one of 'superscript', 'subscript' or 'text',
        and the second item of the tuple is the corresponding part of the expression
    """
    pattern = (
        r"(?P<subscript>_(?P<sub_single>[^{_\\\^}])|_{(?P<sub_braces>.+?)}|_(?P<sub_command>\\[a-zA-Z]+))|"
        r"(?P<superscript>\^(?P<sup_single>[^{_\\\^}])|\^{(?P<sup_braces>.+?)}|\^(?P<sup_command>\\[a-zA-Z]+))|"
        r"(?P<text>[^{_\^}]+)"
    )

    groups = []
    for match in re.finditer(pattern, latex_expr):
        if match.group("subscript"):
            for group in ["sub_single", "sub_braces", "sub_command"]:
                if match.group(group):
                    groups.append(
                        ("subscript", replace_latex_commands(match.group(group)))
                    )
        elif match.group("superscript"):
            for group in ["sup_single", "sup_braces", "sup_command"]:
                if match.group(group):
                    groups.append(
                        ("superscript", replace_latex_commands(match.group(group)))
                    )
        elif match.group("text"):
            groups.append(("text", replace_latex_commands(match.group("text"))))

    return groups


def strip_latex(input_string):
    """Takes a string, which may have latex maths expressions inside (i.e. surrounded by $),
    and returns a list of tokens determining whether the text should be sub-/superscript,
    and replacing some latex commands such as \\odot

    Parameters
    ----------
    input_string : str
        String to be parsed

    Returns
    -------
    list
        List of tuples, where the first item of the tuple is one of 'superscript', 'subscript' or 'text',
        and the second item of the tuple is the corresponding part of the string
    """
    split = re.findall(r"\$([^$]+)\$|([^$]+)", input_string)

    parsed = []
    for part in split:
        if part[0]:
            parsed.extend(parse_latex_math(part[0]))
        elif part[1]:
            parsed.append(("text", part[1]))
    return parsed


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
        _, dashes = _scale_dashes(
            *_get_dash_pattern(line.get_linestyle()), line.get_linewidth()
        )

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
        "type": "data",
    }


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

    return {
        "meta": ref_line_meta,
        "data": [
            {"x": xs[0], "y": ys[0]},
            {"x": xs[-1], "y": ys[-1]},
        ],
    }


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
        "meta": {
            "label": strip_latex(text.get_text()),
            "colour": get_artist_colour(text),
        },
        "data": {
            "x": text.get_position()[0],
            "y": text.get_position()[1],
        },
    }


def get_line_groups(lines, label):
    """Takes a list of Line2D objects and organises them into groups based on whether or not they have the same
    x data. Each group contains the x data for the group, a list of y axis data for each line in the group, and
    the metadata for each line.


    Parameters
    ----------
    lines : list
        List of pyplot Line2D objects
    label : str
        Label of the plot, used to give unique labels for the data

    Returns
    -------
    list
        List of groups containing the x data of the group and the y data and metadata for all lines in the group
    """
    groups = []
    for i, line in enumerate(lines):
        x_data = line.get_xdata()
        y_data = line.get_ydata()
        meta = get_line_meta(line, f"y{i}-{label}")
        for group in groups:
            if np.array_equal(group["x_data"], x_data):
                group["y_data"].append(y_data)
                group["meta"].append(meta)
                break
        else:
            groups.append({"x_data": x_data, "y_data": [y_data], "meta": [meta]})
    return groups


def get_plot_data(axes_map):
    """Takes a pyplot Figure instance and outputs JSON data to render the plots on the GWLandscape service

    Parameters
    ----------
    axes_map : list
        List up tuples, where first item in each tuple is the label of the plot,
        and the second item is the pyplot.Axes object

    Returns
    -------
    dict
        Dictionary containing necessary data to render the plots on a webpage
    """
    json_data = {"plots": []}

    for label, ax in axes_map:
        x_inverted = ax.xaxis_inverted()
        if x_inverted:
            ax.invert_xaxis()

        y_inverted = ax.yaxis_inverted()
        if y_inverted:
            ax.invert_yaxis()

        artists = ax.get_children()
        ref_lines = [
            get_ref_line_data(ref_line, f"refLine{i}")
            for i, ref_line in enumerate(filter(_artist_is_ref_line, artists))
        ]
        texts = [get_text_data(text) for text in filter(_artist_is_text, artists)]

        line_groups = get_line_groups(filter(_artist_is_line, artists), label)

        groups = []
        for line_group in line_groups:
            group = {"data": [], "meta": line_group["meta"]}
            for i, x in enumerate(line_group["x_data"]):
                row = {
                    line["yKey"]: line_group["y_data"][j][i]
                    for j, line in enumerate(group["meta"])
                }
                row[group["meta"][0]["xKey"]] = x
                group["data"].append(row)
            groups.append(group)

        json_data["plots"].append(
            {
                "meta": {
                    "label": label,
                    "xAxis": {
                        "label": strip_latex(ax.get_xlabel()),
                        "ticks": ax.get_xticks(),
                        "limits": ax.get_xlim(),
                        "scale": ax.get_xscale(),
                        "inverted": x_inverted,
                    },
                    "yAxis": {
                        "label": strip_latex(ax.get_ylabel()),
                        "ticks": ax.get_yticks(),
                        "limits": ax.get_ylim(),
                        "scale": ax.get_yscale(),
                        "inverted": y_inverted,
                    },
                },
                "groups": groups,
                "refLines": ref_lines,
                "texts": texts,
            }
        )

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
                "eventChar": chr(ord("@") + 1 + i),
                "time": event.time,
                "a": [event.aprev, event.a],
                "m1": [event.m1prev, event.m1],
                "m2": [event.m2prev, event.m2],
                "eventString": event.eventString,
                "imageNum": event.image_num,
                "flipImage": event.rotate_image,
            }
            for i, event in enumerate(events)
        ]
    }


def get_plot_json(data_path):
    """Get a JSON string containing the information needed to render line and VDH plots on the
    GWLandscape service for a specific COMPAS output file

    Parameters
    ----------
    data_path : str or Path
        Path to the COMPAS output file

    Returns
    -------
    str
        JSON string containing 
    """
    detailed_fig, _, events = run_main_plotter(data_path, outdir=None, show=False)
    axes = detailed_fig.get_axes()
    plots_data = get_plot_data([('mass_plot', axes[0]), ('length_plot', axes[1]), ('hr_plot', axes[3])])
    events_data = get_events_data(events)
    return json.dumps({**plots_data, **events_data}, cls=NumpyEncoder)
    