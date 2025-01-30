import pytest

import matplotlib.pyplot as plt
import numpy as np
import json

from compas_python_utils.detailed_evolution_plotter import plot_to_json


@pytest.fixture
def setup_fig_ax():
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel("Test x label")
    ax.set_ylabel("Test y label")
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
    ax.set_yscale("log")
    yield fig, ax
    plt.close("all")


@pytest.fixture
def plot_lines(setup_fig_ax):
    _, ax = setup_fig_ax
    return [
        ax.plot([1, 2, 3], [1, 2, 3], "r--", lw=1, label="Test line 1")[0],
        ax.plot([2, 4, 6], [1, 2, 3], "b-.", lw=1, label="Test line 2")[0],
        ax.plot([1, 2, 3], [2, 4, 6], "g", lw=1, label="Test line 3")[0],
    ]


@pytest.fixture
def plot_ref_lines(setup_fig_ax):
    _, ax = setup_fig_ax
    return [
        ax.plot([0.5, 1], [1, 1], "r", dashes=[1, 2], lw=1)[0],
        ax.plot([1, 1], [0.5, 1], "b", dashes=[1, 2, 3], lw=1)[0],
        ax.plot([1, 2], [1, 2], "g", dashes=[], lw=1)[0],
    ]


@pytest.fixture
def plot_texts(setup_fig_ax):
    _, ax = setup_fig_ax
    return [
        ax.text(0.25, 0.25, "Test text 1", color="r"),
        ax.text(0.5, 0.5, "Test text 2", color="b"),
        ax.text(0.75, 0.75, "Test text 3", color="g"),
    ]


def test_replace_latex_commands():
    assert (
        plot_to_json.replace_latex_commands(r"string with \odot in it")
        == "string with ⊙ in it"
    )
    assert (
        plot_to_json.replace_latex_commands(r"string with \; in it")
        == "string with in it"
    )
    assert (
        plot_to_json.replace_latex_commands(r"string    with multiple  spaces in it")
        == "string with multiple spaces in it"
    )
    assert (
        plot_to_json.replace_latex_commands(r"replacing everything \odot \;")
        == "replacing everything ⊙ "
    )


def test_parse_latex_math():
    assert plot_to_json.parse_latex_math(r"normal text") == [("text", "normal text")]

    assert plot_to_json.parse_latex_math(r"^a") == [("superscript", "a")]
    assert plot_to_json.parse_latex_math(r"^{a}") == [("superscript", "a")]
    assert plot_to_json.parse_latex_math(r"^{x = 1 + 2}") == [
        ("superscript", "x = 1 + 2")
    ]
    assert plot_to_json.parse_latex_math(r"^{\odot}") == [("superscript", "⊙")]
    assert plot_to_json.parse_latex_math(r"^\odot") == [("superscript", "⊙")]

    assert plot_to_json.parse_latex_math(r"_a") == [("subscript", "a")]
    assert plot_to_json.parse_latex_math(r"_{a}") == [("subscript", "a")]
    assert plot_to_json.parse_latex_math(r"_{x = 1 + 2}") == [
        ("subscript", "x = 1 + 2")
    ]
    assert plot_to_json.parse_latex_math(r"_{\odot}") == [("subscript", "⊙")]
    assert plot_to_json.parse_latex_math(r"_\odot") == [("subscript", "⊙")]

    assert plot_to_json.parse_latex_math(r"all in a row_a^b") == [
        ("text", "all in a row"),
        ("subscript", "a"),
        ("superscript", "b"),
    ]
    assert plot_to_json.parse_latex_math(r"a_b and b^a") == [
        ("text", "a"),
        ("subscript", "b"),
        ("text", " and b"),
        ("superscript", "a"),
    ]
    assert plot_to_json.parse_latex_math(
        r"multiple_{subscripts} in_a single string_\odot"
    ) == [
        ("text", "multiple"),
        ("subscript", "subscripts"),
        ("text", " in"),
        ("subscript", "a"),
        ("text", " single string"),
        ("subscript", "⊙"),
    ]
    assert plot_to_json.parse_latex_math(
        r"multiple^{superscripts} in^a single string^\odot"
    ) == [
        ("text", "multiple"),
        ("superscript", "superscripts"),
        ("text", " in"),
        ("superscript", "a"),
        ("text", " single string"),
        ("superscript", "⊙"),
    ]

    # A genuine example
    assert plot_to_json.parse_latex_math(r"R_\odot^{+2}") == [
        ("text", "R"),
        ("subscript", "⊙"),
        ("superscript", "+2"),
    ]


def test_strip_latex():
    assert plot_to_json.strip_latex(r"Mass $/ \; M_{\odot}$") == [
        ("text", "Mass "),
        ("text", "/ M"),
        ("subscript", "⊙"),
    ]
    assert plot_to_json.strip_latex(r"Luminosity [log($L/L_\odot$)]") == [
        ("text", "Luminosity [log("),
        ("text", "L/L"),
        ("subscript", "⊙"),
        ("text", ")]"),
    ]


def test_artist_is_line(setup_fig_ax, plot_lines, plot_ref_lines, plot_texts):
    _, ax = setup_fig_ax

    for artist in ax.get_children():
        if artist in plot_lines:
            assert plot_to_json._artist_is_line(artist)
        else:
            assert not plot_to_json._artist_is_line(artist)


def test_artist_is_ref_line(setup_fig_ax, plot_lines, plot_ref_lines, plot_texts):
    _, ax = setup_fig_ax

    for artist in ax.get_children():
        if artist in plot_ref_lines:
            assert plot_to_json._artist_is_ref_line(artist)
        else:
            assert not plot_to_json._artist_is_ref_line(artist)


def test_artist_is_text(setup_fig_ax, plot_lines, plot_ref_lines, plot_texts):
    _, ax = setup_fig_ax

    for artist in ax.get_children():
        if artist in plot_texts:
            assert plot_to_json._artist_is_text(artist)
        else:
            assert not plot_to_json._artist_is_text(artist)


def test_get_artist_color(plot_lines, plot_ref_lines, plot_texts):
    for group in [plot_lines, plot_ref_lines, plot_texts]:
        # Check that colours are correctly converted to hex code
        assert plot_to_json.get_artist_colour(group[0]) == "#ff0000ff"
        assert plot_to_json.get_artist_colour(group[1]) == "#0000ffff"
        assert plot_to_json.get_artist_colour(group[2]) == "#008000ff"

        # Check that changes to alpha also show up
        group[0].set_alpha(0.75)
        group[1].set_alpha(0.5)
        group[2].set_alpha(0.25)
        assert plot_to_json.get_artist_colour(group[0]) == "#ff0000bf"
        assert plot_to_json.get_artist_colour(group[1]) == "#0000ff80"
        assert plot_to_json.get_artist_colour(group[2]) == "#00800040"


def test_get_line_dashes(plot_lines, plot_ref_lines):
    # Check that dash patterns from both the linestyle and the dashes argument are correctly collected
    assert plot_to_json.get_line_dashes(plot_lines[0]) == "3.7 1.6"
    assert plot_to_json.get_line_dashes(plot_lines[1]) == "6.4 1.6 1.0 1.6"
    assert plot_to_json.get_line_dashes(plot_lines[2]) == None

    assert plot_to_json.get_line_dashes(plot_ref_lines[0]) == "1.0 2.0"
    assert plot_to_json.get_line_dashes(plot_ref_lines[1]) == "1.0 2.0 3.0"
    assert plot_to_json.get_line_dashes(plot_ref_lines[2]) == None

    # Check that linewidth also affects the dash pattern
    plot_lines[0].set_linewidth(2)
    plot_ref_lines[0].set_linewidth(2)
    assert plot_to_json.get_line_dashes(plot_lines[0]) == "7.4 3.2"
    assert plot_to_json.get_line_dashes(plot_ref_lines[0]) == "2.0 4.0"


def test_get_line_meta(plot_lines):
    assert plot_to_json.get_line_meta(plot_lines[0], "y") == {
        "colour": "#ff0000ff",
        "dashes": "3.7 1.6",
        "width": 1.0,
        "xKey": "x",
        "yKey": "y",
        "label": "Test line 1",
        "type": "data",
    }
    assert plot_to_json.get_line_meta(plot_lines[1], "y") == {
        "colour": "#0000ffff",
        "dashes": "6.4 1.6 1.0 1.6",
        "width": 1.0,
        "xKey": "x",
        "yKey": "y",
        "label": "Test line 2",
        "type": "data",
    }
    assert plot_to_json.get_line_meta(plot_lines[2], "y") == {
        "colour": "#008000ff",
        "dashes": None,
        "width": 1.0,
        "xKey": "x",
        "yKey": "y",
        "label": "Test line 3",
        "type": "data",
    }


def test_get_ref_line_data(plot_ref_lines):
    assert plot_to_json.get_ref_line_data(plot_ref_lines[0], "Test ref line 1") == {
        "meta": {
            "colour": "#ff0000ff",
            "dashes": "1.0 2.0",
            "width": 1.0,
            "xKey": "x",
            "yKey": "y",
            "label": "Test ref line 1",
            "type": "hline",
        },
        "data": [
            {"x": 0.5, "y": 1.0},
            {"x": 1.0, "y": 1.0},
        ],
    }
    assert plot_to_json.get_ref_line_data(plot_ref_lines[1], "Test ref line 2") == {
        "meta": {
            "colour": "#0000ffff",
            "dashes": "1.0 2.0 3.0",
            "width": 1.0,
            "xKey": "x",
            "yKey": "y",
            "label": "Test ref line 2",
            "type": "vline",
        },
        "data": [
            {"x": 1.0, "y": 0.5},
            {"x": 1.0, "y": 1.0},
        ],
    }
    assert plot_to_json.get_ref_line_data(plot_ref_lines[2], "Test ref line 3") == {
        "meta": {
            "colour": "#008000ff",
            "dashes": None,
            "width": 1.0,
            "xKey": "x",
            "yKey": "y",
            "label": "Test ref line 3",
            "type": "ref",
        },
        "data": [
            {"x": 1.0, "y": 1.0},
            {"x": 2.0, "y": 2.0},
        ],
    }


def test_get_text_data(plot_texts):
    assert plot_to_json.get_text_data(plot_texts[0]) == {
        "meta": {"label": [("text", "Test text 1")], "colour": "#ff0000ff"},
        "data": {"x": 0.25, "y": 0.25},
    }
    assert plot_to_json.get_text_data(plot_texts[1]) == {
        "meta": {"label": [("text", "Test text 2")], "colour": "#0000ffff"},
        "data": {"x": 0.5, "y": 0.5},
    }
    assert plot_to_json.get_text_data(plot_texts[2]) == {
        "meta": {"label": [("text", "Test text 3")], "colour": "#008000ff"},
        "data": {"x": 0.75, "y": 0.75},
    }


def test_get_line_groups(plot_lines):
    # Can't just use bare assert because the data in the plot is in terms of numpy arrays and highly nested
    # Easier to just convert to a JSON string and then compare those, though this may fail if dicts
    # are returned to being unordered in the future
    returned = json.dumps(
        plot_to_json.get_line_groups(plot_lines), cls=plot_to_json.NumpyEncoder
    )
    expected = json.dumps(
        [
            {
                "x_data": np.array([1, 2, 3]),
                "y_data": [np.array([1, 2, 3]), np.array([2, 4, 6])],
                "meta": [
                    {
                        "colour": "#ff0000ff",
                        "dashes": "3.7 1.6",
                        "width": 1.0,
                        "xKey": "x",
                        "yKey": "y0",
                        "label": "Test line 1",
                        "type": "data",
                    },
                    {
                        "colour": "#008000ff",
                        "dashes": None,
                        "width": 1.0,
                        "xKey": "x",
                        "yKey": "y2",
                        "label": "Test line 3",
                        "type": "data",
                    },
                ],
            },
            {
                "x_data": np.array([2, 4, 6]),
                "y_data": [np.array([1, 2, 3])],
                "meta": [
                    {
                        "colour": "#0000ffff",
                        "dashes": "6.4 1.6 1.0 1.6",
                        "width": 1.0,
                        "xKey": "x",
                        "yKey": "y1",
                        "label": "Test line 2",
                        "type": "data",
                    }
                ],
            },
        ],
        cls=plot_to_json.NumpyEncoder,
    )
    assert returned == expected


def test_get_plot_data(setup_fig_ax, plot_lines, plot_ref_lines, plot_texts):
    returned = json.dumps(
        plot_to_json.get_plot_data([("test", setup_fig_ax[1])]),
        cls=plot_to_json.NumpyEncoder,
    )
    expected = json.dumps(
        {
            "plots": [
                {
                    "meta": {
                        "label": "test",
                        "xAxis": {
                            "label": [("text", "Test x label")],
                            "ticks": np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0]),
                            "limits": (np.float64(0.0), np.float64(1.0)),
                            "scale": "linear",
                            "inverted": np.False_,
                        },
                        "yAxis": {
                            "label": [("text", "Test y label")],
                            "ticks": np.array([0.01, 0.1, 1.0, 10.0]),
                            "limits": (np.float64(0.0), np.float64(1.0)),
                            "scale": "log",
                            "inverted": np.False_,
                        },
                    },
                    "groups": [
                        {
                            "data": [
                                {
                                    "y0": np.int64(1),
                                    "y2": np.int64(2),
                                    "x": np.int64(1),
                                },
                                {
                                    "y0": np.int64(2),
                                    "y2": np.int64(4),
                                    "x": np.int64(2),
                                },
                                {
                                    "y0": np.int64(3),
                                    "y2": np.int64(6),
                                    "x": np.int64(3),
                                },
                            ],
                            "meta": [
                                {
                                    "colour": "#ff0000ff",
                                    "dashes": "3.7 1.6",
                                    "width": 1.0,
                                    "xKey": "x",
                                    "yKey": "y0",
                                    "label": "Test line 1",
                                    "type": "data",
                                },
                                {
                                    "colour": "#008000ff",
                                    "dashes": None,
                                    "width": 1.0,
                                    "xKey": "x",
                                    "yKey": "y2",
                                    "label": "Test line 3",
                                    "type": "data",
                                },
                            ],
                        },
                        {
                            "data": [
                                {"y1": np.int64(1), "x": np.int64(2)},
                                {"y1": np.int64(2), "x": np.int64(4)},
                                {"y1": np.int64(3), "x": np.int64(6)},
                            ],
                            "meta": [
                                {
                                    "colour": "#0000ffff",
                                    "dashes": "6.4 1.6 1.0 1.6",
                                    "width": 1.0,
                                    "xKey": "x",
                                    "yKey": "y1",
                                    "label": "Test line 2",
                                    "type": "data",
                                }
                            ],
                        },
                    ],
                    "refLines": [
                        {
                            "meta": {
                                "colour": "#ff0000ff",
                                "dashes": "1.0 2.0",
                                "width": 1.0,
                                "xKey": "x",
                                "yKey": "y",
                                "label": "refLine0",
                                "type": "hline",
                            },
                            "data": [
                                {"x": np.float64(0.5), "y": np.int64(1)},
                                {"x": np.float64(1.0), "y": np.int64(1)},
                            ],
                        },
                        {
                            "meta": {
                                "colour": "#0000ffff",
                                "dashes": "1.0 2.0 3.0",
                                "width": 1.0,
                                "xKey": "x",
                                "yKey": "y",
                                "label": "refLine1",
                                "type": "vline",
                            },
                            "data": [
                                {"x": np.int64(1), "y": np.float64(0.5)},
                                {"x": np.int64(1), "y": np.float64(1.0)},
                            ],
                        },
                        {
                            "meta": {
                                "colour": "#008000ff",
                                "dashes": None,
                                "width": 1.0,
                                "xKey": "x",
                                "yKey": "y",
                                "label": "refLine2",
                                "type": "ref",
                            },
                            "data": [
                                {"x": np.int64(1), "y": np.int64(1)},
                                {"x": np.int64(2), "y": np.int64(2)},
                            ],
                        },
                    ],
                    "texts": [
                        {
                            "meta": {
                                "label": [("text", "Test text 1")],
                                "colour": "#ff0000ff",
                            },
                            "data": {"x": 0.25, "y": 0.25},
                        },
                        {
                            "meta": {
                                "label": [("text", "Test text 2")],
                                "colour": "#0000ffff",
                            },
                            "data": {"x": 0.5, "y": 0.5},
                        },
                        {
                            "meta": {
                                "label": [("text", "Test text 3")],
                                "colour": "#008000ff",
                            },
                            "data": {"x": 0.75, "y": 0.75},
                        },
                    ],
                }
            ]
        },
        cls=plot_to_json.NumpyEncoder,
    )
    assert returned == expected
