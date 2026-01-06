#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import annotations

__author__ = "Chao Yang"
__version__ = "1.0"

"""
A series of functions to plot scientific data.
"""

from pathlib import Path

import matplotlib as mpl

ACS_STYLE = {
    # Font
    "font.family": "sans-serif",
    "font.sans-serif": [
        "Arial",
        "Helvetica",
        "DejaVu Sans",
        "Bitstream Vera Sans",
        "sans-serif",
    ],
    "font.size": 8,
    # Axes
    "axes.labelsize": 8,
    "axes.titlesize": 9,
    "axes.linewidth": 0.8,
    "axes.edgecolor": "black",
    "axes.spines.top": True,
    "axes.spines.right": True,
    "axes.grid": False,
    # Lines and markers
    "lines.linewidth": 1.2,
    "lines.markersize": 4,
    # Ticks
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    # Legend
    "legend.fontsize": 7,
    "legend.frameon": False,
    # Figure
    "figure.dpi": 300,
    "savefig.dpi": 600,
    "figure.figsize": (3.5, 2.5),  # single column JACS figure
    "figure.facecolor": "white",
    # Savefig
    "savefig.transparent": False,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
}


NATURE_STYLE = {
    # Font
    "font.family": "sans-serif",
    "font.sans-serif": [
        "Arial",
        "Helvetica",
        "DejaVu Sans",
        "Bitstream Vera Sans",
        "sans-serif",
    ],
    "font.size": 8,
    # Axes
    "axes.labelsize": 8,
    "axes.titlesize": 9,
    "axes.linewidth": 0.8,
    "axes.edgecolor": "black",
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.grid": False,
    # Lines and markers
    "lines.linewidth": 1.2,
    "lines.markersize": 4,
    # Ticks
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    # Grid
    "grid.linewidth": 0.5,
    "grid.linestyle": "--",
    # Legend
    "legend.fontsize": 7,
    "legend.frameon": False,
    # Figure
    "figure.dpi": 300,
    "savefig.dpi": 600,
    "figure.figsize": (3.5, 2.5),  # roughly Nature single-column size in inches
    "figure.facecolor": "white",
    # Savefig
    "savefig.transparent": False,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
}


def write_style(style_name, style_dict, file_path):
    """
    Write a style dictionary to a .mplstyle file.

    Parameters:
    - style_name: Name of the style.
    - style_dict: Dictionary containing style parameters.
    - file_path: Path to save the .mplstyle file.
    """
    with open(file_path, "w") as f:
        f.write(f"# {style_name} style\n")
        for key, value in style_dict.items():
            f.write(f"{key}: {value}\n")


def acs_style_plot(style="acs"):
    if style in mpl.style.available:
        mpl.style.use(style)
    else:
        config_dir = mpl.get_configdir()
        if not (Path(config_dir) / "stylelib").exists():
            (Path(config_dir) / "stylelib").mkdir(parents=True, exist_ok=True)
        write_style(
            style, ACS_STYLE, Path(config_dir) / "stylelib" / f"{style}.mplstyle"
        )
