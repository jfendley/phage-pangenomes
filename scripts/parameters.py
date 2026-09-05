"""
This script provides matplotlib formatting for consistency across all plots.

Author: Jemma M. Fendley
"""

import matplotlib as mpl, shutil
from cycler import cycler

color_palette = [
    "#332288",
    "#CC6677",
    "#DDCC77",
    "#117733",
    "#88CCEE",
    "#882255",
    "#44AA99",
    "#999933",
    "#AA4499",
    "#BBBBBB",
]
mpl.rcParams["axes.prop_cycle"] = cycler(color=color_palette)
mpl.rcParams["font.size"] = 14
mpl.rcParams["axes.titlesize"] = 14
mpl.rcParams["axes.labelsize"] = 12
mpl.rcParams["legend.fontsize"] = 10
mpl.rcParams["legend.handletextpad"] = 0.1
mpl.rcParams["legend.handlelength"] = 1
mpl.rcParams["legend.title_fontsize"] = 10
mpl.rcParams["xtick.labelsize"] = 10
mpl.rcParams["ytick.labelsize"] = 10
mpl.rcParams["font.family"] = "serif"

# LaTeX and TeX Live need to be installed for the desired figure formatting.
mpl.rcParams["text.usetex"] = True if shutil.which("latex") else False
mpl.rcParams["font.serif"] = (
    "Computer Modern" if shutil.which("latex") else "Liberation Serif"
)

# For consistent math text
mpl.rcParams["mathtext.fontset"] = "cm" if shutil.which("latex") else "Liberation Serif"
mpl.rcParams["text.latex.preamble"] = r"\usepackage{amsmath} \usepackage{amsfonts}"
