"""
This script provides matplotlib formatting for consistency across all plots.

Author: Jemma M. Fendley
"""

import matplotlib as mpl
from cycler import cycler
import shutil

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
    "#DDDDDD",
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
#   If one or both of these is not installed, please comment out the following three lines,
#   and uncomment the next two.
mpl.rcParams["text.usetex"] = True if shutil.which("latex") else False
mpl.rcParams["font.serif"] = "cm" if shutil.which("latex") else "Liberation Serif"
mpl.rcParams["text.latex.preamble"] = r"\usepackage{amsmath} \usepackage{amsfonts}"
# mpl.rcParams["text.usetex"] = False
# mpl.rcParams["font.serif"] = "Liberation Serif"
