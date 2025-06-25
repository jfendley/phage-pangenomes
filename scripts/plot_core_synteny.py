"""
This script creates a figure in the paper that illustrates core synteny.

Author: Jemma M Fendley
"""

import json, argparse
import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.transforms import ScaledTranslation
import matplotlib.pyplot as plt
import parameters  # preset figure formatting


def main():
    parser = argparse.ArgumentParser(description="checks each group for core synteny")
    parser.add_argument("-f", "--figure", help="paper figure", type=str, required=True)
    parser.add_argument(
        "-s", "--figure_svg", help="svg format", type=str, required=True
    )
    parser.add_argument(
        "-i",
        "--input",
        help="the JSON file with synteny information",
        required=True,
        type=str,
    )
    args = parser.parse_args()

    # Load the synteny information
    with open(args.input) as f:
        dict = json.load(f)

    # Extract the data needed for the figure
    count_syntenic = dict["n_syntenic"]
    nonsyntenic_info = [
        (x["percent_syntenic"], x["max_cayley"])
        for x in dict["nonsyntenic_groups"]
        if x["synteny"] == "nonsyntenic"
    ]
    nonsyntenic_percent = [x[0] for x in nonsyntenic_info]
    nonsyntenic_cayley = [x[1] for x in nonsyntenic_info]
    cyclic_percents = [
        x["percent_syntenic"]
        for x in dict["nonsyntenic_groups"]
        if x["synteny"] == "cyclic"
    ]
    unique_cyclic_percents, counts_cyclic_percents = np.unique(
        cyclic_percents, return_counts=True
    )
    cyclic_labels = [
        "cyclic permutations ({0:0.0f} groups)".format(len(cyclic_percents)),
        "_nolegend_",
        "_nolegend_",
    ]

    # Initialize the figure the figure
    fig = plt.figure(
        layout="constrained",
        figsize=(6.2, 3.4),
        dpi=450,
    )
    gs = GridSpec(
        1,
        2,
        width_ratios=[1.3, 2],
        figure=fig,
    )
    ax1 = fig.add_subplot(gs[0, 1])
    ax = fig.add_subplot(gs[0, 0])
    ax.set_visible(False)  # The first panel will have a figure drawn in Inkscape.

    # Add transparent horiztonal lines to help interpret figure
    for x in range(17):
        ax1.axhline(x, color="k", alpha=0.2, linewidth=0.75)

    # Plot the nonsyntenic data
    ax1.scatter(
        nonsyntenic_percent,
        nonsyntenic_cayley,
        s=5,
        rasterized=True,
        label="partially syntenic ({0:0.0f} groups)".format(len(nonsyntenic_percent)),
    )

    # Plot the cyclic data
    for index, cyclic_percent in enumerate(unique_cyclic_percents):
        ax1.scatter(
            cyclic_percent,
            0,
            s=5 * counts_cyclic_percents[index],
            color="C1",
            rasterized=True,
            label=cyclic_labels[index],
        )

    # Plot the syntenic data
    ax1.scatter(
        100,
        0,
        color="C2",
        s=5 * count_syntenic,
        alpha=1,
        rasterized=True,
        label="completely syntenic ({0:0.0f} groups)".format(count_syntenic),
    )

    # Figure formatting (hard-coded)
    ax1.set_ylim([-1.4, 17])
    ax1.set_xlim([65, 102.5])
    ax1.set_yticks([0, 2, 4, 6, 8, 10, 12, 14, 16], [0, 2, 4, 6, 8, 10, 12, 14, 16])

    lgnd = ax1.legend()
    for handle in lgnd.legend_handles:
        handle.set_sizes([6.0])
    ax1.set_ylabel("max n. transpositions to consensus")
    ax1.set_xlabel("\% genomes following consensus core pham order")

    # add panel labels
    ax1.text(
        0.0,
        1.0,
        "b)",
        transform=(
            ax1.transAxes + ScaledTranslation(-31 / 72, -5 / 72, fig.dpi_scale_trans)
        ),
        va="bottom",
    )
    fig.text(0.0, 0.987, "a)", va="top")

    # save the figure
    fig.savefig(args.figure, dpi=450)
    fig.savefig(args.figure_svg, dpi=450)


if __name__ == "__main__":
    main()
