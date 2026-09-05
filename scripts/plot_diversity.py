"""
This script creates paper Figure 1, showing the spread of diversity of the groups.

Author: Jemma M. Fendley
"""

import matplotlib.pyplot as plt, seaborn as sns
import pandas as pd, argparse, numpy as np
import parameters  # preset Matplotlib formatting
import matplotlib as mpl
from cycler import cycler

# modify the order of the color palette for this figure only
color_palette = [
    "#BBBBBB",
    "#332288",
    "#CC6677",
    "#DDCC77",
    "#117733",
    "#88CCEE",
    "#882255",
    "#44AA99",
    "#999933",
    "#AA4499",
]
mpl.rcParams["axes.prop_cycle"] = cycler(color=color_palette)


def main():
    parser = argparse.ArgumentParser(
        description="creates a figure showing group diversity"
    )
    parser.add_argument(
        "-o", "--output", help="figure png file", type=str, required=True
    )
    parser.add_argument(
        "-p", "--output_pdf", help="figure pdf file", type=str, required=True
    )
    parser.add_argument(
        "-i",
        "--input",
        help="TSV file with mean pairwise metrics for all groups",
        type=str,
        required=True,
    )
    args = parser.parse_args()

    # load the TSV file
    df = pd.read_csv(args.input, sep="\t")
    n_groups = len(df)
    assert n_groups == 88, "There are not 88 groups, modify hard-coded label."

    # initialize the figure
    fig, axes = plt.subplots(
        nrows=1, ncols=1, figsize=(6.5, 4.65), layout="constrained"
    )

    # extract the statistics of interest
    df["mean pairwise ANI in core genome (\%)"] = 100 - 100 * df["mean_hamming"]
    df["mean pairwise coverage (\%)"] = df["mean_percent_pairwise_coverage"]

    # separate out groups mentioned in the paper (hard-coded)
    select_groups_list = ["88 total groups", "A11", "E", "EE", "DE1", "F1"]

    def figure_label(name):
        # returns the label for the group in the figure
        if name in select_groups_list:
            return name
        else:
            return "88 total groups"

    # add the figure labels and sort
    df["figure_label"] = df["name"].apply(figure_label)
    df["figure_label"] = pd.Categorical(
        df["figure_label"], categories=select_groups_list, ordered=True
    )
    df = df.sort_values(by="figure_label")

    # plot the data
    sns.scatterplot(
        data=df,
        x="mean pairwise ANI in core genome (\%)",
        y="mean pairwise coverage (\%)",
        alpha=0.8,
        size="figure_label",
        sizes=[50, 90, 90, 90, 90, 90],
        hue="figure_label",
        style="figure_label",
    )

    cutoffs, labels = [9500, 7000, 2500], ["95\% ANI", "70\% ANI", "50\% ANI"]
    rotation, offset, colors = [-26, -20, -20], [106, -160, -360], ["C6", "C7", "C8"]

    # add the species and genus thresholds
    for i, label in enumerate(labels[:2]):
        x = np.linspace(int(cutoffs[i] / 100), 100, 1000)
        y = cutoffs[i] / x
        y2 = np.ones(1000) * 100
        axes.plot(x, y, color=colors[i], alpha=0.9, linewidth=1)
        axes.text(
            x[offset[i]],
            y[offset[i]],
            "approx. " + label,
            fontsize=8,
            rotation=rotation[i],
            rotation_mode="anchor",
        )

        # add shading in the regions
        if i == 0:
            axes.fill_between(x, y, y2, color=colors[i], alpha=0.125)

    x1, x2 = np.linspace(70, 95, 100), np.linspace(95, 100, 100)
    axes.fill_between(x1, 7000 / x1, 9500 / x1, color="C7", alpha=0.05)
    axes.fill_between(x2, 7000 / x2, np.ones(100) * 85, color="C7", alpha=0.05)

    # add vOTU thresholds and shading
    x, y, y2 = np.ones(100) * 95, np.linspace(85, 100, 100), np.linspace(50, 85, 100)
    axes.plot(x, y, color="k", linewidth=1, alpha=0.9)
    axes.plot(x, y2, color="k", linewidth=1, alpha=0.5, linestyle="dashed")
    x, x2, y = np.linspace(95, 100, 100), np.linspace(50, 95, 100), np.ones(100) * 85
    axes.plot(x, y, color="k", linewidth=1, alpha=0.9)
    axes.fill_between(x, y, 9500 / x, color="k", alpha=0.05)
    axes.plot(x2, y, color="k", linewidth=1, alpha=0.5, linestyle="dashed")
    axes.text(95.35, 83.75, "vOTU diversity", fontsize=8.5)

    # formatting
    axes.legend(
        handlelength=1,
        handletextpad=0.5,
        framealpha=0.98,
        loc="upper left",
        ncol=1,
        columnspacing=0.8,
        borderaxespad=0.3,
    )
    axes.set_xlim([69.5, 100])
    axes.set_ylim([57, 100])

    # save the figure
    fig.savefig(args.output, dpi=450)
    fig.savefig(args.output_pdf, dpi=450)


if __name__ == "__main__":
    main()
