"""
This script creates a plot for the SI that shows the extent of synteny in the accesory phams

Author: Jemma M. Fendley
"""

import numpy as np, pandas as pd
import argparse
import matplotlib.pyplot as plt, seaborn as sns
import parameters  # preset matplotlib formatting
from matplotlib.transforms import ScaledTranslation


def main():
    parser = argparse.ArgumentParser(description="plot accessory synteny figures")
    parser.add_argument(
        "-o", "--output", help="png file for plot", type=str, required=True
    )
    parser.add_argument(
        "-p", "--output_pdf", help="pdf file for plot", type=str, required=True
    )
    parser.add_argument(
        "-i", "--input", help="file for synteny dataframe", type=str, required=True
    )
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")

    # print key statistics
    print(
        "N. groups with all single-junction phams: ",
        len(df[df["fraction_single_junction"] == 1]),
    )
    n_syntenic = df["n_nonsyntenic_junctions"].value_counts()[0]
    print("N. groups with all syntenic junctions: ", n_syntenic)

    # convert to percent and change to plot-friendly names
    df["percent of non-empty junctions that are not syntenic"] = (
        100 * df["n_nonsyntenic_junctions"] / df["n_nonempty_junctions"]
    )
    df["percent of non-syntenic pairs in non-syntenic junctions"] = (
        100 * df["n_nonsyntenic_pairs"] / df["n_pham_pairs"]
    )
    df["\% single-junction accessory phams"] = df["fraction_single_junction"] * 100

    # calculate maximum values to create appropriate histogram bins
    max_junctions = np.max(df["percent of non-empty junctions that are not syntenic"])
    max_pairs = np.max(df["percent of non-syntenic pairs in non-syntenic junctions"])

    # count the number in each histogram bin for colorbar labeling
    counts, xedges, yedges = np.histogram2d(
        df["percent of non-empty junctions that are not syntenic"],
        df["percent of non-syntenic pairs in non-syntenic junctions"],
        bins=[
            np.arange(-2.5, max_junctions + 5, 5),
            np.arange(-2.5, max_pairs + 5, 5),
        ],
    )
    # there will be a lot of groups in the histogram bin (0,0), meaning syntenic, and then only a
    #   few phages in each of the nonsyntenic bins
    n_syntenic = df["n_nonsyntenic_junctions"].value_counts()[0]

    # this is the maximum number of groups per bin, excluding the large syntenic bin
    max_small = np.max([int(x) for x in np.unique(counts) if x != n_syntenic])

    # plot the histogram of pairwise synteny
    fig, axes = plt.subplots(
        1, 2, width_ratios=[1, 1], layout="constrained", figsize=(6, 3)
    )
    hist = sns.histplot(
        data=df,
        x="percent of non-empty junctions that are not syntenic",
        y="percent of non-syntenic pairs in non-syntenic junctions",
        ax=axes[1],
        bins=[
            np.arange(-2.5, max_junctions + 5, 5),
            np.arange(-2.5, max_pairs + 5, 5),
        ],
        cmap="viridis",
        norm="log",
        vmin=1,
        vmax=n_syntenic,
    )
    cbar = hist.figure.colorbar(hist.collections[0])
    cbar.set_label("n. groups", rotation=270, labelpad=11)
    tick_list = list(range(1, max_small + 1)) + [
        x * 10 for x in range(1, n_syntenic // 10 + 1)
    ]
    cbar.set_ticks(tick_list, labels=tick_list)

    # plot the distribution of % single-junction accesory phams
    assert (
        np.min(df["\% single-junction accessory phams"]) > 92.5
    ), "Error: need to modify bins for accessory plots"
    sns.histplot(
        data=df,
        x="\% single-junction accessory phams",
        ax=axes[0],
        bins=np.arange(92.5, 101, 1),
    )
    axes[0].set_ylabel("n. groups (/{0:0.0f})".format(len(df)))
    axes[1].set_xlabel("\% non-syntenic non-empty junctions")
    axes[1].set_ylabel("\% non-syntenic pairs \n in non-syntenic junctions")

    # add panel labels and save figure
    label_list = ["a)", "b)"]
    for i in range(2):
        axes[i].text(
            0.0,
            1.0,
            label_list[i],
            transform=(
                axes[i].transAxes
                + ScaledTranslation(-35 / 72, +1 / 72, fig.dpi_scale_trans)
            ),
            va="bottom",
        )
    fig.savefig(args.output, dpi=450)
    fig.savefig(args.output_pdf, dpi=450)


if __name__ == "__main__":
    main()
