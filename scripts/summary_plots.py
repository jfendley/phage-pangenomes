"""
This script plots Figure S1 in the paper, a summary of properties
    of the groups in the databse.

Author: Jemma M. Fendley
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.transforms import ScaledTranslation
import parameters  # preset matplotlib formatting


def main():
    parser = argparse.ArgumentParser(description="create PA basic LD plot")
    parser.add_argument(
        "-o", "--output", help="png output figure file", type=str, required=True
    )
    parser.add_argument(
        "-p", "--output_pdf", help="pdf output figure file", type=str, required=True
    )
    parser.add_argument(
        "-i",
        "--input",
        help="dataframe with group information",
        type=str,
        required=True,
    )

    args = parser.parse_args()

    # load the dataframe
    df = pd.read_csv(args.input, sep="\t")

    # initialize the figure
    fig, axes = plt.subplots(nrows=4, ncols=2, figsize=(6, 8.5), layout="constrained")

    # plot the histogram of number of phages per group
    sns.histplot(data=df, x="n_phages", ax=axes[0, 0])
    axes[0, 0].set_xlabel("n. phages")

    # plot the histogram of mean genome length of each group
    sns.histplot(data=df, x="mean_genome_length", ax=axes[0, 1])
    axes[0, 1].set_xlabel("mean genome length")

    # plot the histogram of Hamming distances
    df["mean p.w. core Hamming distance (\%)"] = df["mean_hamming"] * 100
    sns.histplot(data=df, x="mean p.w. core Hamming distance (\%)", ax=axes[1, 0])
    axes[1, 0].set_yticks(np.arange(0, 25, 5), np.arange(0, 25, 5))

    # plot the histogram of mean percent core
    sns.histplot(data=df, x="mean_percent_core", ax=axes[1, 1])
    axes[1, 1].set_xlabel("mean percent core (of genome length)")

    # plot the life cycle
    sns.countplot(
        data=df,
        x="life_cycle",
        ax=axes[2, 0],
        order=df["life_cycle"].value_counts().index,
    )
    axes[2, 0].set_xlabel("life cycle")

    # plot the most common morphotype
    sns.countplot(
        data=df,
        x="most_common_morphotype",
        ax=axes[2, 1],
        order=df["most_common_morphotype"].value_counts().index,
    )
    axes[2, 1].set_xlabel("most common morphotype")

    # plot the most common host genus
    sns.countplot(
        data=df,
        x="most_common_host_genus",
        ax=axes[3, 0],
        order=df["most_common_host_genus"].value_counts().index,
    )
    axes[3, 0].set_xlabel("isolation host genus")
    x_axis = axes[3, 0].get_xticks()
    x_labels = axes[3, 0].get_xticklabels()
    axes[3, 0].set_xticks(x_axis, x_labels, rotation=90)

    # plot the most common host species
    sns.countplot(
        data=df,
        x="most_common_host_species",
        ax=axes[3, 1],
        order=df["most_common_host_species"].value_counts().index,
        hue="most_common_host_genus",
    )
    x_axis = axes[3, 1].get_xticks()
    x_labels = axes[3, 1].get_xticklabels()
    axes[3, 1].set_xticks(x_axis, x_labels, rotation=90)
    axes[3, 1].set_xlabel("most common isolation host species")
    axes[3, 1].legend(
        ncol=2, columnspacing=0.2, borderaxespad=0.3, fontsize=8, handlelength=0.5
    )
    axes[3, 1].set_ylim([0, 32.5])  # removes tick label

    # tidy y axis labels
    for i in range(4):
        axes[i, 0].set_ylabel("n. groups")
        y_axis = axes[i, 1].axes.get_yaxis()
        y_label = y_axis.get_label()
        y_label.set_visible(False)

    # add panel labels
    list_of_labels = ["a)", "b)", "c)", "d)", "e)", "f)", "g)", "h)"]
    count = 0
    for i in range(4):
        for j in range(2):
            axes[i, j].text(
                0.0,
                1.0,
                list_of_labels[count],
                transform=(
                    axes[i, j].transAxes
                    + ScaledTranslation(-18 / 72, -2 / 72, fig.dpi_scale_trans)
                ),
                va="bottom",
                # fontfamily="serif",
            )
            count += 1

    # save the figure
    fig.savefig(args.output, dpi=450)
    fig.savefig(args.output_pdf, dpi=450)


if __name__ == "__main__":
    main()
