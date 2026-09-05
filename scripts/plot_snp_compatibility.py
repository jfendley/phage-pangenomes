"""
This script plots a figure showing summary statistics of the snp compatibility analysis
    for all of the groups.

Author: Jemma M. Fendley
"""

import numpy as np, pandas as pd
import matplotlib.pyplot as plt, seaborn as sns
import matplotlib as mpl, argparse
import parameters  # preset matplotlib formatting
from matplotlib.transforms import ScaledTranslation

mpl.rcParams["legend.handletextpad"] = 0.4


def main():
    parser = argparse.ArgumentParser(description="plot SNP compatibility statistics")
    parser.add_argument(
        "-o", "--output", help="png figure file", type=str, required=True
    )
    parser.add_argument(
        "-p", "--output_pdf", help="pdf figure file", type=str, required=True
    )
    parser.add_argument(
        "-i", "--input", help="snp compatibility TSV", type=str, required=True
    )
    args = parser.parse_args()

    # initalize figure and load dataframe
    fig, axes = plt.subplots(
        nrows=1,
        ncols=2,
        figsize=(6, 3),
        gridspec_kw={"width_ratios": [1, 1]},
        layout="constrained",
    )
    df = pd.read_csv(args.input, sep="\t")

    # print some summary statistics
    print("Mean interval length: ", np.mean(df["data_biallelic_mean_interval_size"]))
    df["snp_enrichment"] = (
        df["expected_n_snps_biallelic"] / df["data_biallelic_mean_n_snps"]
    )
    mean_enrichment_snps = df["snp_enrichment"].mean()
    df["length_enrichment"] = (
        df["expected_distance_biallelic"] / df["data_biallelic_mean_interval_size"]
    )
    mean_enrichment_length = df["length_enrichment"].mean()
    print("Mean ratio of expected length to length: ", mean_enrichment_length)
    print("Mean ratio of expected n. snps to n. snps: ", mean_enrichment_snps)
    print(
        "N. groups with smaller expected length: ",
        len(df[df["length_enrichment"] <= 1]),
    )
    print(
        "N. groups with smaller expected n. snps interval size: ",
        len(df[df["snp_enrichment"] <= 1]),
    )

    # plot the data
    sns.scatterplot(
        data=df,
        x="data_biallelic_mean_n_snps",
        y="expected_n_snps_biallelic",
        ax=axes[0],
        alpha=0.8,
        s=8,
        label="data",
    )
    sns.scatterplot(
        data=df,
        x="data_biallelic_mean_interval_size",
        y="expected_distance_biallelic",
        ax=axes[1],
        alpha=0.8,
        s=8,
        label="data",
    )
    # plot line y=x for comparison
    axes[0].set_xlim([2, 100])
    axes[1].set_xlim([9, 600])
    axes[0].plot(
        range(2, 100), range(2, 100), color="k", alpha=0.5, label=r"$y=x$", linewidth=1
    )
    axes[1].plot(
        range(9, 600), range(9, 600), color="k", alpha=0.5, label=r"$y=x$", linewidth=1
    )

    # plot also y=2x, y=5x, and y=15x for comparison
    colors = ["C1", "C2", "C3", "C4"]
    linestyle_list = ["--", "-.", ":"]
    for j, i in enumerate([2, 5, 15]):
        axes[0].plot(
            np.arange(2, 100),
            i * np.arange(2, 100),
            color=colors[j],
            linewidth=1,
            linestyle=linestyle_list[j],
            alpha=0.5,
            label=r"$y={0:0.0f}x$".format(i),
        )
        axes[1].plot(
            np.arange(9, 600),
            i * np.arange(9, 600),
            color=colors[j],
            linewidth=1,
            linestyle=linestyle_list[j],
            alpha=0.5,
            label=r"$y={0:0.0f}x$".format(i),
        )

    # figure formatting
    letters = ["a)", "b)"]
    for i in range(2):
        axes[i].set_yscale("log")
        axes[i].set_xscale("log")
        axes[i].text(
            0.0,
            1.0,
            letters[i],
            transform=(
                axes[i].transAxes
                + ScaledTranslation(-20 / 72, +7 / 72, fig.dpi_scale_trans)
            ),
            va="bottom",
            # fontfamily="serif",
        )
    axes[0].legend(loc="upper left", labelspacing=0)
    axes[1].legend().set_visible(False)
    axes[0].set_ylabel(
        r"$\mathbb{E}\left[\text{n. snps between recurrent mutations}\right]$"
    )
    axes[1].set_ylabel(
        r"$\mathbb{E}\left[\text{distance between recurrent mutations}\right]$"
    )
    axes[0].set_xlabel("mean n. compatible snps")
    axes[1].set_xlabel("mean compatible interval size")

    # save figure
    fig.savefig(args.output, dpi=450)
    fig.savefig(args.output_pdf, dpi=450)


if __name__ == "__main__":
    main()
