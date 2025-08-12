"""
This script plots a figure showing summary statistics of the snp compatibility analysis
    for all of the groups.

Author: Jemma M. Fendley
"""

import numpy as np, pandas as pd
import matplotlib.pyplot as plt, seaborn as sns
import matplotlib as mpl
import argparse, json
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
    print("Mean: ", np.mean(df["data_biallelic_mean_interval_size"]))

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
    axes[0].plot(range(2, 100), range(2, 100), color="k", label=r"$y=x$")
    axes[1].plot(range(9, 600), range(9, 600), color="k", label=r"$y=x$")

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
    axes[0].legend()
    axes[1].legend().set_visible(False)
    axes[0].set_ylabel(r"$\mathbb{E}\left[\text{n. snps between recurrent mutations}\right]$")
    axes[1].set_ylabel(r"$\mathbb{E}\left[\text{distance between recurrent mutations}\right]$")
    axes[0].set_xlabel("mean n. compatible snps")
    axes[1].set_xlabel("mean compatible interval size")

    # save figure
    fig.savefig(args.output, dpi=450)
    fig.savefig(args.output_pdf, dpi=450)


if __name__ == "__main__":
    main()
