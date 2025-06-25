"""
This script plots the distribution of regions of high linkage across the core genome.

Author: Jemma M. Fendley
"""

import pandas as pd, numpy as np
import argparse, matplotlib.pyplot as plt
from Bio import AlignIO
from datetime import datetime
import matplotlib as mpl
import parameters  # preset matplotlib formatting

# find_snp_positions finds the positions of SNPs given a core genome alignment
from utils import find_snp_positions


def main():
    parser = argparse.ArgumentParser(description="plot distribution of linkage")
    parser.add_argument(
        "-o", "--output", help="output figure file", type=str, required=True
    )
    parser.add_argument(
        "-l", "--linkage", help="linkage disequilibrium file", type=str, required=True
    )
    parser.add_argument(
        "-c",
        "--core_genome",
        help="core genome alignment file",
        type=str,
        required=True,
    )

    args = parser.parse_args()

    # load the core genome alignment and extract the snp poistions
    core_genome_alignment = np.array(AlignIO.read(args.core_genome, "fasta"))
    n_positions = core_genome_alignment.shape[1]
    snp_positions = find_snp_positions(core_genome_alignment)[0]

    # load the linkage disequilibrium file (slow) and drop columns to reduce memory
    df = pd.read_feather(args.linkage)
    df = df.drop(columns=["p1", "p2", "bg_ld"])

    # choose the size of the bins
    binsize = 100

    def get_bin(position):
        """
        Returns a bin identifier given a position
        """
        return position // binsize

    # add bin columns to the dataframe
    df["bin1"] = df["position1"].apply(get_bin)
    df["bin2"] = df["position2"].apply(get_bin)

    # initalize matrix
    max_bin = (n_positions - 1) // binsize
    enrichment_matrix = np.full((max_bin + 1, max_bin + 1), np.nan)

    # choose a threshold for what is considered "high" linkage
    threshold = np.percentile(df["ld"], 90)
    data = df.groupby(["bin1", "bin2"])["ld"].count()
    data_high = df[df["ld"] >= threshold].groupby(["bin1", "bin2"])["ld"].count()

    # each bin will have fraction of that bin with high linkage
    for i in data_high.index:
        if data[i] >= 10:
            enrichment = data_high[i] / data[i]
            enrichment_matrix[min(i[0], i[1]), max(i[0], i[1])] = enrichment
    for i in [x for x in data.index if x not in data_high.index and data[x] >= 10]:
        enrichment_matrix[min(i[0], i[1]), max(i[0], i[1])] = 0

    # initialize figure
    fig, axes = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=(4, 4.24),
        gridspec_kw={"height_ratios": [5, 1]},
    )

    # plot the linkage enrichment matrix
    mat = axes[0].matshow(enrichment_matrix, vmin=0, cmap="YlGnBu")
    cb = fig.colorbar(mat, ax=axes, shrink=0.4)
    cb.ax.get_yaxis().labelpad = 15
    cb.ax.set_ylabel("fraction high LD", rotation=270)
    tick_range = range(0, n_positions // binsize, binsize)
    axes[0].set_xticks(
        ticks=tick_range,
        labels=binsize * np.array(tick_range, dtype=np.int64),
    )
    axes[0].set_yticks(
        ticks=tick_range,
        labels=binsize * np.array(tick_range, dtype=np.int64),
        rotation=90,
        ha="center",
        rotation_mode="anchor",
        va="baseline",
    )

    axes[0].set_title(
        args.linkage.split("/")[-1].split("_")[0]
        + r" enrichment: (n. ld $> {0:0.2f})$/ (n. ld)".format(threshold)
    )
    axes[0].set_xlabel("binsize = {0:0.0f}".format(binsize))
    axes[0].set_ylabel("core genome position")

    # plot the snp distribution across the core genome
    axes[1].hist(
        snp_positions,
        bins=range(0, n_positions + binsize, binsize),
        density=False,
        cumulative=False,
    )
    axes[1].set_xlabel("core genome position")
    axes[1].set_ylabel("snp count")
    axes[1].set_xlim([0, n_positions - 1])
    axes[1].set_ylim([0, binsize])
    axes[1].set_xticks(
        ticks=binsize * np.array(tick_range, dtype=np.int64),
        labels=binsize * np.array(tick_range, dtype=np.int64),
    )
    axes[1].set_yticks(
        ticks=[0, 40, 80],
        labels=[0, 40, 80],
        rotation=90,
        ha="center",
        rotation_mode="anchor",
        va="baseline",
    )

    # save the figure
    fig.savefig(args.output, dpi=450)


if __name__ == "__main__":
    main()
