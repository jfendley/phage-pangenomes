"""
This script looks plots an pairwise snp statistics for one pair of phages.

Author: Jemma M. Fendley
"""

from Bio import AlignIO
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import json, argparse, matplotlib as mpl
import parameters  # preset matplotlib formatting
from utils import calculate_hamming  # calculates hamming distance
from utils import convert_to_list  # converts string list of ints to list of ints
from matplotlib.transforms import ScaledTranslation


mpl.rcParams["legend.handletextpad"] = 0.4


def main():
    parser = argparse.ArgumentParser(description="pairwise SNP analysis for one group")
    parser.add_argument(
        "-o", "--output", help="figure png file", type=str, required=True
    )
    parser.add_argument(
        "-p", "--output_pdf", help="figure pdf file", type=str, required=True
    )
    parser.add_argument(
        "-s",
        "--core_pham_starts",
        help="name of core pham starts file",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-c",
        "--core_genome_file",
        help="name of core genome file",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-t",
        "--snp_analysis",
        help="pairwise snp analysis table",
        type=str,
        required=True,
    )
    args = parser.parse_args()

    # load the pairwise snp analysis dataframe
    df = pd.read_csv(args.snp_analysis, sep="\t")
    df["distance_list"] = df["distances"].apply(convert_to_list)

    # extract the phages and estimated maximum distance from the first row
    row = df.iloc[0]
    chosen_phages, group = [row["phage_1"], row["phage_2"]], row["group"]
    estimated_recombination_length = np.max(row["distance_list"])
    print("Estimated recombination length: ", estimated_recombination_length, "bp")

    # load the core phams starts dictionary
    with open(args.core_pham_starts) as f:
        core_starts_dict = json.load(f)

    # load the core genome alignment and extract the phage indices
    alignment = AlignIO.read(args.core_genome_file, "fasta")
    phage_IDs = [record.id for record in alignment]
    ind1, ind2 = [phage_IDs.index(x) for x in chosen_phages]
    core_genome_alignment = np.array(alignment)
    n_positions = core_genome_alignment.shape[1]

    binsize = 100

    # Calculate Hamming distance
    hamming = calculate_hamming(core_genome_alignment, ind1, ind2)
    assert hamming < 0.02 and hamming > 0.001, "Need to chose different pair of phages"

    # find the snp positions and initiliaze the figure
    snp_positions = np.where(
        core_genome_alignment[ind1] != core_genome_alignment[ind2]
    )[0]
    fig, axes = plt.subplots(
        nrows=1,
        ncols=2,
        figsize=(6, 2.8),
        gridspec_kw={"width_ratios": [1, 1]},
        layout="constrained",
    )

    # plot the distribution and cumulative distribution of snps across the core genome
    axes[0].hist(
        snp_positions,
        bins=np.arange(0, n_positions + binsize, binsize),
        label="n. SNPs per bin",
    )
    n_snps_cumulative, bins_cumulative, patches = axes[1].hist(
        snp_positions,
        bins=np.arange(0, n_positions + binsize, binsize),
        cumulative=True,
        histtype="step",
        label="cumulative n. SNPs",
    )

    # plot the diagonal for reference
    axes[1].plot((0, n_positions - 1), (0, n_snps_cumulative[-1]), label="diagonal")

    # also plot the core pham boundaries for reference
    for ind in range(2):
        axes[ind].axvline(
            list(core_starts_dict.values())[0],
            color="k",
            linewidth=0.25,
            alpha=0.6,
            label="core pham boundaries",
        )
        for x in list(core_starts_dict.values())[1:]:
            axes[ind].axvline(
                x, color="k", linewidth=0.25, alpha=0.6, label="_nolegend_"
            )

        # figure formatting
        axes[ind].set_xlabel("core genome position")
        axes[ind].set_xlim([-1, n_positions + 1])
        axes[ind].legend()
    axes[0].set_ylabel("n. SNPs (binsize = {0:0.0f})".format(binsize))
    axes[1].set_ylabel("cumulative n. SNPs")
    fig.suptitle(
        group
        + ": "
        + chosen_phages[0]
        + " and "
        + chosen_phages[1]
        + " (distance: {0:0.2f}\%)".format(100 * hamming)
    )
    label_list, x_loc = ["a)", "b)"], [-35, -37]
    for i in range(2):
        axes[i].text(
            0.0,
            1.0,
            label_list[i],
            transform=(
                axes[i].transAxes
                + ScaledTranslation(x_loc[i] / 72, +8 / 72, fig.dpi_scale_trans)
            ),
            va="bottom",
        )

    # save the figure
    fig.savefig(args.output, dpi=450)
    fig.savefig(args.output_pdf, dpi=450)


if __name__ == "__main__":
    main()
