"""
This script plots the distribution of SNPs across the core genome alignment.

Author: Jemma M. Fendley
"""

import matplotlib.pyplot as plt, numpy as np
import argparse, matplotlib as mpl
from Bio import AlignIO
import parameters  # preset matplotlib formatting
from matplotlib.gridspec import GridSpec
from matplotlib.transforms import ScaledTranslation

# find_snp_positions finds the positions of SNPs given a core genome alignment
from utils import find_snp_positions

mpl.rcParams["font.size"] = 12


def main():
    parser = argparse.ArgumentParser(
        description="plot distribution of SNPs across core genome"
    )
    parser.add_argument("-o", "--output", help="figure file", type=str, required=True)
    parser.add_argument(
        "-i", "--core_genome", help="core genome alignment", type=str, required=True
    )
    args = parser.parse_args()

    # find the snp positions, and the n. of each alllele at each position
    dna_alphabet = ["G", "T", "A", "C"]  # might need to be modified for lowercase
    core_genome_alignment = np.array(AlignIO.read(args.core_genome, "fasta"))
    n_phages, n_positions = core_genome_alignment.shape
    allele_counts = np.array(
        [np.sum(core_genome_alignment == nuc, axis=0) for nuc in dna_alphabet]
    )
    snp_positions, snp_positions_biallelic = find_snp_positions(core_genome_alignment)

    # extract the frequency of the majority allele at each snp position
    max_frequency = np.max(allele_counts, axis=0) / n_phages
    snp_frequencies = max_frequency[snp_positions]
    snp_frequencies_biallelic = max_frequency[snp_positions_biallelic]

    # initialize figure
    fig = plt.figure(layout="constrained", figsize=(7, 6), dpi=450)
    gs = GridSpec(
        2,
        2,
        width_ratios=[1, 1],
        height_ratios=[1, 1],
        figure=fig,
    )
    ax1, ax2 = fig.add_subplot(gs[0, :]), fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[1, 1])
    axes = [ax1, ax2, ax3]

    # plot the distribution across the core genome of snp locations
    ax1.hist(
        snp_positions,
        bins=range(0, n_positions + 100, 100),
        density=False,
        cumulative=False,
    )

    # plot the distribution of the frequency of the majority allele
    binsize = 1 / n_phages
    ax2.hist(
        snp_frequencies,
        bins=np.arange(
            np.min(snp_frequencies) - binsize / 2,
            np.max(snp_frequencies) + binsize,
            binsize,
        ),
    )
    ax2.set_xlim([np.min(snp_frequencies) - binsize, 1])
    ax3.hist(
        snp_frequencies_biallelic,
        bins=np.arange(
            np.min(snp_frequencies_biallelic) - binsize / 2,
            np.max(snp_frequencies_biallelic) + binsize,
            binsize,
        ),
    )
    ax3.set_xlim([np.min(snp_frequencies_biallelic) - binsize, 1])

    # format and save figure
    ax1.set_xlabel("core genome position")
    ax1.set_ylabel("SNP count")
    ax1.set_xlim([0, n_positions])
    ax1.set_ylim([0, 100])
    ax1.set_yticks(
        ticks=[0, 40, 80],
        labels=[0, 40, 80],
        rotation=90,
        ha="center",
        rotation_mode="anchor",
        va="baseline",
    )
    ax2.set_xlabel("snp frequencies (frac. majority allele)")
    ax3.set_xlabel("biallelic snp frequencies")
    for i in range(1, 3):
        axes[i].set_ylabel("n. snps")
        axes[i].set_yscale("log")
        axes[i]
    list_of_labels = ["a)", "b)", "c)"]
    for i in range(3):
        axes[i].text(
            0.0,
            1.0,
            list_of_labels[i],
            transform=(
                axes[i].transAxes
                + ScaledTranslation(-38 / 72, -5 / 72, fig.dpi_scale_trans)
            ),
            va="bottom",
        )
    fig.suptitle(
        "Group "
        + args.core_genome.split("/")[-1].split("_")[0]
        + ": {0:0.0f} phages".format(n_phages)
    )
    fig.savefig(args.output)


if __name__ == "__main__":
    main()
