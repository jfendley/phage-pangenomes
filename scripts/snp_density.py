"""
This script calculates the snp density (n. of snps in rolling windows) across the core genome, used
    to weight the linkage disequilbrium analysis. It also plots the distribution of the weights.

Author: Jemma M. Fendley
"""

import numpy as np
import argparse, json
from Bio import AlignIO
import matplotlib.pyplot as plt
import parameters  # preset matplotlib formatting

# find_snp_positions finds the positions of SNPs given a core genome alignment
from utils import find_snp_positions


def main():
    parser = argparse.ArgumentParser(description="caculate rolling SNP density")
    parser.add_argument(
        "-o", "--output", help="SNP density JSON file", type=str, required=True
    )
    parser.add_argument(
        "-f", "--output_fig", help="output figure file", type=str, required=True
    )
    parser.add_argument(
        "-i", "--core_genome", help="core genome alignment", type=str, required=True
    )
    args = parser.parse_args()

    # find the snp positions and the length of the core genome
    core_genome_alignment = np.array(AlignIO.read(args.core_genome, "fasta"))
    n_positions = core_genome_alignment.shape[1]
    snp_positions = find_snp_positions(core_genome_alignment)[0]

    # calculate weights that are invesely proportional to the density
    half_window_size, weights_dict = 199, {}
    for i in snp_positions:
        if i < half_window_size:
            j = half_window_size
        elif i > n_positions - half_window_size:
            j = n_positions - half_window_size
        else:
            j = i
        weight = (half_window_size * 2 + 1) / (
            np.sum(np.abs(snp_positions - j) <= half_window_size)
        )
        weights_dict[int(i)] = float(weight)

    # save the dictionary
    with open(args.output, "w") as f:
        json.dump(weights_dict, f, indent=4)

    # plot the distribution of the weights
    fig, axes = plt.subplots(1, 1, figsize=(3, 3), layout="constrained")
    axes.hist(list(weights_dict.values()))
    axes.set_yscale("log")
    axes.set_xlabel("weight")
    axes.set_ylabel("n. snps")
    fig.savefig(args.output_fig, dpi=450)


if __name__ == "__main__":
    main()
