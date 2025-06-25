"""
This script determines the frequency of biallelic sites in the core genome.

Author: Jemma M. Fendley
"""

import numpy as np, pandas as pd
import argparse, json
import matplotlib.pyplot as plt, seaborn as sns
import parameters  # preset matplotlib formatting
from Bio import AlignIO
from tqdm import tqdm

# find_snp_positions finds the positions of SNPs given a core genome alignment
from utils import find_snp_positions


def main():
    parser = argparse.ArgumentParser(
        description="determines the frequency of biallelic sites"
    )
    parser.add_argument(
        "-o", "--output", help="biallelic frequency figure", type=str, required=True
    )
    parser.add_argument(
        "-g",
        "--groups",
        nargs="+",
        help="list of group JSON files for all groups",
        required=True,
    )
    args = parser.parse_args()

    dict_list, dna_alphabet = [], ["G", "T", "A", "C"]  # might need to be lowercase
    for file in tqdm(args.groups):
        # load the core genome alignment and extract the number of phages
        core_genome_alignment = np.array(AlignIO.read(file, "fasta"))

        # find positions with snps
        snp_positions, snp_positions_biallelic = find_snp_positions(
            core_genome_alignment
        )

        # also investigate the effect of codon position
        snp_positions_third_codon = [x for x in snp_positions if x % 3 == 2]
        snp_positions_biallelic_third_codon = [
            x for x in snp_positions_biallelic if x % 3 == 2
        ]

        # save the information into a dataframe
        row = {
            "name": file.split("/")[-1].split("_")[0],
            "fraction_biallelic": len(snp_positions_biallelic) / len(snp_positions),
            "fraction_biallelic_third_position": len(
                snp_positions_biallelic_third_codon
            )
            / len(snp_positions_third_codon),
            "fraction_biallelic_non_third_position": (
                len(snp_positions_biallelic) - len(snp_positions_biallelic_third_codon)
            )
            / (len(snp_positions) - len(snp_positions_third_codon)),
            "n_snps": len(snp_positions),
            "n_snps_third_position": len(snp_positions_third_codon),
            "n_snps_non_third_position": len(snp_positions)
            - len(snp_positions_biallelic_third_codon),
        }
        dict_list.append(row)
    group_df = pd.DataFrame(dict_list)

    # plot the distribution across the groups
    fig, axes = plt.subplots(1, 1, layout="constrained", figsize=(3.5, 3.5))
    bins = np.arange(0.225, 1 + 0.05, 0.05)
    sns.histplot(
        data=group_df,
        x="fraction_biallelic",
        ax=axes,
        bins=bins,
        legend=False,
    )
    axes.set_xlabel("fraction of snps that are biallelic")
    axes.set_ylabel("n. groups")
    fig.savefig(args.output, dpi=450)

    # print some key statistics to quote in paper
    print(
        "Mean % biallelic per group: {0:0.1%} "
        "\nMean % biallelic across all groups: {1:0.1%} "
        "\nMean % biallelic across all groups third codon position: {2:0.1%} "
        "\nMean % biallelic across all groups non-third codon position: {3:0.2%}".format(
            np.mean(group_df["fraction_biallelic"]),
            np.average(group_df["fraction_biallelic"], weights=group_df["n_snps"]),
            np.average(
                group_df["fraction_biallelic_third_position"],
                weights=group_df["n_snps_third_position"],
            ),
            np.average(
                group_df["fraction_biallelic_non_third_position"],
                weights=group_df["n_snps_non_third_position"],
            ),
        )
    )


if __name__ == "__main__":
    main()
