"""
This script calculates the linkage disequilibrium for each pair of SNPs in the core genome.

This script is identical to linkage_disequilibrium.py except it uses find_snp_positions_processed instead of
    find_snp_positions.

Author: Jemma M. Fendley
"""

import numpy as np, pandas as pd
import argparse
from tqdm import tqdm
from Bio import AlignIO
from math import comb
from functools import cache

# find_snp_positions finds the positions of SNPs given a core genome alignment
from utils import find_snp_positions_processed


@cache
def calculate_linkage(p1, p2, p12):
    """
    Calculate the normalized linkage disequilibrium (r)

    Input:
        p1 (float):
        p2 (float):
        p12 (float):
    """
    linkage, denominator = p12 - p1 * p2, p1 * (1 - p1) * p2 * (1 - p2)
    return abs(linkage) / np.sqrt(denominator)


@cache
def background_linkage(N1, p1, N2, p2, N):
    """
    Calculate the background linkage - the expected value if the SNPs were shuffled randomly.

    Input:
        N1 (integer): the number of phages with the relevant allele at site 1
        N2 (integer): the number of phages with the relevant allele at site 2
        N (integer): the total number of phages
    """
    lower_bound = max(0, N1 + N2 - N)
    expected_linkage = [
        abs(j / N - p1 * p2) * comb(N2, j) * comb(N - N2, N1 - j) / comb(N, N1)
        for j in range(lower_bound, min(N1, N2) + 1)
    ]
    return np.sum(expected_linkage) / np.sqrt((p1 * p2 * (1 - p1) * (1 - p2)))


def main():
    parser = argparse.ArgumentParser(description="create linkage disequilibrium file")
    parser.add_argument(
        "-o", "--output", help=".feather linkage file", type=str, required=True
    )
    parser.add_argument(
        "-i",
        "--core_genome",
        help="core genome alignment file",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-p",
        "--positions",
        help="core genome positions to pham file",
        type=str,
        required=True,
    )
    args = parser.parse_args()

    dna_alphabet = ["G", "T", "A", "C", "P"]  # might need to be modified for lowercase

    # load the core genome alignment and find the snp positions
    core_genome_alignment = np.array(AlignIO.read(args.core_genome, "fasta"))
    n_phages = core_genome_alignment.shape[0]
    all_snp_positions = find_snp_positions_processed(core_genome_alignment)[0]

    # find the overlap positions and remove from analysis
    position_df = pd.read_csv(args.positions, sep="\t")
    overlap_positions = set(position_df[position_df["overlap"]]["position"])
    snp_positions = np.sort(list(set(all_snp_positions) - overlap_positions))

    # if there are a lot of SNPs, running this script would require too much memory, and
    #   then we will only consider SNP pairs within a certain distance
    large_memory = True if len(snp_positions) > 20000 else False

    # for each snp positions, find the phages with the majority allele, and those with gaps
    allele_counts = np.array(
        [np.sum(core_genome_alignment == nuc, axis=0) for nuc in dna_alphabet]
    )
    max_array = np.argmax(allele_counts, axis=0)
    snp_phages = {
        i: set(np.where(core_genome_alignment[:, i] == dna_alphabet[max_array[i]])[0])
        for i in snp_positions
    }
    snp_gaps = {
        i: set(np.where(core_genome_alignment[:, i] == "-")[0]) for i in snp_positions
    }

    # initialize dictionary
    dict = {
        "position1": [],
        "position2": [],
        "p1": [],
        "p2": [],
        "ld": [],
        "bg_ld": [],
    }
    # if statement outside of loop for computational efficiency
    if not large_memory:
        for i, ind1 in enumerate(tqdm(snp_positions[:-1])):
            # load the gaps and majority-allele phages for this snp position
            gaps1, all_phages1 = snp_gaps[ind1], snp_phages[ind1]
            for ind2 in snp_positions[i + 1 :]:
                # ignore rows (phages) where either snp position has a gap
                gap_phages = gaps1.union(snp_gaps[ind2])
                phage_count = n_phages - len(gap_phages)
                if len(gap_phages) == 0:
                    phages2 = snp_phages[ind2]
                    phages1 = all_phages1
                else:
                    phages1 = all_phages1.difference(gap_phages)
                    phages2 = snp_phages[ind2].difference(gap_phages)
                count1, count2 = len(phages1), len(phages2)

                # ensure thresholds still met with removed gap rows
                if phage_count - count1 < 2:
                    continue
                elif phage_count - count2 < 2:
                    continue
                elif count1 < 2 or count2 < 2:
                    continue

                # record the frequency of the dominant alleles
                p1 = count1 / phage_count
                p2 = count2 / phage_count

                # calculate linkage disequilibrium and the random expectation
                count12 = len(phages1.intersection(phages2))
                ld = calculate_linkage(p1, p2, count12 / phage_count)
                b_ld = background_linkage(count1, p1, count2, p2, phage_count)

                # record to dictionary
                dict["position1"].append(ind1)
                dict["position2"].append(ind2)
                dict["p1"].append(round(p1, 2))
                dict["p2"].append(round(p2, 2))
                dict["ld"].append(round(ld, 4))
                dict["bg_ld"].append(round(b_ld, 4))
    else:
        # same code, just with an additional if statement for performance purposes
        for i, ind1 in enumerate(tqdm(snp_positions[:-1])):
            gaps, phages1 = snp_gaps[ind1], snp_phages[ind1]
            for ind2 in snp_positions[i + 1 :]:
                # for large memory, only consider pairs of SNPs closer than 10kbp
                if ind2 - ind1 > 10000:
                    continue
                gap_phages = gaps.union(snp_gaps[ind2])
                phage_count = n_phages - len(gap_phages)
                if len(gap_phages) == 0:
                    phages2 = snp_phages[ind2]
                else:
                    phages1 = phages1.difference(gap_phages)
                    phages2 = snp_phages[ind2].difference(gap_phages)
                count1, count2 = len(phages1), len(phages2)
                if phage_count - count1 < 2:
                    continue
                elif phage_count - count2 < 2:
                    continue
                elif count1 < 2 or count2 < 2:
                    continue
                count12 = len(phages1.intersection(phages2))
                p1 = count1 / phage_count
                p2 = count2 / phage_count
                ld = calculate_linkage(p1, p2, count12 / phage_count)
                b_ld = background_linkage(count1, p1, count2, p2, phage_count)
                dict["position1"].append(ind1)
                dict["position2"].append(ind2)
                dict["p1"].append(round(p1, 2))
                dict["p2"].append(round(p2, 2))
                dict["ld"].append(round(ld, 4))
                dict["bg_ld"].append(round(b_ld, 4))

    # save the dataframe
    background_linkage.cache_clear()
    calculate_linkage.cache_clear()
    df = pd.DataFrame(dict)
    df.to_feather(args.output)


if __name__ == "__main__":
    main()
