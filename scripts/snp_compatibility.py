"""
This scripts investigates the compatibility of SNPs (whether or not they pass the 4-allele test)
    as a function of distance along the core genome alignment.

Author: Jemma M. Fendley
"""

import numpy as np, pandas as pd
from Bio import AlignIO
import argparse
from tqdm import tqdm
import matplotlib.pyplot as plt
import parameters  # preset matplotlib formatting
import matplotlib as mpl

# find_snp_positions finds the positions of SNPs given a core genome alignment
from utils import find_snp_positions

mpl.rcParams["legend.handletextpad"] = 0.4


def calculate_expected_distance(n_snps, N):
    """
    This function calculates the expected distance between expected recurrent mutation.

    N: (integer) the number of positions in the core genome alignment
    n_snps: (integer) the number of positions with SNPs in the core genome alignment

    Given N, n_snps, this function estimates the mutation number K and then the expected
    number of sites with more than one mutation (n_reccurrent) and then returns the expected
    distance between recurrent mutations.

    It also returns the expected number of SNPs between reccurent mutations.

    For more detail on this estimation, please see the SI F SNP compatibility
    """
    K = np.log((N - n_snps) / N) / np.log((N - 1) / N)
    n_reccurent = N - N * ((N - 1) / N) ** K - K * ((N - 1) / N) ** (K - 1)
    recurrent_distance = N / n_reccurent
    return round(recurrent_distance), round(n_snps / n_reccurent)


def four_allele_test(first_phages_1, second_phages_1, first_phages_2, second_phages_2):
    """
    This function performs the four allele test on two biallelic SNPs, labelled 1 and 2.

    The "first phages" are those that have the majority allele, and the "second phages" are
    those that have the minority allele.

    The input are four sets, and the output is True if it FAILS the test (has all four possible
    allele combinations), and False if it PASSES (does not have all four possible allele combinations).
    """
    num1_1 = len(first_phages_1.intersection(first_phages_2))
    num2_2 = len(second_phages_1.intersection(second_phages_2))
    num1_2 = len(first_phages_1.intersection(second_phages_2))
    num2_1 = len(second_phages_1.intersection(first_phages_2))
    if num1_1 > 0 and num1_2 > 0 and num2_1 > 0 and num2_2 > 0:
        return True
    else:
        return False


def main():
    parser = argparse.ArgumentParser(description="test SNP compatibility")
    parser.add_argument(
        "-o",
        "--output",
        help="TSV file with summary statistics",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-f", "--output_fig", help="figure file", type=str, required=True
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
        "--position_to_pham",
        help="core genome position to pham TSV",
        type=str,
        required=True,
    )
    args = parser.parse_args()

    np.random.seed(0)  # for reproducibility
    group_name = args.core_genome.split("/")[-1].split("_")[0]

    # find the positions in the core genome alignments which have SNPs, and record also
    #   the counts of each alelle at each site
    dna_alphabet = ["G", "T", "A", "C"]
    core_genome_alignment = np.array(AlignIO.read(args.core_genome, "fasta"))
    n_positions = core_genome_alignment.shape[1]
    allele_counts = np.array(
        [np.sum(core_genome_alignment == nuc, axis=0) for nuc in dna_alphabet]
    )
    all_snp_positions, all_biallelic_snp_positions = find_snp_positions(
        core_genome_alignment
    )

    # find the overlap positions and remove from analysis
    position_df = pd.read_csv(args.position_to_pham, sep="\t")
    overlap_positions = set(position_df[position_df["overlap"]]["position"])
    snp_positions = np.sort(list(set(all_snp_positions) - overlap_positions))
    biallelic_snp_positions = np.sort(
        list(set(all_biallelic_snp_positions) - overlap_positions)
    )

    # record the expected n. snps and distance between recurrent mutations
    dict = {"group": group_name}
    dict["expected_distance"], dict["expected_n_snps"] = calculate_expected_distance(
        len(snp_positions), n_positions
    )
    dict["expected_distance_biallelic"], dict["expected_n_snps_biallelic"] = (
        calculate_expected_distance(len(biallelic_snp_positions), n_positions)
    )

    # create dictionaries that record the phages at each site with the top two alleles
    second_allele, max_allele = np.argsort(allele_counts, axis=0)[-2:]
    first_phage_dict, second_phage_dict = {}, {}
    for i in snp_positions:
        first_phage_dict[i] = set(
            np.where(core_genome_alignment[:, i] == dna_alphabet[max_allele[i]])[0]
        )
        second_phage_dict[i] = set(
            np.where(core_genome_alignment[:, i] == dna_alphabet[second_allele[i]])[0]
        )

    # shuffle the order of the core genome alignment for comparison
    new_ordering = np.random.permutation(n_positions)
    random_dict = {i: x for i, x in enumerate(new_ordering)}
    random_reverse_dict = {x: i for i, x in enumerate(new_ordering)}

    # initialize figure
    fig, axes = plt.subplots(
        nrows=2,
        ncols=2,
        figsize=(7, 7),
        layout="constrained",
        gridspec_kw={"width_ratios": [1, 1]},
    )
    binsize = 5
    label_list = ["", "_biallelic"]

    # index 1 represents if using data (0) or random (1)
    for index1, label in enumerate(["data", "random"]):
        if index1 == 0:
            snp_positions_sample = snp_positions
            biallelic_snp_positions_sample = biallelic_snp_positions
        else:
            snp_positions_sample = np.sort([random_dict[i] for i in snp_positions])
            biallelic_snp_positions_sample = np.sort(
                [random_dict[i] for i in biallelic_snp_positions]
            )
        # index 2 represents if using all snps (0) or only biallelic snps (1)
        for index2, snps in enumerate(
            [snp_positions_sample, biallelic_snp_positions_sample]
        ):
            # initialize, start with the first snp
            label2, start = label_list[index2], snps[0]
            start_list, n_snps, list_n_snps = [start], 1, []
            reset_variable = True

            # Iterate through all of the SNPs. Start with a "start" position (originally 0).
            #   check if all pairwise SNPs from start to i are compatible. Keep iterating until
            #   you reach incompatibilty. Then restart and change the "start" position
            for i in tqdm(snps[1:]):
                if reset_variable == True:
                    first_phages, second_phages = [], []
                    sites = [x for x in snps if x <= i and x >= start]
                else:
                    sites = [i]
                    n_snps += 1
                # add the new sites to the list of phages
                for new_site in sites:
                    if index1 == 0:
                        first_phages.append(first_phage_dict[new_site])
                        second_phages.append(second_phage_dict[new_site])
                    else:
                        first_phages.append(
                            first_phage_dict[random_reverse_dict[new_site]]
                        )
                        second_phages.append(
                            second_phage_dict[random_reverse_dict[new_site]]
                        )
                current_length = len(first_phages)
                reset_variable = False
                # check all pairwise compatiblity. If any fail 4-allele test, then reset start variable
                #   and continue onwards.
                for ind in range(current_length - 1):
                    if reset_variable == True:
                        break
                    for j in range(ind + 1, current_length):
                        if four_allele_test(
                            first_phages[ind],
                            second_phages[ind],
                            first_phages[j],
                            second_phages[j],
                        ):
                            start = i
                            start_list.append(i)
                            reset_variable = True
                            list_n_snps.append(n_snps)
                            n_snps = 1
                            break
                        else:
                            continue
            # once iteration complete, just need to add the last to list
            if reset_variable == False:
                list_n_snps.append(n_snps + 1)
            else:
                list_n_snps.append(n_snps)
            start_list.append(n_positions)

            interval_sizes = [
                start_list[i + 1] - start_list[i] for i in range(len(start_list) - 1)
            ]
            # quick sanity check
            list_n_snps_check = [
                np.count_nonzero((snps >= start_list[i]) & (snps < start_list[i + 1]))
                for i in range(len(start_list) - 1)
            ]
            assert list_n_snps == list_n_snps_check, "Error in script"

            # plot the distributions and save the means and maxes
            label3 = ["n_snps", "interval_size"]
            for index3, plot_data in enumerate([list_n_snps, interval_sizes]):
                my_bins = np.arange(-0.5, np.max(plot_data) + binsize, binsize)
                dict[label + label2 + "_mean_" + label3[index3]] = round(
                    float(np.mean(plot_data)), 5
                )
                dict[label + label2 + "_max_" + label3[index3]] = int(np.max(plot_data))
                axes[index2, index3].hist(
                    plot_data, bins=my_bins, alpha=0.5, label=label
                )
                if index1 == 0:
                    axes[index2, index3].axvline(
                        np.mean(plot_data),
                        linestyle="--",
                        color="k",
                        label=label + " mean",
                    )
            axes[index1, index2].set_yscale("log")

        # addd some labels
        axes[0, index1].set_ylabel("count")
        axes[1, index1].set_ylabel("count (bi-allelic only)")
        axes[index1, 0].set_xlabel("compatible interval size (n. snps)")
        axes[index1, 1].set_xlabel("compatible interval size (nucleotides)")
        axes[0, index1].set_ylabel("count")
        axes[1, index1].set_ylabel("count (bi-allelic only)")

    # add recurrent mutations information and legend
    distance_list = [dict["expected_distance"], dict["expected_distance_biallelic"]]
    n_snps_list = [dict["expected_n_snps"], dict["expected_n_snps_biallelic"]]
    for i in range(2):
        axes[i, 0].axvline(
            n_snps_list[i], color="k", label="expected recurrent mutation"
        )
        axes[i, 1].axvline(
            distance_list[i], color="k", label="expected recurrent mutation"
        )
    axes[0, 0].legend()

    fig.suptitle("group " + group_name, fontsize=24)
    fig.savefig(args.output_fig, dpi=300)
    df = pd.DataFrame([dict])
    df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
