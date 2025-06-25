"""
This script calculates pairwise distance metrics for the phages in the group.

Author: Jemma M. Fendley
"""

import numpy as np
import argparse, json
import pandas as pd
from Bio import AlignIO
from utils import calculate_hamming  # calculates hamming distance


def get_phage_info(phage_file):
    """
    Returns the phage ID and the list of genes from the phage genes file.
    """
    phage_ID = phage_file.split("/")[-1].split("_")[0]
    with open(phage_file) as f:
        phage_dict = json.load(f)
    return phage_ID, phage_dict["gene_list"]


def pair_ID(phage1, phage2):
    """
    Returns a unique ID for a pair of phages
    """
    if phage1 == phage2:
        return phage1 + "_" + phage2
    else:
        first_phage, second_phage = min(phage1, phage2), max(phage1, phage2)
        assert first_phage != second_phage
    return first_phage + "_" + second_phage


def main():
    parser = argparse.ArgumentParser(
        description="calculates pairwise distance metrics for all pairs of phages in a group"
    )
    parser.add_argument(
        "-o", "--output", help="name of output TSV file", type=str, required=True
    )
    parser.add_argument(
        "-s", "--summary", help="group summary JSON file", type=str, required=True
    )
    parser.add_argument(
        "-g", "--group_json", help="group JSON file", type=str, required=True
    )
    parser.add_argument(
        "-p",
        "--phages",
        nargs="+",
        help="list of phage gene files",
        required=True,
    )
    parser.add_argument(
        "-c", "--core", nargs="+", help="dna alignments of core phams", required=True
    )
    args = parser.parse_args()

    # retrieve the genome lengths from the summary JSON
    with open(args.summary) as f:
        summary_dict = json.load(f)
    genome_length_dict = summary_dict["genome_lengths"]

    # create an array with concatenated alignments of all of the core phams
    n_phages = len(args.phages)
    name_list, core_genome_alignment = [], np.empty([n_phages, 0])  # initialize
    for core_file_path in args.core:
        # extract the alignment and add to the core genome alignment array
        alignment = AlignIO.read(core_file_path, "fasta")
        alignment_array = np.array(alignment)
        core_genome_alignment = np.concatenate(
            (core_genome_alignment, alignment_array), axis=1
        )

        # record the order of the phages in the alignment
        name_list.append([record.id for record in alignment])

    # double check all of the alignments had the phages in the same order
    assert all(
        name_test == name_list[0] for name_test in name_list
    ), "Error with order in alignment file"
    # create dictionary of phage name to index in the alignment array
    phage_to_index = {name_list[0][i]: i for i in range(len(name_list[0]))}
    phage_list = list(phage_to_index.keys())

    # create a phage to phams dictionary and a phage to duplicate phams dictionary
    with open(args.group_json) as f:
        G = json.load(f)
    group_name = G["name"]
    pham_dict = {
        phage_ID: [
            x["pham_ID"]
            for x in [y for y in G["paths"] if y["ID"] == phage_ID][0]["path"]
        ]
        for phage_ID in phage_list
    }
    duplicate_phams = {
        phage_ID: set(
            [z for z in pham_dict[phage_ID] if pham_dict[phage_ID].count(z) > 1]
        )
        for phage_ID in phage_list
    }

    pairwise_dict_list = []  # initialize

    # iterate over all pairs of phages
    for i in range(n_phages - 1):
        phage_ID_1, gene_list_1 = get_phage_info(args.phages[i])

        for j in range(i + 1, n_phages):
            phage_ID_2, gene_list_2 = get_phage_info(args.phages[j])

            # set up for easy referencing later
            genome_lengths = [
                genome_length_dict[phage_ID_1],
                genome_length_dict[phage_ID_2],
            ]
            gene_lists = [gene_list_1, gene_list_2]

            # record the IDs of the phams shared
            shared_phams_all = set(pham_dict[phage_ID_1]).intersection(
                set(pham_dict[phage_ID_2])
            )

            # calculate the Jaccard distance
            jaccard = 1 - len(shared_phams_all) / len(
                set(pham_dict[phage_ID_1]).union(set(pham_dict[phage_ID_2]))
            )

            # calculate the Hamming distance
            hamming = calculate_hamming(
                core_genome_alignment,
                phage_to_index[phage_ID_1],
                phage_to_index[phage_ID_2],
            )

            # The next section will calculate the "coverage" of each pair, i.e. what percentage
            #   of the phages are covered by pairwise shared phams.

            # phams that are duplicate in either phage
            shared_phams_duplicates = set(
                [
                    x
                    for x in shared_phams_all
                    if x
                    in duplicate_phams[phage_ID_1].union(duplicate_phams[phage_ID_2])
                ]
            )

            # phams that are not duplicate in either
            shared_phams = [
                x for x in shared_phams_all if x not in shared_phams_duplicates
            ]

            # the next section processes the duplicate phams
            add_indices_list = [[], []]
            if len(shared_phams_duplicates) > 0:
                phage_ID_list = [phage_ID_1, phage_ID_2]
                for duplicate_pham in shared_phams_duplicates:
                    # if they have the same number of copies of the pham, then can just add it to the list
                    if pham_dict[phage_ID_1].count(duplicate_pham) == pham_dict[
                        phage_ID_2
                    ].count(duplicate_pham):
                        shared_phams.append(duplicate_pham)

                    # otherwise we only want to add the number that they share
                    else:
                        shared_count = min(
                            pham_dict[phage_ID_1].count(duplicate_pham),
                            pham_dict[phage_ID_2].count(duplicate_pham),
                        )
                        assert shared_count >= 1, "Error with shared phams"

                        for relevant_index, relevant_gene_list in enumerate(gene_lists):
                            # add the first X occurences to the list, where X is the number shared
                            add_indices = [
                                x
                                for x, y in enumerate(
                                    pham_dict[phage_ID_list[relevant_index]]
                                )
                                if y == duplicate_pham
                            ][:shared_count]

                            # double check that the indices added correspond to the correct pham
                            assert all(
                                [
                                    relevant_gene_list[x]["phams"][0] == duplicate_pham
                                    for x in add_indices
                                ]
                            ), "Genes in incorrect order"

                            add_indices_list[relevant_index].extend(add_indices)

            row = {
                "group": group_name,
                "pair_ID": pair_ID(phage_ID_1, phage_ID_2),
                "phage_1": phage_ID_1,
                "phage_2": phage_ID_2,
                "length_1": genome_lengths[0],
                "length_2": genome_lengths[1],
                "Jaccard": jaccard,
                "Hamming": hamming,
            }

            relevant_keys = ["percent_shared_1", "percent_shared_2"]
            for i in range(2):
                # the list of genes (including duplicates) that is shared
                shared_list = [
                    x
                    for y, x in enumerate(gene_lists[i])
                    if x["phams"][0] in shared_phams or y in add_indices_list[i]
                ]

                # the follow calculates what percentage of the genome is covered by the shared genes
                starts = [x["Start"] for x in shared_list]
                ends = [x["Stop"] for x in shared_list]

                tot_length, initial_start = 0, starts[0]
                for j in range(len(ends) - 1):
                    if starts[j + 1] > ends[j]:
                        tot_length += ends[j] - initial_start
                        initial_start = starts[j + 1]
                tot_length += ends[-1] - initial_start
                row[relevant_keys[i]] = 100 * tot_length / genome_lengths[i]
            pairwise_dict_list.append(row)

    df = pd.DataFrame(pairwise_dict_list)

    # a few checks to ensure everything is the correct size
    assert len(df) == len(set(df["pair_ID"]))
    assert len(df) == int((n_phages * (n_phages - 1)) / 2)

    # save the dataframe
    df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
