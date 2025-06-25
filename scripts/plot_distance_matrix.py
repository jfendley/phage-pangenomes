"""
This script calculates the pairwise Hamming distance of the core genome, as well as the pairwise Jaccard
    distance (gene content similarity) of all pairs of phages in a group. It plots them both on a
    heircharically-clustered distance matrix. It also suggets subgroups of phages based on the distance matrix.

Author: Jemma M. Fendley
"""

import numpy as np, pandas as pd
import argparse, json
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
from Bio import AlignIO
import parameters  # preset matplotlib formatting

# sort_matrix sorts a distance matrix by heirarchically clustering
# calculate_hamming calculates the hamming distance of two rows in an alignment
from utils import sort_matrix, calculate_hamming


def main():
    parser = argparse.ArgumentParser(
        description="create a pairwise distance matrix of a group of phages"
    )
    parser.add_argument(
        "-o",
        "--output",
        help="TSV file with summary pairwise statistics",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-f", "--output_fig", help="distance matrix file", type=str, required=True
    )
    parser.add_argument(
        "-s", "--output_subgroups", help="subgroups JSON file", type=str, required=True
    )
    parser.add_argument(
        "-c",
        "--core_genome",
        help="core genome alignment fasta",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-g", "--group_json", help="group JSON file", type=str, required=True
    )
    args = parser.parse_args()

    group_name = args.core_genome.split("/")[-1].split("_")[0]

    # create a phage to set of phams dictionary
    with open(args.group_json) as f:
        group_info = json.load(f)
    phage_phams_dictionary = {
        x["ID"]: set([y["pham_ID"] for y in x["path"]]) for x in group_info["paths"]
    }
    n_phages = len(group_info["paths"])

    # load the core genome alignment and record the order of the phages in it
    core_genome_alignment = AlignIO.read(args.core_genome, "fasta")
    phage_order = [record.id for record in core_genome_alignment]
    core_genome_array = np.array(core_genome_alignment)

    # calculate the hamming distance of core genome and jaccard distance of shared phams
    #   for each pair of phages
    hamming, jaccard = np.zeros((n_phages, n_phages)), np.zeros((n_phages, n_phages))
    for i in range(n_phages - 1):
        for j in range(i + 1, n_phages):
            hamming[i, j] = calculate_hamming(core_genome_array, i, j)
            jaccard[i, j] = 1 - len(
                phage_phams_dictionary[phage_order[i]].intersection(
                    phage_phams_dictionary[phage_order[j]]
                )
            ) / len(
                phage_phams_dictionary[phage_order[i]].union(
                    phage_phams_dictionary[phage_order[j]]
                )
            )

    # find the indicies of the upper and lower triangle of the matrix
    i_lower, i_upper = np.tril_indices(n_phages, -1), np.triu_indices(n_phages, 1)
    i_upper_1, i_upper_2 = i_upper

    # populate the rest of the matrix using symmetry
    hamming[i_lower], jaccard[i_lower] = hamming.T[i_lower], jaccard.T[i_lower]

    # use the distance matrix to create a dendrogram, and calculate the branch lengths
    dn = dendrogram(linkage(squareform(hamming), "complete"), no_plot=True)
    hamming_branch_lengths = [
        y for z in [[x[1] - x[0], x[2] - x[3]] for x in dn["dcoord"]] for y in z
    ]

    # flatten the matrix to calculate statistics, then save them in a dataframe
    flattened_hamming, flattened_jaccard = hamming[i_upper], jaccard[i_upper]
    df = pd.DataFrame(
        [
            {
                "group": group_name,
                "mean_hamming": np.mean(flattened_hamming),
                "std_hamming": np.std(flattened_hamming, ddof=1),
                "max_hamming": np.max(flattened_hamming),
                "mean_hamming_branch_length": np.mean(hamming_branch_lengths),
                "std_hamming_branch_length": np.std(hamming_branch_lengths, ddof=1),
                "max_hamming_branch_length": np.max(hamming_branch_lengths),
                "mean_jaccard": np.mean(flattened_jaccard),
                "std_jaccard": np.std(flattened_jaccard, ddof=1),
                "max_jaccard": np.max(flattened_jaccard),
            }
        ]
    )
    df.to_csv(args.output, sep="\t", index=False)

    # sort the hamming matrix by heirarchically clustering
    sorted_hamming, res_order, res_linkage = sort_matrix(hamming, "complete")

    # find and save heirachical clusters
    max_distance = np.mean(hamming) + np.std(hamming) / 2
    grouping = fcluster(res_linkage, max_distance, criterion="distance")
    subgroups_dictionary = {
        int(x): [phage_order[i] for i in np.where(grouping == x)[0]]
        for x in set(grouping)
    }
    with open(args.output_subgroups, "w") as f:
        json.dump(subgroups_dictionary, f, indent=4)

    # create a matrix with bottom half hamming, top half jaccard
    distance_matrix = np.zeros((n_phages, n_phages))
    distance_matrix[i_lower] = sorted_hamming[i_lower]
    distance_matrix[i_upper_1, i_upper_2] = jaccard[
        [res_order[i] for i in i_upper_1], [res_order[j] for j in i_upper_2]
    ]

    # plot the distance matrix
    fig, axes = plt.subplots(layout="constrained", figsize=(4, 4))
    im = axes.matshow(distance_matrix, cmap="rainbow")
    axes.set_xlabel("Jaccard distance (gene content difference)")
    axes.xaxis.set_label_position("top")
    axes.set_ylabel("Hamming distance (ANI distance)")
    axes.set_title(str(group_name) + ": distance matrix")
    plt.colorbar(im)
    fig.savefig(args.output_fig)


if __name__ == "__main__":
    main()
