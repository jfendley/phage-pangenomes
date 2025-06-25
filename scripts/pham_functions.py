"""
This script creates a table with the functions of all of the phams in the group.

Author: Jemma M. Fendley
"""

import json, argparse
import numpy as np
import pandas as pd
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(description="create a dataframe of pham functions")
    parser.add_argument(
        "-o", "--output", help="output pham functions TSV file", type=str, required=True
    )
    parser.add_argument(
        "-a",
        "--accessory",
        help="accessory pham locations JSON",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-d", "--defense", help="pham defense info TSV file", type=str, required=True
    )
    parser.add_argument(
        "-c", "--core", help="list of core phams", nargs="+", required=True
    )
    parser.add_argument(
        "-p", "--phages", help="phage gene files", nargs="+", required=True
    )
    args = parser.parse_args()

    # load the defense dataframe and convert to dictionary
    defense_df = pd.read_csv(args.defense, sep="\t")
    defense_dict = pd.Series(
        defense_df.activity.values, index=defense_df.pham_ID
    ).to_dict()

    # load the accessory locations json and record the group name
    group_name = args.accessory.split("/")[-1].split("_")[0]
    with open(args.accessory) as f:
        all_locations = json.load(f)

    # simplify the dictionary to a set of values instead of a list (i.e. avoid listing the same junction twice)
    accessory_locations = {
        key: [int(x) for x in set(value)] for key, value in all_locations.items()
    }

    # unpack all possible locations in order to count how many are in each junction
    all_locations = [int(x) for y in list(accessory_locations.values()) for x in y]
    n_core = len(args.core)

    # count how many accessory phams are in each junction and take the top 4 with the most accessory phams
    junctions, junction_counts = np.unique(all_locations, return_counts=True)
    top_four = np.argsort(junction_counts)[::-1][:4]
    sorted_junctions, sorted_counts = junctions[top_four], junction_counts[top_four]

    # threshold of at least 5 accessory phams chosen so that each group has at least 1 hotspot.
    hotspot_indices = sorted_junctions[np.where(sorted_counts[:4] >= 5)]

    hotspots_gene_list = [
        x
        for x, y in accessory_locations.items()
        if any([z in y for z in hotspot_indices])
    ]
    # if a pham is in multiple locations and one of those locations is a hotspot, we consider it to be a hotspot

    # add the flanking core genes to the list of hotspot genes
    for i in hotspot_indices:
        if i == 0:
            hotspots_gene_list.append(args.core[0])
        elif i == n_core:
            hotspots_gene_list.append(args.core[-1])
        else:
            hotspots_gene_list.append(args.core[i - 1])
            hotspots_gene_list.append(args.core[i])

    hotspots_gene_list = list(set(hotspots_gene_list))  # remove any duplicates

    # create a dictionary of phams toa  list of their function annotations
    function_dict = defaultdict(list)
    for phage in args.phages:
        with open(phage) as f:
            gene_info = json.load(f)
        for gene in gene_info["gene_list"]:
            function_dict[gene["phams"][0]].append(gene["Notes"])

    dict_list = []
    for pham_ID, all_functions in function_dict.items():
        # find the list of pham functions and their counts
        functions, counts = np.unique(all_functions, return_counts=True)

        # Find the majority function
        majority_function = functions[np.argsort(counts)[::-1]][0]

        # remove quatation marks in the annotation
        majority_function = majority_function[2:-1]

        # tidy the annotations
        if "DNA-binding" in majority_function:
            majority_function = majority_function.replace("DNA-binding", "DNA binding")
        elif majority_function == "":
            majority_function = "no majority annotation"

        # if pham is in the DefenseFinder results, override any previous function and record it as Defense/AntiDefense
        if int(pham_ID) in list(defense_dict.keys()):
            majority_function = defense_dict[int(pham_ID)]

        # save to the list, recording if the pham is core and in/near a hostpot
        row = {
            "group": group_name,
            "pham_ID": pham_ID,
            "majority_function": majority_function,
            "core": pham_ID in args.core,
            "hotspot": pham_ID in hotspots_gene_list,
        }
        dict_list.append(row)

    # save the dataframe
    df = pd.DataFrame(dict_list)
    df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
