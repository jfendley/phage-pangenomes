"""
This script modifies a JSON file for a group of phages that contains information to create a pangenome graph.
This script specifically cyclically permutes a few of the phages to create a group with
    a perfectly syntenic core genome.

Author: Jemma M. Fendley
"""

import argparse, json
from statistics import mode
import numpy as np

# create_edge takes sequential phams (nodes) and converts them to an edge object
from utils import create_edge


def main():
    parser = argparse.ArgumentParser(description="build syntenic group JSON")
    parser.add_argument(
        "-o", "--output", help="syntenic group JSON file", required=True, type=str
    )
    parser.add_argument(
        "-i", "--input", help="nonsyntenic group JSON file", type=str, required=True
    )
    parser.add_argument(
        "-n",
        "--nonsyntenic_info",
        help="nonsyntenic information dictionary",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-s", "--synteny", help="synteny type of the group", type=str, required=True
    )
    args = parser.parse_args()

    if args.synteny != "cyclic":
        raise ValueError("This script only works with cyclic groups")

    # load the JSON files
    with open(args.input) as f:
        G = json.load(f)
    with open(args.nonsyntenic_info) as f:
        nonsyntenic_dict = json.load(f)

    # extract the nonsyntenic phages
    group_nonsyntenic_info = [
        x for x in nonsyntenic_dict["nonsyntenic_groups"] if x["name"] == G["name"]
    ]
    assert (
        len(group_nonsyntenic_info) == 1
    ), "Error: nonsyntenic group not found, or found more than once"
    group_nonsyntenic_info = group_nonsyntenic_info[0]
    nonsyntenic_phages = group_nonsyntenic_info["nonsyntenic_phages"]

    # split the phage paths into those that follow the consensus ordering and those that don't
    syntenic_paths = [x for x in G["paths"] if x["ID"] not in nonsyntenic_phages]
    nonsyntenic_paths = [x for x in G["paths"] if x["ID"] in nonsyntenic_phages]

    # remove the strandedness to get the consenus core pham ordering
    core_pham_ordering = [x[:-1] for x in group_nonsyntenic_info["orderings"][0]]
    first_core = core_pham_ordering[0]

    # find the most commmon position that the first core pham is located
    # This will be the anchor of the cyclic permutations
    first_core_start_index = mode(
        [
            [x["pham_ID"] for x in path["path"]].index(first_core)
            for path in syntenic_paths
        ]
    )
    for path in nonsyntenic_paths:
        new_path_dict = {"ID": path["ID"], "name": path["name"]}  # initialize new path

        # find the location of the first core pham in the nonsyntenic ordering
        reduced_path = [x["pham_ID"] for x in path["path"]]
        first_core_loc = reduced_path.index(first_core)

        # cylically permute the order of the phams
        new_start_index = first_core_loc - first_core_start_index
        new_index = [
            x % len(path["path"])
            for x in range(new_start_index, new_start_index + len(path["path"]))
        ]
        new_path = [path["path"][i] for i in new_index]
        new_path_dict["path"] = new_path

        # check that the new order was correctly cyclically permuted
        new_reduced_path = [x["pham_ID"] for x in new_path_dict["path"]]
        new_first_core_loc = new_reduced_path.index(first_core)
        assert new_first_core_loc == first_core_start_index, "Error in permutation"

        syntenic_paths.append(new_path_dict)

    assert len(syntenic_paths) == len(G["paths"]), "Error in number of phages"

    # reorder the new paths to match the original ordering
    original_order = [x["ID"] for x in G["paths"]]
    new_order = [x["ID"] for x in syntenic_paths]
    ordering = [new_order.index(x) for x in original_order]
    paths = [syntenic_paths[i] for i in ordering]

    # recalculate the edges since there will be some descrepancy at the ends
    path_lists = [x["path"] for x in paths]
    edges_per_path = [
        [
            create_edge(path_list[x], path_list[x + 1])
            for x in np.arange(len(path_list) - 1)
        ]
        for path_list in path_lists
    ]

    # consolidate all of the edges
    all_edges = [x for y in edges_per_path for x in y]
    unique_edges = np.unique(all_edges)
    edges = [
        {
            "ID": x,
            "pham1": x.split("_")[0],
            "strand1": x.split("_")[1],
            "pham2": x.split("_")[2],
            "strand2": x.split("_")[3],
            "total_count": all_edges.count(x),
            "phage_count": len(paths) - [y.count(x) for y in edges_per_path].count(0),
        }
        for x in unique_edges
    ]

    # save the dictionary
    dict = {
        "name": G["name"],
        "type": G["type"],
        "phams": G["phams"],
        "edges": edges,
        "paths": paths,
    }
    with open(args.output, "w") as f:
        json.dump(dict, f, indent=4)


if __name__ == "__main__":
    main()
