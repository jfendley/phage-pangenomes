"""
This script builds a JSON file for a group of phages that contains information to create a pangenome graph.

Author: Jemma M. Fendley
"""

import numpy as np
import argparse, json
from collections import defaultdict

# create_edge takes sequential phams (nodes) and converts them to an edge object
# convert_strand converts the strandedness (forward or reverse) of a gene to the needed format
from utils import convert_strand, create_edge


def main():
    """
    The output of this script is a json containing phams, edges, and paths. This script is inspired by
    Pangraph (see https://docs.pangraph.org/).

    Phams: the "nodes" of the graph, phage gene families
    Edges: the "edges" or "links" of the graph, connections between nodes
    Paths: the "paths" through the graph. Each path represents a genome as an ordered list of phams
    """
    parser = argparse.ArgumentParser(
        description="creates a JSON file for the basis of pangenome graph"
    )
    parser.add_argument(
        "-o",
        "--output",
        help="output json file name",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-t",
        "--type",
        help="type of group, e.g. subgroup",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-d",
        "--dict",
        help="file with phage ID to phage name dictionary",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-p",
        "--phages",
        nargs="+",
        help="Phage gene files",
        required=True,
    )

    args = parser.parse_args()
    with open(args.dict) as f:
        name_dict = json.load(f)

    paths, phams, edges_per_path, phams_per_path = [], [], [], []  # initialize
    pham_length_dict = defaultdict(list)

    for phage_file in args.phages:
        phage_ID = phage_file.split("/")[-1].split("_")[0]
        with open(phage_file) as f:
            gene_list_dict = json.load(f)

        # Convert the list of genes into a simple path object
        path_list = [
            {"pham_ID": x["phams"][0], "strand": convert_strand(x["Orientation"])}
            for x in gene_list_dict["gene_list"]
        ]
        path = {"ID": phage_ID, "name": name_dict[phage_ID], "path": path_list}
        paths.append(path)

        # extract the edges/links
        edges_per_path.append(
            [
                create_edge(path_list[x], path_list[x + 1])
                for x in np.arange(len(path_list) - 1)
            ]
        )
        # extract the phams and record their lengths
        phams_per_path.append([x["pham_ID"] for x in path_list])
        for x in gene_list_dict["gene_list"]:
            pham_length_dict[x["phams"][0]].append(x["Length"])

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

    # consolidate all of the phams
    all_phams = [x for y in phams_per_path for x in y]
    unique_phams = np.unique(all_phams)
    phams = [
        {
            "pham_ID": x,
            "total_count": all_phams.count(x),
            "phage_count": len(paths) - [y.count(x) for y in phams_per_path].count(0),
            "mean_length": round(np.mean(pham_length_dict[x])),
            "read_length": round(np.mean(pham_length_dict[x])) * all_phams.count(x),
        }
        for x in unique_phams
    ]

    # save final dictionary
    dict = {
        "name": args.output.split("/")[-1].split("_")[0],
        "type": args.type,
        "phams": phams,
        "edges": edges,
        "paths": paths,
    }
    with open(args.output, "w") as f:
        json.dump(dict, f, indent=4)


if __name__ == "__main__":
    main()
