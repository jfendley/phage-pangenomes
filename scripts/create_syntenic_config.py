"""
This script creates the config files necessary for the next snakemake workflow.

Author: Jemma M. Fendley
"""

import argparse, json
import numpy as np


def main():
    parser = argparse.ArgumentParser(
        description="create config files for syntenic groups"
    )
    parser.add_argument(
        "-c",
        "--core",
        help="output group to core phams dictionary",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-p",
        "--phages",
        help="output group to phages dictionary",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-n", "--nonsyntenic", help="nonsyntenic group info", type=str, required=True
    )
    parser.add_argument(
        "-g", "--groups", nargs="+", required=True, help="group JSON files"
    )
    args = parser.parse_args()

    phage_dict, core_dict = {}, {}  # initialize

    with open(args.nonsyntenic) as f:
        nonsyntenic_dict = json.load(f)
    # record the groups that are not syntenic (not including cyclic permutations)
    nonsyntenic_groups = {
        x["name"]: x
        for x in nonsyntenic_dict["nonsyntenic_groups"]
        if x["synteny"] == "nonsyntenic"
    }

    for group in args.groups:
        with open(group) as f:
            G = json.load(f)  # load group JSON
        name = G["name"]
        if name in list(nonsyntenic_groups.keys()):
            old_pham_list = [x["pham_ID"] for x in G["phams"]]  # all original phams

            nonsyntenic_info = nonsyntenic_groups[name]
            if nonsyntenic_info["percent_syntenic"] < 80:
                continue  # exclude groups with not a large majority syntenic
            if np.max(nonsyntenic_info["counts"]) < 10:
                continue  # exclude groups with fewer than 10 phages syntenic

            new_paths = [
                x
                for x in G["paths"]
                if x["ID"] not in nonsyntenic_info["nonsyntenic_phages"]
            ]  # exclude the nonsyntenic phages

            phage_dict[name] = [x["ID"] for x in new_paths]

            pham_only_paths = [[x["pham_ID"] for x in y["path"]] for y in new_paths]
            core_phams = [
                x
                for x in old_pham_list
                if all(
                    [pham_only_path.count(x) == 1 for pham_only_path in pham_only_paths]
                )
            ]  # core phams are in every phage exactly once
            core_dict[name] = core_phams
        else:
            # can save directly with no modifications, even for the cyclic groups
            core_dict[name] = [
                x["pham_ID"]
                for x in G["phams"]
                if x["phage_count"] == len(G["paths"])
                and x["total_count"] == len(G["paths"])
            ]
            phage_dict[name] = [x["ID"] for x in G["paths"]]

    # save the dictionaries
    with open(args.core, "w") as f:
        json.dump(core_dict, f, indent=4)
    with open(args.phages, "w") as f:
        json.dump(phage_dict, f, indent=4)


if __name__ == "__main__":
    main()
