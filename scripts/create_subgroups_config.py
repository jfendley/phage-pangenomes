"""
This script creates the config files necessary for the next snakemake workflow.

Author: Jemma M. Fendley
"""

import argparse, json
import numpy as np


def main():
    parser = argparse.ArgumentParser(description="create config files for subgroups")
    parser.add_argument(
        "-o",
        "--output_group",
        help="output subgroup to group dictionary",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-p",
        "--output_phages",
        help="output subgroup to phages dictionary",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-g",
        "--groups",
        nargs="+",
        required=True,
        help="groups subgroup info JSON files",
    )
    args = parser.parse_args()

    phages_dict, group_dict = {}, {}  # initialize

    # iterate through all of the groups' subgroup information JSON files
    for group in args.groups:
        group_name = group.split("/")[-1].split("_")[0]
        with open(group) as f:
            subgroup_list = json.load(f)
        # for each subgroup record the name and the phages
        for subgroup, phage_list in subgroup_list.items():
            subgroup_name = group_name + "-" + subgroup
            phages_dict[subgroup_name] = phage_list
            group_dict[subgroup_name] = group_name

    # save the dictionaries
    with open(args.output_group, "w") as f:
        json.dump(group_dict, f, indent=4)
    with open(args.output_phages, "w") as f:
        json.dump(phages_dict, f, indent=4)


if __name__ == "__main__":
    main()
