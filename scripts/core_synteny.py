"""
This script investigates the synteny of the core genome for each group.

Author: Jemma M. Fendley
"""

import json, argparse
import numpy as np


def main():
    parser = argparse.ArgumentParser(description="checks each group for core synteny")
    parser.add_argument(
        "-o",
        "--output",
        help="JSON file that describes synteny",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-g", "--groups", nargs="+", help="the JSON files for all groups", required=True
    )
    args = parser.parse_args()

    count_syntenic, dict_list = 0, []  # initialize

    for group_json in args.groups:
        with open(group_json) as f:
            G = json.load(f)
        group_name, paths = G["name"], G["paths"]

        # find the core phams
        core_phams = [
            x["pham_ID"]
            for x in G["phams"]
            if x["total_count"] == len(paths) and x["phage_count"] == len(paths)
        ]

        # extract only the order (including standedness) of the core phams for each genome/path
        core_paths = [
            [
                y["pham_ID"] + y["strand"]
                for y in x["path"]
                if y["pham_ID"] in core_phams
            ]
            for x in paths
        ]

        # find all of the possible core pham orderings
        unique_unsorted, counts_unsorted = np.unique(
            core_paths, axis=0, return_counts=True
        )

        if len(unique_unsorted) > 1:  # nonsynentic
            row = {"name": group_name}  # intialize dictionary

            # sort the possible core paths in order of most common to least common
            ordering = np.argsort(counts_unsorted)[::-1]
            unique, counts = unique_unsorted[ordering], counts_unsorted[ordering]
            assert np.sum(counts) == len(paths), "Error: bug with number of phages"

            row["counts"] = [int(x) for x in counts]  # necessary to save to JSON
            row["percent_syntenic"] = 100 - 100 * np.sum(counts[1:]) / np.sum(counts)
            row["orderings"] = [list(x) for x in unique]

            # check if the orderings are cyclic permutations of each other
            if np.array(
                [
                    " ".join(map(str, x)) in " ".join(map(str, list(unique[0]) * 2))
                    for x in unique
                ]
            ).all():
                row["synteny"] = "cyclic"

                # record the phages that do not follow the consensus orderings
                row["nonsyntenic_phages"] = []
                for unique_list in unique[1:]:
                    nonsyntenic_phages = [
                        paths[i]["ID"]
                        for i in range(len(core_paths))
                        if all(core_paths[i] == unique_list)
                    ]
                    row["nonsyntenic_phages"].extend(nonsyntenic_phages)
            else:
                row["synteny"] = "nonsyntenic"
                row["hamming"], row["cayley"], row["nonsyntenic_phages"] = [], [], []

                consensus = unique[0]
                for unique_list in unique[1:]:
                    different = np.where(unique_list != consensus)[0]

                    # Hamming distance between the ordering and the consensus ordering
                    row["hamming"].append(float(len(different) / len(consensus)))

                    # The calculation of the Cayley distance is done by hand (hard-coded)
                    # To better use this code for other applications, one would need to properly
                    # write an algorithm to calculate the Cayley distance properly
                    if group_name == "EC":
                        cayley = 9
                    else:
                        cayley = len(different) - 1
                        assert cayley != 0
                    row["cayley"].append(int(cayley))

                    # record the phages that do not follow the consensus orderings
                    nonsyntenic_phages = [
                        paths[i]["ID"]
                        for i in range(len(core_paths))
                        if all(core_paths[i] == unique_list)
                    ]
                    row["nonsyntenic_phages"].extend(nonsyntenic_phages)
                row["max_cayley"] = int(np.max(row["cayley"]))
            # consistency check
            assert len(row["nonsyntenic_phages"]) == np.sum(
                row["counts"][1:]
            ), "Error in number of nonsyntenic phages"
            dict_list.append(row)
        else:
            count_syntenic += 1

    # Save the information into a JSON file
    dict = {"n_syntenic": count_syntenic, "nonsyntenic_groups": dict_list}
    with open(args.output, "w") as f:
        json.dump(dict, f, indent=4)


if __name__ == "__main__":
    main()
