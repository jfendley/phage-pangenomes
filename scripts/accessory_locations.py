"""
This scripts create a dictionary of all of the locations of the accessory phams,
    with respect to the syntenic core genome.

Author: Jemma M. Fendley
"""

import json, argparse
import pandas as pd
import numpy as np
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(
        description="create dictionary for which core backbone junction each pham belongs to"
    )
    parser.add_argument(
        "-o", "--output", help="accessory locations JSON", type=str, required=True
    )
    parser.add_argument(
        "-i", "--input", help="group JSON file", type=str, required=True
    )
    parser.add_argument(
        "-c", "--core_phams", help="list of core phams", nargs="+", required=True
    )
    parser.add_argument(
        "-d", "--duplicate", help="include duplicate phams?", type=str, required=True
    )
    args = parser.parse_args()

    # load the group JSON file, convert the phams to a dataframe
    with open(args.input) as f:
        G = json.load(f)
    pham_df = pd.DataFrame(G["phams"]).set_index("pham_ID")

    # extract just the phams from the paths (neglect strandedness for this analysis)
    reduced_paths = {
        phage["ID"]: [x["pham_ID"] for x in phage["path"]] for phage in G["paths"]
    }
    num_phages = len(reduced_paths)

    # find the locations (indices) of the core phams for each of the phages
    core_indices = {
        phage: np.array(
            [i for i, pham in enumerate(reduced_path) if pham in args.core_phams]
        )
        for phage, reduced_path in reduced_paths.items()
    }

    if args.duplicate == "No":
        # do not include phams that appear multiple times in the same phage
        accessory_phams = pham_df[
            (pham_df["phage_count"] < num_phages)
            & (pham_df["total_count"] == pham_df["phage_count"])
        ].index.tolist()

        # find the location of each of the phams in each of the phages it is in,
        #   with respect to the core genome backbone. Each pham is located in a "junction".
        location_dict = defaultdict(list)
        for pham_ID in accessory_phams:
            for phage, reduced_path in reduced_paths.items():
                if pham_ID in reduced_path:
                    location = reduced_path.index(pham_ID)
                    if location < np.min(core_indices[phage]):
                        junction = 0
                    else:
                        junction = (
                            np.max(np.where(location > core_indices[phage])[0]) + 1
                        )
                    location_dict[pham_ID].append(junction)
    elif args.duplicate == "Yes":
        # include all phams that are not core (which are in every phage exactly once)
        accessory_phams = pham_df[
            (pham_df["phage_count"] < num_phages)
            | (
                (num_phages == pham_df["phage_count"])
                & (pham_df["total_count"] > pham_df["phage_count"])
            )
        ].index.tolist()

        # find the location(s) of each of the phams in each of the phages it is in,
        #   with respect to the core genome backbone. Each pham is located in a "junction".
        location_dict = defaultdict(list)
        for pham_ID in accessory_phams:
            for phage, reduced_path in reduced_paths.items():
                if pham_ID in reduced_path:
                    pham_locations = [
                        i for i, pham in enumerate(reduced_path) if pham == pham_ID
                    ]
                    for location in pham_locations:
                        if location < np.min(core_indices[phage]):
                            junction = 0
                        else:
                            junction = (
                                np.max(np.where(location > core_indices[phage])[0]) + 1
                            )
                        location_dict[pham_ID].append(junction)

    # save the dictionary
    dict = {
        str(x): [str(y) for y in list(location_dict[x])]
        for x in list(location_dict.keys())
    }
    with open(args.output, "w") as f:
        json.dump(dict, f, indent=4)


if __name__ == "__main__":
    main()
