"""
This script investigates the synteny of the accessory phams.

Author: Jemma M. Fendley
"""

import argparse, json
import numpy as np, pandas as pd
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(
        description="investigate the synteny within the accessory phams"
    )
    parser.add_argument(
        "-o", "--output", help="file for synteny dataframe", type=str, required=True
    )
    parser.add_argument(
        "-g", "--group_json", help="group JSON file", type=str, required=True
    )
    parser.add_argument(
        "-a",
        "--accessory",
        help="accessory pham locations JSON file",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-c", "--core", help="list of core phams", nargs="+", required=True
    )
    args = parser.parse_args()

    with open(args.accessory) as f:
        initial_dict = json.load(f)  # load the acccessory locations

    # Count the number of different junctions each pham is located in
    # Phams that only exist in one phage trivially only exist in one junction,
    #   so those are excluded from the statistic.
    n_junctions = [len(set(y)) for y in initial_dict.values() if len(y) > 1]

    # record the fraction of accessory phams that exist in only 1 junction
    fraction_single_junction = n_junctions.count(1) / len(n_junctions)

    with open(args.group_json) as f:
        G = json.load(f)  # load the group JSON

    paths, n_phages = G["paths"], len(G["paths"])
    pham_df = pd.DataFrame(G["phams"]).set_index("pham_ID")

    # we disregard accessory phams that appear multiple times in the same phage,
    #   in order for synteny to be well-defined.
    accessory_phams = pham_df[
        (pham_df["phage_count"] < n_phages)
        & (pham_df["total_count"] == pham_df["phage_count"])
    ].index.tolist()

    # record the indices of the location of core phams for each phage
    all_core_locations = {
        phage["ID"]: np.array(
            [i for i, x in enumerate(phage["path"]) if x["pham_ID"] in args.core]
        )
        for phage in paths
    }

    # disregard the strandedness and focus solely on pham order
    simple_paths = {
        phage["ID"]: [x["pham_ID"] for x in phage["path"]] for phage in paths
    }

    # This initializes a matrix containing strings of an appropriate length
    synteny_matrix = np.full((n_phages, len(accessory_phams)), "__________")

    # initialize: indices of all the phages of the pham in a specific junction
    pham_phages = defaultdict(lambda: defaultdict(set))
    pham_junctions = defaultdict(list)  # junctions per accessory pham
    junction_phams = defaultdict(set)  # accessory phams per junction
    nonempty_junctions = set()
    for j, pham in enumerate(accessory_phams):
        for i, (phage_ID, core_locations) in enumerate(all_core_locations.items()):
            if pham in simple_paths[phage_ID]:
                pham_location = simple_paths[phage_ID].index(pham)
                # find the junction that the accessory pham is located in and its
                #   relative location within that junction
                if pham_location < np.min(core_locations):
                    junction = 0
                    relative_location = pham_location
                else:
                    junction = np.max(np.where(pham_location > core_locations)[0]) + 1
                    relative_location = pham_location - core_locations[junction - 1]
                # record the information
                synteny_matrix[i, j] = str(junction) + "_" + str(relative_location)
                pham_junctions[pham].append(junction)
                nonempty_junctions.add(junction)
                junction_phams[junction].add(j)
                pham_phages[j][junction].add(i)
            else:
                synteny_matrix[i, j] = "-"

    # quick check that this matches the input file
    assert pham_junctions == {
        x: [int(i) for i in y] for x, y in initial_dict.items()
    }, "Error in input accessory locations JSON file"

    dict_list = []
    for i in nonempty_junctions:  # iterate through the junctions
        phams = list(junction_phams[i])  # phams in the junction
        if len(phams) > 1:
            for j, pham1 in enumerate(phams[:-1]):
                for pham2 in phams[j + 1 :]:
                    # find the phages that contain both pham1 and pham2 in junction i
                    phages = pham_phages[pham1][i].intersection(pham_phages[pham2][i])
                    if len(phages) > 1:
                        row = {
                            "junction": i,
                            "pham1": accessory_phams[pham1],
                            "pham2": accessory_phams[pham2],
                            "n_phages": len(phages),
                        }
                        pair_order = []
                        for k in phages:
                            position1 = synteny_matrix[k, pham1].split("_")[1]
                            position2 = synteny_matrix[k, pham2].split("_")[1]
                            if (
                                position1 > position2
                            ):  # pham1 comes after pham2 in phage k
                                pair_order.append(1)
                            elif (
                                position1 < position2
                            ):  # pham1 comes before pham2 in phage k
                                pair_order.append(0)
                        # find the number of phages that do not follow the consensus pairwise ordering
                        counts = [pair_order.count(x) for x in [0, 1]]
                        row["n_nonsyntenic"] = np.min(counts)
                        dict_list.append(row)

    # save all of the information
    df = pd.DataFrame(dict_list)
    df_nonsyntenic = df[df["n_nonsyntenic"] > 0]
    save_dict = [
        {
            "group": G["name"],
            "n_nonempty_junctions": len(nonempty_junctions),
            "n_nonsyntenic_junctions": df_nonsyntenic["junction"].nunique(),
            "n_pham_pairs": len(df),
            "n_nonsyntenic_pairs": len(df_nonsyntenic),
            "fraction_single_junction": fraction_single_junction,
        }
    ]
    pd.DataFrame(save_dict).to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
