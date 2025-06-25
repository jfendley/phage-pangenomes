"""
This script counts the number of duplicate phams: phams that appear
    more than once in the same phage.

Author: Jemma M. Fendley
"""

import argparse, json
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="counts the duplicate phams")
    parser.add_argument(
        "-o", "--output", help="output TSV file", type=str, required=True
    )
    parser.add_argument(
        "-g", "--groups", nargs="+", help="group JSON files", required=True
    )
    args = parser.parse_args()

    all_phams, all_duplicate_phams, all_duplicate_core = set(), set(), set()
    dict_list = []  # initialize
    for group in args.groups:
        with open(group) as f:
            G = json.load(f)

        # record some basic statistics
        row = {
            "name": G["name"],
            "n_phages": len(G["paths"]),
            "n_phams": len(G["phams"]),
        }

        # add to the list of all phams
        all_phams.update([x["pham_ID"] for x in G["phams"]])

        # phams that appear multiple times in at least one phage
        duplicate_phams_list = [
            x for x in G["phams"] if x["total_count"] > x["phage_count"]
        ]
        duplicate_phams = [x["pham_ID"] for x in duplicate_phams_list]
        row["n_duplicate"] = len(duplicate_phams)
        all_duplicate_phams.update(duplicate_phams)

        # potential core are phams that appear in every phage in the group but not exactly once
        potential_core = [
            x["pham_ID"]
            for x in duplicate_phams_list
            if x["phage_count"] == len(G["paths"])
        ]
        row["n_duplicate_core"] = len(potential_core)
        all_duplicate_core.update(potential_core)

        # count which phages have phams that appear more than once
        row["n_phages_with_duplicates"] = len(
            [
                x
                for x in G["paths"]
                if any(
                    [y in [i["pham_ID"] for i in x["path"]] for y in duplicate_phams]
                )
            ]
        )
        dict_list.append(row)

    df = pd.DataFrame(dict_list)
    df = df.sort_values(by=["n_duplicate"], ascending=False)

    # add summary statistics for all of the groups
    df.loc[-1] = {
        "name": "all_groups",
        "n_phages": df["n_phages"].sum(),
        "n_phams": len(all_phams),
        "n_duplicate": len(all_duplicate_phams),
        "n_duplicate_core": len(all_duplicate_core),
        "n_phages_with_duplicates": df["n_phages_with_duplicates"].sum(),
    }

    # calculate some statistics
    df["percent_duplicate_phams"] = (100 * df["n_duplicate"] / df["n_phams"]).round(2)
    df["percent_duplicate_phages"] = (
        100 * df["n_phages_with_duplicates"] / df["n_phages"]
    ).round(2)

    # save the dataframe
    df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
