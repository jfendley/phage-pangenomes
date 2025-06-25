"""
This script creates a dictionary of groups to the list of their core phams.

Author: Jemma M. Fendley
"""

import argparse, json


def main():
    parser = argparse.ArgumentParser(description="creates group to core phams JSON")
    parser.add_argument(
        "-o", "--output", required=True, type=str, help="output file name"
    )
    parser.add_argument(
        "-g", "--groups", nargs="+", required=True, help="list of group JSON files"
    )

    args = parser.parse_args()

    core_dict = {}  # initialize dictionary

    for group in args.groups:
        with open(group) as f:
            G = json.load(f)

        # list of core phams, present exactly once in every phage in the group
        core_dict[G["name"]] = [
            x["pham_ID"]
            for x in G["phams"]
            if x["phage_count"] == len(G["paths"])
            and x["total_count"] == x["phage_count"]
        ]

    # save dictionary
    with open(args.output, "w") as f:
        json.dump(core_dict, f, indent=4)


if __name__ == "__main__":
    main()
