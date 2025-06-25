"""
This script calculates the mean of two different pairwise metrics for all of the groups,
    used for a paper figure.

Author: Jemma M. Fendley
"""

import argparse
import pandas as pd
import numpy as np


def get_shorter_coverage(length_1, length_2, coverage_1, coverage_2):
    """
    Returns the coverage of the shorter genome, given two genome lengths and two
        respective coverages.
    """
    if length_1 < length_2:
        return coverage_1
    elif length_1 > length_2:
        return coverage_2
    else:
        return max(coverage_1, coverage_2)


def main():
    parser = argparse.ArgumentParser(
        description="calculates the mean of pairwise metrics for all groups"
    )
    parser.add_argument(
        "-o", "--output", help="name of output TSV file", type=str, required=True
    )
    parser.add_argument(
        "-g",
        "--groups",
        nargs="+",
        help="list of group TSV files with pairwise metrics",
        required=True,
    )
    args = parser.parse_args()

    dict_list = []  # initialize
    for group in args.groups:
        # extract the group name from the file name
        group_name = group.split("/")[-1].split("_")[0]

        # load the dataframe
        df_group = pd.read_csv(group, sep="\t")

        # calculate the coverage of the shorter genome
        df_group["min_coverage"] = df_group.apply(
            lambda x: get_shorter_coverage(
                x.length_1, x.length_2, x.percent_shared_1, x.percent_shared_2
            ),
            axis=1,
        )

        # record the mean coverage and hamming distance (to be used for ANI downstream)
        dict_list.append(
            {
                "name": group_name,
                "mean_percent_pairwise_coverage": np.mean(df_group["min_coverage"]),
                "mean_hamming": np.mean(df_group["Hamming"]),
                "mean_jaccard": np.mean(df_group["Jaccard"]),
                "std_hamming": np.std(
                    df_group["Hamming"], ddof=1
                ),  # sample standard deviation
                "std_jaccard": np.std(df_group["Jaccard"], ddof=1),
            }
        )

    # save the dataframe
    df = pd.DataFrame(dict_list)
    df.to_csv(args.output, index=False, sep="\t")


if __name__ == "__main__":
    main()
