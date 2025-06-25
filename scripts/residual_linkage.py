"""
This script plots a figure which shows the continuum of residual linkage across all the groups.

Author: Jemma M. Fendley
"""

import pandas as pd, argparse, json
from matplotlib.transforms import ScaledTranslation
import matplotlib.pyplot as plt, seaborn as sns
import parameters  # preset matplotlib formatting


def new_classify(non_distinct, multiple_haplotypes):
    """
    Returns a description of the group
    """
    if multiple_haplotypes and non_distinct:
        return "small haplotype region, \n not distinct subgroups"
    elif non_distinct:
        return "not distinct subgroups"
    elif multiple_haplotypes:
        return "small haplotype region"
    else:
        return "all remaining groups"


def main():
    parser = argparse.ArgumentParser(description="plot residual linkage of all groups")
    parser.add_argument(
        "-o", "--output", help="figure png file", type=str, required=True
    )
    parser.add_argument(
        "-p", "--output_pdf", help="figure pdf file", type=str, required=True
    )
    parser.add_argument(
        "-s", "--statistics", help="pairwise statistics file", type=str, required=True
    )
    parser.add_argument(
        "-g",
        "--haplotype_groups",
        help="list of groups with multiple haplotypes",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-l",
        "--linkage_information",
        help="linkage statistics file",
        type=str,
        required=True,
    )
    args = parser.parse_args()

    # extract the list of groups with a multiple-haplotype region
    with open(args.haplotype_groups) as f:
        haplotype_information = json.load(f)
    multiple_haplotype_groups = [x["name"] for x in haplotype_information]

    # load the data frames, combine, and add relevant columns
    df_linkage = pd.read_csv(args.linkage_information, sep="\t")
    df_distance = pd.read_csv(args.statistics, sep="\t")
    df = pd.merge(df_linkage, df_distance, on="group")
    df["multiple_haplotypes"] = df["group"].isin(multiple_haplotype_groups)
    df["residual"] = df["asymptote"] - df["background"]
    df["ratio_residual"] = (df["first"] - df["asymptote"]) / (
        df["first"] - df["background"]
    )
    df["percent_residual"] = 100 - 100 * df["ratio_residual"]

    # record the number of groups with low residual linkage
    low_residual_groups = df[(df["ratio_residual"] > 0.8) & (df["residual"] < 0.1)][
        "group"
    ]
    print(
        "N. groups with low residual linkage: {0:0.0f} ({1:0.01%})".format(
            len(low_residual_groups),
            len(low_residual_groups) / len(df),
        )
    )

    # This is a hard-coded list of groups that appear (by eye) to not have distinct subgroup structure
    not_distinct_subgroups = ["A1", "A2", "A6", "A8", "A11", "AM", "AU1", "AW", "AY"]
    not_distinct_subgroups.extend(["AZ1", "B1", "B2", "B3", "BU", "C1", "CA", "CQ1"])
    not_distinct_subgroups.extend(["CS3", "CZ2", "CZ4", "D1", "DC1", "DN1", "DV"])
    not_distinct_subgroups.extend(["E", "EA2", "EF", "F1", "K4", "L1", "M1", "O", "P1"])
    not_distinct_subgroups.extend(["Q", "S"])

    # record the number of groups with no distinct subgroups
    print(
        "N. groups with not distinct subgroups: {0:0.0f} ({1:0.01%})".format(
            len(not_distinct_subgroups),
            len(not_distinct_subgroups) / len(df),
        )
    )
    print(
        "N. groups with low residual linkage that HAVE distinct subgroups: {0:0.0f}".format(
            len([x for x in low_residual_groups if x not in not_distinct_subgroups])
        )
    )

    # add the group classification to the dataframe
    df["non_distinct"] = df["group"].isin(not_distinct_subgroups)
    df["Classification"] = df.apply(
        lambda x: new_classify(x.non_distinct, x.multiple_haplotypes), axis=1
    )

    # initialize the figure
    fig, axes = plt.subplots(
        nrows=1,
        ncols=2,
        figsize=(7, 3.5),
        gridspec_kw={"width_ratios": [1, 1]},
        layout="compressed",
    )

    # plot the continuum of residual linkage
    sns.scatterplot(
        data=df,
        x="residual",
        y="percent_residual",
        hue="Classification",
        hue_order=[
            "small haplotype region, \n not distinct subgroups",
            "small haplotype region",
            "not distinct subgroups",
            "all remaining groups",
        ],
        ax=axes[0],
    )

    # plot the standard deviation of hamming statistics
    sns.scatterplot(
        data=df,
        x="std_hamming_branch_length",
        y="std_hamming",
        hue="Classification",
        ax=axes[1],
        hue_order=[
            "small haplotype region, \n not distinct subgroups",
            "small haplotype region",
            "not distinct subgroups",
            "all remaining groups",
        ],
    )

    # format and save figure
    handles, labels = axes[0].get_legend_handles_labels()
    axes[0].legend(handles=handles, labels=labels)
    axes[1].legend().set_visible(False)
    axes[0].set_xlabel("residual linkage (difference)")
    axes[0].set_ylabel("residual linkage (\% difference)")
    axes[1].set_xlabel("std. dev. of branch lengths")
    axes[1].set_ylabel("std. dev. of p.w. Hamming distances")
    label_list, x_loc = ["a)", "b)"], [-30, -38]
    for i in range(2):
        axes[i].text(
            0.0,
            1.0,
            label_list[i],
            transform=(
                axes[i].transAxes
                + ScaledTranslation(x_loc[i] / 72, +1 / 72, fig.dpi_scale_trans)
            ),
            va="bottom",
        )
    fig.savefig(args.output, dpi=500)
    fig.savefig(args.output_pdf, dpi=500)


if __name__ == "__main__":
    main()
