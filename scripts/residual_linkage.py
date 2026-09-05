"""
This script plots a figure which shows the continuum of residual linkage across all the groups for the SI.

Author: Jemma M. Fendley
"""

import pandas as pd, argparse
from matplotlib.transforms import ScaledTranslation
import matplotlib.pyplot as plt, seaborn as sns
import parameters  # preset matplotlib formatting


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
        "-l",
        "--linkage_information",
        help="linkage statistics file",
        type=str,
        required=True,
    )
    args = parser.parse_args()

    # load the data frames, combine, and add relevant columns
    df_linkage = pd.read_csv(args.linkage_information, sep="\t")
    df_distance = pd.read_csv(args.statistics, sep="\t")
    df = pd.merge(df_linkage, df_distance, on="group")
    df["residual linkage"] = df["asymptote"] - df["background"]

    # record the number of groups with low residual linkage
    low_residual_groups = df[df["residual linkage"] < 0.1]["group"]
    print(
        "N. groups with low residual linkage: {0:0.0f} ({1:0.01%})".format(
            len(low_residual_groups),
            len(low_residual_groups) / len(df),
        )
    )
    print("Low residual linkage gruops: ", list(low_residual_groups))

    # initialize the figure
    fig, axes = plt.subplots(
        nrows=1,
        ncols=2,
        figsize=(6.3, 3.15),
        gridspec_kw={"width_ratios": [1, 1]},
        layout="compressed",
    )

    # plot the continuum of residual linkage
    sns.scatterplot(
        data=df,
        x="asymptote",
        y="residual linkage",
        ax=axes[0],
        label="{0:0.0f} total groups".format(len(df)),
    )
    axes[0].set_ylabel("residual: asymptote - background")

    # plot the standard deviation of hamming distances and residual linkage
    sns.scatterplot(
        data=df,
        x="std_hamming",
        y="residual linkage",
        ax=axes[1],
    )

    # plot the level of detectability
    axes[1].axhline(
        0.1,
        color="C2",
        linestyle="--",
        label="level of detectability,\nlow residual linkage",
    )
    axes[0].axhline(
        0.1,
        color="C2",
        linestyle="--",
        label="level of detectability,\nlow residual linkage",
    )

    # format and save figure
    axes[1].legend().set_visible(False)
    axes[0].legend()
    axes[1].set_xlabel("std. dev. of p.w. Hamming distances")

    label_list, x_loc = ["a)", "b)", "c)"], [-32, -28, -32]
    y_loc = [-2, -2, -5]
    for i in range(2):
        axes[i].text(
            0.0,
            1.0,
            label_list[i],
            transform=(
                axes[i].transAxes
                + ScaledTranslation(x_loc[i] / 72, y_loc[i] / 72, fig.dpi_scale_trans)
            ),
            va="bottom",
        )
    fig.savefig(args.output, dpi=500)
    fig.savefig(args.output_pdf, dpi=500)


if __name__ == "__main__":
    main()
