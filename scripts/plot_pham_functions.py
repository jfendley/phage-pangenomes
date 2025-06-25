"""
This script creates a plot that shows the functions that are enriched in/near hotspots.

Author: Jemma M. Fendley
"""

import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import parameters  # preset matplotlib formatting
import matplotlib.colors as mcolors
import matplotlib as mpl, seaborn as sns
from matplotlib.transforms import ScaledTranslation
import argparse

mpl.rcParams["ytick.labelsize"] = 8


def make_label(x, y):
    """
    Returns a string of x and y used as labels in the figure
    """
    return str(x) + "/" + str(y)


def main():
    parser = argparse.ArgumentParser(
        description="creates a plot that shows pham functions enriched in hotspots"
    )
    parser.add_argument(
        "-i", "--input", help="pham functions TSV", type=str, required=True
    )
    parser.add_argument(
        "-o", "--output", help="output file png", type=str, required=True
    )
    parser.add_argument("-p", "--pdf", help="output file pdf", type=str, required=True)
    args = parser.parse_args()

    # load the dataframe
    df = pd.read_csv(args.input, sep="\t")
    n_phams = df["pham_ID"].nunique()

    # find the phams with more than one annotation across different groups, and
    #   by default, we will consider these to be annotated.
    count_functions = df.groupby("pham_ID")["majority_function"].nunique()
    different_function_phams = count_functions.index[count_functions > 1]

    # count the number of phams not annotated in any group
    n_not_annotated = (
        df[~df["pham_ID"].isin(different_function_phams)]
        .drop_duplicates(subset=["pham_ID"])["majority_function"]
        .value_counts()["no majority annotation"]
    )
    print("Fraction annotated: ", 1 - n_not_annotated / n_phams)

    # record all possible functions, and the number that are core and accessory
    all_functions = [
        x for x in df["majority_function"].unique() if x != "no majority annotation"
    ]
    n_list = [df["core"].sum(), (~df["core"]).sum()]

    # record the number that are core and flanking a hotspot, and the number that are
    #   accessory and in a hotspot
    n_hotspot_list = [
        (df["core"] & df["hotspot"]).sum(),  # core
        (~df["core"] & df["hotspot"]).sum(),  # accessory
    ]

    dict_list = []
    # for each function, record its statistics
    for function in all_functions:
        this_df = df[df["majority_function"] == function]
        row = {
            "function": function,
            "n": len(this_df),
            "n_core": this_df["core"].sum(),
            "n_accessory": (~this_df["core"]).sum(),
            "n_hotspot": this_df["hotspot"].sum(),
            "n_core_hotspot": (this_df["core"] & this_df["hotspot"]).sum(),
            "n_accessory_hotspot": (~this_df["core"] & this_df["hotspot"]).sum(),
        }
        dict_list.append(row)
    function_df = pd.DataFrame(dict_list)

    # initalize the figure
    list_of_labels = ["a)", "b)"]
    fig, axes = plt.subplots(
        nrows=1,
        ncols=2,
        figsize=(6.5, 7.5),
        gridspec_kw={"width_ratios": [1, 1]},
        layout="compressed",
    )
    # plot the core and then the accessory
    for i, type in enumerate(["core", "accessory"]):
        # restrict analysis to phams that appear at least 12 times (for presentation)
        #   and those that exist at least twice in/next to a hotspot
        plot_df = function_df[
            (function_df["n_" + type] > 11)
            & (function_df["n_" + type + "_hotspot"] > 1)
        ]
        if i == 0:
            plot_df["label"] = plot_df.apply(
                lambda x: make_label(x.n_core_hotspot, x.n_core), axis=1
            )
        elif i == 1:
            plot_df["label"] = plot_df.apply(
                lambda x: make_label(x.n_accessory_hotspot, x.n_accessory), axis=1
            )
        plot_df["fraction_" + type] = plot_df["n_" + type] / n_list[i]
        plot_df["fraction_" + type + "_hotspot"] = (
            plot_df["n_" + type + "_hotspot"] / n_hotspot_list[i]
        )

        # core/accessory enrichment, given the name "core"/"accessory" for plotting purposes
        plot_df[type] = (
            plot_df["fraction_" + type + "_hotspot"] / plot_df["fraction_" + type] - 1
        )

        # sort by decreasing enrichment
        plot_data = plot_df[["function", type]].set_index("function")
        plot_data_sorted = plot_data.sort_values(by=type, ascending=False)
        plot_labels = plot_df[["function", "label"]].set_index("function")
        plot_labels_sorted = plot_labels.reindex(plot_data_sorted.index)

        # calculate the maximum possible enrichment and set colorbar scale accordingly
        max_stat = n_list[i] / n_hotspot_list[i] - 1
        offset = mcolors.TwoSlopeNorm(vmin=-1, vcenter=0, vmax=max_stat)

        sns.heatmap(
            plot_data_sorted,
            annot=plot_labels_sorted,
            fmt="",
            cmap="RdBu",
            norm=offset,
            cbar_kws=dict(
                use_gridspec=False, label="(O-E)/E", location="bottom", pad=0.02
            ),
            ax=axes[i],
            annot_kws={"fontsize": 11},
        )

        # add panel label
        axes[i].text(
            0.0,
            1.0,
            list_of_labels[i],
            transform=(
                axes[i].transAxes
                + ScaledTranslation(-190 / 72, -12 / 72, fig.dpi_scale_trans)
            ),
            va="bottom",
        )

    # save the figure
    fig.savefig(args.output, dpi=450)
    fig.savefig(args.pdf, dpi=450)


if __name__ == "__main__":
    main()
