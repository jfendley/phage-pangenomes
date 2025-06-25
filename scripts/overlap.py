"""
This script investigates core genes that "overlap", i.e. share DNA sequence. The output is
    a figure that is in the SI.

Author: Jemma M. Fendley
"""

import argparse, json
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns, matplotlib as mpl
import parameters  # preset matplotlib formatting
from matplotlib.transforms import ScaledTranslation

mpl.rcParams["legend.handletextpad"] = 0.4


def second_pham(string):
    """
    This function analyzes the output of core_position_pham.py. If a position is
        located in one pham, return that pham. If a position is located in two phams,
        return the second (greater index) pham.

    This is used simply to count the number of core phams per group.
    """
    if "_" in str(string):
        return np.max([int(x) for x in string.split("_")])
    else:
        return int(string)


def main():
    parser = argparse.ArgumentParser(description="create PA basic LD plot")
    parser.add_argument(
        "-o", "--output", help="figure png file", type=str, required=True
    )
    parser.add_argument(
        "-p", "--output_pdf", help="figure pdf file", type=str, required=True
    )
    parser.add_argument(
        "-c", "--cyclic", help="list of cyclic groups", nargs="+", required=True
    )
    parser.add_argument(
        "-i",
        "--core_positions",
        help="core position to pham files",
        nargs="+",
        required=True,
    )
    parser.add_argument(
        "-a",
        "--accessory",
        help="accessory pham junction locations file",
        nargs="+",
        required=True,
    )

    args = parser.parse_args()

    df_list = []
    for core_file in args.core_positions:
        # record the group name and find the corresponding accessory locations file
        row = {"name": core_file.split("/")[-1].split("_")[0]}
        accessory_file = [x for x in args.accessory if "/" + row["name"] + "_" in x]
        assert len(accessory_file) == 1, "Error with input files"
        with open(accessory_file[0]) as f:
            accessory_locations = json.load(f)

        # load the dataframe with core position to pham information
        core_df = pd.read_csv(core_file, sep="\t")

        # record the number of core phams, and the number that are completely contained in another pham
        core_df["second_pham"] = core_df["pham"].apply(second_pham)
        n_core = np.max(core_df["second_pham"]) + 1
        n_complete_overlap = n_core - core_df[~core_df["overlap"]]["pham"].nunique()

        # find the list of empty junctions
        nonempty_junctions = set(
            [int(x) for y in list(accessory_locations.values()) for x in y]
        )
        empty_junctions = set(range(n_core + 1)) - nonempty_junctions

        # find the list of overlapping junctions
        overlap = core_df[core_df["overlap"]]["pham"]
        overlap_junctions = set([int(y.split("_")[1]) for y in overlap])

        # record the percentage of core that are involved in overlap (not currently in figure)
        core_overlap = set([int(x) for y in overlap for x in y.split("_")])
        row["percent_core_overlap"] = len(core_overlap) / n_core * 100
        row["percent_core_contained"] = n_complete_overlap / n_core * 100

        if row["name"] not in args.cyclic:
            # the junctions on the edge of the pangenome can not have overlapping core phams,
            #   since they only have a pham on one side. There are n_core - 1 non-edge junctions.
            row["percent_junction_overlap"] = (
                len(overlap_junctions) / (n_core - 1) * 100
            )
            empty_non_edge_junctions = empty_junctions - set([0, n_core])

            # record the junctions that are empty and have overlapping flanking core genes
            overlap_empty_junctions = overlap_junctions.intersection(empty_junctions)
        else:
            row["percent_junction_overlap"] = len(overlap_junctions) / n_core * 100
            # if cyclic then junction 0 = junction n_core. Hence if they are both empty,
            #   it only needs to be considered once
            if 0 in empty_junctions and n_core in empty_junctions:
                empty_non_edge_junctions = empty_junctions - set([n_core])

            # If at least one of 0 or n_core is not empty, then we can consider neither to be empty
            #   and can remove from both lists.
            else:
                empty_non_edge_junctions = empty_junctions - set([0, n_core])
            overlap_empty_junctions = overlap_junctions.intersection(
                empty_non_edge_junctions
            )

        # percentage of empty junctions that do not have any overlapping flanking core genes
        row["percent_empty_not_overlap"] = (
            (len(empty_non_edge_junctions) - len(overlap_empty_junctions))
            / len(empty_non_edge_junctions)
            * 100
        )

        # percentage of junctions with overlapping flanking core genes that are not empty
        #   This is possible since in one phage there might be overlapping core gene, while
        #   in another phage, there is an accessory pham
        if len(overlap_junctions) > 0:
            row["percent_overlap_not_empty"] = (
                (len(overlap_junctions) - len(overlap_empty_junctions))
                / len(overlap_junctions)
                * 100
            )
        else:
            row["percent_overlap_not_empty"] = 0
        df_list.append(row)

    df = pd.DataFrame(df_list)

    # plot the figure
    options = ["junction_overlap", "empty_not_overlap", "overlap_not_empty"]
    column_list = ["percent_" + x for x in options]
    fig, axes = plt.subplots(
        nrows=1,
        ncols=3,
        figsize=(7.6, 3),
        gridspec_kw={"width_ratios": [1, 1, 1.025]},
        layout="constrained",
        dpi=450,
    )
    list_of_labels = ["a)", "b)", "c)"]
    for i, column_name in enumerate(column_list):
        sns.histplot(
            data=df,
            x=column_name,
            ax=axes[i],
            bins=np.arange(-2.5, np.max(df[column_name]) + 5, 5),
            label="data",
        )
        axes[i].axvline(
            np.mean(df[column_name]),
            color="k",
            linestyle="--",
            linewidth=3,
            label="mean",
        )
        if i != 0:
            axes[i].set_ylabel("")
        # the yticks are hard-coded, this may need modification if the data changes
        if i == 2:
            new_ticks = [0, 5, 10, 15]
            new_labels = [0, 5, 10, 15]
        else:
            new_ticks = [0, 5, 10, 15, 20]
            new_labels = [0, 5, 10, 15, 20]
        axes[i].set_yticks(
            ticks=new_ticks,
            labels=new_labels,
            rotation=90,
            ha="center",
            rotation_mode="anchor",
            va="baseline",
        )
        axes[i].text(
            0.0,
            1.0,
            list_of_labels[i],
            transform=(
                axes[i].transAxes
                + ScaledTranslation(-15 / 72, 5 / 72, fig.dpi_scale_trans)
            ),
            va="bottom",
        )

    # add labels and save the figure
    axes[0].set_xlabel("junctions with overlap (\%)")
    axes[0].legend(loc=2)
    axes[1].set_xlabel("non-overlap empty junctions (\%)")
    axes[2].set_xlabel("non-empty overlap junctions (\%)")
    axes[0].set_ylabel("n. groups")
    fig.savefig(args.output, dpi=450)
    fig.savefig(args.output_pdf, dpi=450)


if __name__ == "__main__":
    main()
