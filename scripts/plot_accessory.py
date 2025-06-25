"""
This script plots the main accessory localization figure in the main text for input into Inkscape.

Author: Jemma M. Fendley
"""

import json, argparse
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from itertools import combinations
from tqdm import tqdm
from matplotlib.gridspec import GridSpec
from matplotlib.transforms import ScaledTranslation
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import parameters  # preset matplotlib formatting
import matplotlib as mpl
from collections import defaultdict

mpl.rcParams["legend.handletextpad"] = 0.4


def main():
    parser = argparse.ArgumentParser(
        description="create accessory phams localization plots"
    )
    parser.add_argument(
        "-o", "--output", help="output figure png", type=str, required=True
    )
    parser.add_argument(
        "-s", "--output_svg", help="output figure svg", type=str, required=True
    )
    parser.add_argument(
        "-e", "--EE_json", help="group EE JSON", type=str, required=True
    )
    parser.add_argument("-c", "--core", help="core phams JSON", type=str, required=True)
    parser.add_argument(
        "-g",
        "--groups",
        help="list of all groups accessory locations JSONs",
        nargs="+",
        required=True,
    )
    args = parser.parse_args()

    # load the group to core phams JSON
    with open(args.core) as f:
        core_phams = json.load(f)

    # load the group JSON for group EE
    with open(args.EE_json) as f:
        G = json.load(f)

    # extract the lengths of core phams in the correct core pham ordering for group EE
    core_phams_lengths = {
        x["pham_ID"]: x["mean_length"]
        for x in G["phams"]
        if x["pham_ID"] in core_phams["EE"]
    }
    core_phams_order = [
        x["pham_ID"] for x in G["paths"][0]["path"] if x["pham_ID"] in core_phams["EE"]
    ]
    core_lengths = [core_phams_lengths[pham_ID] for pham_ID in core_phams_order]

    # record the cumulative lenghts (or starting positions) to input into figure
    core_starts = [int(np.sum(core_lengths[:i])) for i in range(len(core_lengths) + 1)]

    threshold_frac = 0.75  # extract what fraction of junctions contain this fraction of accessory phams
    np.random.seed(seed=0)

    # initialize the figure formatting
    fig = plt.figure(layout="compressed", figsize=(4.5, 6), dpi=450)
    gs = GridSpec(3, 1, height_ratios=[1.15, 0.35, 2], figure=fig)
    ax1, ax2 = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[2, 0])
    height, width, x_offset, y_offset = 1.2, 2.45, -0.007, -0.01
    ax_ins = inset_axes(
        ax2,
        width=width,
        height=height,
        bbox_to_anchor=(0 + x_offset, 0 - y_offset, 1, 0.68),
        bbox_transform=ax2.transAxes,
    )
    axes, label_list = [ax1, ax2], ["a)", "b)", "c)"]
    for i in range(2):
        axes[i].text(
            0.0,
            1.0,
            label_list[i],
            transform=(
                axes[i].transAxes
                + ScaledTranslation(-35 / 72, +3 / 72, fig.dpi_scale_trans)
            ),
            va="bottom",
        )

    dict_list = []
    for group in tqdm(args.groups):
        row = {"name": group.split("/")[-1].split("_")[0]}
        with open(group) as f:
            initial_dict = json.load(f)  # load the acccessory locations

        # extract the unique list of junction locations for each pham
        dict = {str(x): list(set(y)) for x, y in initial_dict.items()}

        # flatten the dictionary into an array
        all_junctions = [int(x) for y in list(dict.values()) for x in y]

        # save some important numbers for later use
        n_junctions, n_accessory = len(core_phams[row["name"]]) + 1, len(dict)
        threshold = threshold_frac * n_accessory  # threshold as a number of phams
        n_nonempty = len(set(all_junctions))

        # plot panel a) with localization of accessory phams in group EE
        if row["name"] == "EE":
            width = 250  # scale for figure x-axis
            # plot the number of accessory phams in each junction
            counts = [all_junctions.count(x) for x in range(n_junctions)]
            axes[0].bar(core_starts, counts, width=np.ones(n_junctions) * width)

            # remove some of the x-axis labels to save space and format figure
            hide_labels = [9, 5, 7, 9, 13, 16]
            labels = ["" if x in hide_labels else x for x in range(n_junctions)]
            axes[0].set_xticks(core_starts, labels, rotation=90)
            axes[0].set_yticks(np.arange(0, 11, 2), np.arange(0, 11, 2))
            axes[0].set_xlim(left=0 - width, right=core_starts[-1] + width)
            axes[0].set_xlabel("core junction index")
            axes[0].set_ylabel("n. accessory phams")
            axes[0].set_title("group EE")

        # consider only accessory phams that are located in one single junction
        single_locations = [
            int(x) for y in list(dict.values()) for x in y if len(y) == 1
        ]

        # count the nunmber of phams in each junction
        single_counts = [single_locations.count(x) for x in range(n_junctions)]
        n_counted_phams, n_counted_junctions, n_phams_list = 0, 0, [0]  # initialize
        if n_accessory == np.sum(single_counts):  # if all phams single-location

            sorted_single_counts = np.sort(single_counts)[::-1]
            while n_counted_phams < n_accessory:  # stop when all phams counted
                # add a junction and record the number of accessory phams contained
                n_counted_junctions += 1
                n_counted_phams = np.sum(sorted_single_counts[:n_counted_junctions])
                n_phams_list.append(n_counted_phams)

            # calculate number of junctions needed to contain the threshold number of phams
            n_junctions_threshold = next(
                x for x, y in enumerate(n_phams_list) if y >= threshold
            )
        else:
            # extract the phams that are located in multiple junctions
            actual_multi_lists = [x for x in list(dict.values()) if len(x) > 1]

            # extract the junctions that contain these phams
            multi_junctions = list(set([x for y in actual_multi_lists for x in y]))

            # extract the junctions that only contain single-location accessory phams
            single_junctions = [
                x for x in range(n_junctions) if str(x) not in multi_junctions
            ]

            # sort the counts of phams in these junctions
            single_junctions_counts = np.sort(
                [single_counts[x] for x in single_junctions]
            )[::-1]

            multi_dict = defaultdict(list)  # initialize

            # calculate the number of phams that are contained in all possible combinations
            #   of the junctions that contain multi-location phams
            # Note: this is computationally slow
            for i in range(1, len(multi_junctions) + 1):
                for comb in combinations(multi_junctions, i):
                    multi_bonus = len(
                        [x for x in actual_multi_lists if all([y in comb for y in x])]
                    )
                    single = np.sum([single_counts[int(j)] for j in comb])
                    multi_dict[i].append(single + multi_bonus)

            # record the maximum number of phams given i junctions that contain multi-location phams
            max_dict = {x: np.max(y) for x, y in multi_dict.items()}

            while n_counted_phams < n_accessory:
                # add a junction and record the number of accessory phams contained
                n_counted_junctions += 1

                # number of phams contained if only considering junctions with single-location phams
                all_single_count = np.sum(single_junctions_counts[:n_counted_junctions])

                max_multi = []
                for i in range(1, len(multi_junctions) + 1)[::-1]:
                    if n_counted_junctions >= i:
                        multi_contribution = max_dict[i]
                        single_contribution = np.sum(
                            single_junctions_counts[: n_counted_junctions - i]
                        )
                        max_multi.append(multi_contribution + single_contribution)
                        # record the maximum number of phams contained n_counted_junctions, where i of the junctions
                        #   contain phams with multiple locations, and (n_counted_junctions - i) of the junctions only
                        #   contain single-location phams.
                max_multi = np.max(max_multi)

                # record the maximum of the two possible values
                n_counted_phams = max(max_multi, all_single_count)
                n_phams_list.append(n_counted_phams)
            n_junctions_threshold = next(
                x[0] for x in enumerate(n_phams_list) if x[1] >= threshold
            )

        # save the relevant statistics
        row["percent_junctions_threshold"] = n_junctions_threshold / n_junctions
        row["percent_junctions_all"] = n_nonempty / n_junctions
        row["list_percent_junctions"] = np.arange(n_nonempty + 1) * 100 / n_junctions
        row["list_percent_accessory"] = [x * 100 / n_accessory for x in n_phams_list]
        dict_list.append(row)
    df = pd.DataFrame(dict_list)

    # plot the inset of the distribution of number of junctions needed to contain 100%
    #   and then threshold % (currently 75%)
    binstart, binsize, numbers = 0.025, 0.05, [100, threshold_frac * 100]
    colors, label = ["C0", "C1"], ["{0:0.0f}\%".format(x) for x in numbers]
    data_list = [df["percent_junctions_all"], df["percent_junctions_threshold"]]
    means = [x.mean() for x in data_list]
    for i, data in enumerate(data_list):
        ax_ins.hist(
            data,
            alpha=0.5,
            label=label[i],
            bins=np.arange(
                binstart,
                np.max(data) + binsize,
                binsize,
            ),
        )
        ax_ins.axvline(
            means[i],
            color=colors[i],
            linestyle="dashed",
            linewidth=3,
            label="mean",
            alpha=0.75,
        )
        axes[1].axhline(
            y=numbers[i], color=colors[i], linestyle="--", alpha=0.75, linewidth=3
        )

    # in the main figure for each group plot the % of accessory phams contained in x%
    #   of junctions
    for i in range(len(args.groups)):
        axes[1].plot(
            df["list_percent_junctions"][i],
            df["list_percent_accessory"][i],
            color="k",
            alpha=0.1,
        )

    # figure formatting
    ax_ins.set_ylabel("n. groups")
    ax_ins.set_xlabel("core junctions (\%)")
    ax_ins.set_yticks([0, 10, 20, 30], [0, 10, 20, 30])
    ax_ins.legend()
    ax_ins.set_xticks(
        [0, 0.25, means[1], 0.5, means[0], 0.75, 1],
        [0, 25, round(means[1], 3) * 100, 50, round(means[0], 3) * 100, 75, 100],
    )
    ax_ins.set_ylim([0, 37])
    tick_array = [0, 25, 50, 75, 100]
    axes[1].set_xticks(tick_array, tick_array)
    axes[1].set_yticks(tick_array, tick_array)
    axes[1].set_xlabel("core junctions (\%)")
    axes[1].set_ylabel("cumulative accessory phams (\%)")

    # save figure
    fig.savefig(args.output, dpi=450)
    fig.savefig(args.output_svg, dpi=450)


if __name__ == "__main__":
    main()
