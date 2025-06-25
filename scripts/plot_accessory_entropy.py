"""
This script plots the accessory localization figures in the SI that show
    entropy and participation ratio.

Author: Jemma M. Fendley
"""

import json, argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.transforms import ScaledTranslation
import parameters  # preset matplotlib formating
import matplotlib as mpl
import seaborn as sns

mpl.rcParams["legend.handletextpad"] = 0.4


def max_stats(N, k):
    """
    calculates the maximum entropy and participation ratio given N accessory phams and k hotspots

    N = number of accessory phams (integer)
    k = number of hostpots (integer)
    """
    if N > k:
        quotient, remainder = N // k, N % k
        max_entropy = (
            -(k - remainder) * (quotient / N) * np.log(quotient / N)
            - remainder * ((quotient + 1) / N) * np.log((quotient + 1) / N)
        ) / np.log(k)
        max_pr = -np.log(
            (
                (k - remainder) * (quotient / N) ** 2
                + remainder * ((quotient + 1) / N) ** 2
            )
        ) / np.log(k)

    # if there are fewer accessory phams than junctions, the maximum entropy is when each
    #   accessory pham is in a different junction
    if N <= k:
        max_entropy = np.log(N) / np.log(k)
        max_pr = np.log(N) / np.log(k)

    return max_entropy, max_pr


def main():
    parser = argparse.ArgumentParser(
        description="create SI accessory phams localization plots"
    )
    parser.add_argument(
        "-o", "--output", help="output figure png", type=str, required=True
    )
    parser.add_argument(
        "-p", "--output_pdf", help="output figure pdf", type=str, required=True
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

    np.random.seed(seed=0)  # for reproducibility
    dict_list = []
    for group in args.groups:
        row = {"name": group.split("/")[-1].split("_")[0]}
        n_junctions = len(core_phams[row["name"]]) + 1

        with open(group) as f:
            initial_dict = json.load(f)  # load the acccessory locations

        # extract the unique list of junction locations for each pham
        dict = {str(x): list(set(y)) for x, y in initial_dict.items()}

        # consider only accessory phams that are located in one single junction
        single_locations = [
            int(x) for y in list(dict.values()) for x in y if len(y) == 1
        ]

        # count the nunmber of phams in each junction
        single_counts = [single_locations.count(x) for x in range(n_junctions)]
        n_single_phams = np.sum(single_counts)

        # calculate the entropy and participation ratio
        row["entropy"] = np.sum(
            [
                -(i / n_single_phams) * np.log(i / n_single_phams)
                for i in single_counts
                if i > 0
            ]
        ) / np.log(n_junctions)
        row["pr"] = np.sum([(i / n_single_phams) ** 2 for i in single_counts])
        row["norm_pr"] = -np.log(row["pr"]) / np.log(n_junctions)

        # calculate the maximum entropy and participation ratio possible
        row["max_entropy"], row["max_pr"] = max_stats(n_single_phams, n_junctions)

        # create a random equally-distributed sample for comparison
        random_options = np.random.choice(
            n_junctions, size=n_single_phams, replace=True
        )
        random_counts = [(random_options == x).sum() for x in range(n_junctions)]

        # calculate the entropy and participation ratio for this random sample
        row["random_entropy"] = np.sum(
            [
                -(i / n_single_phams) * np.log(i / n_single_phams)
                for i in random_counts
                if i > 0
            ]
        ) / np.log(n_junctions)
        row["random_pr"] = -np.log(
            np.sum([(i / n_single_phams) ** 2 for i in random_counts])
        ) / np.log(n_junctions)

        dict_list.append(row)
    df = pd.DataFrame(dict_list)

    # initialize the figure
    fig, axes = plt.subplots(
        nrows=1,
        ncols=3,
        figsize=(7.8, 3),
        gridspec_kw={"width_ratios": [1, 1, 1]},
        layout="constrained",
    )

    legend_list = ["data", "random", "maximum"]

    # calculate bin sizes for particpation ratio plot
    pr_list = ["norm_pr", "random_pr", "max_pr"]
    min_pr, max_pr = np.min(df[pr_list]), np.max(df[pr_list])
    pr_binsize = (max_pr - min_pr) / 1000
    pr_bins = np.arange(min_pr, max_pr + pr_binsize, pr_binsize)

    # plot the normalized participation ratios
    for i, x in enumerate(pr_list):
        hist = sns.histplot(
            ax=axes[2],
            data=df,
            x=x,
            bins=pr_bins,
            stat="density",
            cumulative=True,
            element="step",
            fill=False,
            label=legend_list[i],
        )
        hist.set(ylabel=None)
    axes[2].legend()
    axes[2].set_xlabel("normalized (log inverse) p.r.")

    # plot the non-normalized participation ratio, along with fractions for comparison
    hist = sns.histplot(
        ax=axes[1],
        data=df,
        x="pr",
        bins=3000,
        stat="density",
        cumulative=True,
        element="step",
        fill=False,
        label=legend_list[0],
    )
    hist.set(ylabel=None)
    axes[1].axvline(
        np.mean(df["pr"]),
        color="C1",
        linewidth=1.5,
        label="data mean: {0:0.2f}".format(np.mean(df["pr"])),
    )
    axes[1].axvline(1 / 5, color="C2", linestyle="dashed", linewidth=1, label="1/5")
    axes[1].axvline(1 / 6, color="C3", linestyle="dashed", linewidth=1, label="1/6")
    axes[1].axvline(1 / 7, color="C4", linestyle="dashed", linewidth=1, label="1/7")
    axes[1].legend()
    axes[1].set_xlabel("participation ratio (p.r.)")

    # plot the entropies
    entropy_list = ["entropy", "random_entropy", "max_entropy"]
    max_entropy, min_entropyy = np.max(df[entropy_list]), np.min(df[entropy_list])
    entropy_binsize = (max_entropy - min_entropyy) / 1000
    bins = np.arange(min_entropyy, max_entropy + entropy_binsize, entropy_binsize)
    for i, x in enumerate(entropy_list):
        sns.histplot(
            ax=axes[0],
            data=df,
            x=x,
            bins=bins,
            stat="density",
            cumulative=True,
            element="step",
            fill=False,
            label=legend_list[i],
        )
    axes[0].legend()
    axes[0].set_ylabel("cumulative n. groups (density)")
    axes[0].set_xlabel("normalized entropy")

    # add panel labels
    label_list = ["a)", "b)", "c)"]
    for i in range(3):
        axes[i].text(
            0.0,
            1.0,
            label_list[i],
            transform=(
                axes[i].transAxes
                + ScaledTranslation(-20 / 72, +5 / 72, fig.dpi_scale_trans)
            ),
            va="bottom",
        )

    # simplify x axis labels
    for i in [0, 2]:
        axes[i].set_xticks(
            np.arange(0.25, 1.05, 0.25), ["0.25", "0.50", "0.75", "1.00"]
        )

    # save the figure
    fig.savefig(args.output, dpi=450)
    fig.savefig(args.output_pdf, dpi=450)


if __name__ == "__main__":
    main()
