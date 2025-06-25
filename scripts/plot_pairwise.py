"""
This script plots summary statistics for the pairwise snp analysis for all groups.

Author: Jemma M. Fendley
"""

import numpy as np, pandas as pd
import matplotlib.pyplot as plt, seaborn as sns
import parameters  # preset matplotlib formatting
import argparse, json
from utils import convert_to_list  # converts string list of ints to list of ints
import matplotlib as mpl

mpl.rcParams["legend.handletextpad"] = 0.4


def main():
    parser = argparse.ArgumentParser(
        description="plot summary statistics for pairwise snp analyses"
    )
    parser.add_argument(
        "-s", "--SNP_info", help="name of SNP info file", type=str, required=True
    )
    parser.add_argument(
        "-g", "--groups", help="all group JSON files", nargs="+", required=True
    )
    parser.add_argument("-o", "--output", help="output file", type=str, required=True)
    parser.add_argument("-p", "--pdf", help="output pdfd file", type=str, required=True)

    args = parser.parse_args()

    # create a dictionary of group to mean core pham length
    core_length_dict = {}
    for group in args.groups:
        with open(group) as f:
            G = json.load(f)
        core_pham_lengths = [
            x["mean_length"]
            for x in G["phams"]
            if x["total_count"] == len(G["paths"])
            and x["phage_count"] == len(G["paths"])
        ]
        core_length_dict[G["name"]] = np.mean(core_pham_lengths)

    # load the dataframe
    df = pd.read_csv(args.SNP_info, sep="\t")

    # find the maximum estimated recombination length per group
    df["distance_list"] = df["distances"].apply(convert_to_list)
    df["max_distance"] = df["distance_list"].apply(np.max)
    data = df.groupby("group")["max_distance"].max().reset_index()
    print("Mean: ", np.mean(data["max_distance"]))

    # add the mean core pham lengths to the dictionary
    data["mean_core_pham_length"] = data["group"].map(core_length_dict)

    # plot the maximum estimated recombination length v.s. mean core pham length
    fig, axes = plt.subplots(1, 1, layout="constrained", figsize=(5, 4))

    # plot the lines y=x and y=2x for comparison
    xlims = (
        np.min(data["mean_core_pham_length"]) - 10,
        np.max(data["mean_core_pham_length"]) + 10,
    )
    xlims2 = (
        2 * (np.min(data["mean_core_pham_length"]) - 10),
        2 * (np.max(data["mean_core_pham_length"]) + 10),
    )

    # format figure and save
    axes.set_xlim([xlims[0], xlims[1]])
    axes.plot(xlims, xlims, label=r"$y=x$", color="C1")
    axes.plot(xlims, xlims2, label=r"$y=2x$", color="C2")

    sns.scatterplot(
        data=data,
        x="mean_core_pham_length",
        y="max_distance",
        label="{0:0.0f} groups".format(len(data)),
        ax=axes,
        alpha=0.8,
        color="C0",
    )
    axes.set_xlabel("mean core pham length")
    axes.set_ylabel("maximum estimated recombination length")
    axes.legend()
    fig.savefig(args.output, dpi=450)
    fig.savefig(args.pdf, dpi=450)


if __name__ == "__main__":
    main()
