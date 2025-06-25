"""
This script plots a figure that shows the effect of processing the multiple-haplotype regions on
    linkage disequilibrium. This script is tailored for group E, but could be modified slightly
    for any of the groups with a multiple-haplotype region.

Author: Jemma M. Fendley
"""

# %%
import pandas as pd, numpy as np
import argparse, json
import matplotlib.pyplot as plt
import parameters  # preset matplotlib formatting
from utils import rolling_average  # calculates the rolling average


# %%
def plot_processed(output, output_pdf, linkage_list, weights_list):
    # initalize the figure and lists to iterate over
    group_name = linkage_list[0].split("/")[-1].split("_")[0]
    label_list = [group_name, group_name + " processed"]
    main_colors = ["C0", "C5"]
    background_colors = ["C4", "C1"]

    fig, axes = plt.subplots(1, 1, figsize=(4.5, 3.5), layout="constrained")

    # iterate through the original and then processed LD files
    for i, linkage_file in enumerate(linkage_list):
        # load the weights dictionary
        with open(weights_list[i]) as f:
            weights_dictionary = json.load(f)
        weights = {int(x): y for x, y in weights_dictionary.items()}

        # load the dataframe and add the appropriate columns
        df = pd.read_feather(linkage_file)
        df["weight"] = df["position1"].map(weights) * df["position2"].map(weights)
        df["distance"] = df["position2"] - df["position1"]
        df["codon_distance"] = df["distance"] // 3 + 1

        # extract the relevant data
        all_data = df.groupby("codon_distance").apply(
            lambda x: np.average(x.ld, weights=x.weight)
        )
        all_data_background = df.groupby("codon_distance").apply(
            lambda x: np.average(x.bg_ld, weights=x.weight)
        )
        weight_sums = df.groupby("codon_distance")["weight"].sum()
        counts = df.groupby("codon_distance")["ld"].count()

        # filter the data to only distances with sufficient data
        weight_threshold = np.mean(weight_sums.values[:1000]) / 2
        sufficient_data = counts.index[
            np.where((counts.values >= 100) & (weight_sums.values >= weight_threshold))
        ]
        data = all_data[sufficient_data]
        data_background = all_data_background[sufficient_data]

        # calculate the rolling average
        roll_large = 300
        data_rolling = rolling_average(data, roll_large)
        # plot the random expectation
        axes.scatter(
            data_background.index,
            data_background.values,
            label=label_list[i] + " random expectation",
            s=2,
            color=background_colors[i],
            rasterized=True,
        )

        # plot the linkage data and the rolling averages
        axes.scatter(
            data.index,
            data.values,
            label="_nolegend_",
            s=1,
            alpha=0.1,
            color=main_colors[i],
            rasterized=True,
        )
        axes.scatter(
            data[::-1].index[roll_large - 1 :],
            data_rolling,
            label=label_list[i] + " data",
            s=2,
            color=main_colors[i],
            rasterized=True,
        )

    # format and save figure
    axes.set_ylabel(r"linkage disequilibrium ($r$)")
    axes.set_xlabel("core genome codon distance")
    axes.set_xscale("log")
    axes.set_yscale("log")
    axes.set_xlim([0.9, 11000])  # hard-coded to match other paper figures
    axes.set_ylim([0.065, 1])
    axes.legend(
        loc="lower left",
        markerscale=2,
    )
    axes.set_yticks(
        ticks=[0.1, 1],
        labels=["$10^{-1}$", "$10^0$"],
        rotation=90,
        ha="center",
        rotation_mode="anchor",
        va="baseline",
    )
    fig.savefig(output, dpi=450)
    fig.savefig(output_pdf, dpi=450)


# %%


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="plot linkage of subgroups of DE1")
    parser.add_argument(
        "-o", "--output", help="figure png file", type=str, required=True
    )
    parser.add_argument(
        "-p", "--output_pdf", help="figure pdf file", type=str, required=True
    )
    parser.add_argument(
        "-l", "--linkage", help="group linkage file", type=str, required=True
    )
    parser.add_argument(
        "-w", "--weights", help="group weights file", type=str, required=True
    )
    parser.add_argument(
        "-d",
        "--linkage_processed",
        help="processed linkage file",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-s",
        "--weights_processed",
        help="processed weights file",
        type=str,
        required=True,
    )
    args = parser.parse_args()
    plot_processed(
        args.output,
        args.output_pdf,
        [args.linkage, args.linkage_processed],
        [args.weights, args.weights_processed],
    )
