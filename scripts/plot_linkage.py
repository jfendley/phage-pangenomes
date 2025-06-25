"""
This script plots linkage disequilibrium as a function of distanceance along the core genome. It
    also records statistics including the asymptote and maximum.

Author: Jemma M. Fendley
"""

import argparse, json
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import parameters  # preset matplotlib formatting


def main():
    parser = argparse.ArgumentParser(description="plot normalized for SNP density")
    parser.add_argument("-o", "--output", help="figure file", type=str, required=True)
    parser.add_argument(
        "-s", "--output_stats", help="statistics file", type=str, required=True
    )
    parser.add_argument(
        "-l", "--linkage", help="linkage disequilibrium file", type=str, required=True
    )
    parser.add_argument(
        "-r", "--reference", help="position to reference file", type=str, required=True
    )
    parser.add_argument(
        "-w", "--weights", help="position to weights file", type=str, required=True
    )
    args = parser.parse_args()

    group_name = args.output.split("/")[-1].split("_")[0]

    # load the weights dictionary
    with open(args.weights) as f:
        weights_dictionary = json.load(f)
    weights = {int(x): y for x, y in weights_dictionary.items()}

    # load the reference dictionary
    with open(args.reference) as f:
        reference_dictionary = json.load(f)
    reference = {int(x): y for x, y in reference_dictionary.items()}
    phage = args.reference.split("_")[-1].split(".json")[0]

    # load the dataframe (slow) and drop some columns to reduce memory
    df = pd.read_feather(args.linkage)
    df = df.drop(columns=["p1", "p2"])

    # calculate the distance and weights
    df["distance"] = df["position2"] - df["position1"]
    df["weight"] = df["position1"].map(weights) * df["position2"].map(weights)
    df["codon_distance"] = df["distance"] // 3 + 1

    # calculate the distance in the reference genome
    position1, position2 = phage + "_position1", phage + "_position2"
    df[position1] = df["position1"].map(reference)
    df[position2] = df["position2"].map(reference)
    df["reference_distance"] = df[position2] - df[position1]
    df["reference_codon_distance"] = df["reference_distance"] // 3 + 1

    # calculate the asymptotic value and the mean background random expectation
    end_df = df[df["codon_distance"] > 3000]
    if np.max(df["codon_distance"]) < 3000:
        end_df = df[df["codon_distance"] > 2500]
        print(group_name)
        print("THIS GROUP HAS NOT 3000 CODONS!!!!")
    asymptote = np.average(end_df["ld"], weights=end_df["weight"])
    background = np.average(df["bg_ld"], weights=df["weight"])

    # prepare the data for plotting
    data = df.groupby("codon_distance").apply(
        lambda x: np.average(x.ld, weights=x.weight)
    )
    counts = df.groupby("codon_distance")["weight"].sum()
    background_data = df.groupby("codon_distance").apply(
        lambda x: np.average(x.bg_ld, weights=x.weight)
    )
    reference_df = df[df["reference_distance"] > 0]
    reference_data = reference_df.groupby("reference_codon_distance").apply(
        lambda x: np.average(x.ld, weights=x.weight)
    )
    reference_counts = reference_df.groupby("reference_codon_distance")["weight"].sum()
    reference_background = reference_df.groupby("reference_codon_distance").apply(
        lambda x: np.average(x.bg_ld, weights=x.weight)
    )

    # initialize the plots
    fig, axes = plt.subplots(
        nrows=2,
        ncols=2,
        figsize=(6, 6),
        gridspec_kw={"width_ratios": [1, 1]},
        layout="constrained",
    )
    backround_plot = [background_data, reference_background]
    counts_plot = [counts, reference_counts]
    for i, plot_data in enumerate([data, reference_data]):
        # plot the data, background data, and asymptote
        axes[0, i].scatter(
            plot_data.index, plot_data.values, s=1, alpha=0.5, label="data"
        )
        axes[0, i].scatter(
            backround_plot[i].index,
            backround_plot[i].values,
            s=1,
            alpha=0.5,
            label="random expectation",
        )
        axes[0, i].axhline(
            asymptote,
            color="k",
            alpha=0.5,
            label="asymptote",
            linewidth=1,
        )

        # plot the sum of the weights
        axes[1, i].scatter(
            counts_plot[i].index, counts_plot[i].values, s=1, alpha=0.5, label="data"
        )

        # figure formatting
        axes[i, 0].legend(loc=3)
        axes[0, i].set_ylim(bottom=0.04, top=1)
        for j in range(2):
            axes[i, j].set_xscale("log")
            axes[i, j].set_yscale("log")
    axes[0, 0].set_ylabel("weighted linkage disequilibrium")
    axes[1, 0].set_ylabel("sum of weights")
    axes[1, 0].set_xlabel("core genome codon distance")
    axes[1, 1].set_xlabel("reference codon distance")
    fig.suptitle("group " + group_name)

    # save the figure and some statistics
    dict = {
        "group": group_name,
        "asymptote": asymptote,
        "background": background,
        "first": np.max(data.values[0:10]),
    }
    save_df = pd.DataFrame([dict])
    save_df.to_csv(args.output_stats, sep="\t", index=False)
    fig.savefig(args.output, dpi=450)


if __name__ == "__main__":
    main()
