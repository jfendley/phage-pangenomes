"""
This script plots a figure that shows the effect on the linkage disequilbirium results
    from weighting and using a reference phage for distances.


Author: Jemma M. Fendley
"""

# %%
import pandas as pd, numpy as np
import argparse, json, matplotlib.pyplot as plt
from matplotlib.transforms import ScaledTranslation
import parameters  # preset matplotlib formatting
from utils import rolling_average  # calculates a rolling average of the data


def plot_linkage_weights(linkage, reference, weights_file, output, output_pdf):
    # load the weights dictionary
    with open(weights_file) as f:
        weights_dictionary = json.load(f)
    weights = {int(x): y for x, y in weights_dictionary.items()}

    # load the dataframe and add the appropriate columns
    df = pd.read_feather(linkage)
    df["weight"] = df["position1"].map(weights) * df["position2"].map(weights)
    df["distance"] = df["position2"] - df["position1"]
    df["codon_distance"] = df["distance"] // 3 + 1

    # initialize figures and colors
    fig, axes = plt.subplots(
        2, 2, figsize=(7, 7), layout="constrained", gridspec_kw={"width_ratios": [1, 1]}
    )
    main_colors, secondary_colors = ["C0", "C5"], ["C4", "C1"]
    list_of_labels = ["a)", "b)", "c)", "d)"]

    # extract the data by codon distance to plot
    data = df.groupby("codon_distance").apply(
        lambda x: np.average(x.ld, weights=x.weight)
    )
    data_random = df.groupby("codon_distance").apply(
        lambda x: np.average(x.bg_ld, weights=x.weight)
    )
    counts = df.groupby("codon_distance")["ld"].count()
    sum_weights = df.groupby("codon_distance")["weight"].sum()
    data_no_weight = df.groupby("codon_distance")["ld"].mean()

    data_random_no_weight = df.groupby("codon_distance")["bg_ld"].mean()

    # calculate the log-scaled rolling averages
    roll_large = 300
    data_all = rolling_average(data, roll_large)
    data_no_weight_all = rolling_average(data_no_weight, roll_large)

    # load the reference positions dictionary
    with open(reference) as f:
        reference_dictionary = json.load(f)
    reference_map = {int(x): y for x, y in reference_dictionary.items()}

    # add the reference information to the dataframe
    position1 = reference.split("_")[-1].split(".json")[0] + "_position1"
    position2 = reference.split("_")[-1].split(".json")[0] + "_position2"
    df[position1] = df["position1"].map(reference_map)
    df[position2] = df["position2"].map(reference_map)
    df["reference_distance"] = df[position2] - df[position1]
    assert all(df["reference_distance"] >= 0)
    df["codon_reference_distance"] = df["reference_distance"] // 3 + 1

    # extract the reference data
    data_reference = df.groupby("codon_reference_distance").apply(
        lambda x: np.average(x.ld, weights=x.weight)
    )
    data_reference_random = df.groupby("codon_reference_distance").apply(
        lambda x: np.average(x.bg_ld, weights=x.weight)
    )
    data_reference_no_weight = df.groupby("codon_reference_distance")["ld"].mean()
    data_reference_no_weight_random = df.groupby("codon_reference_distance")[
        "bg_ld"
    ].mean()
    counts_reference = df.groupby("codon_reference_distance")["ld"].count()
    sum_weights_reference = df.groupby("codon_reference_distance")["weight"].sum()
    data_reference_all = rolling_average(data_reference, roll_large)
    data_reference_no_weight_all = rolling_average(data_reference_no_weight, roll_large)

    # organize the data for iterative plotting
    main_data = [[data_no_weight, data], [data_reference_no_weight, data_reference]]
    data_background = [
        [data_random_no_weight, data_random],
        [data_reference_no_weight_random, data_reference_random],
    ]
    all_data = [
        [data_no_weight_all, data_all],
        [data_reference_no_weight_all, data_reference_all],
    ]
    count_data = [[counts, sum_weights], [counts_reference, sum_weights_reference]]
    data_labels = ["data", "data SNP weighted"]
    count_labels = ["n. pairs of SNPs", "sum of weights"]
    background_labels = ["random expectation", "random expectation SNP weighted"]

    for i in range(2):
        for j in range(2):
            # plot the random expectations
            axes[0, i].scatter(
                data_background[i][j].index,
                data_background[i][j].values,
                label=background_labels[j],
                s=2,
                rasterized=True,
                color=secondary_colors[j],
                alpha=0.8,
            )
            # plot the main data
            axes[0, i].scatter(
                main_data[i][j].index,
                main_data[i][j].values,
                rasterized=True,
                color=main_colors[j],
                label="_nolegend_",
                s=0.7,
                alpha=0.1,
            )
            # plot the rolling averages
            axes[0, i].scatter(
                main_data[i][j][::-1].index[roll_large - 1 :],
                all_data[i][j],
                label=data_labels[j],
                s=2,
                alpha=0.8,
                color=main_colors[j],
                rasterized=True,
            )
            # plot the amount of data
            axes[1, i].scatter(
                count_data[i][j].index,
                count_data[i][j].values,
                label=count_labels[j],
                s=2,
                alpha=0.8,
                color=main_colors[j],
                rasterized=True,
            )

            # figure formatting
            label = list_of_labels[j] if i == 0 else list_of_labels[j + 2]
            axes[i, j].text(
                0.0,
                1.0,
                label,
                transform=(
                    axes[i, j].transAxes
                    + ScaledTranslation(-30 / 72, +5 / 72, fig.dpi_scale_trans)
                ),
                va="bottom",
            )
            axes[0, i].set_ylim([0.0248, 1])
            axes[1, i].set_ylim([0.8, 4.5e5])
            axes[i, j].set_yscale("log")
            axes[i, j].set_xscale("log")
            axes[i, j].set_xlim([1, 2e4])

    # more figure formatting and save figuerw
    for ind1 in range(2):
        axes[ind1, 0].legend(markerscale=4)
    axes[0, 0].set_ylabel(r"linkage disequilibrium ($r$)")
    axes[1, 0].set_ylabel("")
    axes[1, 0].set_xlabel("core genome codon distance")
    axes[1, 1].set_xlabel("reference genome codon distance")
    axes[1, 0].set_ylabel("amount of data")
    fig.savefig(output, dpi=450)
    fig.savefig(output_pdf, dpi=450)


# %%
if __name__ == "__main__":
    base = "/home/jemma/Projects/phagephams/analyzeData/results/"
    linkage = base + "groups/A11/A11_core/A11_linkage.feather"
    reference = base + "groups/A11/A11_core/A11_positions_MN703405.json"
    weights = base + "groups/A11/A11_core/A11_snp_density.json"

    output = base + "paper_figures/linkage_weights_A11.png"
    output_pdf = base + "paper_figures/linkage_weights_A11.pdf"

    plot_linkage_weights(linkage, reference, weights, output, output_pdf)
