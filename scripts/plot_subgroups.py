"""
This script plots the linkage of DE1 and its subgroups. It could be modified to handle
    a different group, but the formatting is tailored to group DE1.

Author: Jemma M. Fendley
"""

# %%
import pandas as pd, numpy as np
import argparse, json, matplotlib.pyplot as plt
from Bio import AlignIO
from matplotlib.gridspec import GridSpec
from matplotlib.transforms import ScaledTranslation
import parameters  # preset matplotlib formatting
from utils import sort_matrix  # sorts a matrix by hierarchically clustering
from utils import calculate_hamming  # calculates hamming distance
from utils import rolling_average  # calculates a rolling average of data


# %%
def plot_subgroups(output, output_pdf, core_genome, linkage_list, weights_list):

    # load the core genome alignment and extract distance matrix
    core_genome_alignment = np.array(AlignIO.read(core_genome, "fasta"))
    n_phages = core_genome_alignment.shape[0]
    distance_matrix = np.array(
        [
            [
                (calculate_hamming(core_genome_alignment, i, j) if i < j else 0)
                for j in range(n_phages)
            ]
            for i in range(n_phages)
        ]
    )
    i_lower = np.tril_indices(n_phages, -1)
    distance_matrix[i_lower] = distance_matrix.T[i_lower]
    sorted_distance_matrix = sort_matrix(distance_matrix, "complete")[0]

    # initialize figure and relevant lists and add panel labels
    fig = plt.figure(layout="constrained", figsize=(6.5, 3.5))
    gs = GridSpec(2, 2, height_ratios=[5, 1], width_ratios=[1, 2.8], figure=fig)
    ax1, ax2 = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0:2, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.set_visible(False)
    axes, x_loc, y_loc = [ax1, ax2], [-26, -32], [+25.8, -8]
    cmap_max = 0.26  # hard-corded maximum distance to match the other plot
    panel_labels, cmap = ["a)", "b)"], "rainbow"
    main_colors = ["C3", "C0", "C5"]
    background_colors = ["C7", "C4", "C1"]
    label_list = ["DE1", "subgroup 1", "subgroup 2"]
    for i in range(2):
        axes[i].text(
            0.0,
            1.0,
            panel_labels[i],
            transform=(
                axes[i].transAxes
                + ScaledTranslation(x_loc[i] / 72, y_loc[i] / 72, fig.dpi_scale_trans)
            ),
            va="bottom",
        )

    # plot the distance matrix
    im = axes[0].matshow(sorted_distance_matrix, vmin=0, vmax=cmap_max, cmap=cmap)
    fig.colorbar(im, ax=axes[0], shrink=0.4, location="right")
    axes[0].set_yticks(
        ticks=np.arange(0, 44, 20),
        labels=np.arange(0, 44, 20),
        rotation=90,
        ha="center",
        rotation_mode="anchor",
        va="baseline",
    )
    axes[0].set_xticks(
        ticks=np.arange(0, 44, 20),
        labels=np.arange(0, 44, 20),
    )
    axes[0].set_xlabel("phage")
    axes[0].set_ylabel("phage")
    axes[0].xaxis.set_label_position("top")
    axes[0].xaxis.set_ticks_position("top")

    # add arrows to label the differents subgroups in the distance matrix
    x = axes[0].annotate(
        "subgroup 1",
        xy=(0.14, -0.03),
        xytext=(0.14, -0.47),
        fontsize=10,
        ha="center",
        va="bottom",
        xycoords="axes fraction",
        bbox=dict(boxstyle="square", fc="0.8"),
        arrowprops=dict(arrowstyle="-[, widthB=1.1, lengthB=0.3", lw=2.0),
    )
    y = axes[0].annotate(
        "subgroup 2",
        xy=(0.645, -0.03),
        xytext=(0.645, -0.27),
        fontsize=10,
        ha="center",
        va="bottom",
        xycoords="axes fraction",
        bbox=dict(boxstyle="square", fc="0.8"),
        arrowprops=dict(arrowstyle="-[, widthB=3.3, lengthB=0.3", lw=2.0),
    )
    x.set_in_layout(False)
    y.set_in_layout(False)

    random_label = [" random expectation", " random exp.", " random exp."]
    # iterate through DE1, subgroup 1, and subgroup 2
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
        axes[1].scatter(
            data_background.index,
            data_background.values,
            label=label_list[i] + random_label[i],
            s=2,
            color=background_colors[i],
            rasterized=True,
        )

        # plot the linkage data and the rolling averages
        axes[1].scatter(
            data.index,
            data.values,
            label="_nolegend_",
            s=1,
            alpha=0.1,
            color=main_colors[i],
            rasterized=True,
        )
        axes[1].scatter(
            data[::-1].index[roll_large - 1 :],
            data_rolling,
            label=label_list[i] + " data",
            s=2,
            color=main_colors[i],
            rasterized=True,
        )

    # format figure and save
    axes[1].set_ylabel("linkage disequilibrium (r)")
    axes[1].set_xlabel("core genome codon distance")
    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    axes[1].set_xlim([0.9, 8100])  # hard-coded
    axes[1].set_ylim([0.065, 1])  # hard-coded to match the other figure
    axes[1].set_yticks(
        ticks=[0.1, 1],
        labels=["$10^{-1}$", "$10^0$"],
        rotation=90,
        ha="center",
        rotation_mode="anchor",
        va="baseline",
    )
    handles, labels = axes[1].get_legend_handles_labels()
    legend_order = [1, 0, 3, 2, 5, 4]
    new_handles, new_labels = [handles[j] for j in legend_order], [
        labels[j] for j in legend_order
    ]
    axes[1].legend(new_handles, new_labels, loc="lower left", markerscale=2)
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
        "-g", "--group_linkage", help="group linkage file", type=str, required=True
    )
    parser.add_argument(
        "-s", "--group_weights", help="group weights file", type=str, required=True
    )
    parser.add_argument(
        "-c", "--core_genome", help="core genome file", type=str, required=True
    )
    parser.add_argument(
        "-l",
        "--linkage",
        help="list of subgroup linkage files",
        nargs="+",
        required=True,
    )
    parser.add_argument(
        "-w",
        "--weights",
        help="list of subgroup weight files",
        nargs="+",
        required=True,
    )

    args = parser.parse_args()

    plot_subgroups(
        args.output,
        args.output_pdf,
        args.core_genome,
        [args.group_linkage] + args.linkage,
        [args.group_weights] + args.weights,
    )
