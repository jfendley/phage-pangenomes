"""
This script plots Figure 7 in the main text. It shows the linkage disequilibrium of group DE1 and its
    subgroups, along with the spread of residual linkage as a function of distance statistics.

Author: Jemma M. Fendley
"""

import pandas as pd, numpy as np, seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import argparse, json, matplotlib.pyplot as plt, matplotlib as mpl
from Bio import AlignIO, Phylo
from matplotlib.gridspec import GridSpec
from matplotlib.transforms import ScaledTranslation
import parameters  # preset matplotlib formatting
from utils import sort_matrix  # sorts a matrix by hierarchically clustering
from utils import calculate_hamming  # calculates hamming distance
from utils import rolling_average  # calculates a rolling average of data
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, _DistanceMatrix
from matplotlib.lines import Line2D
from cycler import cycler

# change the color cycle to better match other paper figures
color_palette = [
    "#BBBBBB",
    "#332288",
    "#CC6677",
    "#DDCC77",
    "#117733",
    "#88CCEE",
    "#882255",
    "#44AA99",
    "#999933",
    "#AA4499",
]
mpl.rcParams["axes.prop_cycle"] = cycler(color=color_palette)


def plot_subgroups(
    output,
    output_pdf,
    linkage_information,
    statistics,
    core_genome,
    linkage_list,
    weights_list,
):
    """
    Plots the desired figure given the input arguments.
    """

    # load the data frames, combine, and add relevant columns
    df_linkage = pd.read_csv(linkage_information, sep="\t")
    df_distance = pd.read_csv(statistics, sep="\t")
    df_residual = pd.merge(df_linkage, df_distance, on="group")

    df_residual["residual"] = df_residual["asymptote"] - df_residual["background"]

    # load the core genome alignment for DE1 and extract distance matrix
    core_genome_alignment_list = AlignIO.read(core_genome, "fasta")
    phage_order = [record.id for record in core_genome_alignment_list]
    core_genome_alignment = np.array(core_genome_alignment_list)
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

    # construct a neighbour-joining tree from the Hamming distance matrix
    constructor = DistanceTreeConstructor()
    dist_matrix = [
        [distance_matrix[i, j] for j in range(i + 1)]
        for i in range(len(distance_matrix))
    ]
    dm = _DistanceMatrix(phage_order, dist_matrix)
    nj_tree = constructor.nj(dm)

    # find all of the internal branches
    internal_clades = nj_tree.get_nonterminals()

    # create a dictionary of branch to length
    branch_lengths = [x.branch_length for x in internal_clades]

    # find the longest internal branch
    longest_internal_branch = internal_clades[np.argmax(branch_lengths)]

    nj_tree.ladderize(reverse=True)

    # initialize figure and relevant lists and add panel labels
    fig = plt.figure(layout="constrained", figsize=(4.5, 6.25))
    gs = GridSpec(
        3,
        4,
        height_ratios=[0.85, 2.25, 1.4],
        width_ratios=[2, 0.01, 0.01, 2],
        figure=fig,
    )
    ax1, ax2 = fig.add_subplot(gs[1, 0:4]), fig.add_subplot(gs[2, 0:2])
    ax3 = fig.add_subplot(gs[2, 2:4])

    axins = fig.add_subplot(gs[0, 0])
    ax4 = fig.add_subplot(gs[0, 2:4])
    axes, x_loc, y_loc = (
        [axins, ax1, ax2, ax3, ax4],
        [-37, -30, -28, -28, -30],
        [-12, -9, -6, -6, -9],
    )

    # plot tree with subgroups and longest branch colored
    nj_tree.root.color = "#1965b0"
    if longest_internal_branch:
        longest_internal_branch.color = "#BBBBBB"
        for child in longest_internal_branch.clades:
            child.color = "#762A83"
    mpl.rcParams["lines.linewidth"] = 1
    Phylo.draw(nj_tree, axes=axes[4], do_show=False, label_func=lambda x: "")

    # add legend and lines
    custom_lines = [
        Line2D([0], [0], color="#762A83", lw=1),
        Line2D([0], [0], color="#1965b0", lw=1),
        Line2D([0], [0], color="#BBBBBB", lw=1),
    ]
    axes[4].legend(
        custom_lines,
        ["subgroup 1", "subgroup 2", "longest branch"],
        fontsize=8,
        loc="lower right",
        handletextpad=0.5,
        handlelength=1,
        borderaxespad=0.2,
        labelspacing=0.2,
    )

    # set colors for rest of figures
    cmap_max = 0.26  # hard-corded maximum distance to match the other plot
    panel_labels, cmap = ["a)", "b)", "c)", "d)", "b)"], "rainbow"
    main_colors = ["C4", "#762A83", "#1965b0"]
    background_colors = ["C8", "#9970AB", "#7bafde"]
    label_list = ["DE1", "subgroup 1", "subgroup 2"]

    # add panel labels
    for i in range(4):
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
    edge_size, x_offset, y_offset = 0.65, 0.88, -0.95

    # format colorbar
    caxins = inset_axes(
        axins,
        height=edge_size / 14,
        width=edge_size - 0.08,
        borderpad=0,
        bbox_to_anchor=(0 + x_offset, 0 + y_offset, 1, 1),
        bbox_transform=axins.transAxes,
    )
    cbar = fig.colorbar(im, cax=caxins, shrink=0.2, orientation="horizontal")
    cbar.ax.set_xticks(ticks=[0, 0.1, 0.2], labels=[0, 0.1, 0.2], fontsize=8)
    cbar.ax.set_xlabel("core distance", labelpad=1, fontsize=8)

    # format figure
    axes[0].set_yticks(
        ticks=[0, 19, 39],
        labels=[1, 20, 40],
        rotation=90,
        ha="center",
        rotation_mode="anchor",
        va="baseline",
    )
    axes[4].set_yticks(
        ticks=[1, 20, 40],
        labels=[1, 20, 40],
        rotation=90,
        ha="center",
        rotation_mode="anchor",
        va="baseline",
    )
    axes[4].set_xlim([-0.01, 0.21])
    axes[4].set_ylabel("phage")
    axes[0].set_xticks(
        ticks=[0, 19, 39],
        labels=[1, 20, 40],
    )
    axes[0].set_xlabel("phage")
    axes[0].set_ylabel("phage")
    axes[0].xaxis.set_label_position("bottom")
    axes[0].xaxis.set_ticks_position("bottom")

    # add arrows to label the differents subgroups in the distance matrix
    x = axes[0].annotate(
        "subgroup 1",
        xy=(1.04, 0.86),
        xytext=(1.22, 0.86),
        fontsize=8,
        ha="left",
        va="center",
        xycoords="axes fraction",
        bbox=dict(boxstyle="square", fc="0.8"),
        arrowprops=dict(arrowstyle="-[, widthB=0.9, lengthB=0.35", lw=1.5),
    )
    y = axes[0].annotate(
        "subgroup 2",
        xy=(1.04, 0.355),
        xytext=(1.22, 0.355),
        fontsize=8,
        ha="left",
        va="center",
        xycoords="axes fraction",
        bbox=dict(boxstyle="square", fc="0.8"),
        arrowprops=dict(arrowstyle="-[, widthB=2.6, lengthB=0.35", lw=1.5),
    )
    x.set_in_layout(False)
    y.set_in_layout(False)

    x = axes[0].annotate(
        "",
        xy=(2.06, 0.86),
        xytext=(1.86, 0.86),
        ha="right",
        va="center",
        xycoords="axes fraction",
        bbox=dict(boxstyle="square", fc="0.8"),
        arrowprops=dict(arrowstyle="-[, widthB=0.45, lengthB=0.25", lw=1.4),
    )
    y = axes[0].annotate(
        "",
        xy=(2.06, 0.355),
        xytext=(1.86, 0.355),
        ha="right",
        va="center",
        xycoords="axes fraction",
        bbox=dict(boxstyle="square", fc="0.8"),
        arrowprops=dict(arrowstyle="-[, widthB=1.53, lengthB=0.25", lw=1.4),
    )
    x.set_in_layout(False)
    y.set_in_layout(False)

    # iterate through DE1, subgroup 1, and subgroup 2 to plot linkage

    random_label = [" random expectation", " random exp.", " random exp."]
    shape_list = ["o", "^", "s"]
    background_shape_list = ["*", "D", "x"]
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
            marker=background_shape_list[i],
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
            marker=shape_list[i],
        )
        axes[1].scatter(
            data[::-1].index[roll_large - 1 :],
            data_rolling,
            label=label_list[i] + " data",
            s=2,
            color=main_colors[i],
            rasterized=True,
            marker=shape_list[i],
        )

    # format figure
    axes[1].set_ylabel(r"linkage disequilibrium ($r$)")
    axes[1].set_xlabel("core genome codon distance")
    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    axes[1].set_xlim([0.9, 7500])  # hard-coded
    axes[1].set_ylim([0.099, 1])  # hard-coded to match the other figure
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
    axes[1].legend(
        new_handles,
        new_labels,
        loc="lower left",
        markerscale=2,
        borderaxespad=0.5,
        fontsize=10,
        labelspacing=0.2,
        framealpha=0.6,
        borderpad=0.3,
    )

    # now plot the residual linkage

    # separate out the named groups
    select_groups_list = ["all", "A11", "E", "EE", "DE1", "F1"]
    df_residual["mean_hamming_percent"] = 100 - df_residual["mean_hamming"] * 100

    def is_special(name):
        if name in select_groups_list:
            return name
        else:
            return "all"

    df_residual["is_special"] = df_residual["group"].apply(is_special)
    df_residual["is_special"] = pd.Categorical(
        df_residual["is_special"], categories=select_groups_list, ordered=True
    )
    df_residual = df_residual.sort_values(by="is_special")

    # plot the mean hamming and residual linkage
    sns.scatterplot(
        data=df_residual,
        x="mean_hamming_percent",
        y="residual",
        ax=axes[2],
        size="is_special",
        sizes=[20, 50, 50, 50, 50, 50],
        hue="is_special",
        style="is_special",
    )

    # format figure
    axes[2].set_xlim([69, 101.5])
    for i in [2, 3]:
        axes[i].set_ylim([-0.01, 0.53])
        axes[i].set_yticks(
            ticks=[0.0, 0.2, 0.4],
            labels=[0.0, 0.2, 0.4],
            rotation=90,
            ha="center",
            rotation_mode="anchor",
            va="baseline",
        )
    axes[2].get_legend().set_visible(False)
    axes[2].set_ylabel("residual linkage")
    axes[2].set_xlabel("mean p.w. core ANI (\%)")

    # plot the max internal branch length and residual linkage
    sns.scatterplot(
        data=df_residual,
        x="nj_max_non_terminal_branch_length",
        y="residual",
        ax=axes[3],
        size="is_special",
        sizes=[20, 50, 50, 50, 50, 50],
        hue="is_special",
        style="is_special",
    )

    # format figure
    handles, labels = axes[3].get_legend_handles_labels()
    axes[3].legend(
        handles[1:],
        labels[1:],
        loc="lower right",
        markerscale=0.75,
        borderaxespad=0.2,
        borderpad=0.3,
        labelspacing=0.2,
        fontsize=8,
    ).set_title(None)
    axes[3].set_ylabel("residual linkage")
    axes[3].set_xlabel("max. branch length")

    # save figure
    fig.savefig(output, dpi=450)
    fig.savefig(output_pdf, dpi=450)


if __name__ == "__main__":

    # load the input
    parser = argparse.ArgumentParser(
        description="plot residual linkage and linkage of subgroups of DE1"
    )
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
        "-t", "--group_weights", help="group weights file", type=str, required=True
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

    parser.add_argument(
        "-s", "--statistics", help="pairwise statistics file", type=str, required=True
    )
    parser.add_argument(
        "-i",
        "--linkage_information",
        help="linkage statistics file",
        type=str,
        required=True,
    )
    args = parser.parse_args()

    # run function to plot the figure
    plot_subgroups(
        args.output,
        args.output_pdf,
        args.linkage_information,
        args.statistics,
        args.core_genome,
        [args.group_linkage] + args.linkage,
        [args.group_weights] + args.weights,
    )
