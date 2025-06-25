"""
This script plots the two main linkage diseqiulibrium figures in the main text.

Author: Jemma M. Fendley
"""

# %%
import pandas as pd, numpy as np
import argparse, json, matplotlib.pyplot as plt
from Bio import AlignIO
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.transforms import ScaledTranslation
from tqdm import tqdm
import parameters  # preset matplotlib formatting
from utils import sort_matrix, calculate_hamming, find_snp_positions, rolling_average


def get_bin(position):
    """
    Returns a bin identifier given a position
    """
    binsize = 100
    return position // binsize


# %%
def plot_figures(
    groups,
    linkage_files,
    weights_files,
    core_genome_files,
    output,
    output_pdf,
    output_structure,
    output_structure_pdf,
):

    # initalizing colors
    main_colors = ["C0", "C5", "C3"]
    secondary_colors = ["C4", "C1", "C7"]
    binsize, first_cmap, second_cmap = 100, "PuBuGn", "rainbow"
    inset_max = 0.26  # hard-coded, this is the largest pairwise hamming distanc
    list_of_labels = ["a)", "b)", "c)", "d)", "e)", "f)"]

    # initializing the figures
    fig2a, axes_main = plt.subplots(
        1,
        3,
        figsize=(7, 2.6),
        gridspec_kw={"width_ratios": [1, 1, 1]},
        dpi=450,
        layout="compressed",
    )
    fig = plt.figure(layout="compressed", figsize=(7, 3.4))
    gs = GridSpec(2, 3, width_ratios=[1, 1, 1], height_ratios=[5, 1.2], figure=fig)
    ax1, ax4 = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0])
    ax2, ax5 = fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 1])
    ax3, ax6 = fig.add_subplot(gs[0, 2]), fig.add_subplot(gs[1, 2])
    axes_structure = [ax1, ax2, ax3, ax4, ax5, ax6]
    edge_size, x_offset, y_offset, inset_fontsize = 0.7, -0.02, -0.02, 10
    axins = [
        inset_axes(
            ax,
            width=edge_size,
            height=edge_size,
            borderpad=0,
            bbox_to_anchor=(0 + x_offset, 0 + y_offset, 1, 1),
            bbox_transform=ax.transAxes,
        )
        for ax in [ax1, ax2, ax3]
    ]

    random_label = ["random expectation", "random exp.", "random exp."]
    for i, group in enumerate(tqdm(groups)):
        # find the relevant files
        linkage_file, weights_file = linkage_files[i], weights_files[i]
        core_genome_file = core_genome_files[i]

        # load the weights dictionary
        with open(weights_file) as f:
            weights_dictionary = json.load(f)
        weights = {int(x): y for x, y in weights_dictionary.items()}

        # load the dataframe and add the appropriate columns
        df = pd.read_feather(linkage_file)
        df["bin1"] = df["position1"].apply(get_bin)
        df["bin2"] = df["position2"].apply(get_bin)
        df["weight"] = df["position1"].map(weights) * df["position2"].map(weights)
        df["distance"] = df["position2"] - df["position1"]
        df["codon_distance"] = df["distance"] // 3 + 1

        # load the core genome alignment and extract distance matrix and snp positions
        core_genome_alignment = np.array(AlignIO.read(core_genome_file, "fasta"))
        n_phages, n_positions = core_genome_alignment.shape
        snp_positions = find_snp_positions(core_genome_alignment)[0]
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

        # plot the core hamming distance matrix
        im = axins[i].matshow(
            sorted_distance_matrix, vmin=0, vmax=inset_max, cmap=second_cmap
        )
        axins[i].xaxis.set_tick_params(labelsize=inset_fontsize - 2)
        axins[i].yaxis.set_tick_params(labelsize=inset_fontsize - 2)
        axins[i].xaxis.set_ticks_position("bottom")
        axins[i].set_xlabel("phage", fontsize=inset_fontsize, labelpad=1)
        axins[i].set_ylabel("phage", fontsize=inset_fontsize, labelpad=1)

        # initalize enrichment matrix
        max_bin = (n_positions - 1) // binsize
        enrichment_matrix = np.full((max_bin + 1, max_bin + 1), np.nan)

        # choose a threshold for what is considered "high" linkage
        threshold = np.percentile(df["ld"], 90)
        data_bins_high = df.groupby(["bin1", "bin2"])["ld"].count()
        data_bins_high_high = (
            df[df["ld"] >= threshold].groupby(["bin1", "bin2"])["ld"].count()
        )

        # each bin will have fraction of that bin with high linkage
        for j in data_bins_high_high.index:
            if data_bins_high[j] >= 10:
                enrichment = data_bins_high_high[j] / data_bins_high[j]
                enrichment_matrix[max(j[0], j[1]), min(j[0], j[1])] = enrichment
        for j in [
            x
            for x in data_bins_high.index
            if x not in data_bins_high_high.index and data_bins_high[x] >= 10
        ]:
            enrichment_matrix[max(j[0], j[1]), min(j[0], j[1])] = 0

        # plot the enrichment matrix and add labels
        mat = axes_structure[i].matshow(
            enrichment_matrix, vmin=0, vmax=1, cmap=first_cmap
        )
        axes_structure[i].xaxis.set_ticks_position("bottom")
        ticks = [np.arange(0, 410, 100), np.arange(0, 410, 100), np.arange(0, 310, 100)]
        snp_ticks = [
            np.arange(0, 40010, 10000),
            np.arange(0, 40010, 10000),
            np.arange(0, 30010, 10000),
        ]
        tick_labels = [
            ["0", "10000", "20000", "30000", ""],
            ["0", "", "20000", "", "40000"],
            ["0", "10000", "20000", "30000"],
        ]
        axes_structure[i].set_xticks(ticks=ticks[i], labels=tick_labels[i])
        axes_structure[i].set_yticks(
            ticks=ticks[i],
            labels=tick_labels[i],
            rotation=90,
            ha="center",
            rotation_mode="anchor",
            va="baseline",
        )
        axes_structure[i + 3].set_xticks(ticks=snp_ticks[i], labels=tick_labels[i])
        axes_structure[i].set_title(group + " high-linkage regions", fontsize=12)
        axes_structure[i + 3].hist(
            snp_positions,
            bins=range(0, n_positions, binsize),
            density=False,
            cumulative=False,
        )
        axes_structure[i + 3].set_xlabel("core genome position")
        axes_structure[i + 3].set_xlim([0, n_positions])
        axes_structure[i + 3].set_ylim([0, 100])
        axes_structure[i + 3].set_yticks(
            ticks=[0, 40, 80],
            labels=[0, 40, 80],
            rotation=90,
            ha="center",
            rotation_mode="anchor",
            va="baseline",
        )

        # sort the linkage data by codon distance
        all_data = df.groupby("codon_distance").apply(
            lambda x: np.average(x.ld, weights=x.weight)
        )
        all_data_background = df.groupby("codon_distance").apply(
            lambda x: np.average(x.bg_ld, weights=x.weight)
        )

        # filter the data only to distances with sufficient data
        weight_sums = df.groupby("codon_distance")["weight"].sum()
        weight_threshold = np.mean(weight_sums.values[:1000]) / 2
        counts = df.groupby("codon_distance")["ld"].count()
        sufficient_data = counts.index[
            np.where((counts.values >= 100) & (weight_sums.values >= weight_threshold))
        ]
        data = all_data[sufficient_data]
        data_background = all_data_background[sufficient_data]

        # calculate rolling averages
        roll_large = 300
        data_all = rolling_average(data, roll_large)

        # plot the data and random expectation
        axes_main[i].scatter(
            data[::-1].index[roll_large - 1 :],
            data_all,
            label="data",
            s=2,
            color=main_colors[i],
            rasterized=True,
        )
        axes_main[i].scatter(
            data_background[::-1].index,
            data_background[::-1].values,
            label=random_label[i],
            s=2,
            color=secondary_colors[i],
            rasterized=True,
        )
        axes_main[i].scatter(
            data[::-1].index,
            data[::-1].values,
            label="_nolegend_",
            s=0.7,
            alpha=0.07,
            color=main_colors[i],
            rasterized=True,
        )

        # figure formatting
        axes_main[i].set_title("group " + group, pad=10)
        axes_structure[i].text(
            0.0,
            1.0,
            list_of_labels[i],
            transform=(
                axes_structure[i].transAxes
                + ScaledTranslation(-15 / 72, +4 / 72, fig.dpi_scale_trans)
            ),
            va="bottom",
        )
        axes_main[i].set_xlabel("core genome codon distance")
        axes_main[i].set_xscale("log")
        axes_main[i].set_yscale("log")
        axes_main[i].set_xlim([0.9, 11000])  # hard-coded
        axes_main[i].set_ylim([0.065, 1])  # hard-coded
        axes_main[i].set_yticks(
            ticks=[0.1, 1],
            labels=["$10^{-1}$", "$10^0$"],
            rotation=90,
            ha="center",
            rotation_mode="anchor",
            va="baseline",
        )
        axes_main[i].text(
            0.0,
            1.0,
            list_of_labels[i],
            transform=(
                axes_main[i].transAxes
                + ScaledTranslation(-18 / 72, +8 / 72, fig2a.dpi_scale_trans)
            ),
            va="bottom",
        )

    # figure formatting
    axins[0].set_yticks(
        ticks=np.arange(0, 24, 10),
        labels=np.arange(0, 24, 10),
        rotation=90,
        ha="center",
        rotation_mode="anchor",
        va="baseline",
    )
    axins[1].set_yticks(
        ticks=np.arange(0, 104, 50),
        labels=np.arange(0, 104, 50),
        rotation=90,
        ha="center",
        rotation_mode="anchor",
        va="baseline",
    )
    axins[2].set_yticks(
        ticks=np.arange(0, 44, 20),
        labels=np.arange(0, 44, 20),
        rotation=90,
        ha="center",
        rotation_mode="anchor",
        va="baseline",
    )
    axins[0].set_xticks(
        ticks=np.arange(0, 24, 10),
        labels=np.arange(0, 24, 10),
    )
    axins[1].set_xticks(
        ticks=np.arange(0, 104, 50),
        labels=np.arange(0, 104, 50),
    )
    axins[2].set_xticks(
        ticks=np.arange(0, 44, 20),
        labels=np.arange(0, 44, 20),
    )
    axes_main[0].set_ylabel(r"linkage disequilibrium ($r$)")
    axes_main[0].legend(loc="lower left", markerscale=2)
    axes_structure[3].text(25.0, 70.0, "binsize = 100bp", fontsize=10)
    axes_structure[3].set_ylabel("n. SNPs")
    axes_structure[0].set_ylabel("core genome position")

    # add colorbars
    cbar1 = fig.colorbar(mat, ax=axes_structure[0], shrink=0.6, location="left")
    cbar2 = fig.colorbar(im, ax=axes_structure[2], shrink=0.6, location="right")
    cbar2.ax.set_yticks(
        ticks=[0, 0.1, 0.2],
        labels=[0, 0.1, 0.2],
        rotation=270,
        ha="center",
        rotation_mode="anchor",
        va="baseline",
    )
    cbar2.ax.set_ylabel(
        "core Hamming distance", rotation=270, fontsize=10, labelpad=10.0
    )
    cbar1.ax.set_yticks(
        ticks=[0, 0.5, 1],
        labels=[0, 0.5, 1],
        rotation=90,
        ha="center",
        rotation_mode="anchor",
        va="baseline",
    )
    cbar1.ax.set_ylabel("fraction high linkage", fontsize=10, labelpad=0)

    # save figures
    fig.savefig(output_structure, dpi=450)
    fig.savefig(output_structure_pdf, dpi=450)
    fig2a.savefig(output, dpi=450)
    fig2a.savefig(output_pdf, dpi=450)


# %%
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="create LD file")
    parser.add_argument(
        "-o", "--output_png", help="main figure png", type=str, required=True
    )
    parser.add_argument(
        "-p", "--output_pdf", help="main figure pdf", type=str, required=True
    )
    parser.add_argument(
        "-s",
        "--output_structure_png",
        help="structure figure png",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-f",
        "--output_structure_pdf",
        help="structure figure pdf",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-g", "--groups", help="list of groups", nargs="+", required=True
    )
    parser.add_argument(
        "-l", "--linkage", help="list of linkage files", nargs="+", required=True
    )
    parser.add_argument(
        "-w", "--weights", help="list of weights files", nargs="+", required=True
    )
    parser.add_argument(
        "-c", "--core", help="list of core genome files", nargs="+", required=True
    )
    args = parser.parse_args()

    plot_figures(
        args.groups,
        args.linkage,
        args.weights,
        args.core,
        args.output_png,
        args.output_pdf,
        args.output_structure_png,
        args.output_structure_pdf,
    )
