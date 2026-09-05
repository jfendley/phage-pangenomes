"""
This script plots a figure for the SI which shows the length scales on which linkage decays.

Author: Jemma M. Fendley
"""

import pandas as pd, argparse, json
from matplotlib.transforms import ScaledTranslation
import matplotlib.pyplot as plt, seaborn as sns, numpy as np
import parameters  # preset matplotlib formatting
from scipy.optimize import curve_fit


def main():
    parser = argparse.ArgumentParser(
        description="plot distance of linkage decay for all groups"
    )
    parser.add_argument(
        "-o", "--output", help="figure png file", type=str, required=True
    )
    parser.add_argument(
        "-p", "--output_pdf", help="figure pdf file", type=str, required=True
    )
    parser.add_argument(
        "-s",
        "--linkage_information",
        help="linkage statistics file",
        type=str,
        required=True,
    )
    # this script contains two different example groups
    parser.add_argument(
        "-l1", "--linkage1", help="linkage file 1", type=str, required=True
    )
    parser.add_argument(
        "-w1", "--weights1", help="position to weights file 1 ", type=str, required=True
    )
    parser.add_argument(
        "-l2", "--linkage2", help="linkage file 2", type=str, required=True
    )
    parser.add_argument(
        "-w2", "--weights2", help="position to weights file 2", type=str, required=True
    )
    args = parser.parse_args()

    # load the data frames, combine, and add relevant columns
    df = pd.read_csv(args.linkage_information, sep="\t").set_index("group")
    df["residual"] = df["asymptote"] - df["background"]

    # filter by those that have sufficiently low residual linkage and asymptotic values
    df = df[(df["residual"] <= 0.2) & (df["asymptote"] <= 0.32)]

    # create a dictionary of group to half-way value (the distance at which half the linkage is lost)
    halfway_dict = pd.Series(df.halfway, index=df.index).to_dict()

    # initialize the figure
    fig, axes = plt.subplots(
        nrows=1,
        ncols=3,
        figsize=(7, 2.6),
        gridspec_kw={"width_ratios": [1, 1, 1]},
        layout="compressed",
    )

    # plot the continuum of distance statistics
    sns.scatterplot(data=df, x="distance", y="halfway", ax=axes[2], alpha=0.8)

    # format figure and add means
    axes[2].set_xlabel("distance scale (codons)")
    axes[2].set_ylabel("half-distance (codons)")
    axes[2].axvline(
        np.mean(df["distance"]),
        color="C1",
        label="mean: {0:0.0f}".format(np.mean(df["distance"])),
    )
    axes[2].axhline(
        np.mean(df["halfway"]),
        color="C2",
        linestyle="--",
        label="mean: {0:0.0f}".format(np.mean(df["halfway"])),
    )

    # now plot a figure for each of the two example groups
    weights_list = [args.weights1, args.weights2]
    data_list = [args.linkage1, args.linkage2]
    group = ["AY", "BU"]
    for i, linkage_file in enumerate(data_list):
        # double check that it is the correct group
        assert group[i] in linkage_file, "Error in file"
        assert group[i] in weights_list[i], "Error in file"

        # load the weights dictionary
        with open(weights_list[i]) as f:
            weights_dictionary = json.load(f)
        weights = {int(x): y for x, y in weights_dictionary.items()}

        # load the dataframe (slow) and drop some columns to reduce memory
        df = pd.read_feather(linkage_file)
        df = df.drop(columns=["p1", "p2"])

        # calculate the distance and weights
        df["distance"] = df["position2"] - df["position1"]
        df["weight"] = df["position1"].map(weights) * df["position2"].map(weights)
        df["codon_distance"] = df["distance"] // 3 + 1
        df = df.drop(columns=["distance"])  # reduce memory

        # calculate the asymptotic value
        end_df = df[df["codon_distance"] > 3000]
        if np.max(df["codon_distance"]) < 3000:
            end_df = df[df["codon_distance"] > 2500]
        asymptote = np.average(end_df["ld"], weights=end_df["weight"])
        del end_df  # reduce memory

        # prepare the data for plotting
        data = df.groupby("codon_distance").apply(
            lambda x: np.average(x.ld, weights=x.weight)
        )

        # find the linear (on a semilog-scale) fit
        def func(x, a, b):
            return a * np.exp(-x / b)

        # only look at the initial decay
        linear_start, linear_end = 20, 200
        inds = np.where((data.index >= linear_start) & (data.index <= linear_end))[0]
        x = data.index[inds]
        y = data.values[inds] - asymptote
        popt, pcov = curve_fit(func, x, y)

        # plot the data
        axes[i].scatter(
            data.index,
            data.values - asymptote,
            s=1,
            alpha=0.5,
            label="group " + group[i],
        )

        # plot the linear fit
        axes[i].plot(
            data.index,
            func(data.index, *popt),
            color="C1",
            # "k-",
            alpha=0.8,
            label=r"$y = %5.2f e^{-x/%5.0f}$" % tuple(popt),
            linewidth=2,
        )

        # figure formatting
        axes[i].set_ylim(bottom=0.01, top=1)
        axes[i].set_xlim([-30, 900])
        axes[i].set_yscale("log")
        axes[i].axvline(
            halfway_dict[group[i]],
            linestyle="--",
            linewidth=1,
            color="C2",
            label="half-distance = {0:0.0f}".format(halfway_dict[group[i]]),
        )
        axes[i].legend(handletextpad=0.2, handlelength=0.5, markerscale=2)
        axes[i].tick_params(axis="y", labelrotation=90)
        axes[i].set_ylabel(r"linkage ($r$) - asymptote")
        axes[i].set_xlabel("core genome codon distance")
    axes[2].tick_params(axis="y", labelrotation=90)
    axes[2].legend()
    axes[2].set_xlim([100, 475])
    label_list, x_loc = ["a)", "b)", "c)"], [-32, -32, -32]
    y_loc = [-2, -2, -2]
    for i in range(3):
        axes[i].text(
            0.0,
            1.0,
            label_list[i],
            transform=(
                axes[i].transAxes
                + ScaledTranslation(x_loc[i] / 72, y_loc[i] / 72, fig.dpi_scale_trans)
            ),
            va="bottom",
        )

    # save figure
    fig.savefig(args.output, dpi=500)
    fig.savefig(args.output_pdf, dpi=500)


if __name__ == "__main__":
    main()
