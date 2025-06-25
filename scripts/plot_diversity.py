"""
This script creates paper Figure 1, showing the spread of diversity of the groups.

Author: Jemma M. Fendley
"""

import matplotlib.pyplot as plt, seaborn as sns
import pandas as pd, argparse, numpy as np
import parameters  # preset Matplotlib formatting


def main():
    parser = argparse.ArgumentParser(
        description="creates a figure showing group diversity"
    )
    parser.add_argument(
        "-o", "--output", help="figure png file", type=str, required=True
    )
    parser.add_argument(
        "-p", "--output_pdf", help="figure pdf file", type=str, required=True
    )
    parser.add_argument(
        "-i",
        "--input",
        help="TSV file with mean pairwise metrics for all groups",
        type=str,
        required=True,
    )
    args = parser.parse_args()

    # load the TSV file
    df = pd.read_csv(args.input, sep="\t")
    n_groups = len(df)

    # initialize the figure
    fig, axes = plt.subplots(
        nrows=1, ncols=1, figsize=(6.5, 4.65), layout="constrained"
    )

    # extract the statistics of interest
    df["mean pairwise ANI in core genome (\%)"] = 100 - 100 * df["mean_hamming"]
    df["mean pairwise coverage (\%)"] = df["mean_percent_pairwise_coverage"]

    # separate out groups mentioned in the paper (hard-coded)
    select_groups_list = ["A11", "E", "EE", "DE1", "F1"]
    df_select_groups = df[df["name"].isin(select_groups_list)]
    df_other_groups = df[
        df["name"].isin([x for x in df["name"] if x not in select_groups_list])
    ]

    # plot the data separately
    sns.scatterplot(
        data=df_other_groups,
        x="mean pairwise ANI in core genome (\%)",
        y="mean pairwise coverage (\%)",
        alpha=1,
        label="{0:0.0f} total groups".format(n_groups),
        color="C9",
        s=50,
    )

    sns.scatterplot(
        data=df_select_groups,
        x="mean pairwise ANI in core genome (\%)",
        y="mean pairwise coverage (\%)",
        alpha=0.8,
        s=50,
        hue="name",
    )

    # add the species threshold
    cutoffs, labels = [9500, 7000], ["species", "genus"]
    rotation, offset, colors = [-26, -20], [119, -150], ["C5", "C6"]
    # add the species and genus thresholds
    for i, label in enumerate(labels):
        x = np.linspace(int(cutoffs[i] / 100), 100, 1000)
        y = cutoffs[i] / x
        y2 = np.ones(1000) * 100
        axes.plot(x, y, color=colors[i], alpha=0.9, linewidth=1)
        axes.text(
            x[offset[i]],
            y[offset[i]],
            label + " diversity",
            fontsize=8.5,
            rotation=rotation[i],
            rotation_mode="anchor",
        )
        if i == 0:
            axes.fill_between(x, y, y2, color=colors[i], alpha=0.125)
        elif i == 1:
            x1, x2 = np.linspace(70, 95, 100), np.linspace(95, 100, 100)
            axes.fill_between(x1, 7000 / x1, 9500 / x1, color="C6", alpha=0.05)
            axes.fill_between(x2, 7000 / x2, np.ones(100) * 85, color="C6", alpha=0.05)

    x, y, y2 = np.ones(100) * 95, np.linspace(85, 100, 100), np.linspace(50, 85, 100)
    axes.plot(x, y, color="k", linewidth=1, alpha=0.9)
    axes.plot(x, y2, color="k", linewidth=1, alpha=0.5, linestyle="dashed")
    x, x2, y = np.linspace(95, 100, 100), np.linspace(50, 95, 100), np.ones(100) * 85
    axes.plot(x, y, color="k", linewidth=1, alpha=0.9)
    axes.fill_between(x, y, 9500 / x, color="k", alpha=0.05)
    axes.plot(x2, y, color="k", linewidth=1, alpha=0.5, linestyle="dashed")
    axes.text(95.35, 83.75, "vOTU diversity", fontsize=8.5)

    # formating
    axes.legend(
        handlelength=1,
        handletextpad=0.5,
        framealpha=0.98,
        loc="upper left",
        ncol=1,
        columnspacing=0.8,
        borderaxespad=0.3,
    )
    axes.set_xlim([69.5, 100])
    axes.set_ylim([57, 100])

    # save the figure
    fig.savefig(args.output, dpi=450)
    fig.savefig(args.output_pdf, dpi=450)


if __name__ == "__main__":
    main()
