"""
This script creates a plot for the SI that shows how diversity in groups affects synteny.

Note: the jitter in the seabron striplot is random, and so the figure created may not
    perfectly match the one in the paper.

Author: Jemma M. Fendley
"""

import pandas as pd, argparse, json
import matplotlib.pyplot as plt, seaborn as sns
import parameters  # preset matplotlib formatting


def main():
    parser = argparse.ArgumentParser(description="plot synteny and diversity figure")
    parser.add_argument(
        "-o", "--output", help="png file for plot", type=str, required=True
    )
    parser.add_argument(
        "-p", "--output_pdf", help="pdf file for plot", type=str, required=True
    )
    parser.add_argument(
        "-i", "--input", help="basic info table", type=str, required=True
    )
    parser.add_argument(
        "-n", "--nonsyntenic", help="nonsyntenic group info", type=str, required=True
    )
    parser.add_argument(
        "-d",
        "--distances",
        help="pairwise distances for each group",
        nargs="+",
        required=True,
    )
    args = parser.parse_args()

    # load basic info dataframe (df)
    df = pd.read_csv(args.input, sep="\t")

    # record a list of all the groups
    all_groups = df["name"].unique()

    # create dictionaries for group to synteny
    df = df.set_index("name")
    group_to_synteny = pd.Series(df.synteny, index=df.index).to_dict()

    # add a hyphen to nonsyntenic in dataframe and dictionary
    group_to_synteny = {
        x: y if y != "nonsyntenic" else "non-syntenic"
        for x, y in group_to_synteny.items()
    }
    df["synteny"] = df["synteny"].apply(
        lambda x: x if x != "nonsyntenic" else "non-syntenic"
    )

    # load the file that contains the non-syntenic group information
    with open(args.nonsyntenic) as f:
        non_syntenic_groups = json.load(f)

    # initialize and list cyclic groups
    all_pairs_dict_list, non_syntenic_groups_list = [], []
    cyclic_groups = [
        x["name"]
        for x in non_syntenic_groups["nonsyntenic_groups"]
        if x["synteny"] == "cyclic"
    ]

    # cycle through nonsyntenic (incl. cyclic) groups
    for group_info in non_syntenic_groups["nonsyntenic_groups"]:
        group_name = group_info["name"]
        non_syntenic_groups_list.append(group_name)
        non_syntenic_phages = group_info["nonsyntenic_phages"]

        # load the pairwise distance metrics
        pw_data_file = [
            x for x in args.distances if f"/{group_name}_pairwise_metrics.tsv" in x
        ]
        assert len(pw_data_file) == 1
        pw_df = pd.read_csv(pw_data_file[0], sep="\t")

        for phage in list(set(pw_df["phage_1"]).union(set(pw_df["phage_2"]))):

            # calculate the mean distance between this phage and everyone else in the group
            mean_hamming = pw_df[
                (pw_df["phage_1"] == phage) | (pw_df["phage_2"] == phage)
            ]["Hamming"].mean()
            row = {"group": group_name, "phage": phage, "mean_hamming": mean_hamming}

            # record if the phage's core gene order does not reflect the consensus
            if group_name in cyclic_groups:
                row["non_syntenic"] = (
                    "cyclic" if phage in non_syntenic_phages else "syntenic"
                )
            else:
                row["non_syntenic"] = (
                    "non-syntenic" if phage in non_syntenic_phages else "syntenic"
                )
            all_pairs_dict_list.append(row)

    # cycle through the rest of the groups
    for group in [x for x in list(all_groups) if x not in non_syntenic_groups_list]:
        # load the pairwise distance metrics
        pw_data_file = [
            x for x in args.distances if f"/{group}_pairwise_metrics.tsv" in x
        ]
        assert len(pw_data_file) == 1
        pw_df = pd.read_csv(pw_data_file[0], sep="\t")

        # record the mean distance between each phage and all other phages
        for phage in list(set(pw_df["phage_1"]).union(set(pw_df["phage_2"]))):
            mean_hamming = pw_df[
                (pw_df["phage_1"] == phage) | (pw_df["phage_2"] == phage)
            ]["Hamming"].mean()
            row = {
                "group": "syntenic \n groups",
                "phage": phage,
                "mean_hamming": mean_hamming,
            }
            row["non_syntenic"] = "syntenic"
            all_pairs_dict_list.append(row)

    # save to dataframe, sort and separate
    pairs_df = pd.DataFrame(all_pairs_dict_list)
    pairs_df = pairs_df.sort_values(["group", "non_syntenic"], ascending=[True, False])
    syntenic_df = pairs_df[pairs_df["group"] == "syntenic \n groups"]
    non_syntenic_df = pairs_df[pairs_df["group"] != "syntenic \n groups"]

    # initalize figure
    fig, axes = plt.subplots(
        nrows=1,
        ncols=1,
        figsize=(5, 4),
        layout="constrained",
    )
    axes = [axes]
    palette = {
        "syntenic": "C0",
        "non-syntenic": "C1",
        "cyclic": "C2",
    }

    # plot each of the nonsyntenic groups separately
    sns.stripplot(
        data=non_syntenic_df[non_syntenic_df["non_syntenic"] == "syntenic"],
        x="group",
        y="mean_hamming",
        hue="non_syntenic",
        hue_order=["syntenic"],
        alpha=0.75,
        ax=axes[0],
        palette=palette,
    )
    sns.stripplot(
        data=non_syntenic_df[non_syntenic_df["non_syntenic"] == "non-syntenic"],
        x="group",
        y="mean_hamming",
        hue="non_syntenic",
        hue_order=["non-syntenic"],
        alpha=0.75,
        ax=axes[0],
        palette=palette,
        marker="s",
    )
    sns.stripplot(
        data=non_syntenic_df[non_syntenic_df["non_syntenic"] == "cyclic"],
        x="group",
        y="mean_hamming",
        hue="non_syntenic",
        hue_order=["cyclic"],
        alpha=0.75,
        ax=axes[0],
        palette=palette,
        marker="^",
        s=6,
    )

    # plot the syntenic groups
    sns.stripplot(
        data=syntenic_df,
        x="group",
        y="mean_hamming",
        alpha=0.4,
        ax=axes[0],
        s=2,
    )
    sns.boxplot(
        data=syntenic_df,
        x="group",
        y="mean_hamming",
        boxprops=dict(alpha=0.2),
        color="C0",
        ax=axes[0],
        width=0.9,
        showfliers=False,
    )

    # figure formatting
    axes[0].tick_params("x", rotation=90)
    axes[0].legend(title="")
    axes[0].set_ylabel("per phage: mean distance to other phages")
    axes[0].set_ylim([0, 0.46])
    axes[0].set_ylabel("mean p.w. core distance per phage")
    axes[0].set_xlabel("group", labelpad=-12)

    # save figure
    fig.savefig(args.output, dpi=500)
    fig.savefig(args.output_pdf, dpi=500)


# %%

if __name__ == "__main__":
    main()
