"""
Creates and saves a dataframe (table) with basic information about phage groups

Author: Jemma M. Fendley
"""

import argparse, json
import pandas as pd


# It might be worth updating the config file created from the download-phage pipeline,
#   then this function would not be necessary.
def get_type(type):
    """
    Returns a human-readable "type" for each group
    """
    if type == "irred":
        return "cluster"
    elif type == "subcluster":
        return type
    else:
        raise ValueError("Group is not of correct type.")


def main():
    parser = argparse.ArgumentParser(
        description="builds a table with information about all groups"
    )
    parser.add_argument(
        "-o", "--output", help="TSV file for dataframe", type=str, required=True
    )
    parser.add_argument(
        "-i",
        "--group_info",
        help="group info dictionary (data folder)",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-l",
        "--lifestyle",
        help="TSV file with cluster lifestyles",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-c",
        "--percent_core",
        help="TSV file with percent core",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-m",
        "--mean_pairwise_metrics",
        help="TSV with mean pariwise metrics",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-g", "--groups", nargs="+", required=True, help="summary file for each group"
    )
    parser.add_argument(
        "-n",
        "--nonsyntenic",
        required=True,
        type=str,
        help="nonsynstenic information JSON",
    )
    args = parser.parse_args()

    dict_list = []  # initialize

    # record the key names for looping over
    dictionary_keys = ["morphotypes", "host_genus", "host_species"]
    singular_keys = [
        "most_common_morphotype",
        "most_common_host_genus",
        "most_common_host_species",
    ]
    plural_keys = ["all_morphotypes", "all_host_genera", "all_host_species"]

    for group in args.groups:
        with open(group) as f:
            summary_dict = json.load(f)

        # record basic values from the dictionary
        row = {
            "name": summary_dict["name"],
            "type": get_type(summary_dict["type"]),
            "cluster": "".join(i for i in summary_dict["name"] if not i.isdigit()),
            "n_phages": summary_dict["n_phages"],
            "n_phams": summary_dict["n_phams"],
            "n_core": summary_dict["n_core"],
            "mean_genome_length": round(summary_dict["mean_genome_length"]),
            "mean_n_phams": round(summary_dict["mean_n_phams"], 1),
        }

        # for morphotypes, host genera, and host species, record the most common as
        #   well as a list of all of them
        for i, dictionary_key in enumerate(dictionary_keys):
            ordered_list = [
                y["value"]
                for y in sorted(
                    summary_dict[dictionary_key], key=lambda x: x["n_phage"]
                )
                if y["value"] != ""
            ]  # sort by prevalence
            row[singular_keys[i]] = ordered_list[0]  # most common
            row[plural_keys[i]] = ", ".join(ordered_list)  # list of all
        dict_list.append(row)
    info_df = pd.DataFrame(dict_list)

    # load the lifestyle information and add to the dataframe
    df_lifestyle = pd.read_csv(args.lifestyle, sep="\t")
    df_lifestyle = df_lifestyle.rename(columns={"Life Cycle ": "life_cycle"})
    df_lifestyle["cluster"] = [x.strip() for x in df_lifestyle["Cluster "]]
    lifestyle_dict = pd.Series(
        df_lifestyle.life_cycle.values, index=df_lifestyle.cluster
    ).to_dict()
    info_df["life_cycle"] = info_df["cluster"].map(lifestyle_dict)

    # load the percent core and percent genes info and add to the dataframe
    df_percent_core = pd.read_csv(args.percent_core, sep="\t")
    for x in ["percent_core", "percent_genes"]:
        df_mean = df_percent_core.groupby("name")[x].mean().round(2)
        mean_dict = pd.Series(df_mean.values, index=df_mean.index).to_dict()
        info_df["mean_" + x] = info_df["name"].map(mean_dict)

    # load the nonsyntenic information and record in the dataframe
    with open(args.nonsyntenic) as f:
        nonsyntenic_information = json.load(f)

    def return_synteny(group):
        """
        parses the nonsyntenic information dictionary to return the synteny of a groups
        """
        nonsyntenic_groups = [
            x["name"] for x in nonsyntenic_information["nonsyntenic_groups"]
        ]
        if group in nonsyntenic_groups:
            group_dict = [
                x
                for x in nonsyntenic_information["nonsyntenic_groups"]
                if x["name"] == group
            ]
            assert len(group_dict) == 1, "Error in nonsyntenic JSON file"
            return group_dict[0]["synteny"]
        else:
            return "syntenic"

    info_df["synteny"] = info_df["name"].apply(return_synteny)

    # load the mean pairwise metric dataframe and add it to the main dataframe
    pairwise_metrics_df = pd.read_csv(args.mean_pairwise_metrics, sep="\t").round(4)
    df = pd.merge(info_df, pairwise_metrics_df, on="name")
    df["mean_percent_pairwise_coverage"] = df["mean_percent_pairwise_coverage"].round(2)

    # reorder the columns for display purposes
    column_order = [
        "name",
        "type",
        "cluster",
        "life_cycle",
        "n_phages",
        "n_phams",
        "n_core",
        "synteny",
        "mean_genome_length",
        "mean_n_phams",
        "mean_percent_core",
        "mean_percent_genes",
        "mean_percent_pairwise_coverage",
        "mean_hamming",
        "mean_jaccard",
        "std_hamming",
        "std_jaccard",
        "most_common_morphotype",
        "all_morphotypes",
        "most_common_host_genus",
        "all_host_genera",
        "most_common_host_species",
        "all_host_species",
    ]
    df = df[column_order]

    # save the dataframe
    df.to_csv(args.output, index=False, sep="\t")


if __name__ == "__main__":
    main()
