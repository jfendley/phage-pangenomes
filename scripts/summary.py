"""
This script saves a JSON file with basic information about the group needed for downstream analysis.

Author: Jemma M. Fendley
"""

import numpy as np
import argparse, json


def main():
    parser = argparse.ArgumentParser(
        description="create JSON file with basic group information"
    )
    parser.add_argument(
        "-o", "--output_file", help="name of JSON file", type=str, required=True
    )
    parser.add_argument(
        "-g,", "--group_json", help="group JSON file", type=str, required=True
    )
    parser.add_argument(
        "-p",
        "--phages",
        nargs="+",
        help="phage info files",
        required=True,
    )

    args = parser.parse_args()

    final_dict = {}  # Initialize dictionary to save

    with open(args.group_json) as f:
        G = json.load(f)

    final_dict["name"], final_dict["type"] = G["name"], G["type"]
    final_dict["n_phages"], final_dict["n_phams"] = len(G["paths"]), len(G["phams"])

    # count the number of phams per phage
    n_phams = [len(x["path"]) for x in G["paths"]]

    # extract the core phams
    core_phams = [
        x["pham_ID"]
        for x in G["phams"]
        if x["total_count"] == len(G["paths"]) and x["phage_count"] == len(G["paths"])
    ]
    final_dict["n_core"], final_dict["core"] = len(core_phams), core_phams

    # initializing lists and dictionary
    morphotype_list, isolation_host_genus, genome_lengths = [], [], {}
    isolation_host_species, lysogeny_notes = [], []

    for phage_file in args.phages:
        phage_ID = phage_file.split("/")[-1].split("_")[0]
        assert phage_ID in [x["ID"] for x in G["paths"]], "Error: phage not found"

        with open(phage_file) as f:
            phage_dict = json.load(f)
            # record information about the phage
            lysogeny_notes.append(phage_dict["lysogeny_notes"])
            morphotype_list.append(phage_dict["morphotype"])
            isolation_host_genus.append(phage_dict["isolation_host"]["genus"])
            isolation_host_species.append(phage_dict["isolation_host"]["species"])
            genome_lengths[phage_ID] = int(phage_dict["genome_length"])

    # extract all morphotypes and respective counts
    morphotypes, morphotype_counts = np.unique(morphotype_list, return_counts=True)
    morphotype_dict_list = [
        {"value": morphotypes[i], "n_phage": int(morphotype_counts[i])}
        for i in range(len(morphotypes))
    ]

    # extract all possible lysogeny notes and respective counts
    lysogeny, lysogeny_counts = np.unique(lysogeny_notes, return_counts=True)
    lysogeny_dict_list = [
        {"lysogeny_notes": lysogeny[i], "n_phage": int(lysogeny_counts[i])}
        for i in range(len(lysogeny))
    ]

    # extract all genera and respective counts
    genus, genus_counts = np.unique(isolation_host_genus, return_counts=True)
    genus_dict_list = [
        {"value": genus[i], "n_phage": int(genus_counts[i])} for i in range(len(genus))
    ]

    # extract all species and respective counts
    species, species_counts = np.unique(isolation_host_species, return_counts=True)
    species_dict_list = [
        {"value": species[i], "n_phage": int(species_counts[i])}
        for i in range(len(species))
    ]

    # load the information into the dictionary
    final_dict["morphotypes"] = morphotype_dict_list
    final_dict["genome_lengths"] = genome_lengths
    final_dict["mean_genome_length"] = np.mean(list(genome_lengths.values()))
    final_dict["lysogeny"] = lysogeny_dict_list
    final_dict["host_species"] = species_dict_list
    final_dict["host_genus"] = genus_dict_list
    final_dict["mean_n_phams"] = np.mean(n_phams)

    # save the dictionary
    with open(args.output_file, "w") as f:
        json.dump(final_dict, f, indent=4)


if __name__ == "__main__":
    main()
