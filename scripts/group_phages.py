"""
This script creates a dictionary of groups to the list of phages they contain.

Author: Jemma M. Fendley
"""

import argparse, json, subprocess, time
from Bio import SeqIO
from tqdm import tqdm
import numpy as np
from utils import (
    get_cds,
)  # This function returns the CDS feature (gene) of a genbank file corresponding to an amino acid sequence
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(
        description="creates the group to phages dictionary"
    )
    parser.add_argument(
        "-i",
        "--input",
        help="original config file (in data folder)",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="output group to phages dictionary",
        type=str,
        required=True,
    )
    args = parser.parse_args()

    with open(args.input) as f:
        group_info = json.load(f)

    data_dir = args.input.split("/")[0]  # input needs to be in the data directory

    inconsistent_phage_count = 0  # Initialize
    group_phages_dict = defaultdict(list)

    for group in tqdm(group_info):
        # Exclude singletons (phages not assigned to a cluster) as well as clusters if they
        #   have been broken into subclusters (named a "super cluster")
        if group["type"] == "super_cluster" or group["name"] == "Singleton":
            continue

        # Exclude groups with fewer than 10 phages
        if len(group["phage_list"]) < 10:
            continue

        for phage in group["phage_list"]:
            all_genes_found = True  # Initialize

            # This is the gene file from the Actinobacteriophage Database (phagesdb.org)
            gene_file = data_dir + "/phages/" + phage + "/" + phage + "_genes.json"
            with open(gene_file) as f:
                phage_genes = json.load(f)

            # This is the file from GenBank
            gbk_file = data_dir + "/phages/" + phage + "/" + phage + ".gbk"
            try:
                genome_record = SeqIO.read(gbk_file, "genbank")
            except ValueError:
                print("Attempting to download GenBank file again")
                time.sleep(1)  # to not accidentally exceed GenBank API limit
                return_code = subprocess.call(
                    "ncbi-acc-download --out " + gbk_file + " " + phage,
                    shell=True,
                )
                genome_record = SeqIO.read(gbk_file, "genbank")

            all_translations = []
            for gene in phage_genes["gene_list"]:
                # The index records which number appearance the gene is, in case there are multiple genes with the same amino acid sequence
                index = all_translations.count(gene["translation"])
                cds_feature = get_cds(genome_record, gene["translation"], index)

                if cds_feature is None:
                    all_genes_found = (
                        False  # This gene was not found in the GenBank file
                    )
                    inconsistent_phage_count += 1
                    print("Phage files do not match for: ", phage)
                    print("Removing phage from analysis")
                    # If the genes in the database gene file don't match the GenBank file, we remove the phage from our analysis
                    break
                all_translations.append(gene["translation"])

            if all_genes_found:
                group_phages_dict[group["name"]].append(phage)

    group_phages = {
        x: group_phages_dict[x]
        for x in list(group_phages_dict.keys())
        if len(group_phages_dict[x]) >= 10
    }  # Remove any groups that now contain fewer than 10 phages

    print(inconsistent_phage_count, " phages removed from analysis")

    # Double check that no phage is in multiple groups
    all_phages = np.array([x for y in list(group_phages.values()) for x in y])
    assert len(all_phages) == len(set(all_phages))

    print("Number of phages to analyze: ", len(all_phages))
    print("Number of groups to analyze: ", len(group_phages))

    with open(args.output, "w") as f:
        json.dump(group_phages, f, indent=4)


if __name__ == "__main__":
    main()
