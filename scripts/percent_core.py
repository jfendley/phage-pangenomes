"""
Calculate the percent of each phage (by length) that is genes (coding regions)
    and that is core with respect to the group

Author: Jemma M. Fendley
"""

import argparse, json
import pandas as pd
import json


def main():
    parser = argparse.ArgumentParser(
        description="calculates perecent (length) core and genes of each genome "
    )
    parser.add_argument(
        "-o", "--output", help="name of output file", type=str, required=True
    )
    parser.add_argument(
        "-s", "--summary", help="group summary JSON file", type=str, required=True
    )
    parser.add_argument(
        "-p",
        "--phages",
        nargs="+",
        help="list of phage gene files",
        required=True,
    )

    args = parser.parse_args()

    with open(args.summary) as f:
        summary_dict = json.load(f)
    genome_length_dict = summary_dict["genome_lengths"]
    core_phams = summary_dict["core"]

    dict_list = []  # initialize
    for phagefile in args.phages:
        phage_ID = phagefile.split("/")[-1].split("_")[0]
        row = {"name": summary_dict["name"], "phage_ID": phage_ID}  # start dictionary

        with open(phagefile) as f:
            phagesdict = json.load(f)
        # extract the list of all genes and the list of core genes
        gene_list = phagesdict["gene_list"]
        core_list = [x for x in gene_list if x["phams"][0] in core_phams]

        # will iterate over both the gene list and the core list
        lists = [core_list, gene_list]
        keys = ["percent_core", "percent_genes"]
        for index, current_list in enumerate(lists):
            starts = [x["Start"] for x in current_list]  # start positions
            ends = [x["Stop"] for x in current_list]  # end positions
            total_length, initial_start = 0, starts[0]  # initialize

            # The following calculates the total length in core or all genes,
            #   making sure not to overcount overlapping genes
            for i in range(len(ends) - 1):
                if starts[i + 1] > ends[i]:  # if no overlap, otherwise go to next one
                    total_length += ends[i] - initial_start
                    initial_start = starts[i + 1]
            total_length += ends[-1] - initial_start
            row[keys[index]] = 100 * total_length / genome_length_dict[phage_ID]
        dict_list.append(row)

    # save the dataframe
    df = pd.DataFrame(dict_list)
    df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
