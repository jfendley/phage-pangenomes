"""
This script processes the output of DefenseFinder to record which phams have defense/antidefense function.

DefenseFinder: https://github.com/mdmparis/defense-finder

Author: Jemma M. Fendley
"""

import argparse
from Bio import SeqIO
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="process DefenseFinder output")
    parser.add_argument(
        "-o",
        "--output",
        help="TSV file with pham and defense/antidefense function",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-g", "--genes", help="all genes FASTA", type=str, required=True
    )
    parser.add_argument(
        "-d",
        "--defense",
        help="DefenseFinder output, ends in genes.tsv",
        type=str,
        required=True,
    )

    args = parser.parse_args()

    # extract the pham for each gene from the all genes fasta
    pham_dict = {
        x.id: x.description.split(" pham_")[1] for x in SeqIO.parse(args.genes, "fasta")
    }

    # load the DefenseFinder results and add a column for phams
    df = pd.read_csv(args.defense, sep="\t")
    df["pham_ID"] = df["hit_id"].map(pham_dict)

    # Check that each pham is mapped to either a Defense or AntiDefense system, not both.
    # The "activity" column options are "Defense" or "AntiDefense"
    counts = df.groupby("pham_ID")["activity"].nunique()
    assert (
        len(set(counts)) == 1
    ), "Modify code to handle phams with Defense and AntiDefense hits"

    # Record the type of system, as well if it's Defense or AntiDefense, then save the dataframe
    pham_to_type = df.groupby("pham_ID").agg(
        {"type": lambda x: list(set(x))[0], "activity": lambda y: list(set(y))[0]}
    )
    pham_to_type.to_csv(args.output, index=True, sep="\t")


if __name__ == "__main__":
    main()
