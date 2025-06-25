"""
This scripts compiles a list of all the genes across all of the phages,
    to be used as input into DefenseFinder.

Author: Jemma M. Fendley
"""

from Bio import SeqIO
import json, argparse
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def main():
    parser = argparse.ArgumentParser(description="create FASTA containing all genes")
    parser.add_argument(
        "-o", "--output", help="FASTA with all genes", type=str, required=True
    )
    parser.add_argument(
        "-p", "--phages", help="list of phage gene files", nargs="+", required=True
    )
    args = parser.parse_args()

    list_of_sequences = []  # initialize

    # loop through all of the phage gene files
    for gene_file in args.phages:
        phage_ID = gene_file.split("/")[-1].split("_")[0]  # extract the phage ID
        with open(gene_file) as f:
            phage_json = json.load(f)

        # loop through all of the genes and add to the list
        for i, gene in enumerate(phage_json["gene_list"]):
            list_of_sequences.append(
                SeqRecord(
                    Seq(gene["translation"]),
                    id=phage_ID + "_" + str(i),
                    description=gene["GeneID"] + " pham_" + gene["phams"][0],
                )
            )

    # save to the file
    SeqIO.write(list_of_sequences, args.output, "fasta")


if __name__ == "__main__":
    main()
