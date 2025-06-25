"""
This script creates a FASTA file of all of the occurences of a given pham in a group.

Author: Jemma M. Fendley
"""

import argparse, json


def main():
    parser = argparse.ArgumentParser(description="create FASTA file from pham data")
    parser.add_argument(
        "-o", "--output_file", help="output FASTA file", type=str, required=True
    )
    parser.add_argument("-p", "--pham", help="pham ID", type=str, required=True)
    parser.add_argument(
        "-l", "--list", nargs="+", help="list of phage gene files", required=True
    )

    args = parser.parse_args()

    for phage_file in args.list:
        phage_ID = phage_file.split("/")[-1].split("_")[0]
        with open(phage_file) as f:
            phage = json.load(f)

        # find the genes that match the input pham ID
        subset = [y for y in phage["gene_list"] if y["phams"][0] == args.pham]

        # write to the fasta file
        for z in subset:
            with open(args.output_file, "a") as fa:
                fa.write(
                    ">" + phage_ID + " " + z["GeneID"] + "\n" + z["translation"] + "\n"
                )
            # Note that the previous step adds to an existing file. In the workflow, before running this
            #   script, any file with the same name is deleted to prevent problems


if __name__ == "__main__":
    main()
