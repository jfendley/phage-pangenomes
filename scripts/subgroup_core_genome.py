"""
This script extracts the core genome alignment for only the phages that belong to a subgroup

Author: Jemma M. Fendley
"""

import argparse
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment


def main():
    parser = argparse.ArgumentParser(
        description="restrict core genome alignment to subgroup"
    )
    parser.add_argument(
        "-o",
        "--output",
        help="subgroup core genome alignment file",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-i",
        "--core_genome",
        help="core genome alignment file",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-p", "--phages", help="list of phages in subgroup", nargs="+", required=True
    )
    args = parser.parse_args()

    # load the core genome alignment
    core_genome_alignment = AlignIO.read(args.core_genome, "fasta")

    # extract the rows for phages that are in the subgroup
    subgroup_alignment = [
        record for record in core_genome_alignment if record.id in args.phages
    ]

    # save the new alignment
    AlignIO.write(MultipleSeqAlignment(subgroup_alignment), args.output, "fasta")


if __name__ == "__main__":
    main()
