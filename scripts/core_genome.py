"""
This script concatenates core pham alignments into one core genome alignment file.

Author: Jemma M. Fendley
"""

import json, argparse
import numpy as np
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord


def main():
    parser = argparse.ArgumentParser(description="concatenate core genome alignment")
    parser.add_argument(
        "-g", "--group_info", help="JSON file with group info", type=str, required=True
    )
    parser.add_argument(
        "-o", "--output", help="core genome fasta", type=str, required=True
    )
    parser.add_argument(
        "-s", "--starts", help="core pham start positions JSON", type=str, required=True
    )
    parser.add_argument(
        "-r", "--strands", help="core phams strandedness JSON", type=str, required=True
    )
    parser.add_argument(
        "-c", "--core", nargs="+", help="list of core phams", required=True
    )
    parser.add_argument(
        "-a", "--aligns", nargs="+", help=".fna files of each core pham", required=True
    )
    args = parser.parse_args()

    with open(args.group_info) as f:
        G = json.load(f)

    # put the core phams in the correct order
    core_phams = [
        x["pham_ID"] for x in G["paths"][0]["path"] if x["pham_ID"] in args.core
    ]
    phage_IDs = [x["ID"] for x in G["paths"]]

    # record and save the core phams strandedness dictionary
    core_pham_strand = {
        x["pham_ID"]: x["strand"]
        for x in G["paths"][0]["path"]
        if x["pham_ID"] in args.core
    }
    with open(args.strands, "w") as f:
        json.dump(core_pham_strand, f, indent=4)

    # create an array with concatenated alignments of all of the core phams
    n_phages = len(phage_IDs)
    name_list, core_genome_alignment = [], np.empty([n_phages, 0])
    core_pham_starts = {}
    for core in core_phams:
        # find the alignment file and do some quick checks
        core_path = [x for x in args.aligns if "/" + core + "_" in x]
        assert len(core_path) == 1, "Error in core pham files input"
        alignment = AlignIO.read(core_path[0], "fasta")
        assert len(alignment) == n_phages, "Error in alignment files"

        # record the starting position in the alignment of the pham
        core_pham_starts[core] = core_genome_alignment.shape[1]

        # extract the alignment and add to the core genome alignment array
        alignment_array = np.array(alignment)
        core_genome_alignment = np.concatenate(
            (core_genome_alignment, alignment_array), axis=1
        )

        # record the order of the phages in the alignment
        name_list.append([record.id for record in alignment])

    # double check all of the alignments had the phages in the same order
    assert all(
        name_test == name_list[0] for name_test in name_list
    ), "Error with order in alignment files"

    # save the alignment and the core pham starts dictionary
    list_of_alignment = [
        SeqRecord(Seq("".join(row)), id=name_list[0][i])
        for i, row in enumerate(core_genome_alignment)
    ]
    AlignIO.write(MultipleSeqAlignment(list_of_alignment), args.output, "fasta")
    with open(args.starts, "w") as f:
        json.dump(core_pham_starts, f, indent=4)


if __name__ == "__main__":
    main()
