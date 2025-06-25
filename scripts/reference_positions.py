"""
This script maps a position in the core genome alignment to the actual position in a
    reference phage genome.

Author: Jemma M. Fendley
"""

import json, argparse
import numpy as np
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq

# get_cds returns the CDS feature (gene) of a genbank file corresponding to an aa sequence
from utils import get_cds


def main():
    parser = argparse.ArgumentParser(
        description="find positions of core alignment in reference phage"
    )
    parser.add_argument(
        "-o", "--output", help="output JSON file", type=str, required=True
    )
    parser.add_argument(
        "-r", "--reference_phage", help="phage gbk file", type=str, required=True
    )
    parser.add_argument(
        "-c",
        "--core_genome",
        help="core genome alignment file",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-s", "--pham_starts", help="pham start positions JSON", type=str, required=True
    )
    parser.add_argument(
        "-d", "--strands", help="core phams strandedness JSON", type=str, required=True
    )
    args = parser.parse_args()

    # load the starting position of every pham in the alignment
    with open(args.pham_starts) as f:
        pham_starts_dict = json.load(f)
    pham_starts = np.array(list(pham_starts_dict.values()))
    core_phams = list(pham_starts_dict.keys())

    # find the alignment corresponding to the phage
    accession = args.reference_phage.split("/")[-1].split(".")[0]
    genome_record = SeqIO.read(args.reference_phage, "genbank")
    alignment = AlignIO.read(args.core_genome, "fasta")
    phage_alignment = np.array(
        [record for record in alignment if record.id == accession][0]
    )

    # load the core phams standedness dictionary
    with open(args.strands) as f:
        pham_strands_dict = json.load(f)

    # iterate through every core pham
    pham_positions_map, reference_pham_starts = [], []
    for core_index, pham_start in enumerate(pham_starts):
        if core_index < len(pham_starts) - 1:
            pham_end = pham_starts[core_index + 1]
        else:
            pham_end = len(phage_alignment)

        # extract the part of the core genome alignment of the pham
        pham_alignment = phage_alignment[pham_start:pham_end]

        # find the dna sequence to translate
        if pham_strands_dict[core_phams[core_index]] == "+":
            dna_sequence = Seq(
                "".join(
                    np.delete(pham_alignment, np.where(pham_alignment == "-"))[3:-3]
                )
            )
        elif pham_strands_dict[core_phams[core_index]] == "-":
            dna_sequence = Seq(
                "".join(
                    np.delete(pham_alignment, np.where(pham_alignment == "-"))[3:-3]
                )
            ).reverse_complement()

        # find the starting position of the pham in the reference phage
        aa_sequence = "M" + str(dna_sequence.translate(table=11))
        cds_feature = get_cds(genome_record, aa_sequence)
        reference_pham_starts.append(cds_feature.location.start)

        # make a dictionary of pham position in alignment to reference position
        # handle gaps by mapping them to the last nongap position
        position_dict, last_nongap = {}, 0
        for i in range(len(pham_alignment)):
            if pham_alignment[i] != "-":
                position_dict[i] = i - np.count_nonzero(pham_alignment[:i] == "-")
                last_nongap = position_dict[i]
            else:
                position_dict[i] = last_nongap
        pham_positions_map.append(position_dict)

    # put it all together and create a map of core genome alignment position to reference pham position
    final_dict = {}
    for i in range(len(phage_alignment)):
        pham_index = np.max(np.where(pham_starts <= i))
        relative_position = i - pham_starts[pham_index]
        reference_position = pham_positions_map[pham_index][relative_position]
        reference_pham_start = reference_pham_starts[pham_index]
        final_dict[i] = reference_pham_start + reference_position

    # save the dictionary
    with open(args.output, "w") as f:
        json.dump(final_dict, f, indent=4)


if __name__ == "__main__":
    main()
