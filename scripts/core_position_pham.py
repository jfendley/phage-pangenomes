"""
This script creates a table with core genome position to core pham index in the core genome alignment.
    It also records which positions correspond to overlapping phams.

Author: Jemma M. Fendley
"""

import json, argparse
import numpy as np, pandas as pd
from Bio import AlignIO, SeqIO
from tqdm import tqdm

# get_cds returns the CDS feature (gene) of a genbank file corresponding to an aa sequence
from utils import get_cds


def main():
    parser = argparse.ArgumentParser(
        description="creates dictionary of core position to gene"
    )
    parser.add_argument("-o", "--output", help="output table", type=str, required=True)
    parser.add_argument(
        "-i",
        "--core_genome",
        help="core genome alignment file",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-s", "--starts", help="core pham start JSON", type=str, required=True
    )
    parser.add_argument(
        "-n",
        "--nonsyntenic",
        help="nonsyntenic information file",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-r", "--strands", help="core phams strandedness JSON", type=str, required=True
    )
    parser.add_argument(
        "-a", "--alignments", help="core pham alignment files", nargs="+", required=True
    )
    parser.add_argument(
        "-p", "--phages", help="phage .gbk files", nargs="+", required=True
    )
    args = parser.parse_args()

    # record the nonsyntenic phages in the group for the cyclic phages
    with open(args.nonsyntenic) as f:
        nonsyntenic_info = json.load(f)
    group = args.output.split("/")[-1].split("_")[0]
    nonsyntenic = [
        x
        for x in nonsyntenic_info["nonsyntenic_groups"]
        if x["name"] == group and x["synteny"] == "cyclic"
    ]
    if len(nonsyntenic) == 1:
        nonsyntenic_phages = nonsyntenic[0]["nonsyntenic_phages"]
    elif len(nonsyntenic) == 0:
        nonsyntenic_phages = []

    # load the core phams and their starting position in the core genome alignment
    with open(args.starts) as f:
        pham_starts_dict = json.load(f)
    core_phams = list(pham_starts_dict.keys())
    pham_starts = list(pham_starts_dict.values())

    # load the core phams standedness dictionary
    with open(args.strands) as f:
        pham_strands_dict = json.load(f)

    # load the core genome alignment
    core_genome_alignment = np.array(AlignIO.read(args.core_genome, "fasta"))

    # each matrix entry i,j represents the distance between the end of core pham j-1 and
    #   the beginning of core pham j in phage i
    overlap_matrix = np.full((len(args.phages), len(core_phams)), np.nan)
    for i, file in enumerate(tqdm(args.phages)):
        phage_ID = file.split("/")[-1].split(".gbk")[0]
        genome_record = SeqIO.read(file, "genbank")
        gene_end, gene_starts = 0, []
        for j, core in enumerate(core_phams):
            # find the alignment file and load the alignment
            path = [x for x in args.alignments if "/" + core + "_dna" in x]
            assert len(path) == 1, "Error in input alignment files"
            pham_alignment = AlignIO.read(path[0], "fasta")

            # find the sequence that corresponds to the gene in the phage
            sequence = [
                record.seq for record in pham_alignment if record.id.strip() == phage_ID
            ]
            assert len(sequence) == 1, "Error in finding gene in pham alignment file"

            # find the dna sequence to translate
            if pham_strands_dict[core] == "+":
                dna_sequence = sequence[0].replace("-", "")
            elif pham_strands_dict[core] == "-":
                dna_sequence = sequence[0].replace("-", "").reverse_complement()

            # find the gene in the genbank file
            cds_feature = get_cds(
                genome_record,
                dna_sequence.translate(cds=True, table=11),
            )

            # record distance between start of this gene and end of previous gene
            overlap_matrix[i, j] = cds_feature.location.start - gene_end

            # record the gene start in case of cyclic phage
            gene_starts.append(cds_feature.location.start)

            # update the end for the next gene
            gene_end = cds_feature.location.end
        if phage_ID not in nonsyntenic_phages:
            assert all(np.sort(gene_starts) == gene_starts), "Synteny error"
        else:  # need to recalculate the matrix if the phage has a different cylic core pham ordering
            # find the correct core pham ordering and go through the phams in that order
            ordering = np.argsort(gene_starts)
            gene_end, gene_starts2 = 0, []
            for j in ordering:
                # find the alignment file and load the alignment
                core = core_phams[j]
                path = [x for x in args.alignments if "/" + core + "_dna" in x]
                assert len(path) == 1, "Error in input alignment files"
                pham_alignment = AlignIO.read(path[0], "fasta")

                # find the sequence that corresponds to the gene in the phage
                sequence = [
                    record.seq
                    for record in pham_alignment
                    if record.id.strip() == phage_ID
                ]
                assert (
                    len(sequence) == 1
                ), "Error in finding gene in pham alignment file"

                # find the dna sequence to translate
                if pham_strands_dict[core] == "+":
                    dna_sequence = sequence[0].replace("-", "")
                elif pham_strands_dict[core] == "-":
                    dna_sequence = sequence[0].replace("-", "").reverse_complement()

                # find the gene in the genbank file
                cds_feature = get_cds(
                    genome_record,
                    dna_sequence.translate(cds=True, table=11),
                )

                # record distance between start of this gene and end of previous gene
                overlap_matrix[i, j] = cds_feature.location.start - gene_end

                # record the gene starts to double check new order is correct
                gene_starts2.append(cds_feature.location.start)

                # update the end for the next gene
                gene_end = cds_feature.location.end
            assert all(np.sort(gene_starts2) == gene_starts2)

    # find the junctions and phages that contain overlap
    overlap_junctions = np.unique(np.where(overlap_matrix < 0)[1])
    overlap_phages = {
        i: np.where(overlap_matrix[:, i] < 0)[0] for i in overlap_junctions
    }

    # find the phams that contain overlap, splitting into those that overlapping starts and ends
    start_overlap_phams = overlap_junctions
    end_overlap_phams = [
        x - 1 if x > 0 else len(core_phams) - 1 for x in overlap_junctions
    ]
    all_overlap_phams = set(start_overlap_phams).union(set(end_overlap_phams))

    # lengths of the core phams in the alignment
    n_positions = core_genome_alignment.shape[1]
    pham_lengths = [
        pham_starts[i + 1] - pham_starts[i] for i in range(len(pham_starts) - 1)
    ] + [n_positions - pham_starts[-1]]

    # Identify which positions in the core genome alignment correspond to overlap in which phages
    start_overlap_dict, end_overlap_dict = {}, {}
    for i in all_overlap_phams:
        if i in start_overlap_phams:
            positions_overlap = []
            for phage_ind in overlap_phages[i]:
                assert overlap_matrix[phage_ind, i] < 0
                # find the section of pham alignment that contains overlap in this phage
                overlap_amount = int(np.absolute(overlap_matrix[phage_ind, i]))
                alignment = core_genome_alignment[
                    phage_ind, pham_starts[i] : pham_starts[i] + pham_lengths[i]
                ]
                overlap = alignment[:overlap_amount]

                # if there are gaps in the section, need to go further into the alignment
                #   to reach the correct overlap amount
                if "-" in overlap:
                    length_overlap = np.count_nonzero(np.array(overlap) != "-")
                    j = 0
                    while length_overlap < overlap_amount:
                        j += 1
                        overlap = alignment[: overlap_amount + j]
                        length_overlap = np.count_nonzero(np.array(overlap) != "-")
                # record the overlap positions
                positions = np.where(overlap != "-")[0]
                assert len(positions) == overlap_amount
                positions_overlap.append(positions)
            # record all the overlap positions and how many phages have overlap in that position
            all_positions = [x for y in positions_overlap for x in y]
            start_overlap_dict[i] = {
                x: all_positions.count(x) for x in set(all_positions)
            }
        if i in end_overlap_phams:
            positions_overlap = []
            for phage_ind in overlap_phages[i + 1]:
                # find the part of pham alignment that contains overlap in this phage
                overlap_amount = int(np.absolute(overlap_matrix[phage_ind, i + 1]))
                alignment = core_genome_alignment[
                    phage_ind, pham_starts[i] : pham_starts[i] + pham_lengths[i]
                ]
                overlap = alignment[-overlap_amount:]

                # if there are gaps in the section, need to go further into the alignment
                #   to reach the correct overlap amount
                j = 0
                if "-" in overlap:
                    length_overlap = np.count_nonzero(np.array(overlap) != "-")
                    while length_overlap < overlap_amount:
                        j += 1
                        overlap = alignment[-overlap_amount - j :]
                        length_overlap = np.count_nonzero(np.array(overlap) != "-")
                # record the overlap positions
                positions = [
                    x + len(alignment) - overlap_amount - j
                    for x in np.where(overlap != "-")[0]
                ]
                assert len(positions) == overlap_amount
                positions_overlap.append(positions)
            # record all the overlap positions and how many phages have overlap in that position
            all_positions = [x for y in positions_overlap for x in y]
            end_overlap_dict[i] = {
                x: all_positions.count(x) for x in set(all_positions)
            }

    dict_list = []
    # save a row in the table for each position in the core genome file
    for position in tqdm(range(n_positions)):
        row = {"position": position}

        # find the pham in which the position is located in the alignment
        pham = np.max(np.where(np.array(pham_starts) <= position)[0])

        # count the number of phages in which there is overlap in the position
        start, end, overlap_count = False, False, 0
        if pham in all_overlap_phams:
            relative_position = position - pham_starts[pham]
            assert relative_position >= 0, "Error in relative position"
            if pham in start_overlap_phams:
                if relative_position in list(start_overlap_dict[pham].keys()):
                    overlap_count += start_overlap_dict[pham][relative_position]
                    start = True
            if pham in end_overlap_phams:
                if relative_position in list(end_overlap_dict[pham].keys()):
                    overlap_count += end_overlap_dict[pham][relative_position]
                    end = True
        # add to the table
        if start and end:
            row["pham"] = str(pham - 1) + "_" + str(pham) + "_" + str(pham + 1)
        elif start:
            row["pham"] = str(pham - 1) + "_" + str(pham)
        elif end:
            row["pham"] = str(pham) + "_" + str(pham + 1)
        else:
            row["pham"] = str(pham)
        row["overlap"], row["n_overlap"] = start or end, overlap_count
        dict_list.append(row)

    # save the table
    df = pd.DataFrame(dict_list)
    df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
