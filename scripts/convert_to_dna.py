"""
This script creates a DNA alignment file (FASTA) from an amino acid (aa) alignment (FASTA) and GenBank files.

Author: Jemma M. Fendley
"""

from Bio import SeqIO, AlignIO
import argparse, numpy as np

# get_cds returns the CDS feature (gene) of a genbank file corresponding to an aa sequence
from utils import get_cds


def main():
    parser = argparse.ArgumentParser(
        description="converts aa FASTA file to dna FASTA file"
    )
    parser.add_argument("-i", "--input", help="aa fasta file", type=str, required=True)
    parser.add_argument(
        "-o", "--output", help="output dna fasta file", type=str, required=True
    )
    parser.add_argument(
        "-g",
        "--gbk",
        nargs="+",
        required=True,
        help="GenBank files for all phages in the group",
    )
    args = parser.parse_args()

    aa_alignment = AlignIO.read(args.input, "fasta")

    for alignment_row in aa_alignment:
        accession = alignment_row.id  # phage accession number (ID)
        alignment_row_array = np.array(alignment_row)
        description = alignment_row.description

        file = [i for i in args.gbk if accession in i]
        assert len(file) == 1, "Multiple GenBank files with the same accession"

        genome_record = SeqIO.read(file[0], "genbank")

        # This function returns the cds feature that matches the sequence of the gene of interest
        cds_feature = get_cds(genome_record, alignment_row.seq.replace("-", ""))

        gene_dna_sequence = cds_feature.extract(genome_record.seq)

        # Double check that the dna sequence translates to the aa sequence
        aa_sequence = gene_dna_sequence.translate(cds=True, table=11)
        assert aa_sequence == alignment_row.seq.replace("-", ""), "Incorrect sequence"
        assert len(gene_dna_sequence) % 3 == 0, "Gene length is not a multiple of 3"

        strand = cds_feature.location.strand
        if strand == 1:
            # Split the DNA seq into codons
            codons = [
                gene_dna_sequence[3 * y : 3 * y + 3]
                for y in range(int(len(gene_dna_sequence) / 3))
            ]

            # Convert the aa alignment into a dna alignment
            dna_align = [
                (
                    "---"
                    if alignment_row_array[y] == "-"
                    else str(
                        codons[y - len(np.where(alignment_row_array[:y] == "-")[0])]
                    )
                )
                for y in range(len(alignment_row_array))
            ]
            dna_align.append(str(codons[-1]))  # add the stop codon
        elif strand == -1:
            # take the reverse complement of the dna sequence
            reverse_complement = gene_dna_sequence.reverse_complement()
            assert (
                genome_record.seq[cds_feature.location.start : cds_feature.location.end]
                == reverse_complement
            ), "Error: not on reverse strand as expected"

            # Split the DNA seq into codons
            codons = [
                reverse_complement[3 * y : 3 * y + 3]
                for y in range(int(len(reverse_complement) / 3))
            ]

            # Convert the aa alignment into a dna alignment
            alignment_row_reverse = alignment_row_array[::-1]
            dna_align = [
                (
                    "---"
                    if alignment_row_reverse[y] == "-"
                    else str(
                        codons[
                            1 + y - len(np.where(alignment_row_reverse[:y] == "-")[0])
                        ]
                    )
                )
                for y in range(len(alignment_row_reverse))
            ]
            dna_align = [str(codons[0])] + dna_align  # add the stop codon

        dna_align = "".join(dna_align)  # convert list to string

        # save the dna alignment
        with open(args.output, "a") as fa:
            fa.write(">" + description + "\n" + dna_align + "\n")
        # Note that the previous step adds to an existing file. In the workflow, before running this
        #   script, any file with the same name is deleted to prevent problems


if __name__ == "__main__":
    main()
