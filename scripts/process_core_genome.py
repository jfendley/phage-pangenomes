"""
This script modifies the core genome alignment in regions with multiple haplotypes.
    The modified core genome file is then used in a linkage disequilibrium analysis which
    only considers SNPs which are deviations from the haplotype consensus sequences.

Author: Jemma M. Fendley
"""

import numpy as np, argparse
from Bio import AlignIO
from scipy.cluster.hierarchy import fcluster
from collections import defaultdict
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from utils import sort_matrix  # sorts a matrix by hierarchically clustering
from utils import find_snp_positions  # finds snp positions
from utils import calculate_hamming  # calculates hamming distance


def main():
    parser = argparse.ArgumentParser(
        description="modify the core genome alignment in haplotype region"
    )
    parser.add_argument(
        "-o",
        "--output",
        help="modified core genome alignment file",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-i", "--input", help="core genome alignment file", type=str, required=True
    )
    parser.add_argument(
        "-s", "--start", help="region start positions", nargs="+", required=True
    )
    parser.add_argument(
        "-e", "--end", help="region end positions", nargs="+", required=True
    )
    parser.add_argument(
        "-n", "--n_haplotypes", help="n. haplotypes in region", nargs="+", required=True
    )

    args = parser.parse_args()

    # load the alignment and extract the list of phages
    alignment = AlignIO.read(args.input, "fasta")
    name_list, core_genome_alignment = [x.id for x in alignment], np.array(alignment)
    n_phages = core_genome_alignment.shape[0]

    # find the snp positions and convert the input into lists of integer
    snp_positions = find_snp_positions(core_genome_alignment)[0]
    starts, ends = [int(x) for x in args.start], [int(x) for x in args.end]
    n_haplotypes_list = [int(x) for x in args.n_haplotypes]

    # iterate through the haplotype regions
    i_lower = np.tril_indices(n_phages, -1)
    for index, n_haplotypes in enumerate(n_haplotypes_list):
        # use the distance matrix in the haplotype region to find the haplotypes
        start, end = starts[index], ends[index]
        relevant_alignment = core_genome_alignment[:, start:end]
        distance_matrix = np.array(
            [
                [
                    calculate_hamming(relevant_alignment, i, j) if j > i else 0
                    for j in range(n_phages)
                ]
                for i in range(n_phages)
            ]
        )
        distance_matrix[i_lower] = distance_matrix.T[i_lower]
        res_linkage = sort_matrix(distance_matrix, "complete")[2]
        grouping = fcluster(res_linkage, n_haplotypes, criterion="maxclust")

        # find the snps in the haplotype region to process
        snps_to_process = [x for x in snp_positions if x >= start and x < end]

        # iterate through all of the snp positions and update the alignment
        consensus_phages_dict = defaultdict(list)
        for site in snps_to_process:
            alignment_column = core_genome_alignment[:, site]
            for i in range(1, n_haplotypes + 1):
                # find the phages in the haplotype and the consensus allele
                phages_in_haplotype = np.where(grouping == i)[0]
                alleles, counts = np.unique(
                    alignment_column[phages_in_haplotype], return_counts=True
                )
                haplotype_consensus_allele = alleles[np.argmax(counts)]

                # only modify the allele if not a gap. If a gap, leave it as is and likely the site will not
                #   be considered in the linkage analysis.
                if haplotype_consensus_allele != "-":
                    assert haplotype_consensus_allele != "-"
                    # record the phages in this haplotype that have the consensus allele
                    consensus_phages = [
                        x
                        for x in phages_in_haplotype
                        if core_genome_alignment[x, site] == haplotype_consensus_allele
                    ]
                    consensus_phages_dict[site].extend(consensus_phages)

            # for the phages with the consensus allele, change it to "P" (for "processed")
            #   Then, when linkage disequilibrium is calculated, it will only consider SNPS
            #   that are deviations from the consensus sequence (i.e. not "P")
            for phage in consensus_phages_dict[site]:
                core_genome_alignment[phage, site] = "P"

    # save the new core genome alignment
    list_of_alignment = [
        SeqRecord(Seq("".join(row)), id=name_list[i])
        for i, row in enumerate(core_genome_alignment)
    ]
    AlignIO.write(MultipleSeqAlignment(list_of_alignment), args.output, "fasta")


if __name__ == "__main__":
    main()
