"""
This script contains any function that is used in multiple scripts.

Author: Jemma M. Fendley, with some code adapted from https://gmarti.gitlab.io/ml/2017/09/07/how-to-sort-distance-matrix.html
"""

from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
import numpy as np


def seriation(Z, N, cur_index):
    """
    from: https://gmarti.gitlab.io/ml/2017/09/07/how-to-sort-distance-matrix.html
        input:
            - Z is a hierarchical tree (dendrogram)
            - N is the number of points given to the clustering process
            - cur_index is the position in the tree for the recursive traversal
        output:
            - order implied by the hierarchical tree Z

        seriation computes the order implied by a hierarchical tree (dendrogram)
    """
    if cur_index < N:
        return [cur_index]
    else:
        left = int(Z[cur_index - N, 0])
        right = int(Z[cur_index - N, 1])
        return seriation(Z, N, left) + seriation(Z, N, right)


def sort_matrix(distmatrix, method="complete"):
    """
    from: https://gmarti.gitlab.io/ml/2017/09/07/how-to-sort-distance-matrix.html

    input:
        - distmatrix is a distance matrix
        - method = ["ward","single","average","complete"]
    output:
        - seriated_dist is the input dist_mat,
        but with re-ordered rows and columns
        according to the seriation, i.e. the
        order implied by the hierarchical tree
        - res_order is the order implied by
        the hierarhical tree
        - res_linkage is the hierarhical tree (dendrogram)

    compute_serial_matrix transforms a distance matrix into
        a sorted distance matrix according to the order implied
        by the hierarchical tree (dendrogram)
    """
    condensed_dist_matrix = squareform(distmatrix)

    res_linkage = linkage(condensed_dist_matrix, method)
    N = len(distmatrix)
    res_order = seriation(res_linkage, N, N + N - 2)
    seriated_dist = np.zeros((N, N))
    a, b = np.triu_indices(N, k=1)
    seriated_dist[a, b] = distmatrix[
        [res_order[i] for i in a], [res_order[j] for j in b]
    ]
    seriated_dist[b, a] = seriated_dist[a, b]
    return seriated_dist, res_order, res_linkage


def convert_strand(string):
    """Converts strandedness from json to GFA format"""
    assert string in ["F", "R"], "Error: strandedness in incorrect format"
    return "+" if string == "F" else "-"


def invert_strand(strand):
    """Inverts the strandedness of a pham"""
    assert strand in ["+", "-"], "Error: strandedness in incorrect format"
    return "+" if strand == "-" else "-"


def create_edge(pham1, pham2):
    """
    input: two pham objects, in the order in which they appear in the genome (i.e. pham1 comes just before pham2)
    output: a string edge identifier such that reverse complements have the same identifier,

    create_edge creates a string from a pham pair where the lower pham comes first.
    """
    if (
        int(pham1["pham_ID"]) == int(pham2["pham_ID"])
        and pham1["strand"] == pham2["strand"]
    ):
        edge = (
            pham1["pham_ID"] + "_+_" + pham2["pham_ID"] + "_+"
        )  # Since A+A+ is equivalent to A-A-
    elif int(pham1["pham_ID"]) <= int(pham2["pham_ID"]):
        edge = (
            pham1["pham_ID"]
            + "_"
            + pham1["strand"]
            + "_"
            + pham2["pham_ID"]
            + "_"
            + pham2["strand"]
        )
    else:
        edge = (
            pham2["pham_ID"]
            + "_"
            + invert_strand(pham2["strand"])
            + "_"
            + pham1["pham_ID"]
            + "_"
            + invert_strand(pham1["strand"])
        )
    return edge


def get_cds(genome_record, translation, index=0):
    """
    Extracts the CDS feature of an GenBank file given an amino acid sequence.
    """
    features = [
        x
        for x in genome_record.features
        if x.type == "CDS"
        and "translation" in list(x.qualifiers.keys())
        and translation in x.qualifiers["translation"]
    ]
    if len(features) >= 1:
        return features[index]
    else:
        return None


def calculate_hamming(alignment, row_1, row_2):
    """
    Calculates the hamming distance of two rows in an alignment
    """
    return 1 - (
        np.sum(alignment[row_1] == alignment[row_2])
        - np.count_nonzero((alignment[row_1] == "-") & (alignment[row_2] == "-"))
    ) / np.count_nonzero((alignment[row_1] != "-") | (alignment[row_2] != "-"))


def find_snp_positions(core_genome_alignment):
    """
    input - core genome alignment in np.array format
    """
    dna_alphabet = ["G", "T", "A", "C"]
    n_phages = core_genome_alignment.shape[0]

    # count the number of alleles (and gaps) at each position in the core genome
    allele_counts = np.array(
        [np.sum(core_genome_alignment == nuc, axis=0) for nuc in dna_alphabet]
    )
    gap_counts = np.array(np.sum(core_genome_alignment == "-", axis=0))

    # find positions with snps
    snp_positions = np.where(
        (np.max(allele_counts, axis=0) > 1)
        & (np.max(allele_counts, axis=0) < n_phages - 1)
        & (np.max(allele_counts, axis=0) > gap_counts)
        & (np.sort(allele_counts, axis=0)[-2] > gap_counts)
    )[0]

    # find the biallelic snps
    snp_positions_biallelic = np.where(
        (np.max(allele_counts, axis=0) > 1)
        & (np.max(allele_counts, axis=0) < n_phages - 1)
        & (gap_counts == 0)
        & (np.count_nonzero(allele_counts, axis=0) == 2)
    )[0]

    assert set(snp_positions_biallelic).issubset(snp_positions)

    return snp_positions, snp_positions_biallelic


def find_snp_positions_processed(core_genome_alignment):
    """
    input - core genome alignment in np.array format

    This is identical to find_snp_positions, except with the added "P" in the dna alphabet.
    """
    dna_alphabet = ["G", "T", "A", "C", "P"]
    n_phages = core_genome_alignment.shape[0]

    # count the number of alleles (and gaps) at each position in the core genome
    allele_counts = np.array(
        [np.sum(core_genome_alignment == nuc, axis=0) for nuc in dna_alphabet]
    )
    gap_counts = np.array(np.sum(core_genome_alignment == "-", axis=0))

    # find positions with snps
    snp_positions = np.where(
        (np.max(allele_counts, axis=0) > 1)
        & (np.max(allele_counts, axis=0) < n_phages - 1)
        & (np.max(allele_counts, axis=0) > gap_counts)
        & (np.sort(allele_counts, axis=0)[-2] > gap_counts)
    )[0]

    # find the biallelic snps
    snp_positions_biallelic = np.where(
        (np.max(allele_counts, axis=0) > 1)
        & (np.max(allele_counts, axis=0) < n_phages - 1)
        & (gap_counts == 0)
        & (np.count_nonzero(allele_counts, axis=0) == 2)
    )[0]

    assert set(snp_positions_biallelic).issubset(snp_positions)

    return snp_positions, snp_positions_biallelic


def convert_to_list(string_list):
    """
    Returns a list of integers from a string-type list of integers

    This is needed as the lists saved into a TSV file are converted to a string, so this
    function converts them back to a list
    """
    new_string_list = string_list.replace("[", "").replace("]", "")
    new_list = new_string_list.split(",")
    return [int(x) for x in new_list]


def rolling_average(data, roll_large):
    """
    Returns a somewhat log-scaled rolling average of the data
    """
    assert roll_large == 300
    roll_small, roll = 3, 30  # how much to rolling average
    hundred = np.where(data[::-1].index < 100)[0][0]
    ten = np.where(data[::-1].index < 10)[0][0]
    thousand = np.where(data[::-1].index < 1000)[0][0]
    if thousand > roll_large - 1:
        rolling_data = (
            list(
                data[::-1].rolling(roll_large).mean().values[roll_large - 1 : thousand]
            )
            + list(data[::-1].rolling(roll).mean().values[thousand:hundred])
            + list(data[::-1].rolling(roll_small).mean().values[hundred:ten])
            + list(data[::-1].values[ten:])
        )
    elif hundred > roll_large - 1:
        rolling_data = (
            list(data[::-1].rolling(roll).mean().values[roll_large - 1 : hundred])
            + list(data[::-1].rolling(roll_small).mean().values[hundred:ten])
            + list(data[::-1].values[ten:])
        )
    assert len(rolling_data) == len(data) - roll_large + 1
    return rolling_data
