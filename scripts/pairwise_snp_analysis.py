"""
This script looks at the SNPs between closely related phages to estimate the
    recombination block size.

Author: Jemma M. Fendley
"""

from Bio import AlignIO
import numpy as np, pandas as pd, argparse
from utils import calculate_hamming  # calculates hamming distance


def main():
    parser = argparse.ArgumentParser(description="pairwise SNP analysis for one group")
    parser.add_argument(
        "-o", "--output", help="output TSV file", type=str, required=True
    )
    parser.add_argument(
        "-c",
        "--core_genome_file",
        help="name of core genome file",
        type=str,
        required=True,
    )
    args = parser.parse_args()

    # load the core genome alignment
    group = args.core_genome_file.split("/")[-1].split("_")[0]
    alignment = AlignIO.read(args.core_genome_file, "fasta")
    phages = [record.id for record in alignment]
    core_genome_alignment = np.array(alignment)
    n_positions = core_genome_alignment.shape[1]

    binsize, dict_list = 100, []
    for i, phage1 in enumerate(phages[:-1]):
        for j in range(i + 1, len(phages)):
            # Calculate Hamming distance
            hamming = calculate_hamming(core_genome_alignment, i, j)

            # Need threshold to find pairs of phages that are closely related (so that recombination
            #   events do not overlap) but distant enough to potentially have a recombination event
            if hamming < 0.02 and hamming > 0.001:
                row = {"group": group, "phage_1": phage1, "phage_2": phages[j]}

                # find the snp positions and their distribubtion across the core genome
                snp_positions = np.where(
                    core_genome_alignment[i] != core_genome_alignment[j]
                )[0]
                counts = np.histogram(
                    snp_positions, bins=np.arange(0, n_positions + binsize, binsize)
                )[0]

                # calculate the squared distance to the diagonal (larger squared distances likely
                #   indicate one recombination event rather than scattered mutations)
                cumulative_sum = [np.sum(counts[0:i]) for i in range(1, len(counts))]
                diagonal = np.linspace(
                    cumulative_sum[0], cumulative_sum[-1], len(cumulative_sum)
                )
                squared_distance = (cumulative_sum - diagonal) ** 2
                row["squared_distance"] = sum(squared_distance)

                # find rolling bins with at least 5 snps
                rolling_counts = np.array(
                    [
                        np.count_nonzero(
                            (snp_positions >= i) & (snp_positions < i + binsize)
                        )
                        for i in range(n_positions - binsize)
                    ]
                )
                high_snp_bins = np.where(rolling_counts > 5)[0]
                # now calculate how long consecutive runs of high snp bins are
                if len(high_snp_bins) < 1:
                    continue
                elif len(high_snp_bins) == 1:
                    consecutive_lengths = [binsize]
                else:
                    consecutive_lengths = []
                consecutives = [high_snp_bins[0]]
                for k in range(1, len(high_snp_bins)):
                    next_position = high_snp_bins[k]
                    if next_position == consecutives[-1] + 1:
                        consecutives.append(next_position)
                    else:
                        consecutive_lengths.append(len(consecutives))
                        consecutives = [next_position]
                    if k == len(high_snp_bins) - 1:
                        consecutive_lengths.append(len(consecutives))
                row["distances"] = consecutive_lengths
                dict_list.append(row)

    # save the pairs with the highest 10 squared distance
    if len(dict_list) > 0:
        df = pd.DataFrame(dict_list)
        save_df = df.nlargest(10, "squared_distance")
    else:
        save_df = pd.DataFrame(
            columns=["group", "phage_1", "phage_2", "squared_distance", "distances"]
        )
    save_df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":  #
    main()
