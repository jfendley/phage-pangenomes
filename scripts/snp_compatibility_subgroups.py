"""
This script considers groups that have two subgroups, and it counts the number of SNPs
    that are incompatible with the subgroup split.

Author: Jemma M. Fendley
"""

from Bio import AlignIO
import numpy as np, json
import matplotlib.pyplot as plt
import pandas as pd, seaborn as sns
import parameters, argparse


def classify_snp(n_alleles, n_alleles_1, n_alleles_2):
    """
    Classifies a SNP based on:

    n_alleles (int): number of different alleles at this position
    n_alleles_1 (int): number of different alleles within subgroup 1
    n_alleles_2 (int): number of different alleles within subgroup 2
    """
    assert n_alleles > 1, "Not a SNP site"

    if n_alleles == 2:  # biallelic
        if n_alleles_1 == 2 and n_alleles_2 == 2:
            return "incompatible"
        elif n_alleles_2 == 2:
            return "single mutation 2"
        elif n_alleles_1 == 2:
            return "single mutation 1"
        elif n_alleles_1 == 1 and n_alleles_2 == 1:
            return "split"  # perfectly corresponds to subgroup split
        else:
            assert (
                1 == 0
            ), "Error somewhere"  # all biallelic SNPs should fall into one of these four categories
    elif n_alleles == 3:  # tri-allelic
        if n_alleles_1 == 2 and n_alleles_2 == 2:
            return "triallelic across"
        elif n_alleles_1 == 2:
            return "single mutation 1"
        elif n_alleles_2 == 2:
            return "single mutation 2"
        else:
            return "other tri-allelic"
    elif n_alleles == 4:
        if n_alleles_1 == 2 and n_alleles_2 == 2:
            return "single mutation both"
        elif n_alleles_2 == 2:
            return "single mutation 2"
        elif n_alleles_1 == 2:
            return "single mutation 1"
        else:
            return "other four-allelic"


def main():
    parser = argparse.ArgumentParser(
        description="creates a figure showing SNPs incompatible with subgroups"
    )
    parser.add_argument(
        "-o", "--output", help="figure png file", type=str, required=True
    )
    parser.add_argument(
        "-p", "--output_pdf", help="figure pdf file", type=str, required=True
    )
    parser.add_argument(
        "-s",
        "--subgroup",
        help="JSON file with subgroup information for each group",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-c",
        "--core_genomes",
        help="core genome files for groups with two subgroups",
        nargs="+",
        required=True,
    )
    args = parser.parse_args()

    # load the subgroup information, phages per subgroup
    with open(args.subgroup) as f:
        subgroup_to_phages = json.load(f)

    # translate this into a nested dictionary of group to subgroup to phages
    all_groups = np.unique([x.split("-")[0] for x in list(subgroup_to_phages.keys())])
    group_to_subgroups = {
        x: [y for y in list(subgroup_to_phages.keys()) if y.split("-")[0] == x]
        for x in all_groups
    }
    group_to_sizes = {
        x: [len(subgroup_to_phages[z]) for z in y]
        for x, y in group_to_subgroups.items()
    }
    group_to_subgroup_phages = {
        x: {
            y.split("-")[1]: z
            for y, z in subgroup_to_phages.items()
            if y.split("-")[0] == x
        }
        for x in all_groups
    }

    # find the groups with exactly two subgroups with at least 10 phages
    groups_to_analyze = [
        x
        for x, y in group_to_sizes.items()
        if np.count_nonzero(np.array(y) >= 10) == 2 and len(y) == 2
    ]

    second_dict_list = []  # initialize

    # loop through all of the groups
    for group in groups_to_analyze:

        # find the core_genome file that corresponds to this group
        core_genome_files = [
            x for x in args.core_genomes if "results/groups/" + group + "/" in x
        ]
        assert len(core_genome_files) == 1
        core_genome_file = core_genome_files[0]

        # load the subgroup information, the phages in each subgroup
        subgroup_initial = group_to_subgroup_phages[group]
        n_phage_dict = {x: len(y) for x, y in subgroup_initial.items()}

        # reorder so that by default the smaller subgroup is subgroup 1
        reorder = np.argsort(list(n_phage_dict.values()))
        subgroup = {
            str(i + 1): y
            for i, y in enumerate(
                [
                    x
                    for x in [list(subgroup_initial.values())[z] for z in reorder]
                    if len(x) >= 10
                ]
            )
        }
        all_phages = [x for y in list(subgroup.values()) for x in y]

        # load the core genome alignment
        core_genome_initial = AlignIO.read(core_genome_file, "fasta")
        core_genome = [
            record for record in core_genome_initial if record.id in all_phages
        ]

        # record the order of phages in the core genome alignment
        phage_order = [record.id for record in core_genome]

        # convert to array and count the number of positions
        core_genome_alignment = np.array(core_genome)
        n_positions = core_genome_alignment.shape[1]

        # consider the subgroups separately
        subgroup_indices = {
            x: [phage_order.index(z) for z in y] for x, y in subgroup.items()
        }
        dna_alphabet = ["G", "T", "A", "C"]

        # count the number of alleles (and gaps) at each position in the core genome
        allele_counts = np.array(
            [np.sum(core_genome_alignment == nuc, axis=0) for nuc in dna_alphabet]
        )
        gap_counts = np.array(np.sum(core_genome_alignment == "-", axis=0))

        # count the number of positions that have no gaps
        gap_sites = np.where(gap_counts > 0)[0]
        n_positions_no_gaps = n_positions - len(gap_sites)

        # find the positions that are multi-allelic with no gaps
        any_snps = np.where(
            (gap_counts == 0) & (np.count_nonzero(allele_counts, axis=0) > 1)
        )[0]

        # quick check to ensure that everything adds up correctly
        n_single_allele = np.count_nonzero(
            (np.count_nonzero(allele_counts, axis=0) == 1) & (gap_counts == 0)
        )
        assert n_positions_no_gaps == n_single_allele + len(any_snps), "Problem"

        # compile the SNPs into a dataframe and classsify each SNP
        dict_list = [
            {
                "position": x,
                "n_alleles": len(set(column)),
                "n_alleles_1": len(set(column[subgroup_indices["1"]])),
                "n_alleles_2": len(set(column[subgroup_indices["2"]])),
            }
            for x, column in enumerate(core_genome_alignment.T)
            if x in any_snps
        ]
        df = pd.DataFrame(dict_list)
        df["type"] = df.apply(
            lambda x: classify_snp(x.n_alleles, x.n_alleles_1, x.n_alleles_2), axis=1
        )

        # compile the SNP types that represent single mutations in either subgroup 1 or 2
        both_list = ["triallelic across", "single mutation both", "incompatible"]
        list_1 = both_list + ["single mutation 1"]
        list_2 = both_list + ["single mutation 2"]

        # double check there were no typos and all of them are in fact in the list
        assert all([x in df["type"].unique() for x in list_1])
        assert all([x in df["type"].unique() for x in list_2])

        n_single_mutation_1 = len(df[df["type"].isin(list_1)])
        n_single_mutation_2 = len(df[df["type"].isin(list_2)])
        print(n_single_mutation_1, n_single_mutation_2)
        n_triallelic_across = len(df[df["type"] == "triallelic across"])

        # save the information for plotting
        second_dict_list.append(
            {
                "group": group,
                "n. positions": n_positions_no_gaps,
                "estimate": n_triallelic_across / 2,
                "type": "ratio with tri-allelic estimate (eq. 7)",
                "n. incompatible SNPs": len(df[df["type"] == "incompatible"]),
                "n. SNPs": len(df),
                "n. bialellic": len(df[df["n_alleles"] == 2]),
            }
        )
        second_dict_list.append(
            {
                "group": group,
                "n. positions": n_positions_no_gaps,
                "estimate": n_single_mutation_1
                * n_single_mutation_2
                / (3 * n_positions_no_gaps),
                "type": "ratio with subgroup mutation freq. estimate (eq. 8)",
                "n. incompatible SNPs": len(df[df["type"] == "incompatible"]),
                "n. SNPs": len(df),
                "n. bialellic": len(df[df["n_alleles"] == 2]),
            }
        )
    output_df = pd.DataFrame(second_dict_list)

    # initialize figure
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(5, 4), layout="constrained")

    # calculate and plot the n. incompatible / n. expected
    output_df["ratio of n. incompatible / n. expected"] = (
        output_df["n. incompatible SNPs"] / output_df["estimate"]
    )
    bars = sns.barplot(
        data=output_df,
        x="group",
        y="ratio of n. incompatible / n. expected",
        hue="type",
        ax=axes,
        alpha=0.8,
    )

    # hatch half of them
    for i in range(4):
        bars.patches[i].set_hatch("//")

    # add line y=x for reference
    axes.axhline(
        1,
        color="k",
        label=r"ratio=1 (n. incompatible = n. expected)",
        linewidth=4,
    )

    # format and save the figure
    leg = axes.legend(
        markerscale=0.1, handlelength=0.5, handleheight=0.8, handletextpad=0.5
    )
    leg.get_patches()[0].set_hatch("//")

    fig.savefig(args.output, dpi=450)
    fig.savefig(args.output_pdf, dpi=450)


if __name__ == "__main__":
    main()
