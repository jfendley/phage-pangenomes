# This script splits groups into subgroups (if large enough)

import argparse
import json
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="calculate linkage for one subgroup")
    parser.add_argument("output")
    parser.add_argument("subgroup_info")
    parser.add_argument("group_json")

    args = parser.parse_args()
    group = args.group_json.split("/")[-1].split(".")[0]
    with open(args.group_json) as f:
        G = json.load(f)
    num_phages = len(G["paths"])
    core_phams = set(
        [
            x["pham_ID"]
            for x in G["phams"]
            if x["total_count"] == num_phages and x["phage_count"] == num_phages
        ]
    )
    with open(args.subgroup_info) as g:
        subgroup_info = json.load(g)
    subgroup_indices = [x["subgroup_ind"] for x in subgroup_info]
    dict_list = []
    for index in subgroup_indices:
        subgroup_list = [x for x in subgroup_info if x["subgroup_ind"] == index]
        assert len(subgroup_list) == 1
        subgroup = subgroup_list[0]
        phages = subgroup["phage_ind"]

        phage_paths = [
            set([x["pham_ID"] for x in G["paths"][i]["path"]]) for i in phages
        ]
        sub_core_phams = set.intersection(*phage_paths)
        assert core_phams - sub_core_phams == set(), "Likely error with duplicate phams"

        new_core = sub_core_phams - core_phams
        assert len(sub_core_phams) - len(core_phams) == len(new_core)
        row = {
            "group": group,
            "subgroup_index": index,
            "num_core": len(core_phams),
            "num_sub_core": len(sub_core_phams),
            "num_new_core": len(new_core),
        }
        dict_list.append(row)
    df = pd.DataFrame(dict_list)
    df.to_csv(args.output, index=False, sep="\t")


if __name__ == "__main__":
    main()
