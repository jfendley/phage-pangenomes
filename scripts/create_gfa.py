"""
This script converts the group JSON into a GFA file for visualization with Bandage.

Author: Jemma M. Fendley
"""

import argparse, json


def main():
    parser = argparse.ArgumentParser(description="create GFA for visualization")
    parser.add_argument(
        "-i", "--input", help="name of JSON file", type=str, required=True
    )
    parser.add_argument(
        "-o", "--output", help="name of GFA file", type=str, required=True
    )
    args = parser.parse_args()

    with open(args.input) as f:
        G = json.load(f)

    with open(args.output, "w") as gfa:
        gfa.write("H\tVN:Z:1.0\n# phams")  # Start the GFA file

    # Add the pham information, then the edges, then the paths/sequences
    with open(args.output, "a") as gfa:
        for pham in G["phams"]:
            gfa.write(
                "\nS\t"
                + pham["pham_ID"]
                + "\t*\tLN:i:"
                + str(round(pham["mean_length"]))
                + "\tRC:i:"
                + str(pham["read_length"])
            )

        gfa.write("\n# edges")
        for edge in G["edges"]:
            gfa.write(
                "\nL\t"
                + edge["pham1"]
                + "\t"
                + edge["strand1"]
                + "\t"
                + edge["pham2"]
                + "\t"
                + edge["strand2"]
                + "\t*"
                + "\tRC:i:"
                + str(edge["total_count"])
            )

        gfa.write("\n# sequences")
        for path in G["paths"]:
            gfa.write("\nP\t" + path["ID"] + "\t")
            for pham in path["path"]:
                gfa.write(pham["pham_ID"] + pham["strand"] + ",")
            gfa.seek(gfa.tell() - 1, 0)  # seek to the penultimate char of file
            gfa.truncate()
            gfa.write("\tTP:Z:linear*")


if __name__ == "__main__":
    main()
