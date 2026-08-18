"""Retrieve nullomer statistics from a bit file, snakemake compatible."""

import argparse
import csv

from nulloretriever.analysis.processing import mount_trie_from_bitfile


def setup_argparser() -> argparse.ArgumentParser:
    """Argument passing function for the nullomer statistics script."""
    parser = argparse.ArgumentParser(
        description=(
            "Script for nullomer statistics retrieval and organization in a single csv file."
        )
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output directory for the csv statistics file."
        "Default =workflow/results/custom_statistics_retrieval.csv",
        default="workflow/results/custom_statistics_retrieval.csv",
    )
    parser.add_argument("-n1", "--null1", help="Nullomer file to be analyzed")
    parser.add_argument(
        "-n2",
        "--null2",
        help="Nullomer file of the same"
        "organism in k-1 value. Used to obtain trivial extension.",
    )
    return parser

def main():
    """Given bit nullomer file for k and k-1, extract MAWs
    """
    parser = setup_argparser()
    args = parser.parse_args()
    out_path = args.output

    bigger_null_file = args.null1
    bigger_trie = mount_trie_from_bitfile(bigger_null_file)

    try:
        smaller_null_file = args.null2
    except:
        print("Nullomer file for k-1 not found.")

    counter = bigger_trie.count_kmers()

    if counter > 0 and smaller_null_file != "" and smaller_null_file != None:
        print("Mounting smaller trie")
        smaller_trie = mount_trie_from_bitfile(smaller_null_file)
        print(smaller_trie.count_kmers())
        if smaller_trie.count_kmers() != 0:
            print("Building maw trie")
            bigger_trie.build_maw_trie(smaller_trie,out_path),

    print(counter)
    #with open(out_path, mode="w", newline="") as f:
    #    writer = csv.writer(f)
    #    writer.writerow(final_stats.keys())
    #    writer.writerow(final_stats.values())


if __name__ == "__main__":
    main()
