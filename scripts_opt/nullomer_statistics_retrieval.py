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
    parser.add_argument(
        "-s",
        "--stats",
        help="Stats to be analyzed. Default ="
        "all. Choices = composition, trivial, motifs.",
        choices=["all", "composition", "trivial", "motifs"],
        default="all",
    )
    return parser


def motif_wrapper(trie):
    """Calls all functions related to motif statistics."""
    motifs_results = {
        "cpg": trie.retrieve_nullomers_cpg_stats(),
        "palindromy": trie.retrieve_palindrome_stats(),
        "homopolymers": trie.retrieve_homopolymer_stats(),
    }
    return motifs_results


def dict_flattener(full_dict, final_dict=None, parent_key=""):
    if final_dict is None:
        final_dict = {}

    for key, value in full_dict.items():
        new_key = f"{parent_key}_{key}" if parent_key else key
        if isinstance(value, dict):
            dict_flattener(value, final_dict, parent_key=new_key)
        else:
            final_dict[new_key] = value

    return final_dict


def main():
    """Retrieve nullomer statistics according to user specification.
    Snakemake compatible script aimed to operate on nullomer bit files and
    retrieve statistics of such. The details on each function operation can be
    found in their module of origin.
    """
    parser = setup_argparser()
    args = parser.parse_args()
    out_path = args.output
    stats = args.stats
    if stats == "all":
        stats = ["composition", "trivial", "motifs"]
    else:
        stats = [stats]
    bigger_null_file = args.null1
    if stats == "all" or "trivial" in stats:
        if not args.null2:
            print(
                "Invalid null2 argument given. For trivial analysis the bit"
                "null file for k-1 from the same organism is needed."
            )
            return
        smaller_null_file = args.null2
    bigger_trie = mount_trie_from_bitfile(bigger_null_file)
    counter = bigger_trie.count_kmers()

    if counter > 0 and smaller_null_file != "":
        smaller_trie = mount_trie_from_bitfile(smaller_null_file)
        if smaller_trie.count_kmers() != 0:
            dispatch_table = {
                "composition": bigger_trie.count_gc(),
                "trivial": bigger_trie.find_trivial_ext(smaller_trie),
                "motifs": motif_wrapper(bigger_trie),
            }
        else:
            dispatch_table = {
                "composition": bigger_trie.count_gc(),
                "motifs": motif_wrapper(bigger_trie)
            }
    else:
        dispatch_table = {
            "composition": bigger_trie.count_gc(),
            "motifs": motif_wrapper(bigger_trie)
        }
        print("No nullomers found in the given bit file.")

    retrieved_stats = {}
    retrieved_stats["counter"] = counter

    for stat in stats:
        try:
            print(f"Processing {stat}.")
            retrieved_stats[stat] = dispatch_table[stat]
        except Exception as e:
            print(f"Stat {stat} no available for this organism")
            print(f"Exception= {e}.")

    base_dict = {}
    final_stats = dict_flattener(retrieved_stats, base_dict)
    with open(out_path, mode="w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(final_stats.keys())
        writer.writerow(final_stats.values())


if __name__ == "__main__":
    main()
