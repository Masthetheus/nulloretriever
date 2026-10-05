"""Retrieve nullomer statistics from a bit file, snakemake compatible."""

import csv

from nulloretriever.analysis.processing import mount_trie_from_bitfile
from snakemake.script import Snakemake


def dict_flattener(full_dict, final_dict=None, parent_key=''):
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
    Composition:
        gc_mean -> calculate the mean gc of all existing nullomers
    """
    out_path = snakemake.output[0]
    stats = snakemake.params.stats
    bigger_null_file = snakemake.input.current
    smaller_null_file = str(snakemake.input.get("previous", None))
    organism = snakemake.wildcards.organism
    k_val = snakemake.wildcards.k

    print(f"SMALLER NULL: {smaller_null_file} for k {k_val}")

    bigger_trie = mount_trie_from_bitfile(bigger_null_file)
    counter = bigger_trie.count_kmers()

    retrieved_stats = {"counter": counter}

    needs_composition = "composition" in stats
    needs_motifs = "motifs" in stats

    if counter > 0 and (needs_composition or needs_motifs):
        print("Processing composition and motifs.")
        gc_percent, cpg_stats, palindrome_stats = (
            bigger_trie.retrieve_composition_and_motifs()
        )
        if needs_composition:
            retrieved_stats["composition"] = gc_percent
        if needs_motifs:
            retrieved_stats["motifs"] = {
                "cpg": cpg_stats,
                "palindromy": palindrome_stats,
                "homopolymers": bigger_trie.retrieve_homopolymer_stats(),
            }

    if "trivial" in stats:
        if counter > 0 and smaller_null_file and smaller_null_file != "None":
            smaller_trie = mount_trie_from_bitfile(smaller_null_file)
            smaller_counter = smaller_trie.count_kmers()
            print(f"Smaller trie loaded with {smaller_counter} nulls")
            if smaller_counter != 0:
                retrieved_stats["trivial"] = bigger_trie.find_trivial_ext()
            else:
                print("Stat trivial not available: smaller trie is empty.")
        else:
            print("Stat trivial not available: no previous file provided.")

    base_dict = {"organism": organism, "k": k_val}
    final_stats = dict_flattener(retrieved_stats, base_dict)

    print(retrieved_stats)
    with open(out_path, mode='w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(final_stats.keys())
        writer.writerow(final_stats.values())


if __name__ == "__main__":
    main()
