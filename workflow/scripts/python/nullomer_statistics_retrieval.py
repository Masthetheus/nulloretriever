"""Retrieve nullomer statistics from a bit file, snakemake compatible."""

import csv

from nulloretriever.analysis.processing import mount_trie_from_bitfile
from snakemake.script import Snakemake


def motif_wrapper(trie):
    """Calls all functions related to motif statistics."""
    motifs_results = {
        "cpg": trie.retrieve_nullomers_cpg_stats(),
        "palindromy": trie.retrieve_palindrome_stats(),
        "homopolymers": trie.retrieve_homopolymer_stats()
    }
    return motifs_results


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
    bigger_trie = mount_trie_from_bitfile(bigger_null_file)
    counter = bigger_trie.count_kmers()
    if counter > 0 and smaller_null_file != '':
        smaller_trie = mount_trie_from_bitfile(smaller_null_file)
        if smaller_trie.count_kmers() != 0:
            print("Pegando trivial")
            dispatch_table = {
                "composition": bigger_trie.count_gc(),
                "trivial": bigger_trie.find_trivial_ext(smaller_trie),
                "motifs": motif_wrapper(bigger_trie)
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
    retrieved_stats = {}
    retrieved_stats['counter'] = counter

    for stat in stats:
        try:
            retrieved_stats[stat] = dispatch_table[stat]
            if stat == "trivial":
                print(retrieved_stats[stat])
        except Exception as e:
            print(f"Stat {stat} no available for this organism")
    base_dict = {}
    base_dict['organism'] = organism
    base_dict['k'] = k_val
    print(retrieved_stats)
    final_stats = dict_flattener(retrieved_stats, base_dict)
    with open(out_path, mode='w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(final_stats.keys())
        writer.writerow(final_stats.values())


if __name__ == "__main__":
    main()
