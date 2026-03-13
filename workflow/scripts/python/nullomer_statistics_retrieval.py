"""Retrieve nullomer statistics from a bit file, snakemake compatible."""

import csv

from nulloretriever.analysis.composition import nullomers_gc_mean
from nulloretriever.analysis.motifs import (
    retrieve_nullomers_cpg_stats,
    retrieve_palindrome_stats,
    retrieve_homopolymer_stats
)
from nulloretriever.analysis.counter import quick_nullomer_count
from snakemake.script import snakemake


def motif_wrapper(filename):
    """Calls all functions related to motif statistics."""
    motifs_results = {
        "cpg": retrieve_nullomers_cpg_stats(filename),
        "palindromy": retrieve_palindrome_stats(filename),
        "homopolymers": retrieve_homopolymer_stats(filename)
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
    nullomer_file = snakemake.input[0]
    organism = snakemake.wildcards.organism
    k_val = snakemake.wildcards.k
    dispatch_table = {
        "composition": nullomers_gc_mean,
        "counter": quick_nullomer_count,
        "motifs": motif_wrapper
    }
    retrieved_stats = {}

    for stat in stats:
        try:
            if stat in dispatch_table.keys():
                func = dispatch_table[stat]
                if callable(func):
                    retrieved_stats[stat] = func(nullomer_file)
        except Exception as err:
            print("Unexpected occurence processing the"
                  f"following statistic: {stat}.\n"
                  f"Error: {err=}, {type(err)=}")
            raise

    base_dict = {}
    base_dict['organism'] = organism
    base_dict['k'] = k_val
    final_stats = dict_flattener(retrieved_stats, base_dict)
    with open(out_path, mode='w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(final_stats.keys())
        writer.writerow(final_stats.values())


if __name__ == "__main__":
    main()
