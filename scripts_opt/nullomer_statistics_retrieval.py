"""Retrieve nullomer statistics from a bit file, snakemake compatible."""

import csv

from nulloretriever.analysis.processing import mount_trie_from_bitfile

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
    out_path = "retrieval_test.csv"
    stats = ["composition","counter","motifs"]
    nullomer_file = "workflow/results/k14/GCF_000146045_2/null_bit_format"
    anchor_trie = "othernullfile_result"
    organism = "GCA_000412225_2"
    k_val = 10
    half_k = k_val//2
    trie = mount_trie_from_bitfile(nullomer_file)
    dispatch_table = {
        "composition": trie.count_gc(),
        "counter": trie.count_kmers(),
        "motifs": motif_wrapper(trie)
    }
    retrieved_stats = {}

#    for stat in stats:
#        try:
#            if stat in dispatch_table.keys():
#                func = dispatch_table[stat]
#                if callable(func):
#                    retrieved_stats[stat] = func(nullomer_file)
#        except Exception as err:
#            print("Unexpected occurence processing the"
#                  f"following statistic: {stat}.\n"
#                  f"Error: {err=}, {type(err)=}")
#            raise

    print(trie.kmer_count())
    #for stat in stats:
    #    retrieved_stats[stat] = dispatch_table[stat]
    #base_dict = {}
    #base_dict['organism'] = organism
    #base_dict['k'] = k_val
    #final_stats = dict_flattener(retrieved_stats, base_dict)
    #with open(out_path, mode='w', newline='') as f:
    #    writer = csv.writer(f)
    #    writer.writerow(final_stats.keys())
    #    writer.writerow(final_stats.values())
    #print(retrieved_stats)


if __name__ == "__main__":
    main()
