"""Compares multiple organisms nullomer tries for primes."""
import argparse
import yaml
import csv
from collections import defaultdict

from nulloretriever.analysis.prime_nullomer import *
from nulloretriever.analysis.parsings import *
from nulloretriever.analysis.processing import *


def setup_argparser() -> argparse.ArgumentParser:
    """Parsing function for the prime finder script."""
    parser = argparse.ArgumentParser(description="""
                                     Script aimed at searching different
                                     organisms nullomer information for
                                     primes.""")
    parser.add_argument(
        '-c',
        '--config',
        help="Relative path to the YAML config file containing the organisms"
        " names. Default = workflow/config/config.yaml",
        default="workflow/config/config.yaml"
    )
    parser.add_argument(
        '-o',
        '--output',
        help="Relative path to the CSV final analysis file."
        "Default = workflow/results/prime_nullomers.csv",
        default="workflow/results/prime_nullomers.csv"
    )
    parser.add_argument(
        '-k',
        '--k_range',
        help="Range of k values to search trivial extensions. Default = 10-12",
        nargs=2,
        type=int,
        default=[10, 12]
    )
    return parser

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
    """Search for nullomer primes."""
    parser = setup_argparser()
    args = parser.parse_args()
    output = args.output
    k_range = args.k_range
    config_path = args.config
    all_stats = []
    with open(config_path, 'r') as f:
        try:
            data = yaml.safe_load(f)
            genomes = data.get('organisms')
        except yaml.YAMLError as exc:
            print(f"An error was found parsing the yaml file: {exc}.")

    k = int(k_range[0])
    orgs_analyzed = len(genomes)
    while k <= k_range[1]:
        is_odd = k%2
        half_k = (k//2) + is_odd
        m = 4**(k//2)
        remove_counter = 0
        print(f"======{k}======")
        anchor_trie = mount_trie_from_bitfile(f"workflow/results/k{k}/{genomes[0]}/null_bit_format")
        for genome in genomes[1:]:
            try:
                trie2 = mount_trie_from_bitfile(f"workflow/results/k{k}/{genome}/null_bit_format")
                anchor_trie.retrieve_prime_null(trie2)
                print(anchor_trie.count_kmers())
                if (anchor_trie.count_kmers()) == 0:
                    break
            except Exception as e:
                continue
        count = anchor_trie.count_kmers()
        if count > 0:
            retrieved_stats = {
            'counter':count,
            'org_count': orgs_analyzed,
            'gc': anchor_trie.count_gc(),
            'cpg_stats': anchor_trie.retrieve_nullomers_cpg_stats(),
            'palindrome_stats': anchor_trie.retrieve_palindrome_stats()
            }
            base_dict = {'k': k}
            row_data = dict_flattener(retrieved_stats, base_dict)
            all_stats.append(row_data)
        k += 1
    if all_stats:
        anchor_trie.write_txt_format("prime_null")
        file_name = f"prime_null_k_{k_range[0]}_to_{k}.csv"
        header = list({column: True for line in all_stats for column in line.keys()}.keys())
        with open(file_name, "w", newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=header)
            writer.writeheader()
            writer.writerows(all_stats)

if __name__ == "__main__":
    main()
