import argparse
import yaml
from collections import defaultdict

from nulloretriever.core.triebit_class import TrieBit
from nulloretriever.analysis.prime_nullomer import *
from nulloretriever.analysis.parsings import *


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


def main():
    """Search for nullomer primes."""
    parser = setup_argparser()
    args = parser.parse_args()
    output = args.output
    k_range = args.k_range
    config_path = args.config
    with open(config_path, 'r') as f:
        try:
            data = yaml.safe_load(f)
            genomes = data.get('organisms')
        except yaml.YAMLError as exc:
            print(f"An error was found parsing the yaml file: {exc}.")

    k = int(k_range[0])
    while k <= k_range[1]:
        remove_counter = 0
        print(f"======{k}======")
        anchor_file = f"workflow/results/k{k}/{genomes[0]}/null_bit_format"
        anchor_trie = mount_trie_from_bitfile(anchor_file)
        for genome in genomes[1:]:
            second_file = f"workflow/results/k{k}/{genome}/null_bit_format"
            second_trie = mount_trie_from_bitfile(anchor_trie)
