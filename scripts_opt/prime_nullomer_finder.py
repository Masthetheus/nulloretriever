"""Compares multiple organisms nullomer tries for primes."""
import argparse
import yaml
from collections import defaultdict

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
        anchor_file = f"workflow/results/k{k}/{genomes[0]}/null_bit_format"
        anchor_v1_locations = obtain_v1_positions(anchor_file)
        common_nullomer_set = obtain_nullomer_set(
            anchor_file, anchor_v1_locations)
        common_test = defaultdict(set)
        for genome in genomes[1:]:
            file2 = f"workflow/results/k{k}/{genome}/null_bit_format"
            file2_v1_locations = obtain_v1_positions(file2)
            file2_null_set = obtain_nullomer_set(file2, file2_v1_locations)
            common_v1 = common_nullomer_set.keys() & file2_v1_locations.keys()
            for key in common_v1:
                comp = common_nullomer_set[key] & file2_null_set[key]
                if comp:
                    common_test[key] = comp
                else:
                    continue
            if len(common_v1) == 0:
                print(f"Nothing in similar with organism: {genome}.")
                break
            # common_nullomer_set = search_prime_nullomers(
            #     common_nullomer_set, common_v1, file2_v1_locations, file2)
        k += 1
        # print(len(common_nullomer_set))
        # total_len = sum(len(s) for s in common_nullomer_set.values())
        # print(total_len)
        print(f"Common: {len(common_test)}")
        total_len = sum(len(s) for s in common_test.values())
        print(f"Total len: {total_len}")


if __name__ == "__main__":
    main()
