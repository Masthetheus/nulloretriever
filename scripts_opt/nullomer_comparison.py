"""Receives bit nullomer files for two consecutive k values and compare."""
import argparse
from collections import defaultdict

from nulloretriever.analysis.parsings import (
    obtain_v1_positions, gather_v1_related_data
)
from nulloretriever.analysis.counter import quick_nullomer_count
from nulloretriever.analysis.trivial_extensions import (
    retrieve_trivial_extensions
)


def setup_argparser() -> argparse.ArgumentParser:
    """Receives arguments for use in the trivial extension nullomer search."""
    parser = argparse.ArgumentParser(description="""
                                     Script aimed at searching different
                                     organisms nullomer information for
                                     primes.""")
    parser.add_argument(
        '-g',
        '--genome',
        help="Genome accession code for analysis.",
        required=True
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
    """Search for nullomers in a k value that are a trivial extension."""
    parser = setup_argparser()
    args = parser.parse_args()
    genome = args.genome
    k_range = args.k_range
    k1 = int(k_range[0])
    k2 = k1 + 1
    while k2 <= k_range[1]:
        trivial_extensions_test = defaultdict(set)
        half_k = k1//2
        if k1 % 2 == 0:
            odd = False
        else:
            odd = True
        file1 = f"workflow/results/k{k1}/{genome}/null_bit_format"
        file2 = f"workflow/results/k{k2}/{genome}/null_bit_format"
        total_null_file2 = quick_nullomer_count(file2)
        v1_indexes_file1 = obtain_v1_positions(file1)
        v1_indexes_file2 = obtain_v1_positions(file2)
        extension_counter = 0
        for v1, index in v1_indexes_file1.items():
            v2s = gather_v1_related_data(index, file1)
            extension_counter += retrieve_trivial_extensions(
                v1, v2s, v1_indexes_file2, file2, half_k,
                trivial_extensions_test, odd)
        k1 = k2
        k2 += 1
        total = sum(len(v) for v in trivial_extensions_test.values())
        trivial_new_porc = (total/total_null_file2)*100
        print(f"Novo modo trivial total: {total}.")
        print(f"Porcentagem de trivial extensions nova: {trivial_new_porc}")


if __name__ == "__main__":
    main()
