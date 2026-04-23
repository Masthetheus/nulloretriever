"""Receives bit nullomer files for two consecutive k values and compare."""
import argparse
from collections import defaultdict

from nulloretriever.analysis.parsings import read_per_v1
from nulloretriever.analysis.counter import quick_nullomer_count
from nulloretriever.analysis.trivial_extensions import (
    first_half_extensions,
    second_half_extensions
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
        trivial_extensions = defaultdict(set)
        half_k = k1//2
        skip_counter = 0
        nullomer_count = 0
        new_counter_v1 = 0
        new_counter_v2 = 0
        file1 = f"workflow/results/k{k1}/{genome}/null_bit_format"
        file2 = f"workflow/results/k{k2}/{genome}/null_bit_format"
        count = quick_nullomer_count(file1)
        total_null_file2 = quick_nullomer_count(file2)
        while nullomer_count < count:
            v1, v2s, counter, nullomer_count, skip_counter = read_per_v1(
                file1, skip_counter, nullomer_count)
            new_counter_v1 = first_half_extensions(
                v1, v2s, half_k, file2, new_counter_v1, trivial_extensions)
            new_counter_v2 = second_half_extensions(
                v1, v2s, half_k, file2, new_counter_v2, trivial_extensions)
        k1 = k2
        k2 += 1
        print(f"v1: {new_counter_v1}. v2: {new_counter_v2}. Sum: {
              new_counter_v1+new_counter_v2}")
        trivial_porc = ((new_counter_v1+new_counter_v2)/total_null_file2) * 100
        print(f"Porcentagem de trivial extensions: {trivial_porc}")
        total = sum(len(v) for v in trivial_extensions.values())
        trivial_new_porc = (total/total_null_file2)*100
        print(f"Porcentagem de trivial extensions nova: {trivial_new_porc}")
        print(len(trivial_extensions))
        print(total)


if __name__ == "__main__":
    main()
