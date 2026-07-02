"""K-mer processing and nullomer extraction on the snakemake pipeline."""
import argparse
import subprocess

from nulloretriever.core.triebit_class import TrieBit
from nulloretriever.analysis.processing import *


def setup_argparser() -> argparse.ArgumentParser:
    """Argument passing function for the kmer extraction script."""
    parser = argparse.ArgumentParser(description="""Script for k-mer processing
                                     and nullomer extraction on the snakemake
                                     pipeline""")
    parser.add_argument(
        '-k',
        '--kvalues',
        help="Values of k for k-mer extraction. More than one value can be"
        "informed. Default = 10",
        nargs="*",
        default=[10]
    )
    parser.add_argument(
        '-m',
        '--mode',
        help="Output mode, if in compact txt format or binary format."
        " Default = binary.",
        choices=['binary', 'txt'],
        default='binary'
    )
    parser.add_argument(
        '-o',
        '--output',
        help="Output for the found nullomers file. Default = workflow/data/",
        default='workflow/data'
    )
    parser.add_argument(
        '-g',
        '--genome',
        help="Genome to be analyzed"
    )
    return parser


def main():
    """K-mer extraction and nullomer trie generation."""
    parser = setup_argparser()
    args = parser.parse_args()
    k_values = args.kvalues
    mode = args.mode
    out_path = args.output+'_result'
    genome_path = args.genome
    k = int(k_values[0])
    if k % 2 == 0:
        half_k = int(k/2)
        m = 4**half_k
        k_mask = (2**k) - 1
        v1_size = v2_size = half_k
    else:
        half_k = int(k/2) + 1
        m = 4**(half_k-1)
        k_mask = (2**(k-1)) - 1
        v1_size = half_k
        v2_size = k - half_k
    trie = TrieBit(m, k, half_k)
    proc = subprocess.Popen(
        ["../workflow/scripts/c/bin/kmer_extractor",
         genome_path, k_values[0]],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        bufsize=65536
    )
    sequences_received = 0
    total = 0
    bytes_per_sequence = (k * 2 + 7) // 8
    wasted_space = (bytes_per_sequence * 8) - (k * 2)
    mask = (4**k) - 1
    kmers = set()
    print(f"DEBUG: Variables in use:\n"
          f"bytes_per_sequence: {bytes_per_sequence}"
          f"\nwasted_space:{wasted_space}"
          f"\nk_mask: {k_mask}"
          f"\nmask: {bin(mask)}")
    while True:
        kmer_bytes = proc.stdout.read(bytes_per_sequence)
        if len(kmer_bytes) < bytes_per_sequence:
            break
        kmer_idx = int.from_bytes(kmer_bytes, byteorder='big')
        kmer_idx = (kmer_idx) & (mask)
        kmers.add(kmer_idx)
        sequences_received += 1
        v1 = kmer_idx >> (v2_size*2)
        v2 = kmer_idx & k_mask
        v1_bits = []
        i = half_k - 1
        while i >= 0:
            v1_bits.append(v1 >> (i*2) & 3)
            i -= 1
        trie.insert(tuple(v1_bits), v2)
        total += 1
    proc.wait()
    trie.write_bit_format(out_path)
    trie1 = mount_trie_from_bitfile(out_path)
    print(trie.count_kmers())
    expected = 4**k
    null_count = trie1.count_kmers()
    obtained = null_count + total
    diff = expected - total - null_count
    print(f"{total} k-mkers were inserted, and {trie.count_kmers()}"
          f" kmers were counted.\n {expected} total were expected."
          f"We have total + null equals {obtained}/{expected}.")
    print(f"A total of {null_count} nullomers were found.")
    print(len(kmers))


if __name__ == "__main__":
    main()
