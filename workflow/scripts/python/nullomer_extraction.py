"""Pipeline k-mer processing and nullomer extraction."""
import subprocess

from nulloretriever.core.triebit_class import TrieBit
from snakemake.script import Snakemake


def main():
    """K-mer extraction and nullomer trie generation snakemake compatible."""
    out_path = snakemake.output[0]
    genome_path = snakemake.input.genome
    k_str = str(snakemake.params.k_val)
    k = int(snakemake.params.k_val)
    is_odd = k % 2
    half_k = (k // 2) + is_odd
    m = 4 ** (half_k - is_odd)
    v2_size = k - half_k
    k_mask = (1 << (v2_size * 2)) - 1
    trie = TrieBit(m, k, half_k)
    proc = subprocess.Popen(
        [snakemake.input.bin,
         genome_path, k_str],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        bufsize=65536
    )
    sequences_received = 0
    total = 0
    bytes_per_sequence = (k * 2 + 7) // 8
    wasted_space = (bytes_per_sequence * 8) - (k * 2)
    mask = (4**k) - 1
    print(f"DEBUG: Variables in use:\n"
          f"bytes_per_sequence: {bytes_per_sequence}"
          f"\nwasted_space:{wasted_space}"
          f"\nk_mask: {k_mask}"
          f"\nmask: {mask}")
    while True:
        kmer_bytes = proc.stdout.read(bytes_per_sequence)
        if len(kmer_bytes) < bytes_per_sequence:
            break
        kmer_idx = int.from_bytes(kmer_bytes, byteorder='big')
        kmer_idx = (kmer_idx) & (mask)
        sequences_received += 1
        v1 = kmer_idx >> (v2_size*2)
        v2 = kmer_idx & k_mask
        trie.insert(v1, v2)
        total += 1
    proc.wait()
    kmers_inserted = trie.count_kmers()
    expected = 4**k
    null_count = expected - kmers_inserted
    print(
        f"{trie.count_kmers()} were inserted into the trie."
        f"\n{expected} total were expected."
        f"An total of {null_count} nullomers were found "
        f"for a k of {k}."
    )
    trie.write_bit_format(out_path)


if __name__ == "__main__":
    main()
