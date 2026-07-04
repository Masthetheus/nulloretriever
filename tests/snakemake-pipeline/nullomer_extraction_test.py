"""Teste file for the k-mer processing and nullomer extraction script."""
import argparse
import subprocess
import pytest

from nulloretriever.core.triebit_class import TrieBit

def variaveis_k(k):

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
    return half_k, m, k_mask, v1_size, v2_size

def test_global():
    """K-mer extraction and nullomer trie generation snakemake compatible."""

    k = input("Which k value to test?")
    out_path = f"tests/results/k{k}_extraction_result"
    genome_path = f"tests/data/k{k}_genome_test"
    c_bin = "tests/data/c_kmer_extraction_bin"
    k = int(k)
    half_k, m, k_mask, v1_size, v2_size = variaveis_k(k)
    if k%2!=0:
        assert v1_size > v2_size
    else:
        assert v1_size == v2_size

    trie = TrieBit(m, half_k)
    proc = subprocess.Popen(
        [c_bin,
         genome_path, str(k)],
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
        kmer_idx_init = int.from_bytes(kmer_bytes, byteorder='big')
        kmer_idx = (kmer_idx_init) & (mask)
        assert kmer_idx_init == kmer_idx
        sequences_received += 1
        v1 = kmer_idx >> (v2_size*2)
        v2 = kmer_idx & k_mask
        v1_bits = []
        i = half_k - 1
        while i >= 0:
            v1_bits.append(v1 >> (i*2) & 3)
            i -= 1
        assert len(v1_bits) == half_k
        trie.insert(tuple(v1_bits), v2)
    proc.wait()
    kmers_inserted = trie.count_kmers()
    assert sequences_received == kmers_inserted
    trie.write_bit_format(f"{out_path}_bin")
    print("Trie written in binary.")
    trie.write_txt_format(f"{out_path}_txt")
    print("Trie written in txt.")
    expected = 4**k
    null_count = 4**k - kmers_inserted
    assert expected - null_count == kmers_inserted
    print(f"DEBUG: Null counted: {null_count}.")
    

if __name__ == "__main__":
    main()
