from funcoes_mestrado import *
import re,os,sys,glob
import subprocess
import numpy as np

genoma = "../data/bsubtilis"
k = 12
l = int(k / 2)
m = 4**l
out = "../data/trie_test"

def read_kmers_binary(genome_path, k):
    """Read k-mers in binary format - reading actual bytes"""
    proc = subprocess.Popen(
        ['./fasta_teste_parsing', genome_path, str(k)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    total = 0
    l = k // 2
    trie = inicializar_triebit_teste(l, 4**l)
    
    print(f"DEBUG: Starting to read k-mers, l={l}, m={4**l}")
    
    bits_per_kmer = k * 2  # Each k-mer is k*2 bits = k*2 bytes
    
    # Read all output as bytes
    all_output = proc.stdout.read()
    
    # Process k-mers from the continuous byte stream
    pos = 0
    while pos + bits_per_kmer <= len(all_output):
        # Extract one k-mer (k*2 bytes, each byte is 0 or 1)
        kmer_bytes = all_output[pos:pos + bits_per_kmer]
        
        # Split into v1 and v2 bit sequences
        v1_bits = kmer_bytes[:l * 2]  # First l nucleotides (l*2 bits)
        v2_bits = kmer_bytes[l * 2:]  # Last l nucleotides (l*2 bits)
        
        # Convert v1 bits to base array (for trie navigation)
        v1_bases = []
        for i in range(0, len(v1_bits), 2):
            bit1 = v1_bits[i]      
            bit2 = v1_bits[i+1]    
            nucleotide_code = (bit1 << 1) | bit2
            v1_bases.append(nucleotide_code)
        
        # Convert v2 bits directly back to lexicographic index
        v2_index = 0
        for i, bit in enumerate(v2_bits):
            if bit == 1:
                v2_index |= (1 << (len(v2_bits) - 1 - i))
        
        # Insert into trie
        trie.insert(tuple(v1_bases), v2_index)
        total += 1
        
        # Move to next k-mer
        pos += bits_per_kmer
    
    proc.wait()
    return trie, total

trie, total = read_kmers_binary(genoma, k)
tot = contar_nulomeros_trie_bit_novo(trie,l)
escrever_trie_em_txt_bitarray(trie,out,l)
print(tot)