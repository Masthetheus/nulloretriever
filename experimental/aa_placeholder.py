"""Translates nullomer primes to AA sequences."""
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import numpy as np

file_primes = "prime_null_all"
with open(file_primes, "r") as f:
    primes = [line.strip() for line in f]
k = 13
aa_sequences = []

aa_sequences = []
for prime in primes:
    seq = Seq(prime)
    for start in range(len(prime) - 2):
        fragment = seq[start:]
        if len(fragment) >= 3:
            aa = fragment.translate()
            aa_sequences.append({
                'dna': prime,
                'frame_start': start,
                'aa': aa,
                'aa_no_stop': str(aa).replace('*', '')
            })

for i, prime in enumerate(primes):
    seq = Seq(prime)
    gc = gc_fraction(seq) * 100
    aa = seq.translate()
    stop_count = aa.count('*')
    print(f"Nullomer {i+1}: GC={gc:.1f}% | AA={aa} | stops={stop_count}")

gcs = [gc_fraction(Seq(p))*100 for p in primes]
print(f"\nGC médio dos 24: {np.mean(gcs):.1f}%")
print(f"GC mediano: {np.median(gcs):.1f}%")

unique_peptides = set()
with open('aas_order_all','w') as f:
    for aas in aa_sequences:
        #f.write(f"{aas['dna']}, frame:{aas['frame']}\n")
        if len(aas['aa_no_stop']) >= 4:
            unique_peptides.add(aas['aa_no_stop'])
            f.write(f"{str(aas['aa_no_stop'])}\n")
print(len(unique_peptides))
