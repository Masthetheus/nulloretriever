"""Translates nullomer primes to AA sequences."""
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import numpy as np

file_primes = "prime_null_all_14"

with open(file_primes, "r") as f:
    primes = [line.strip() for line in f]
from Bio.Seq import Seq

aa_sequences = []

for i, prime in enumerate(primes):
    seq = Seq(prime)
    gc = gc_fraction(seq) * 100
    aa = seq.translate()
    stop_count = aa.count('*')
    print(f"Nullomer {i+1}: GC={gc:.1f}% | AA={aa} | stops={stop_count}")
gcs = [gc_fraction(Seq(p))*100 for p in primes]
print(f"\nGC médio dos 24: {np.mean(gcs):.1f}%")
print(f"GC mediano: {np.median(gcs):.1f}%")
for prime in primes:
    seq = Seq(prime)
    rc = seq.reverse_complement()
    for frame in range(3):
        aa_sequences.append({
            'dna': prime,
            'frame': f'+{frame+1}',
            'aa': seq[frame:].translate()
        })
        aa_sequences.append({
            'dna': prime,
            'frame': f'-{frame+1}',
            'aa': rc[frame:].translate()
        })
print(aa_sequences)
