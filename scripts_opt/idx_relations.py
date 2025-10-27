 #!/usr/bin/env python3
import subprocess
import matplotlib.pyplot as plt

# idx = input("What's the sequence index?")
# k = int(input("What's the k value?"))
k = int(12)
pos = k - 1
k_str = str(k)

ranges = {}

for val in range(k):
    if val//2 == 0 :
        continue
    limit = 3*(4**val)
    ranges[limit] = 0

limit = (4**k)
ranges[limit] = 0
print(ranges)
genome_path = "workflow/data/genomes/GCA_001010845_3"
proc = subprocess.Popen(
    ['workflow/scripts/c/bin_fasta_kmer_extraction',
     genome_path, k_str],
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    bufsize=65536
)
sequences_received = 0
total = 0
bytes_per_sequence = (k * 2 + 7) // 8
wasted_space = (bytes_per_sequence * 8) - (k * 2)
k_mask = (2**k) - 1
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
    kmer_idx = (kmer_idx >> wasted_space) & (mask)
    sequences_received += 1
    for key in ranges.keys():
        if kmer_idx <= key:
            ranges[key] += 1
            break
proc.wait()
ranges_percent = {}
i = 0
for key in ranges.keys():
    percentual = (ranges[key]*100)/sequences_received
    ranges_percent[i] = percentual
    i += 1
names = ranges_percent.keys()
values = ranges_percent.values()
plt.bar(names, values, color = "purple")
plt.show()
print(sequences_received)
print(ranges)
print(ranges_percent)
relative_idxs = {}

