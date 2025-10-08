"""Script for k-mer processing and nullomer extraction on the snakemake pipeline"""
from nulloretriever.core.triebit_class import TrieBit, TrieBitNode
import argparse
import subprocess

def setup_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Script for k-mer processing and nullomer extraction on the snakemake pipeline")
    parser.add_argument(
        '-k',
        '--kvalues',
        help="Values of k for k-mer extraction. More than one value can be informed. Default = 10",
        nargs="*",
        default=[10]
    )
    parser.add_argument(
        '-m',
        '--mode',
        help="Output mode, if in compact txt format or binary format. Default = binary.",
        choices=['binary','txt'],
        default='binary'
    )
    parser.add_argument(
        '-o',
        '--output',
        help="Output for the found nullomers file. Default = workflow/data/",
        default='workflow/data/'
    )
    return parser
def main():
    parser = setup_argparser()
    args = parser.parse_args()
    k_values = args.kvalues
    mode = args.mode
    out_path = args.output
    # See later how to sync with snakemake pipeline
    genome_path = "workflow/data/genomes/teste"
    for k in k_values:
        k_str = str(k)
        k = int(k)
        l = int(k/2)
        m = 4**l
        trie=TrieBit(m,l)

        proc = subprocess.Popen(
            ['workflow/scripts/c/bin_fasta_kmer_extraction', genome_path, k_str],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            bufsize=65536
        )
        total = 0
        bytes_per_half = (l * 2 + 7) // 8
        bytes_per_kmer = bytes_per_half * 2
        bytes_per_sequence = (k * 2 + 7) // 8
        k_mask = (2**k) - 1
        while True:
            kmer_bytes = proc.stdout.read(bytes_per_sequence)
            if len(kmer_bytes) < bytes_per_sequence:
                break
            
            print(f"DEBUG PY: Read {len(kmer_bytes)} bytes: {[b for b in kmer_bytes]}")

            kmer_idx = int.from_bytes(kmer_bytes)
            v1 = kmer_idx >> k
            v2 = kmer_idx & k_mask
            v1_bits = []
            i = (k//2) - 1
            while i >= 0:
                v1_bits.append(v1 >> (i*2) & 3)
                i -= 1
            print(f"kmer_idx: {kmer_idx} v1: {v1} v2: {v2}")
            print(v1_bits)
            trie.insert(tuple(v1_bits), v2)
            total += 1
        proc.wait()
        print(trie.count_nullomers())
if __name__ == "__main__":
    main()
