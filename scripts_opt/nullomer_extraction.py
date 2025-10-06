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
    genome_path = "../workflow/data/genomes/teste"
    for k in k_values:
        k_str = str(k)
        k = int(k)
        l = int(k/2)
        m = 4**l
        trie=TrieBit(m,l)

        proc = subprocess.Popen(
            ['../workflow/scripts/c/newversion_kmer_extraction', genome_path, k_str],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            bufsize=65536
        )
        bytes_per_half = (l * 2 + 7) // 8
        bytes_per_kmer = bytes_per_half * 2

        while True:
            kmer_bytes = proc.stdout.read(bytes_per_kmer)
            if len(kmer_bytes) < bytes_per_kmer:
                break
            
            print(f"DEBUG PY: Read {len(kmer_bytes)} bytes: {[b for b in kmer_bytes]}")
            
            v1 = int.from_bytes(kmer_bytes[:bytes_per_half], byteorder='big')
            v2 = int.from_bytes(kmer_bytes[bytes_per_half:], byteorder='big')
            
            full_kmer = (v1 << (l * 2)) | v2
            
            print(f"DEBUG PY: v1={v1}, v2={v2}, full_kmer={full_kmer}")

        proc.wait()

if __name__ == "__main__":
    main()