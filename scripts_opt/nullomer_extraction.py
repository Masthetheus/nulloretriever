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

    for k in k_values:
        k = int(k)
        l = int(k/2)
        m = 4**l

        trie=TrieBit(m,l)

    trie.insert(v1,v2)
print(len(indexes))
expected_null = (4**k) - len(indexes)
print("Checking if all nodes were inserted correctly via nullomer count.")
nullomers = trie.count_nullomers()
print(f"{nullomers} nullomers were found of the {expected_null} expected.")
output_bit = "trie_binary_test"
output_compact_txt = "trie_compact_txt"
print(f"Trying now to write the trie in bit format on the following path: {output_bit}")
trie.write_bit_format(output_bit)
print("Trie in bit format correctly wrriten!")
print(f"Trying now to write the trie in compact txt format on the following path: {output_compact_txt}")
trie.write_compact_txt_format(output_compact_txt)
print("Trie in compact txt format correctly wrriten!")
print("All TrieBit base functionalities tested and passed!")