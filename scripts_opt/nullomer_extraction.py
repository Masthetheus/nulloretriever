"""K-mer processing and nullomer extraction on the snakemake pipeline."""

import argparse
import subprocess

from nulloretriever.core.triebit_class import TrieBit


def setup_argparser() -> argparse.ArgumentParser:
    """Argument passing function for the kmer extraction script."""
    parser = argparse.ArgumentParser(description="""Script for k-mer processing
                                     and nullomer extraction on the snakemake
                                     pipeline""")
    parser.add_argument(
        "-k",
        "--kvalues",
        help="Values of k for k-mer extraction. More than one value can be"
        "informed. Default = 10",
        nargs="*",
        default=[10],
    )
    parser.add_argument(
        "-m",
        "--mode",
        help="Output mode. Available options are compact txt format with"
        "trie indexes, complete txt format with full sequences or binary"
        "format."
        " Default = binary.",
        choices=["binary", "sequence"],
        default="binary",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output for the found nullomers file. Default = workflow/data/",
        default="workflow/data",
    )
    parser.add_argument("-g", "--genome", help="Genome to be analyzed")
    parser.add_argument("-c", "--c_script", help = "Path to c kmer extraction binary.", default = "workflow/scripts/c/bin/kmer_extractor")
    return parser


def main():
    """K-mer extraction and nullomer trie generation."""
    parser = setup_argparser()
    args = parser.parse_args()
    c_bin = args.c_script
    k_values = args.kvalues
    mode = args.mode
    out_path = args.output + "_result"
    genome_path = args.genome
    for k in k_values:
        k = int(k)
        is_odd = k % 2
        half_k = (k // 2) + is_odd
        m = 4 ** (half_k - is_odd)
        v2_size = k - half_k
        k_mask = (1 << (v2_size * 2)) - 1
        trie = TrieBit(m, k, half_k)
        proc = subprocess.Popen(
            [
                c_bin,
                genome_path,
                str(k),
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            bufsize=65536,
        )
        sequences_received = 0
        bytes_per_sequence = (k * 2 + 7) // 8
        wasted_space = (bytes_per_sequence * 8) - (k * 2)
        mask = (4**k) - 1
        print(
            f"DEBUG: Variables in use:\n"
            f"bytes_per_sequence: {bytes_per_sequence}"
            f"\nwasted_space:{wasted_space}"
            f"\nk_mask: {k_mask}"
            f"\nmask: {bin(mask)}"
        )

        v2_buffer = [[] for _ in range(1 << (2 * half_k))]
        buffered_count = 0
        FLUSH_LIMIT = 1_000_000

        while True:
            kmer_bytes = proc.stdout.read(bytes_per_sequence)

            if len(kmer_bytes) < bytes_per_sequence:
                break

            kmer_idx = int.from_bytes(kmer_bytes, byteorder="big")
            kmer_idx = (kmer_idx) & (mask)
            sequences_received += 1

            v1 = kmer_idx >> (v2_size * 2)
            v2 = kmer_idx & k_mask

            v2_buffer[v1].append(v2)
            buffered_count += 1

            if buffered_count >= FLUSH_LIMIT:
                for v1, v2_list in enumerate(v2_buffer):
                    if v2_list:
                        trie.insert_v2_list(v1,v2_list)
                        v2_buffer[v1] = []
                buffered_count = 0
        proc.wait()

        for v1, v2_list in enumerate(v2_buffer):
            if v2_list:
                trie.insert_v2_list(v1, v2_list)

        write_dict = {
            "binary": trie.write_bit_format,
            "sequence": trie.write_sequences,
        }
        if mode == "sequence":
            identifier = 0
            write_dict[mode](out_path, identifier)
        else:
            write_dict[mode](out_path)
        kmers_inserted = trie.count_kmers()
        expected = 4**k
        null_count = expected - kmers_inserted
        print(
            f"{trie.count_kmers()} were inserted into the trie."
            f"\n{expected} total were expected."
            f"An total of {null_count} nullomers were found "
            "for a k of {k}."
        )


if __name__ == "__main__":
    main()
