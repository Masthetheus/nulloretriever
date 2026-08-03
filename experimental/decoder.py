"""Receive bit files with nullomer and decode them into the original sequences"""

import argparse
import yaml

from pathlib import Path
from nulloretriever.core.triebit_class import TrieBit
from nulloretriever.analysis.processing import mount_trie_from_bitfile

def setup_argparser() -> argparse.ArgumentParser:
    """Argument passing function for the decoder script."""
    parser = argparse.ArgumentParser(description="""Script for nullomer decoding
                                     and sequence format output for downstream
                                     analysis.""")
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
        help="Output for the sequence file. Default = workflow/results",
        default="workflow/results",
    )
    parser.add_argument(
        "-c",
        "--config",
        help="Config file with organisms to convert to complete sequence. Default = workflow/config/config.yaml",
        default = "workflow/config/config.yaml",
    )
    parser.add_argument(
        "-d",
        "--directory",
        help="Directory with bit files to converto to complete sequence. Default = workflow/results",
        default = "workflow/results",
    )
    parser.add_argument(
        "-f",
        "--filename",
        help="Name of the nullomer file, if different from default. Default = null_bit_format",
        default="null_bit_format"
    )
    return parser


def main():
    """Convert bit nullomer file into complete sequences."""
    parser = setup_argparser()
    args = parser.parse_args()
    config = args.config
    directory = args.directory
    output = args.output
    filename = args.filename

    if not config and  not directory:
        print("Please inform a filename or a directory with the files to be analyzed")
        sys.exit(1)
    else:
        with open(config, "r") as f:
            config_data = yaml.safe_load(f)
            try:
                org_set = set(config_data["organisms"])
                k_values = set(config_data["k"])
            except Exception as e:
                print(f"Failure reading the config file, please check its integritiy. Error {e}.")

            for k in k_values:
                direct_k = f"{directory}/k{k}"
                for org in org_set:
                    direct_org = f"{direct_k}/{org}/{filename}"
                    trie = mount_trie_from_bitfile(direct_org)
                    final_out=f"{output}/fullseq_k{k}_{org}"
                    print(f"{org} and count {trie.count_kmers()}")
                    trie.write_sequences(final_out, 1) #1 means null are considered as paths with bit set to 1

if __name__ == "__main__":
    main()
