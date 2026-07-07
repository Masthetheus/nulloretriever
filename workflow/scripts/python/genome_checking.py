from nulloretriever.utils.integrity import check_genome_integrity
import sys
import argparse

def setup_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description='Script for genome integrity checking.')
    parser.add_argument(
        '--genome',
        help='Path to the configuration file to check'
    )
    parser.add_argument(
        '--out',
        help='Path to output of genome verification'
    )
    return parser
def main():
    parser = setup_argparser()
    args = parser.parse_args()
    genome = args.genome
    out = args.out
    capslock_file(genome)
    check = check_genome_integrity(genome)
    if check:
        with open(out+'_ok', 'w') as f:
            f.write("Genome integrity checked and approved!")
    else:
        with open(out+'_nok', 'w') as f:
            f.write("Genome integrity checked and not approved!")

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)

