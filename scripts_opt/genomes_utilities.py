"""Optional script for direct genome download via NCBI API."""

from nulloretriever.utils.validation import get_valid_email, get_valid_tool
from nulloretriever.data.ncbidownload import download_genome_bioentrez
from nulloretriever.data.ncbiapidata import read_accession_list
from nulloretriever.utils.integrity import unzip_fasta_file, capslock_file
from Bio import Entrez
import argparse
from pathlib import Path


def setup_argparser() -> argparse.ArgumentParser:
    """Parsing function for the genome utilities script."""
    parser = argparse.ArgumentParser(
        description="Download" "genomes from NCBI by" "assembly accession list."
    )
    parser.add_argument(
        "--accession-list",
        type=str,
        default="workflow/data/ncbi_dataset.tsv",
        help="Path to the file containing assembly accession numbers"
        "(default: ncbi_dataset.tsv)",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="workflow/data/genomes/",
        help="Where to store downloaded genomes."
        "Default = workflow/data/genomes/",
    )
    parser.add_argument(
        "--mode",
        type=str,
        choices=["all", "capitalize"],
        help="Mode for running the script. Options all and capitalize. All downloads"
        "the genome list from NCBI, decompress and capitalize them. Capitalize only"
        "capitalizes all genomes on --output directory."
        "Default = all.",
        default="all",
    )
    return parser


def main():
    """Download the genomes files from given set of accession codes.

    Args:
        Entrez.email(str): e-mail registered on Entrez
        Entrez.tool(str): personal/lab use tool registered on Entrez via e-mail
        accession-list(file): list of NCBI assembly accession numbers
        output: place to store the gathered files
    Returns:
        directory(file): file containing a genome sequence unzipped and in full
                         capslock.
    """
    parser = setup_argparser()
    args = parser.parse_args()
    accession_file = args.accession_list
    output = Path(args.output)
    mode = args.mode
    if mode == "all":
        Entrez.email = get_valid_email()
        Entrez.tool = get_valid_tool()
        print(f"Using accession file: {accession_file}")
        print(
            "Do you want to specify a custom .csv/.tsv column name for the"
            " accession numbers? (y/[n])"
        )
        # remove until the else block to skip verification
        custom_column = input().strip().lower() == "y"
        if custom_column:
            column_name = (
                input(
                    "Enter the column name (NCBI default is 'Assembly"
                    " Accession'): "
                ).strip()
                or "Assembly Accession"
            )
        else:
            column_name = None
        accessions = read_accession_list(accession_file, column=column_name)
        print(f"Found {len(accessions)} accessions in {accession_file}")
        directories = download_genome_bioentrez(accessions, output)
        print("All genomes downloaded successfully.")
        print("Initializing genomes decompression and capitalization!")
        for directory in directories:
            out_path = unzip_fasta_file(directory)
            capslock_file(out_path)
            accession = out_path.stem
            new_path = out_path.parent / f"{accession}.fasta"
            out_path.rename(new_path)
        print("All genomes decompressed and capitalized!")
        print(f"Downloaded genomes can be found in the {output} directory.")
    elif mode == "capitalize":
        target_dir = Path(args.output)
        if not target_dir.is_dir():
            print(f"Directory {target_dir} not found.")
            return
        fasta_files = list(target_dir.glob("*.fasta")) + list(
            target_dir.glob("*.fa")) + list(target_dir.glob("*.fna"))
        if not fasta_files:
            print(f"No FASTA files found in {target_dir}")
            return
        for fasta in fasta_files:
            capslock_file(fasta)
        print(f"All FASTA files in {target_dir} have been capitalized.")
    else:
        print("Mode not recognized, please check the available options.")


if __name__ == "__main__":
    main()
