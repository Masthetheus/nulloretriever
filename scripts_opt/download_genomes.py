from nulloretriever.utils.validation import get_valid_email, get_valid_tool
from nulloretriever.data.ncbidownload import read_accession_list, download_genome_bioentrez
from nulloretriever.utils.integrity import unzip_fasta_file, capslock_file
from Bio import Entrez
import argparse
from pathlib import Path
def main():
    Entrez.email = get_valid_email()
    Entrez.tool = get_valid_tool()
    parser = argparse.ArgumentParser(description="Download genomes from NCBI by assembly accession list.")
    parser.add_argument(
        "--accession-list",
        type=str,
        default="ncbi_dataset.tsv",
        help="Path to the file containing assembly accession numbers (default: ncbi_dataset.tsv)"
    )
    parser.add_argument(
        "--output",
        type=str,
        default="workflow/data/genomes/",
        help="Where to store downloaded genomes"
    )
    args = parser.parse_args()
    accession_file = args.accession_list
    output = Path(args.output)
    print(f"Using accession file: {accession_file}")
    print("Do you want to specify a custom .csv/.tsv column name for the accession numbers? (y/[n])")
    custom_column = input().strip().lower() == 'y'
    if custom_column:
        column_name = input("Enter the column name (NCBI default is 'Assembly Accession'): ").strip() or "Assembly Accession"
    else:
        column_name = None
    acessions = read_accession_list(accession_file, column=column_name)
    print(f"Found {len(acessions)} accessions in {accession_file}")
    directories = download_genome_bioentrez(acessions, output)
    print("All genomes downloaded successfully.")
    print("Initializing genomes decompression and capitalization!")
    for directory in directories:
        out_path = unzip_fasta_file(directory)
        capslock_file(out_path)
    print("All genomes decompressed and capitalized!")
    print(f"Downloaded genomes can be found in the {output} directory.")
if __name__ == "__main__":
    main()