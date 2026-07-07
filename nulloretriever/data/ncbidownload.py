"""Module for automatic genome NCBI download."""

import os
import requests
import time
from nulloretriever.utils.progress_bar import progress_bar
from nulloretriever.data.ncbiapidata import get_genome_download_link
import sys


def download_genome_bioentrez(accessions, out):
    """Download a genome file from NCBI using Bio.Entrez.

    Args:
        assembly_accession (str): Assembly accession number.
        out (str): Path to the output directory where the genome will be saved.
        organism (str): Name of the organism, used for naming the file.
    Returns:
        directiores(arr): donwloaded files path
    """
    # Obtain genome links from Entrez
    links = get_genome_download_link(accessions)
    cont = 1  # For the progress bar
    start = time.time()  # For the progress bar
    if len(links) != len(accessions):
        if not links:
            print("No valid links found for the provided accessions.")
            return
        else:
            print(
                "Some accessions could not be found. Check the log file (data/logs/entrez.log) for details."
            )

    # Makes sure output directory exists
    if not os.path.exists(out):
        print(f"Error: Output directory '{out}' does not exist!")
        print(f"Trying to create {out}")
        try:
            os.makedirs(out)
            print("Directory created!")
        except:
            sys.exit(1)
    directories = []
    try:
        for organism, link in links.items():
            file_name = organism.replace(".", "_")
            # Initiate download
            response = requests.get(link, stream=True)
            if response.status_code == 200:
                file_path = out / f"{file_name}.gz"
                with open(file_path, "wb") as f:
                    for chunk in response.iter_content(chunk_size=1024):
                        f.write(chunk)
                directories.append(file_path)
            else:
                print(f"Error during download of {organism}: {response.status_code}")
            progress_bar(cont, len(links), start=start)
            cont += 1
    except Exception as e:
        print(f"Error during download process: {e}")
        sys.exit(1)
    return directories
