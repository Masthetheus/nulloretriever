"""General utils for NCBI API data gathering."""
import csv
import os
import time
import xml.etree.ElementTree as ET

from Bio import Entrez


def read_accession_txt(filepath):
    """Read a text file containing a list of NCBI accession numbers.

    Each line should contain a valid accession number, starting with "GCF_" or "GCA_".
    Ignores empty lines and comments (lines starting with "#").

    Args:
        filepath (str): Path to the text file.

    Returns:
        list: List of valid accession numbers found in the file.
    """
    accessions = set()
    with open(filepath) as f:
        for line in f:
            acc = line.strip()
            if not acc or acc.startswith("#"):
                continue
            if not acc.startswith(("GCF_", "GCA_")):
                print(f"Warning: strange code found in: '{acc}', ignored.")
                continue
            accessions.add(acc)
    return list(accessions)


def read_accession_list(filepath, column=None):
    """Read a list of NCBI accession numbers from a file.

    Supports .txt, .csv, and .tsv files.
    If the column is not specified, tries to automatically detect the column
    containing the accession numbers.

    Args:
        filepath (str): Path to the file containing accession numbers.
        column (int or str, optional): Index or name of the column containing the accession numbers.
                                       If None, tries to detect automatically.

    Returns:
        list: List of accession numbers found in the file.
    """
    ext = os.path.splitext(filepath)[1].lower()
    accessions = []
    if ext == ".txt":
        accessions = read_accession_txt(filepath)
    elif ext in [".csv", ".tsv"]:
        delimiter = ',' if ext == ".csv" else '\t'
        with open(filepath, newline='') as f:
            reader = csv.reader(f, delimiter=delimiter)
            header = next(reader)
            # Detect column index
            if column is None:
                for i, col in enumerate(header):
                    if "assembly" in col.lower() and ("accession" in col.lower()):
                        column = i
                        break
                else:
                    column = 0
            elif isinstance(column, str):
                column = header.index(column)
            for row in reader:
                if row and len(row) > column:
                    accessions.append(row[column].strip())
    else:
        raise ValueError(
            "Unsupported file format. Use .txt, .csv, or .tsv files.")
    return accessions


def get_accesion_summary_data(acc):
    """Return xml summary data for given accession.

    Args:
        acc (str): Single accession code.
    Returns:
        summary (xml): xml summary of given accesion code ID.
    """
    log = "accession_summary_data_entrez.log"
    try:
        # Search the ID on the DB
        print(f"Searching ID for Assembly Accession: {acc}")
        handle = Entrez.esearch(db="assembly", term=acc, retmode="xml")
        record = Entrez.read(handle)
        handle.close()

        time.sleep(0.35)

        if not record["IdList"]:
            print(f"No result found for {acc}")
            with open(log, 'a') as log:
                log.write(f"ID for Assembly Accession {acc} not found.\n")
            return None

        # Get's the first ID
        assembly_id = record["IdList"][0]
        print(f"Assembly ID found: {assembly_id}")

        # Get assembly summary
        print(f"Searching summary for Assembly ID: {assembly_id}")
        handle = Entrez.esummary(db="assembly", id=assembly_id, retmode="xml")
        summary = Entrez.read(handle)
        handle.close()

        time.sleep(0.35)

        return summary
    except Exception as e:
        print(f"Error during {acc}: {e}")
        return None


def get_genome_download_link(accessions):
    """Return a list of genome Entrez links for given assembly accessions.

    Args:
        accessions (str or list): accession codes of the desired genomes.
    Returns:
        links (list): list of the download links of the desired genomes.
    """
    log = "genome_download.log"
    if isinstance(accessions, str):
        accessions = [accessions]
    links = {}
    try:
        for accession in accessions:
            # Obtain the xml data for given accession code
            summary = get_accesion_summary_data(accession)
            if not summary:
                print("Vai quebrar")
                break
            # Obtain FTP Assembly link
            ftp_path = summary['DocumentSummarySet']['DocumentSummary'][0].get(
                'FtpPath_RefSeq')
            if not ftp_path:
                ftp_path = summary['DocumentSummarySet']['DocumentSummary'][0].get(
                    'FtpPath_GenBank')
            if ftp_path:
                link = ftp_path + "/" + \
                    ftp_path.split("/")[-1] + "_genomic.fna.gz"
                link = link[3:]
                link = 'https' + link
                links[accession] = link
            else:
                print(f"Link FTP not found for {accession}")
                with open(log, 'a') as log:
                    log.write(f"Link FTP not found for {accession}.\n")
                return None
        return links
    except Exception as e:
        print(f"Error during {accession}: {e}")
        return None


def get_genome_metadata(accessions, params=None):
    """Return a list with certain metadata for given accession codes.

    Args:
        accessions (str or list): accessions codes for metadata retrieval
        params (str or list): params to be retrieved
    Returns:
        metadata (list): list with needed metadata of given accession code.
    """
    if isinstance(accessions, str):
        accessions = [accessions]
    if not params:
        params = [
            'Taxid',
            'SpeciesTaxid',
            'SpeciesName',
            'AssemblyStatus',
            'Meta'
        ]
    elif isinstance(params, str):
        params = [params]
    metadata = {}
    try:
        for acc in accessions:
            summary = get_accesion_summary_data(acc)
            if not summary:
                continue
            data = {}
            for param in params:
                data[param] = summary['DocumentSummarySet']['DocumentSummary'][0].get(
                    param)
            metadata[acc] = data
        return metadata
    except Exception as e:
        print(f"Error during {acc}: {e}")
        return None


def get_genome_length(metadata):
    old_metadata = metadata
    for organism in old_metadata:
        xml_content = f"<Root>{old_metadata[organism]['Meta']}</Root>"
        try:
            root = ET.fromstring(xml_content)
            genome_total_length = root.find(
                ".//Stat[@category='total_length']")

            if genome_total_length is not None:
                total_length = genome_total_length.text
                metadata[organism]['Genome length'] = total_length
                metadata[organism].pop("Meta")
        except ET.ParseError as e:
            print(f"Error analyzing Meta string: {e}")
    return metadata


def get_taxonomy_metadata(metadata):
    old_metadata = metadata
    for organism in old_metadata:
        fetch = Entrez.efetch(
            id=metadata[organism]['Taxid'], db='taxonomy', retmode='xml')
        data = Entrez.read(fetch)
        fetch.close()
        taxonomy_data = {d['Rank']: d['TaxId']
                         for d in data[0]['LineageEx']}
        family_id = taxonomy_data.get('family', 'N/A')
        if family_id is not None:
            metadata[organism]['family_id'] = family_id
    return metadata
