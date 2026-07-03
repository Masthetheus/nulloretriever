# Optional Scripts

Inside this directory resides **strictly** auxiliary scripts that are **not** part of the core Snakemake workflow. They are mainly provided as convenience tools for data preparation, format conversion, and post-processing.Stand-alone versions of Snakemake's scripts are also available in this path for manual execution. All scripts are stable and ready for use.

---

## Table of Contents

- [Scripts overview](#scripts-overview)
- [Detailed documentation](#detailed-documentation)
  - [download_genomes.py](#download_genomespy)
  - [nullomer_extraction.py](#nullomer_extractionpy)
  - [nullomer_data_analysis.py](#nullomer_data_analysispy)
  - [nullomer_statistics_retrieval.py](#nullomer_statistics_retrievalpy)
- [Dependencies](#dependencies)
- [Usage notes](#usage-notes)

---

## Scripts overview

- **`genomes_utilities.py`** – Script for direct genome download and decompression via NCBI API. Supports `all` mode (download, decompress, and capitalize) and `capitalize` mode (capitalize existing FASTA files).

- **`nullomer_extraction.py`** – K-mer processing and nullomer extraction used in the Snakemake pipeline. Operates on a given genome and k value, outputting all possible nullomers in binary format.

- **`nullomer_statistics_retrieval.py`** – Retrieve nullomer statistics from a bit file. From the Snakemake pipeline, computes composition, counter, trivial extensions and motif statistics (CPG, palindromy, homopolymers) on top of the previously generated nullomer bit file..

- **`prime_nullomer_finder.py`** – Compares multiple organisms nullomer tries to find primes. Obtain the sequences and write them in a txt file and it's statistics in a centralized csv. 

- **`create_organism_db.py`** – Creates a JSON file with metadata for a given organism list. Fetches genome length and taxonomy metadata from NCBI via API.

- **`snakemake_config_generation.py`** – Script to generate a `config.yaml` file for the Snakemake pipeline. Includes genome integrity checking and configurable k-values.

---

## Detailed documentation

### download_genomes.py

**Description**  
Downloads genome assemblies from the NCBI Entrez API based on a list of accession numbers. After download, it decompresses the archives and converts the sequences to uppercase (required by the main pipeline).

**Usage**

```bash
python download_genomes.py --mode all --accession-list workflow/data/ncbi_dataset.tsv --output workflow/data/genomes/
