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

- **`nullomer_statistics_retrieval.py`** – Retrieve nullomer statistics from a bit file. Snakemake-compatible script that computes composition, counter, and motif statistics (CPG, palindromy, homopolymers).

- **`nullomer_data_analysis.py`** – Multiple nullomer data analysis. Processes CSV data and a JSON organism database to generate graphs and derived metrics.

- **`data_analysis_hub.py`** – Generate graphs and CSVs from nullomer data analysis. Supports grouping, parameter selection, and export of processed data.

- **`prime_nullomer_finder.py`** – Compares multiple organisms' nullomer tries to find primes. Searches for nullomers that are prime across organisms in a given group.

- **`create_organism_db.py`** – Creates a JSON file with metadata for a given organism list. Fetches genome length and taxonomy metadata from NCBI.

- **`snakemake_config_generation.py`** – Script to generate a `config.yaml` file for the Snakemake pipeline. Includes integrity checking and configurable k-values.

---

## Detailed documentation

### download_genomes.py

**Description**  
Downloads genome assemblies from the NCBI Entrez API based on a list of accession numbers. After download, it decompresses the archives and converts the sequences to uppercase (required by the main pipeline).

**Usage**

```bash
python download_genomes.py --mode all --accession-list workflow/data/ncbi_dataset.tsv --output workflow/data/genomes/
