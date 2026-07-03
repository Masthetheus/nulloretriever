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

- **`genome_utilities.py`** – Download genome assemblies from NCBI.
  *Input:* Accession list (`.tsv`/`.csv`).
  *Output:* FASTA files in output directory.

- **`nullomer_extraction.py`** – Extract nullomers from a single genome.
  *Input:* FASTA file, k value.
  *Output:* List of nullomers (CSV).

- **`nullomer_data_analysis.py`** – Analyze nullomer distribution, GC content, and feature overlap.
  *Input:* Nullomer list, genome annotation (optional GFF).
  *Output:* Summary tables and plots.

- **`nullomer_statistics_retrieval.py`** – Aggregate results across multiple genomes into one table.
  *Input:* Multiple nullomer lists.
  *Output:* Combined CSV/TSV.

---

## Detailed documentation

### download_genomes.py

**Description**  
Downloads genome assemblies from the NCBI Entrez API based on a list of accession numbers. After download, it decompresses the archives and converts the sequences to uppercase (required by the main pipeline).

**Usage**

```bash
python download_genomes.py --mode all --accession-list workflow/data/ncbi_dataset.tsv --output workflow/data/genomes/
