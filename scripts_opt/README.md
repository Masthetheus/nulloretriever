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
Detailed documentation

### genomes_utilities.py

Description:
Optional script for direct genome download via NCBI API. Supports 'all' mode (download, decompress, and capitalize) and 'capitalize' mode (capitalize existing FASTA files).

Usage:
python genomes_utilities.py --mode all --accession-list workflow/data/ncbi_dataset.tsv --output workflow/data/genomes/

Arguments:
  --mode            (required) all: download, decompress, and capitalize; capitalize: only capitalize existing FASTA files in the output directory.
  --accession-list  (default: workflow/data/ncbi_dataset.tsv) Path to the file containing assembly accession numbers.
  --output          (default: workflow/data/genomes/) Directory where genomes are stored.
  --column          (optional) Name of the column containing accession numbers. If not provided, the script asks interactively.

Example:
  **Download and prepare genomes from a custom accession list. Accepts txt files, where one accession code is informed per line**
  python genomes_utilities.py --mode all --accession-list my_accessions.txt --output data/genomes/

  **Only capitalize FASTA files already present in the directory**
  python genomes_utilities.py --mode capitalize --output data/genomes/

Output:
Each genome is saved as {accession}.fasta in the specified output directory, with all bases in uppercase.

Notes:
  - The script prompts for an NCBI-registered email address and tool, to comply with Entrez usage policies.
  - Rate limiting is handled by the Bio.Entrez module.

### nullomer_extraction.py

Description:
K-mer processing and nullomer extraction used in the Snakemake pipeline. Operates on a given genome and k value, outputting nullomers in binary trie format.

Usage:
python nullomer_extraction.py -k 15 -m binary -g genome.fasta -o output/

Arguments:
  -k, --kvalues    (default: [10]) Values of k for k-mer extraction. More than one value can be informed.
  -g, --genome     (required) Genome FASTA file to be analyzed.
  -o, --output     (default: workflow/results/) Output directory for the nullomer files.

Example:
  python nullomer_extraction.py -k 15 16 -g data/GCF_000005845_2.fasta -o results/

Output:
  - binary mode: .bit file containing the Trie structure.

Notes:
  - The script assumes the input genome is already in uppercase.
  - For large genomes and high k values, memory usage can be high and the running time as well.


### nullomer_statistics_retrieval.py

Description:
Retrieve nullomer statistics from a trie bit file. Can compute composition, total nullomers, motif statistics (CPG, palindromy, homopolymers) and trivial nullomer occurence.

Usage:
**Executes the complete analysis on nullomers from file n1**
python nullomer_statistics_retrieval.py -n1 nullomers.bit -n2 previous.bit -s all -o stats.csv

Arguments:
  -n1, --null1     (required) Nullomer bit file to be analyzed.
  -n2, --null2     (optional) Nullomer bit file of the same organism for k-1. Used to obtain trivial extensions.
  -s, --stats      (default: all) Stats to be analyzed. Choices: all, composition, trivial, motifs.
  -o, --output     (default: workflow/results/custom_statistics_retrieval.csv) Output CSV file for the statistics.

Example:
  python nullomer_statistics_retrieval.py -n1 results/ecoli_k15.bit -n2 results/ecoli_k14.bit -s composition,motifs -o analysis/ecoli_stats.csv

Output:
CSV file with columns including, depending on selected mode:
  - counter: number of nullomers
  - composition_gc_mean: mean GC content
  - motifs_cpg: CPG-related statistics
  - motifs_palindromy: palindromy-related statistics
  - motifs_homopolymers: homopolymer-related statistics
  - trivial: trivial extension statistics (if --null2 provided)

Notes:
  - The trivial stat requires --null2 to be provided.

### prime_nullomer_finder.py

Description:
Compares multiple organisms' nullomer tries to find primes. Searches for nullomers that are prime across organisms in a given group, obtaining the sequences and writing them in a txt file and their statistics in a centralized CSV.

Usage:
python prime_nullomer_finder.py -c config.yaml -o output/

Arguments:
  -c, --config     (default: workflow/config/config.yaml) Relative path to the YAML config file containing the organism names and nullomer file paths.
  -o, --output     (required) Directory where output files will be saved.

Example:
  python prime_nullomer_finder.py -c workflow/config/config.yaml -o results/primes/

Output:
  - prime_nullomers.txt: list of nullomer sequences that are prime across the specified organisms.
  - prime_statistics.csv: centralized CSV with statistics for each prime nullomer.

Dependencies:
  - nulloretriever package
  - pyyaml
  - json
  - csv

Notes:
  - The config file should follow the same structure as the Snakemake configuration.
  - The script expects nullomer bit files for each organism listed in the config.


### create_organism_db.py

Description:
Creates a JSON file with metadata for a given organism list. Fetches genome length and taxonomy metadata from NCBI via API.

Usage:
python create_organism_db.py -c config.yaml -o organism_db.json

Arguments:
  -c, --config     (required) Relative path to the YAML config file containing the organism list.
  -o, --output     (default: organism_db.json) Output JSON file for the organism metadata.

Example:
  python create_organism_db.py -c workflow/config/config.yaml -o data/organism_db.json

Output:
JSON file with organism metadata including:
  - accession: assembly accession
  - organism_name: scientific name
  - taxonomy_id: NCBI taxonomy ID
  - genome_length: total genome length in bp
  - lineage: taxonomic lineage

Dependencies:
  - Biopython (for NCBI Entrez queries)
  - pyyaml
  - json

Notes:
  - The script prompts for an NCBI-registered email address.
  - Rate limiting is handled by the Bio.Entrez module.


### snakemake_config_generation.py

Description:
Script to generate a config.yaml file for the Snakemake pipeline. Includes genome integrity checking and configurable k-values.

Usage:
python snakemake_config_generation.py --name my_config.yaml --genomes-dir workflow/data/genomes/ --kvalues 15 16 17

Arguments:
  --name           (default: config.yaml) Name of the config file to generate.
  --genomes-dir    (default: workflow/data/genomes/) Directory containing the genome FASTA files.
  --kvalues        (required) List of k-values to be used in the pipeline.
  --output-dir     (default: workflow/results/) Directory for pipeline results.

Example:
  python snakemake_config_generation.py --name config_ecoli.yaml --genomes-dir data/genomes/ --kvalues 15 16 17 --output-dir results/ecoli/

Output:
YAML configuration file with the following structure:
  genomes_dir: "workflow/data/genomes/"
  output_dir: "workflow/results/"
  k_values: [15, 16, 17]
  genomes:
    - GCF_000005845.2
    - GCF_000001405.40

Dependencies:
  - Python standard library (argparse, pathlib, sys)
  - pyyaml
  - nulloretriever package (for integrity checking)

Notes:
  - The script performs integrity checking on the genome files before generating the config.
  - Remember to rename the generated file to config.yaml (or update the Snakefile) before running the pipeline.

