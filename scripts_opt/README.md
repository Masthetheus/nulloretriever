# Optional Scripts

Inside this directory resides **strictly** auxiliary scripts that are **not** part of the core Snakemake workflow. They are mainly provided as convenience tools for data preparation, format conversion, and post-processing.Stand-alone versions of Snakemake's scripts are also available in this path for manual execution. All scripts are stable and ready for use.

---

## Table of Contents

- [Scripts overview](#scripts-overview)
- [Detailed documentation](#detailed-documentation)
  - [create_organism_db.py](#create_organism_dbpy)
  - [genomes_utilities.py](#genomes_utilitiespy)
  - [nullomer_extraction.py](#nullomer_extractionpy)
  - [nullomer_order.c](#nullomer_orderc)
  - [nullomer_statistics_retrieval.py](#nullomer_statistics_retrievalpy)
  - [prime_nullomer_finder.py](#prime_nullomer_finderpy)
  - [snakemake_config_generation.py](#snakemake_config_generationpy)
  - [translate_to_bit.py](#translate_to_bitpy)
- [Dependencies](#dependencies)
- [Usage notes](#usage-notes)

---

## Scripts overview

- **`create_organism_db.py`** – Creates a JSON file with metadata for a given organism list. Fetches genome length and taxonomy metadata from NCBI via API.

- **`genomes_utilities.py`** – Script for direct genome download and decompression via NCBI API. Supports `all` mode (download, decompress, and capitalize) and `capitalize` mode (capitalize existing FASTA files).

- **`nullomer_extraction.py`** – K-mer processing and nullomer extraction used in the Snakemake pipeline. Operates on a given genome and k value, outputting all possible nullomers in binary format.

- **`nullomer_order.c`** - Retrieves nullomeric order for all sequences on a compact bit file, outputting one file per order in sequence format.

- **`nullomer_statistics_retrieval.py`** – Retrieve nullomer statistics from a bit file. From the Snakemake pipeline, computes composition, counter, trivial extensions and motif statistics (CPG, palindromy, homopolymers) on top of the previously generated nullomer bit file..

- **`prime_nullomer_finder.py`** – Compares multiple organisms nullomer tries to find primes. Obtain the sequences and write them in a txt file and it's statistics in a centralized csv. 

- **`snakemake_config_generation.py`** – Script to generate a `config.yaml` file for the Snakemake pipeline. Includes genome integrity checking and configurable k-values.

- **`translate_to_bit.py`** - Translates a nullomer sequence file, with one sequence per line, into a compact bit file. Can be used to compress nullomer data in disk or to allow downstream processing within the pipeline.

---

## Detailed documentation

### create_organism_db.py

Description:
Creates a JSON file with metadata for a given organism list. Fetches genome length and taxonomy metadata from NCBI via API.

Usage:
python create_organism_db.py [-h] [-c CONFIG] [-o OUTPUT]

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

Notes:
  - The script prompts for an NCBI-registered email address.
  - Rate limiting is handled by the Bio.Entrez module.


### genomes_utilities.py

Description:
Optional script for direct genome download via NCBI API. Supports 'all' mode (download, decompress, and capitalize) and 'capitalize' mode (capitalize existing FASTA files).

Usage:
python genomes_utilities.py [-h] [--accession-list ACCESSION_LIST] [--output OUTPUT]
                            [--mode {all,capitalize}]

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
python nullomer_extraction.py [-h] [-k [KVALUES ...]] [-o OUTPUT] [-g GENOME]


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

### nullomer_order.c

Description:
From a compact bit file containing an organism nullomers, obtain each nullomer order up to a limit established by the user.

Usage: ./scripts-opt/bin/null_order <file> <order>

Arguments:
  <file>          (required) Compact bit file containing the nullomers.
  <order>         (required) Up limit of order to calculate.

Example:
  Obtain nullomers pertaining to order 1 and order 2:
  ./scripts-opt/bin/null_order results/ecoli_k15.bit 2

Output:
One TXT file per order with at least one sequence, containing all sequences, one per line, of such order. 

Notes:
  - The script needs a compact bit file for the order extraction work as intended. If only full sequence is available, check [translate_to_bit.py](#translate_to_bitpy).
  - The scripts optimizes the traditional brute-force by operating in a cummulative matter. That being, the search of given order occurs given the found sequences from previous order, until no more candidates are left.
  - All nullomers with order x, are also a nullomer for order x-1, thus, the output reflects this behaviour.

### nullomer_statistics_retrieval.py

Description:
Retrieve nullomer statistics from a trie bit file. Can compute composition, total nullomers, motif statistics (CPG, palindromy, homopolymers) and trivial nullomer occurence.

Usage:
python nullomer_statistics_retrieval.py [-h] [-o OUTPUT] [-n1 NULL1] [-n2 NULL2]
                                        [-s {all,composition,trivial,motifs}]

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
Searches for nullomers that are prime across organisms in a given group, obtaining the sequences and writing them in a txt file. Prime statistics are outputted in a centralized CSV.

Usage:
python prime_nullomer_finder.py [-h] [-c CONFIG] [-o OUTPUT] [-k K_RANGE K_RANGE]
                                [-g {none,phylum_id,order_id,family_id,genus_id}]

Arguments:
  -c, --config     (default: workflow/config/config.yaml) Relative path to the YAML config file containing the organism names and nullomer file paths.
  -o, --output     (required) Directory where output files will be saved.
  -k, --k_range    (default: 10 12) Range of k values to search prime nullomers.
  -g, --group_by   (default: none) Mode of grouping organisms while searching for primes.

Example:
  python prime_nullomer_finder.py -c workflow/config/config.yaml -o results/primes/

Output:
  - prime_nullomers.txt: list of nullomer sequences that are prime across the specified organisms.
  - prime_statistics.csv: centralized CSV with statistics for each prime nullomer.

Notes:
  - The config file should follow the same structure as the Snakemake configuration.
  - The script expects nullomer bit files for each organism listed in the config.


### snakemake_config_generation.py

Description:
Script to generate a config.yaml file for the Snakemake pipeline. Includes genome integrity checking and configurable k-values.

Usage:
python snakemake_config_generation.py [-h] [--name NAME] [--out OUT] [--genomes GENOMES] [--log LOG]
                                      [--cscript CSCRIPT] [-ks KVALUES [KVALUES ...]] [--integrity]

Arguments:
  --name           (default: config.yaml) Name of the config file to generate.
  --genomes-dir    (default: workflow/data/genomes/) Directory containing the genome FASTA files.
  --kvalues        (required) List of k-values to be used in the pipeline.
  --output-dir     (default: workflow/results/) Directory for pipeline results.
  --integrity      (optional) performs integrity checking on the genomes of the config file.

Example:
  python snakemake_config_generation.py --name config_ecoli.yaml --genomes-dir data/genomes/ --kvalues 15 16 17 --output-dir results/ecoli/

Output:
YAML configuration file with the following structure:
  genomes_dir: "workflow/data/genomes/"
  output_dir: "workflow/results/"
  k:
  - 8
  - 9
  - 15
  organisms:
    - GCF_000005845.2
    - GCF_000001405.40
  paths:
    bench: benchmarks/
    checked: results/checked/
    final: runs/
    genomes: data/genomes/
    log: logs/
    results: results/
  statistics:
  - composition
  - counter
  - motifs

Notes:
  - The script performs integrity checking on the genome files before generating the config.
  - Remember that the file snakemake will use shall be named config.yaml.

### translate_to_bit.py

Description:
Script for bit file encoding from nullomeric sequences.

Usage:
python translate_to_bit.py [-h] [-k [KVALUES ...]] [-o OUTPUT] [-i INPUT] [-f FILENAME] [--multiple]

Arguments:
  --kvalues         (default): Values of k from the source files. Default = 10.
  --output          (default): Output directory for the bit files. Default = workflow/results
  --input           (default): Input directory containing the files to be converted. Default = workflow/results
  --filename        (optional): Full path of the nullomer file to be encoded. Used only for single file conversion.
  --multiple        (flag): Allows the processing of multiple input files. All files inside input are read and outputted in bin format.

Example:
  For converting a single file with k = 10:
    python translate_to_bit.py -k 10 -f path/to/sequence_file.txt
  For converting all files inside a folder. For this to work as intended **the files must be grouped by k value, inside folders named kkvalue**:
    python translate_to_bit.py -k 10 11 -i to_convert -o workflow/results/converted --multiple
    In this case, the file structure in the input folder is:
```
    to_convert
            ├── k10
            │   ├── sequence_a.txt
            │   └── sequence_b.txt
            └── k11
                └── sequence_c.txt
```
  Where files a and b are composed by length 10 words and sequence c by length 11 words.

Output:
BINARY file: custom compact bit file containing the Trie structure relative to the input nullomeric sequences.

