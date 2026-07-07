# Nulloretriever

A scalable Snakemake workflow for discovery of nullomers and Minimal Absent Words (MAWs) in complete genomes.

Nulloretriever integrates C, Python, and Snakemake to efficiently process genomes with k > 14, using significantly less memory and time than existing tools.

## Features

- Efficient k-mer extraction with a memory-optimized C core.
- Scalable: handles large genomes (up to 15 Gb) with k >= 15.
- Reproducible: fully containerised workflow with Snakemake and Conda.
- Optional utilities: download genomes via NCBI API, generate configuration files, and retrieve statistics.

## Repository Structure
```
nulloretriever/
├── config/
│   └── config.yaml           # Main configuration file
├── examples/
│   ├── data/                 # Example genome (B. subtilis)
│   └── config.yaml           # Configuration for the example
├── experimental/             # Scripts under development (not stable)
├── nulloretriever/           # Python package (core modules)
│   ├── analysis/
│   ├── core/
│   ├── data/
│   └── utils/
├── scripts_opt/              # Optional utilities (download, etc.) and standalone scripts
├── workflow/
│   ├── scripts/
│   │   ├── c/                # C source code and Makefile
│   │   └── python/           # Python scripts called by Snakemake
│   └── Snakefile             # Main workflow definition
├── .gitignore
├── LICENSE
├── README.md
├── environment.yaml          # Conda environment
├── pyproject.toml            # Python package metadata
└── requirements.txt          # Pip dependencies
```
## Installation

Via Conda or Mamba (recommended):

Clone the repository and create the environment:
```
git clone https://github.com/Masthetheus/nulloretriever.git
cd nulloretriever
conda env create -f environment.yaml
conda activate nulloretriever
pip install -e .
```
This installs all dependencies, including Snakemake, Biopython, and the required Python packages.

For faster installation, use Mamba instead of Conda:
```
mamba env create -f environment.yaml
conda activate nulloretriever
pip install -e .
```
Mamba is a reimplementation of the Conda solver in C++ and resolves dependencies much faster. If you do not have Mamba installed, you can install it with:
```
conda install mamba -n base -c conda-forge
```
Manual installation (pip):

If you prefer bare pip, install the required packages:
```
pip install -r requirements.txt
pip install -e .
```
## Quick Start

A minimal working example using the Bacillus subtilis genome (k=15) is provided in examples/. From the repository root:
```
snakemake -s workflow/Snakefile --configfile examples/config.yaml --cores 2
```
This will:
1. Compile the C binary (if not already built).
2. Extract nullomers for k=15 from the B. subtilis genome.
3. Save results on the results path of the config.yaml. Default = workflow/results/.

For more details, see [examples](examples/README.md).

## Configuration

All pipeline parameters are defined in a YAML configuration file. A template is provided in config/config.yaml:

k: [8, 9, 10, 11, 12, 13, 14, 15, 16]
organisms:
  - GCF_000001405_40
  - GCF_000181295_1
paths:
  bench: workflow/benchmarks
  checked: workflow/results/checked
  final: workflow/runs
  genomes: workflow/data/genomes
  log: workflow/logs
  results: workflow/results
  c: workflow/scripts/c
statistics:
  - composition
  - trivial
  - motifs
fasta_extension: ""
*
**Key parameters:**

- k: List of k-mer sizes to analyze.
- organisms: List of genome identifiers (base names without extension).
- paths: Directory structure for inputs, outputs, logs, and benchmarks.
- statistics: Which statistics to compute (composition, counter, motifs).
- fasta_extension: File extension for genomes (empty if no extension).

## Workflow Rules

The Snakefile defines the following main rules:

- create_base_directories: Creates all required directories.
- generate_c_binaries: Compiles the C k-mer extractor using make.
- extract_nullomers: Runs the C binary on a genome for a given k.
- retrieve_nullomer_statistics: Computes statistics from the extracted nullomers.
- summarize_nullomer_statistics: Aggregates results from all organisms into a single CSV.

## Adding New Genomes

To add a new genome to the pipeline:

1. Place the FASTA file in the directory specified by genomes_dir (default: workflow/data/genomes/).
2. Add the base name (without extension) to the organisms list in config.yaml.
3. Run the workflow.

The pipeline will automatically process all listed organisms for all specified k-values.

## C Compilation

The C binary (kmer_extractor) is compiled automatically by Snakemake. If you need to compile manually:
```
cd workflow/scripts/c/
make
```
The binary will be placed in workflow/scripts/c/bin/kmer_extractor.
The current installation intermediate files can be removed, in case of erros, with:
```
make clean
```
Followed by make to try and solve any found errors.
For further information about the extractor, refer to it's [README.md](workflow/scripts/c/README.md)

## Python scripts

Python scripts in the current project revolves around three maing groups: optional, manual and workflow adherent.

### Optional scripts

Scripts that are entirely optional, aiming to help with file generation and organisms processing.

- **create_organism_db.py**: Creates a JSON file containing metadata for the organisms inside given config.yaml.
- **genomes_utilities.py**: Download and capitalizes a serie of NCBI accession codes. Can also be used only for capitalization, since it's required by the workflow.
- **snakemake_config_generation.py**: Generates a custom config file for usage with the snakemake pipeline.

For further information, refer to it's specific [README](scripts_opt/README.md)

## Manual scripts

Here we have adapted workflow scripts, intended to manual modular execution of some pipeline process, for any given reason. They are usually considered alongisde optional scripts, since aren't necessary for pipeline execution, not adhering to snakemake conventions.

- **nulomer_extraction.py**: Derived from the nullomer extraction rule. Obtain the genome nullomers in bit, sequence or compact txt format for any organism and k value.
- **nullomer_statistics_retrieval.py**: Derived from the nullomer statistics rule. Given a nullomer bit file, retrieve statistics such as GC composition and palindromes occurence.
- **prime_nullomer_finder.py**: Planned to be added to the snakemake workflow in the future. Given a set of nullomer bit files, searches for prime nullomeric sequences and outputs them.

For further information, refer to the optional scripts [README](scripts_opt/README.md)

## Dependencies

All dependencies are managed via Conda (environment.yaml) or pip (requirements.txt). The main packages are:

- Python >= 3.10, < 3.13
- Snakemake >= 8.0
- Biopython
- bitarray
- numpy, pandas 
- pyyaml, requests, psutil

## How to Cite

If you use Nulloretriever in your research, please cite:


## License

This project is licensed under the GNU General Public License v3.0 - see the LICENSE file for details.

## Acknowledgments

This workflow was developed as part of a research project at UFRGS. We thank the bioinformatics community for providing open-source tools and datasets.

Last updated: 2026-07-07
