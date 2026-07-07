# Workflow Directory

This directory contains the Snakemake workflow definition for the Nulloretriever pipeline. It orchestrates the entire process from genome input to nullomer and MAW discovery.

## Structure
```
workflow/
├── Snakefile                 # Main workflow definition
├── config/
│   └── config.yaml           # Configuration file (paths, k-values, organisms)
├── scripts/
│   ├── c/                    # C source code and binaries
│   │   ├── src/              # C source files (.c, .h)
│   │   ├── include/          # Header files
│   │   ├── build/            # Object files (.o) generated during compilation
│   │   ├── bin/              # Compiled binary (kmer_extractor)
│   │   └── Makefile          # Build script for the C binary
│   └── python/               # Python scripts called by Snakemake rules
│       ├── genome_checking.py
│       ├── nullomer_extraction.py
│       ├── nullomer_statistics_retrieval.py
│       ├── nullomer_statistics_summary.py
│       └── path_checking.py
└── (no data/ or results/ directories – these are defined in config.yaml)
```

## Running the Workflow

From the repository root, using the default config.yaml on "workflow/config/config.yaml":
```
snakemake -s workflow/Snakefile --cores N
```
Or using a custom configuration file:
```
snakemake -s workflow/Snakefile --configfile path/to/my_config.yaml --cores N
```
With the minimal example:
```
snakemake -s workflow/Snakefile --configfile examples/data/config.yaml --cores 2
```
## Configuration

The workflow uses a YAML configuration file for further customization, with the following structure (see workflow/config/config.yaml):
```
k: [8,9,15,16]
organisms:
- GCF_000009045_1
- GCF_000005845_2
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
- counter
- motifs
fasta_extension: ""

```
All paths are relative to the repository root. The organisms list should match the base names of your FASTA files (without extension).

Rules Overview

The main rules in the Snakefile are:

- generate_c_binaries: Compiles the C k-mer extractor from source.
- extract_nullomers: Extracts nullomers for a given genome and k-value.
- retrieve_statistics: Aggregates statistics from the extracted nullomers.
- summarize: Combines results into a final summary table.

Each rule uses wildcards for organism and k, making the workflow scalable to multiple genomes and parameter combinations.

Adding New Genomes

To add a new genome to the pipeline:

1. Place the FASTA file (with the correct extension) in the directory specified by genomes_dir.
2. Add the organism identifier (base name without extension) to the organisms list in config.yaml.
3. Run the workflow.

The workflow will automatically process all listed organisms for all specified k-values.

C Binary

The C binary (kmer_extractor) is compiled automatically by Snakemake from the source code in scripts/c/src/. The Makefile handles the build process.

If you need to compile manually:

cd workflow/scripts/c/
make

The binary is placed in scripts/c/bin/kmer_extractor.

Troubleshooting

- Missing input files: Ensure genomes_dir points to the correct directory and that all organism FASTA files exist with the expected extension.
- Compilation errors: Check that gcc is installed and that all source files in scripts/c/src/ are present.
- Snakemake errors: Use snakemake -n to dry-run and check for syntax or missing files.

Dependencies

- Snakemake (>= 7.0)
- Python >= 3.8
- C compiler (gcc)
- Biopython, pandas, pyyaml (installed via environment.yaml)

All dependencies are listed in the root environment.yaml and requirements.txt.

License

This workflow is part of the Nulloretriever project and is licensed under GNU GPL v3.0 (see LICENSE in the repository root).

Last updated: 2026-07-06
