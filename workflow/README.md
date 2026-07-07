# Workflow Directory

This directory contains the Snakemake workflow definition for the Nulloretriever pipeline. It orchestrates the entire process from genome input to nullomer and MAW discovery.

## Structure

workflow/
├── Snakefile                 # Main workflow definition
├── config/
│   └── config.yaml           # Configuration file (paths, k-values, organisms)
├── scripts/
│   ├── c/                    # C source code and binaries
│   │   ├── src/              # Source files (.c, .h)
│   │   ├── include/          # Header files
│   │   ├── build/            # Object files (.o) generated during compilation
│   │   ├── bin/              # Compiled binary (kmer_extractor)
│   │   └── Makefile          # Build script for the C binary
│   └── python/               # Python scripts called by Snakemake rules
│       ├── nullomer_extraction.py
│       ├── nullomer_statistics_retrieval.py
│       └── nullomer_statistics_summary.py
└── (data/ and results/ directories are defined in config.yaml)

## Running the Workflow

From the repository root:
```
snakemake -s workflow/Snakefile --cores N
```
Or using a custom configuration file:
```
snakemake -s workflow/Snakefile --configfile path/to/my_config.yaml --cores N
```
Minimal example:

snakemake -s workflow/Snakefile --configfile examples/config.yaml --cores 2

## Configuration

The workflow uses a YAML configuration file with the following structure (see workflow/config/config.yaml):
```
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
```
All paths are relative to the repository root. The organisms list should match the base names of your FASTA files (without extension). The fasta_extension field can be set to ".fasta" or other variations if your files have an extension; if empty, the pipeline expects files without extension.

## Rules Overview

The main rules in the Snakefile are:

- create_base_directories: Creates all required directories defined in paths.
- generate_c_binaries: Compiles the C k-mer extractor using the Makefile.
- extract_nullomers: Runs the C binary via Python script on a genome for a given k.
- retrieve_nullomer_statistics: Computes statistics from the extracted nullomers.
- summarize_nullomer_statistics: Aggregates results from all nullomers from every organism and k value into a single CSV.

Each rule uses wildcards for organism and k, making the workflow scalable to multiple genomes and parameter combinations.

## Adding New Genomes

To add a new genome to the pipeline:

1. Place the FASTA file in the directory specified by genomes_dir (default: workflow/data/genomes/).
2. Add the organism identifier (base name without extension) to the organisms list in config.yaml.
3. Run the workflow.

The pipeline will automatically process all listed organisms for all specified k-values.

## C Binary Compilation

The C binary (kmer_extractor) is compiled automatically by Snakemake from the source code in workflow/scripts/c/src/. The Makefile handles the build process.

If you need to compile manually:
```
cd workflow/scripts/c/
make
```

If any error is found, first try to clean and re-make the file:

```
make clean
make
```

The binary will be placed in workflow/scripts/c/bin/kmer_extractor.

## Makefile Notes

The Makefile includes -fPIE and -pie flags to ensure compatibility with modern systems that require Position Independent Executables. If you encounter relocation errors, ensure that the bin/ directory exists (the Makefile creates it automatically) and that all dependencies are installed.

## Troubleshooting

- Missing input files: Ensure genomes_dir points to the correct directory and that all organism FASTA files exist with the expected extension.
- Compilation errors: Check that gcc is installed and that all source files in scripts/c/src/ are present.
- Snakemake errors: Use snakemake -n to dry-run and check for syntax or missing files.

## Dependencies

- Snakemake (>= 8.0)
- Python >= 3.10, < 3.13
- C compiler (gcc)
- Biopython, pandas, pyyaml, bitarray, etc. (installed via environment.yaml)

All dependencies are listed in the root environment.yaml and requirements.txt.

## License

This workflow is part of the Nulloretriever project and is licensed under GNU GPL v3.0 (see LICENSE in the repository root).

Last updated: 2026-07-07
