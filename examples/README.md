# Minimal Working Example

This directory contains a minimal working example of the Nulloretriever pipeline
using the Bacillus subtilis ASM904v1 reference genome, under GCF_000009045.1.

**Citation:**
Kunst, F., Ogasawara, N., Moszer, I., et al. (1997). The complete genome sequence of the gram-positive bacterium Bacillus subtilis. *Nature*, 390(6657), 249–256.

## Structure

examples/
├── README.md
├── data
│   ├── GCF_000009045_1        # B. subtilis genome in FASTA format
│   ├── config.yaml            # Configuration file for the current example
├── results                    # Contains the results after running the example
└── workflow                   # Place for .snakemake, that stores example info

## Requirements

- Conda environment with all dependencies installed (see root environment.yaml)
- Snakemake installed (included in the Conda environment)
- C compiler (gcc) for building the C k-mer extractor (done automatically)

## Running the example

1. Make sure the conda environment is active:

```
conda activate nulloretriever
```

2. From the repository root, run:

```
snakemake -s workflow/Snakefile --configfile examples/data/config.yaml --cores 2
```

This will:
- Compile the C binary (if not already built)
- Process the B. subtilis for nullomers on k=12
- Generate output files in examples/results/
## Output

After successful execution, you will find:

- examples/results/k12/GCF_000009045_1/null_bit_format – nullomer data in binary format
- examples/results/nullomer_statistics.csv – retrieved statistics from the null file
- examples/results/nullomer_statistics_summarized.csv - sums different organisms statistics
- examples/results/checked/ - stores non conformant genomes files
- examples/benchmarks/k12/ - stores time and memory benchmarks for the snakemake run

## Customizing

To change k-values or use a different genome, edit config.yaml. Base directories and fasta extension can also be changed. For further information, refer to the README about config files.

## Notes

- The C binary is compiled automatically by Snakemake.
- If you are not at the repository root, path resolution shall fail.
- For more details, refer to the main README in the repository root.
