# Optional Scripts

## Overview

This directory contains supplementary scripts that are not part of the core `nulloretriever` Snakemake workflow. These scripts are intended for ancillary tasks such as data preparation, post-processing, results summarization, or other manual analyses.

## Usage

Each script is designed to be executed as a standalone program from the command line. For detailed usage instructions, each script should implement a command-line help message, accessible via the `-h` or `--help` flag.

Example:
```sh
python opt-scripts/summary_generator.py --input results/file.txt --output summary.csv
```

To view all available parameters for a script:
```sh
python opt-scripts/summary_generator.py --help
```

## Testing

All scripts in this directory should have corresponding tests located in the `tests/opt-scripts/` directory. These tests verify the correctness of the script's logic and command-line interface.

Refer to the `README.md` in the `tests/` directory for instructions on how to run the tests.
