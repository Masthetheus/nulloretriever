# Nullomer extraction pipeline 

The current pipeline aims to extract nullomer data from multiple organisms in an compact data format(triebit and txt files), leadin to low memory usage and time optimization.

## Summary

- [Project Organization](#project-organization)
- [Installation](#installation)
- [Usage](#usage)
- [Workflow](#workflow)
- [Configuration](#configuration)
- [Tests](#tests)
- [License](#license)

## Project Organization

The project is organized to separate raw data, scripts, results, configuration, and documentation. Below is the folder structure and a brief description of each component:

```
nullomer_c_trie_bit/
├── config/           # Configuration files (e.g., config.yaml for pipeline parameters)
├── data/
│   ├── genomas/      # Reference genomes (FASTA files)
│   ├── raw/          # Raw input data (original downloads, unprocessed)
│   └── processed/    # Intermediate processed data (e.g., k-mers, tries)
├── results/
│   ├── k/            # Results organized by k-mer size (e.g., results/8/, results/12/)
│   └── summary/      # Aggregated results, tables, and figures
├── scripts/
│   ├── c/            # C source code for k-mer extraction and processing
│   ├── python/       # Python scripts for trie construction, analysis, and post-processing
├── workflow/
│   ├── Snakefile     # Main Snakemake workflow file
│   └── rules/        # (Optional) Modular Snakemake rules
├── logs/             # Execution and error logs
├── benchmarks/       # Benchmark files for resource usage tracking
├── tests/            # Unit and integration tests
├── README.md         # Project documentation
├── .gitignore        # Files and folders to ignore in version control
├── LICENSE           # Project license
├── requirements.txt  # Python dependencies
└── environment.yml   # Conda environment specification (optional)
```

**Key points:**
- **config/**: Centralizes all configuration files for easy parameter management and snakemake parameter calling.
- **data/**: Stores all input data, both raw and processed, but large files should not be versioned.
- **results/**: Contains all output files, organized for easy retrieval and analysis.
- **scripts/c/**: Contains all C code for efficient k-mer extraction.
- **scripts/python/**: Contains all Python code for trie operations and downstream analysis.
- **workflow/**: Houses the Snakemake workflow, enabling reproducible and automated analyses.
- **logs/** and **benchmarks/**: Facilitate debugging and performance monitoring.
- **tests/**: Ensures code correctness and reproducibility.
- **README.md, LICENSE, requirements.txt, environment.yml**: Provide documentation, licensing, and environment setup for users and collaborators.

## Installation

1. Automated install:

Aiming to mantain workflow integrity it's recommended to install the current pipeline via the setup script.

```sh
./setup.sh
```
If the above doesn't act accordingly, one may check if it's enabled to run. If not, then:
```sh
chmod +x setup.sh
```

2. Manual Installation
    1. Python dependencies
        - It is recommended the use of a conda environment for general package installing and maintenence. A environment.yml file can be found and use to create such.
            ```sh
            conda env create -f environment.yml
            conda activate nullomer-env 
            ```
        - The libraries can be installed also via pip, taking as a base the requirements.txt file.
            ```sh
            pip install -r requirements.txt
            ```
    2. C dependencies:
        - One specific script is used in the snakemake pipeline, thus needing to be compiled accordingly. First, ensures an C compiler is present. For Arch:
            ```sh
            sudo pacman -S base-devel
            ```
        - The path pointed in the compile command must be the same as below, if changed, further pipeline workarounds will be necessary.
            ```sh
            gcc -O2 -o fasta_kmers_novo fasta_kmers_novo.c
            ```

## Usage

1. Set the terminal at the main directory:
    ```sh
    cd ~/nullomer_c_trie_bit/
    ``` 
2. Activate the conda environment
    ```sh
    conda activate nullomer-env
    ```
3. Configuring the config.yaml file.
    1. If you know what you are doing:
        - Open with your editor of preference the config.yaml present in the config/ directory.
        - Insert the organisms to be checked, one per line, preceded by - and one space.
        - Make sure that organisms name are the same as that of the fasta genome files.
        - Inform k values in the same format as organisms.
        - If at least one of the k values is odd, change the odd flag to True.
    2. If you don't know or doesn't want to know what you are doing:
        - Run the config creation python script:
            ```sh
            python3 scripts/python/setup_config.py
            ```
        - Follow the instructions as needed.
        - It is strongly advised to load all genomes that will be part of the analysis on data/genomes and then using the 3rd script option.
        - You can check if the file was created by the command below:
            ```sh
            ls config/
            ```
4. Running the pipeline! (finally)
Run the following on the terminal (Remember to stay at the main directory!!):
```sh
snakemake --cores x
```
Where x stands to the max number of cores that shall be used during the process.
 
## Workflow

## Configuration

## Tests

## License
