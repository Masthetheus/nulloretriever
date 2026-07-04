# Snakemake Pipeline

## Overview

This directory contains all the necessary components to execute the main `nulloretriever` analysis pipeline using Snakemake.

## Directory Structure

-   **`Snakefile`**: The main pipeline definition file that includes the rules from the `rules` directory.
-   **`rules/`**: Contains individual rule files (`.smk`), each defining a specific step in the workflow.
-   **`envs/`**: Contains Conda environment definition files (`.yaml`) for ensuring the reproducibility of each rule.
-   **`scripts/`**: Contains scripts that are executed by the rules in the pipeline. These are distinct from the `opt-scripts` directory.

## How to Run the Pipeline

### Prerequisites

1.  **Install Snakemake and Conda**: Ensure that both are installed and available in your environment.
2.  **Configure the Pipeline**: Copy the template configuration file (`config/config.yaml.template`) to `config/config.yaml` and modify it to specify the input files and parameters for your analysis.

### Execution

To perform a dry run and view the jobs that will be executed:
```sh
snakemake --use-conda -n
```

To execute the full pipeline:
```sh
snakemake --use-conda --cores <number_of_cores>
```
