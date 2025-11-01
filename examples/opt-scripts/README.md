# Optional Scripts

## Overview

This directory contains tests for the optional scripts that are not part of the core `nulloretriever` Snakemake workflow. These scripts are intended for auxiliary tasks, mainly data preparation.

## General Usage

Each script is designed to be executed as a standalone program from the command line. All needed data can be found inside ../examples/data.

Example:
```sh
cd opt-scripts/
python example_create_organism_db.py
```

Some scripts have parameters that can be set. They can be found detailed in **PLACEHOLDER**. To view all available parameters for a script:

```sh
python example_snakemake_config_generation.py --help
```

## Available scripts

In this section you can find each added example script general description and quirks, if any.

### Create Organism DB

From a config.yaml file (default at '../data/config.yaml') retrieves all organisms metadata and outputs them in a JSON file. Needs an e-mail and Entrez tool code. This db is needed for further data manipulation, when used along the snakemake full pipeline execution.

### Snakemake Config Generation

Example of a template config.yaml file generation to be used along the snakemake pipeline or standalone scripts. It searches '../data/' for any existing config.yaml files, keeping it's listed organisms, then generate the config.yaml based on an internal template.

Given the multiple possible options and flags available **PLACEHOLDER**
