# Core Modules

## Overview

This directory contains the core Python modules that form the functional basis of the `nulloretriever` pipeline. These modules are designed to be reusable, modular, and are imported by various scripts and pipeline rules.

## Functionality

Each module encapsulates a specific part of the workflow's logic, such as data processing, file I/O, or algorithmic computations. The code is intended to be documented following standard docstring conventions.

## Usage

These modules are not meant to be executed directly. They should be imported into other scripts or Snakemake rules where their functionality is required.

Example:
```python
from . import utility_module

# Use functions from the module
utility_module.process_data(...)
```

## Testing

All modules in this directory must have corresponding unit tests located in the `tests/modules/` directory. Before committing any changes, ensure that all associated tests pass successfully.

Refer to the `README.md` in the `tests/` directory for instructions on how to run the tests.
