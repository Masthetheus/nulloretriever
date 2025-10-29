# Examples Suite

## Overview

This directory contains base examples for the **NulloRetriever** project. Each is designed to supply a minimal reproducible example of the individual modules, optional scripts and the complete snakemake pipeline.

## Directory Structure

Each example folder is organized to mirror it's actual project usage and implementation. For any module or script, there is a corresponding example file prefixed with `example_`.

-   `examples/modules/example_module_name.py` will supply an example to `modules/module_name.py`.
-   `examples/opt-scripts/test_script_name.py` will supply an example to `opt-scripts/script_name.py`.

## How to Run a Minimal Reproducible Example

### Prerequisites

Ensure that all project dependencies are correctly installed in your system, ideally inside a virtual environment or a conda environment.

**Observation: Given the extension of the project, it is highly advised to fully configure the environment, with all needed dependencies, not only those related to the current example.** 

Please check the [Installation guide](../README.md#installation) for further instructions.

### Execution

To run the entire test suite, execute the following command from the root directory of the repository:

```sh
pytest
```

To run tests for a specific file or directory, provide its path:

```sh
pytest tests/modules/test_module_name.py
```
