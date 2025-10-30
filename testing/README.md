# Examples Suite

## Overview

This directory contains base examples for the **NulloRetriever** project. Each
is designed to supply a minimal reproducible example of the individual modules,
optional scripts and the complete snakemake pipeline.

## Directory Structure

Each example folder is organized to mirror it's actual project usage and
implementation. For any module or script, there is a corresponding example file
prefixed with `example_`.

-   `examples/modules/example_module_name.py` will supply an example to `modules/module_name.py`.
-   `examples/opt-scripts/test_script_name.py` will supply an example to `opt-scripts/script_name.py`.

## Types of Minimal Reproducible Example

In the current structure, we group the examples in two groups: single and suite.
Singles consist on the isolated file of certain module, optional script or
pipeline modular script. In this sense, the output refer directly only to it's
scope. In the other hand, suites aim to integrate the maximum scripts possible,
in their folder ambient, and thus supply a streamlined example of execution.

**Observation: You can differentiate both groups by their file names, since
suites are always preceded by the term itself.**

## How to Run a Minimal Reproducible Example

### Prerequisites

Ensure that all project dependencies are correctly installed in your system,
ideally inside a virtual environment or a conda environment.

**Observation: Given the extension of the project, it is highly advised to**
**fully configure the environment, with all needed dependencies, not only** 
**those related to the current example. Please check the**
**[Installation guide](../README.md#installation) for further instructions.**

- Make sure the NulloRetriever environment is running, or activate it by:

    ```sh
    conda activate nulloretriever
    ```

- Change to the desired examples folder or export it to your PATH before executing the scripts:

    ```sh
    cd desired-folder

    # or

    export PATH=$PATH:/path/of/desired/folder
    ```

**Observation: Detailed information on the execution proccess and available examples can be found in each folder's relative README file.**

### Execution

The majority of examples cases can be executed in two ways:

1. Directly running the script file or calling it with flags in the terminal, for example:

    ```sh
    python example_script.py positional_arg --flag1 --arg1 arg
    ```
    
2. By running the Full Example Suite, if available:

   ```sh
   python full_example_suite --flag1 --arg1 arg
   ```

## Contributing

If some example of your interest don't has an available script, or you would
like a certain suite structure, please submit a issue.
