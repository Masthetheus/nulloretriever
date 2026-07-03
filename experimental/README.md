# Experimental Scripts

This directory contains scripts that are **under development** or intended for **future analyses**. They are not part of the main Snakemake workflow and are not guaranteed to be stable, well-documented, or thoroughly tested.

These scripts are shared here primarily for the author's convenience across different machines and for early collaboration. They may change, break, or be removed without notice.

---

## Scripts overview

- `data_analysis_hub.py` – Generates graphs and CSVs from nullomer data analysis. Supports grouping, parameter selection, and export of processed data. Intended for exploratory visualization.
- `aa_placeholder.py` - Generates aa sequences from prime nullomer file. Prints aa sequences bigger or equal than four. Intended for further aminoacid and prime relationship analysis.
- `degenerate_analysis.py` - From the list of found aminoacids, prints the number of possible sequences that can generate that combination. Intended for further analysis of prime aa behaviour and origin.
- `null_plots.py` - Experimental plots for nullomer data analysis.
- `order_null_plots.py` - Experimental plots for nullomer data analysis grouped by order.

---

## Dependencies

These scripts rely on the same Python environment as the main pipeline. Install dependencies from the root `requirements.txt` or `environment.yaml`:

conda env create -f environment.yaml
conda activate nulloretriever-env

Some scripts may require additional packages not listed in the main environment. They may be installed inside the environment via pip. 

---

## Usage

Each script is designed to be run independently from the command line. Refer to the script's docstring or use the `--help` flag for details on arguments.

Example:

python data_analysis_hub.py --input results/ --output figures/

---

## Status and stability

- `data_analysis_hub.py` – Experimental. Not intended initially to be integrated with workflow.
- `prime_nullomer_finder.py` – Prototype. Requires manual config; use with caution.
- `create_organism_db.py` – Stable. Works, but may be superseded by a future utility.

---

## Contributing

If you are a collaborator and wish to improve or finalize any of these scripts, please move them to `scripts_opt/` once they are stable and document them accordingly. For now, treat this folder as a sandbox.

---

*Last updated: 2026-07-03*
