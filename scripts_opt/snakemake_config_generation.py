"""Script to generate a config.yaml file for the snakemake pipeline"""

import argparse
from pathlib import Path
import yaml
from nulloretriever.utils.paths import gather_files_names, gather_files_paths
from nulloretriever.utils.integrity import check_multiple_genomes_integrity


def setup_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Script for genome integrity checking."
    )
    parser.add_argument(
        "--name",
        help="Name to the config.yaml file. Useful for generating and keeping multiple config files. Remember to rename the one that shall be used to config.yaml before using the pipeline!",
        default="config.yaml",
    )
    parser.add_argument(
        "--out",
        help="Path to output the config.yaml file. Change only if you know what you are doing!",
        default="workflow/data/",
    )
    parser.add_argument(
        "--genomes",
        help="Path to the directory with the genomes files."
        "Default = workflow/data/genomes",
        default="workflow/data/genomes",
    )
    parser.add_argument(
        "--log", help="Path to log directory", default="configlog"
    )
    parser.add_argument(
        "-ks",
        "--kvalues",
        nargs="+",
        type=int,
        help="K values to be analyzed. Default k = 10",
        default=[10],
        required=False,
    )
    parser.add_argument(
        "--integrity",
        action="store_true",
        help="Checks the downloaded genomes integrity and removes in case of non-expected formatting",
    )
    parser.add_argument(
        "--statistics",
        nargs="+",
        type=str,
        default=["composition", "trivial", "motifs"],
        help="""Specifies which statistics must be retrieved from the nullomer files.
        Default = composition, trivial and motifs.""",
    )

    parser.add_argument(
        "--fasta_extension",
        type=str,
        default="" ,
        help="""Specifies which fasta extension is to be expected on genomes files.
        Default ="".""",
    )

    return parser


def main():
    """Creates a config file for the snakemake pipeline."""

    parser = setup_argparser()
    args = parser.parse_args()
    out_path = args.out + args.name
    genomes_path = args.genomes
    kvalues = args.kvalues
    log = args.log
    statistics = args.statistics
    fasta_extension = f".{args.fasta_extension}"

    yaml_dump = {}
    organisms_full = gather_files_paths(genomes_path)
    organisms = [o.removesuffix(fasta_extension) for o in organisms_full]
    already_exists = Path(out_path).exists()
    if already_exists:
        with open(out_path, "r") as f:
            pre_existing_data = yaml.safe_load(f)
            try:
                yaml_dump["organisms"] = pre_existing_data["organisms"]
                org_set = set(yaml_dump["organisms"])
                for org_list in organisms:
                    org_set.update(org_list)
                yaml_dump["organisms"] = list(org_set)
            except Exception as e:
                yaml_dump["organisms"] = ""
                print(f"The following exception was encountered: {e}")
    elif not already_exists and not args.integrity:
        yaml_dump["organisms"] = organisms

    if args.integrity:
        app_organisms, napp_organisms = check_multiple_genomes_integrity(
            organisms, location
        )
        yaml_dump["organisms"] = app_organisms
        if napp_organisms:
            print(
                f"At least one organism didn't pass the integrity check, please verify the log at: {log} for further information."
            )
            try:
                with open(log, "w") as f:
                    f.write(
                        "The following organisms weren't approved on the integrity check:\n"
                    )
                    orgs = "\n".join(napp_organisms)
                    f.write(orgs)
            except:
                print("Error writing the not approved log file!")

    yaml_dump["k"] = kvalues
    yaml_dump["paths"] = {
        "genomes": "workflow/data/genomes",
        "results": "workflow/results",
        "final": "workflow/runs",
        "log": "workflow/logs",
        "bench": "workflow/benchmarks",
        "checked": "workflow/results/checked",
        "c": "workflow/scripts/c",
    }
    yaml_dump["statistics"] = statistics
    yaml_dump["fasta_extension"] = fasta_extension

    with open(out_path, "w") as f:
        yaml.dump(yaml_dump, f)


if __name__ == "__main__":
    main()
