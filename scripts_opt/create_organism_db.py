"""Creates a JSON file with metadata for given organism list."""
import yaml
import json
import argparse

from Bio import Entrez
from nulloretriever.utils.validation import get_valid_email, get_valid_tool
from nulloretriever.data.ncbiapidata import get_genome_metadata, get_genome_length, get_taxonomy_metadata


def setup_argparser() -> argparse.ArgumentParser:
    """Argument parsing function for the JSON creation script."""
    parser = argparse.ArgumentParser(description="""
                                     Script aimed at the
                                     creation of a central JSON file containing
                                     analyzed organisms metadata.""")
    parser.add_argument(
        '-c',
        '--config',
        help="Relative path to the YAML config file containing the organisms"
        " names. Default = workflow/config/config.yaml",
        default="workflow/config/config.yaml"
    )
    parser.add_argument(
        '-o',
        '--output',
        help="Relative path to the JSON output file."
        "Default = workflow/data/organisms.json",
        default="workflow/data/organisms.json"
    )
    return parser


def format_accesion_code(organisms):
    """Adjust a series of accession codes to their original form.

    Args:
        organisms(list): list of accession codes with . replaced by _
    Returns:
        organisms(list): acession codes with the last _ changed back
                             to .
    """
    for index, organism in enumerate(organisms):
        new_org = organism[:-2] + organism[-2:].replace('_', '.')
        organisms[index] = new_org
    return organisms


def main():
    """Generate all organisms central JSON file."""
    parser = setup_argparser()
    args = parser.parse_args()
    config_file = args.config
    Entrez.email = get_valid_email()
    Entrez.tool = get_valid_tool()
    out_path = args.output
    with open(config_file, 'r') as f:
        content = yaml.safe_load(f)
        organisms = content['organisms']
    format_accesion_code(organisms)
    metadata = get_genome_metadata(organisms)
    print(metadata)
    new_metadata = get_genome_length(metadata)
    get_taxonomy_metadata(metadata)
    print(metadata)
    with open(out_path, 'w') as f:
        json.dump(new_metadata, f)


if __name__ == "__main__":
    main()
