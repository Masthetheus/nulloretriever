"""Creates a JSON file with metadata for given organism list."""
import yaml
import json
import argparse


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


def main():
    """Generate all organisms central JSON file."""
    parser = setup_argparser()
    args = parser.parse_args()
    config_file = args.config
    out_path = args.out



if __name__ == "__main__":
    main()
