from nulloretriever.utils.snakemake_path_tools import *
import sys
import argparse

def setup_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description='Script for core path existence checking.')

    parser.add_argument(
        '--config',
        nargs='?',
        help='Path to the configuration file to check',
        default = '../../config/config.yaml'
    )
    parser.add_argument(
        '--output',
        nargs='?',
        help='Path to the output of not found paths',
        default = 'not_found.txt'
    )
    parser.add_argument(
        '--checklist',
        nargs='?',
        help='Paths that need to be checked, if base snakemake or scripts',
        default = 'snakemake'
    )
    return parser
def main():
    parser = setup_argparser()
    args = parser.parse_args()
    config_file = args.config
    not_found_out = args.output
    list_for_checking = args.checklist
    if list_for_checking == 'snakemake':
        yaml_variables = ['paths']
    else:
        yaml_variables = ['py_scripts','c_scripts']
    missing_paths = check_config_paths_existence(config_file, yaml_variables)
    if missing_paths:
        with open(not_found_out, 'w') as f:
            f.write('\n'.join(missing_paths))

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)