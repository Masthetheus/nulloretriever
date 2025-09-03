"""Script to generate a config.yaml file for the snakemake pipeline"""
import sys
import argparse
from pathlib import Path
import yaml
from nulloretriever.utils.paths import gather_files_names, gather_files_paths
from nulloretriever.utils.integrity import check_multiple_genomes_integrity

def setup_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description='Script for genome integrity checking.')
    parser.add_argument(
        '--name',
        help='Name to the config.yaml file. Useful for generating and keeping multiple config files. Remember to rename the one that shall be used to config.yaml before using the pipeline!',
        default='config.yaml'
    )
    parser.add_argument(
        '--out',
        help='Path to output the config.yaml file. Change only if you know what you are doing!',
        default='workflow/config/'
    )
    parser.add_argument(
        '--genomes',
        help='Path to the directory with the genomes files',
        default='workflow/data/genomes/'
    )
    parser.add_argument(
        '--log',
        help='Path to log directory',
        default='configlog'
    )
    parser.add_argument(
        '--cscript',
        help='Path to c script directory',
        default='workflow/scripts/c/fasta_kmer_extraction.c'
    )
    parser.add_argument(
        '-ks',
        '--kvalues',
        nargs="+",
        type=int,
        help='K values to be analyzed. Default k = 10',
        default=[10],
        required=False
    )
    parser.add_argument(
        '--integrity',
        action='store_true',
        help='Checks the downloaded genomes integrity and removes in case of non-expected formatting'
    )
    return parser
def main():
    parser = setup_argparser()
    args = parser.parse_args()
    out_path = args.out + args.name
    genomes_path = args.genomes
    kvalues = args.kvalues
    log = args.log
    c_script = args.cscript
    directories_paths = ['genomes','results','final','log','bench','checked']
    yaml_dump={}
    organisms = gather_files_paths(genomes_path)
    already_exists = Path(out_path).exists()
    if already_exists:
        with open(out_path, 'r') as f:
            pre_existing_data = yaml.safe_load(f)
            yaml_dump['organisms'] = pre_existing_data['organisms']
    elif not already_exists and not args.integrity:
        yaml_dump['organisms'] = gather_files_names(genomes_path)
    if args.integrity:
        app_organisms, napp_organisms = check_multiple_genomes_integrity(organisms)
        yaml_dump['organisms'] = app_organisms
        if napp_organisms:
            print(f"At least one organism didn't pass the integrity check, please verify the log at: {log} for further information.")
            try:
                with open(log,'w') as f:
                    f.write("The following organisms weren't approved on the integrity check:\n")
                    orgs = '\n'.join(napp_organisms)
                    f.write(orgs)
            except:
                print("Error writing the not approved log file!")
    yaml_dump['k'] = kvalues
    yaml_dump['paths']={
        'genomes':"data/genomes",
        'results':"results/",
        'final':"runs/",
        'log':"logs/",
        'bench':"benchmarks/",
        'checked':"results/checked/"
    }
    yaml_dump['genomes']= 'data/genomes/'
    yaml_dump['c_script']= c_script
    try:
        with open(out_path,'w') as f:
            yaml.dump(yaml_dump,f)
    except:
        print("Error in dumping the parameters in the config file!")
    
if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)