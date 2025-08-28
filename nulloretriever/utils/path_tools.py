"""Functions related to path management and manipulation
UNDER DEVELOPMENT:
    testing needed
    compatibility needed
    integration with main pipeline needed
    retrieve_nullomer_files_path: obtain org list needs implementation

"""
import os
import yaml

def retrieve_analyzed_k_values(base_path):
    """Obtains k values analyzed via snakemake pipeline
    Args:
        base_path(str): path that contains all snakemake generated files
    Returns:
        ks(arr): array of analyzed k values
    """
    ks = []
    for name in os.listdir(base_path):
        path = os.path.join(base_path, name)
        print(f"Checking path: {path}")
        if os.path.isdir(path) and name.isdigit():
            ks.append(int(name))
    return sorted(ks)

def retrieve_analyzed_organisms(config_file):
    """Obtains analyzed organisms names as an array
    Args:
        config_file(str): path to the yaml config file used for the nullomer retrieval runs
    Returns:
        organisms(array): array containing all organisms name inside the config.yaml file passed as arg
    """
    organisms = []
    if not os.path.exists(config_file):
        raise FileNotFoundError(f"Config file not found: {config_path}")
    try:
        with open(config_file, 'r') as f:
            config_data = yaml.safe_load(f)
            organisms = config_data.get('organisms',[])
    except Exception as e:
        print("Please check your config file integrity!")
    return organisms

def retrieve_nullomer_files_path(base_path, config_file):
    """Retrieves an array of all nullomer output files path
    Args:
        base_path(str): path that contains all snakemake generated files
        config_file(str): config.yaml file used for the snakemake run
    Returns:
        nullomers_paths(dict): all found paths for given organisms inside base_path paired with k value of the run
        found_errors(arr): all paths to nullomer output files supposed to exist that couldn't be found
    """
    if not os.path.exists(config_file):
        raise FileNotFoundError(f"Config file not found: {config_path}")
    organisms = retrieve_analyzed_organisms(config_file)
    found_errors = []
    k_values = retrieve_analyzed_k_values(base_path)
    nullomers_paths = {}
    try:
        with open(config_file, 'r') as f:
            config_data = yaml.safe_load(f)
            genome_direct = config_data['paths']['genomes']
    except Exception as e:
        print("Please check your config file integrity!")
    for k in k_values:
        l = int(k // 2)
        k_path = base_path + str(k) + '/'
        print(f"Checking directory for k={k}")
        cont = 0
        paths_k = []
        for org in orgs:
            org_path = k_path + org + '/'
            genome_path = genome_direct + org
            nullomer_out_file = org_path + f'nullomers_{org}_{k}'
            if not os.path.exists(nullomer_out_file):
                found_errors.append([org, k, 'File not found!'])
                continue
            if not os.path.exists(genome_path):
                found_errors.append([org, k, 'Genome not found!'])
                continue
            paths_k.append((genome_path, nullomer_out_file))
        if paths_k:
            nullomers_paths[k] = paths_k
    return nullomers_paths, found_errors

