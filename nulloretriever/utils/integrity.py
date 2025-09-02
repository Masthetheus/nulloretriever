"""Module that gather functions aimed at integrity checking and padronization"""
import gzip
from pathlib import Path
import os

def unzip_fasta_file(gzip_path):
    """Unzips files ending with .gz, specifically fasta files
    Args:
        gzip_path(str): path to gzipped file
    Returns:
        out_path(Path): path to uncompressed fasta file
    """
    gzip_path = Path(gzip_path)
    out_path = gzip_path.with_suffix('')
    with gzip.open(gzip_path, 'rb') as fin, open(out_path, 'wb') as fout:
        for line in fin:
            fout.write(line)
    os.remove(gzip_path)
    return out_path

def capslock_file(target_file):
    """Take given file and apply upper to it's whole content
    Args:
        target_file(str): file path
    Returns:
        file: same file as input, but with all it's content in upper case
    """
    with open(target_file, 'r') as f:
        content = f.read()
    content_upper = content.upper()
    with open(target_file, 'w') as f:
        f.write(content_upper)

def check_genome_integrity(genome):
    """Checks given genome file integrity
    Given certain genome, checks:
        composition: if it is composed only by A, T, C and G (further work shall include user guided filter to include also other standard code as R, Y, etc)
    Args:
        genome(str): path to genome file
    Returns:
        log(bool): discloses if the genome passed the composition integrity check 
    """
    bases = set('ATCG')
    with open(genome, 'r') as f:
        for line in f:
            if line.startswith(">"):
                continue
            if not set(line.strip()).issubset(bases):
                return False
    return True

def check_multiple_genomes_integrity(genomes):
    """Checks an array of genomes for file integrity
    Given certain genome, checks:
        composition: if it is composed only by A, T, C and G (further work shall include user guided filter to include also other standard code as R, Y, etc)
    Args:
        genome(dict): dict containing the root folder of genomes location as keys and organisms to process as values
    Returns:
        log(bool): discloses if the genome passed the composition integrity check 
    """
    app_organisms = []
    napp_organisms = []
    for path, organisms in genomes.items():
        path = Path(path)
        for organism in organisms:
            full_path = path / organism
            integrity = check_genome_integrity(full_path)
            if integrity:
                app_organisms.append(organism)
            else:
                napp_organisms.append(organism)
    return app_organisms, napp_organisms