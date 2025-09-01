"""Module that gather functions aimed at integrity checking and padronization"""

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
        padronization: makes sure all base data is in CapsLock
    Args:
        genome(str): path to genome file
    Returns:
        genome(file): same input file but padronized
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