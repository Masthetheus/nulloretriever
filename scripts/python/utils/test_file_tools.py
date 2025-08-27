"""Tools aiming to generate base genomes for testing
UNDER DEVELOPMENT:
    testing needed
"""

from itertools import product
import random

def random_genome(gsize, random_genome): 
    """Creates a random genome of size gsize
    Args:
        gsize(int): total genome lenght in bp
        random_genome(string): output path to random genome
    Returns:
        output: file containing the random genome in a single line
    """
    genome=""
    for count in range(gsize):
        genome+=random.choice("ATCG")
    with open(random_genome, "w") as datafile:
        for i in genome:
            datafile.write(str(i))
        datafile.close()

def generate_kmers(l):
    """Generate all possible combinations of nucleotides of size l
    Args:
        l(int): size of k-mers that shall be generated
    Returns:
        kmers[str]: array of all possible k-mers configurations of size l
    """
    bases = ['A', 'T', 'C', 'G']
    return [''.join(kmer) for kmer in product(bases, repeat=l)]