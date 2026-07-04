"""Functions aimed at motifs analysis of nullomeric sequences from bit files"""

import struct


def calculate_cpg_index(index, k):
    """Calculate cpg occurrence, last and first base of each index sequence
    Args:
        index(int): index of sequence to be analyzed
        k(int): original size of the k-mer sequence
    Returns:
        cpg(int): total count of tuples consisting in C nucleotides directly followed by G nucleotides
    """
    cpg = 0
    for i in range(k - 1):
        shift = (k - i - 2) * 2
        pair = (index >> shift) & 15
        if pair == 6:
            cpg += 1
    return cpg


def generate_cpg_dict(k):
    """Generates a dict with cpg informations by sequence indexes
    Args:
        k(int): original k-mer size
    Returns:
        cpg_dict(dict): dict contaning, for each possible index from sequences sized k"
        "it's total cpg count."
    """
    cpg_dict = []
    for index in range(4**k):
        cpg = calculate_cpg_index(index, k)
        cpg_dict.append(cpg)
    return cpg_dict


def generate_complement_index_dict(half_k):
    """Generates a dict of complementary indexes
    Args:
        half_k(int): original k-mer sequence size
    Returns:
        complement_index_dict(dict): dictionary pairing indexes that represent complimentary k-mer sequences
    """
    complement = {0: 1, 1: 0, 2: 3, 3: 2}
    complement_index_dict = {}
    m = 4**half_k
    for number in range(m):
        temp = number
        bases = []
        for _ in range(half_k):
            bases.append(temp % 4)
            temp //= 4
        comp_bases = [complement[base] for base in bases]
        comp_index = 0
        for i, base in enumerate(comp_bases):
            comp_index += base * (4 ** (half_k - i - 1))
        complement_index_dict[number] = comp_index
    return complement_index_dict


def generate_homopolymer_array(half_k):
    """Generates an array containing all homopolymer indexes for given half_k
    Args:
        half_k(int): original size of k-mer sequences
    Returns:
        homopolymer_array(arr): all possible indexes that represent homopolymeric sequences
    """
    homopolymer_array = []
    for i in range(4):
        idx = 0
        for j in range(half_k):
            idx = (idx << 2) | i
        homopolymer_array.append(idx)
    return homopolymer_array
