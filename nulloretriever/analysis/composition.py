"""Functions aimed at composition analysis of nullomeric sequences from bit files"""


def calculate_gc_index(index, half_k):
    """Calculate nullomeric gc content for each relative index
    Args:
        index(int): index related to a k-mer sequence
        half_k(int): size of the k-mer original sequence
    Returns:
        gc(int): number of occurences of G or C nucleotides on the sequence representend by the given index
    """
    gc = 0
    for i in range(half_k):
        gc += index >> (2 * i) & 1
    return gc


def generate_gc_dict(half_k):
    """Generate a dict pairing each possible index to it's total GC count
    Args:
        half_k(int): size of half k-mer from the original sequence
    Returns:
        gc_dict(dict): dictionary with pairs index:gc_count for all possible k-mers sequence of lenght half_k
    """
    gc_dict = {}
    for index in range(4**half_k):
        gc = calculate_gc_index(index, half_k)
        gc_dict[index] = gc
    return gc_dict
