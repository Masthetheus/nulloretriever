"""Functions aimed at composition analysis of nullomeric sequences from bit files"""


def calculate_gc_index(index, half_k):
    """Calculate nullomeric gc content for each relative index
    Args:
        index(int): index related to a k-mer sequence
        half_k(int): size of the k-mer original sequence
    Returns:
        gc(int): number of occurences of G or C nucleotides on the sequence representend by the given index
    """
    base_gc = 0
    v1_c = 0
    v1_g = 0
    for i in range(half_k):
        base_gc = index >> (2 * i) & 3
        v1_c += (base_gc == 1)
        v1_g += (base_gc == 3)
    return v1_c, v1_g


def generate_gc_dict(half_k):
    """Generate a dict pairing each possible index to it's total GC count
    Args:
        half_k(int): size of half k-mer from the original sequence
    Returns:
        gc_dict(dict): dictionary with pairs index:gc_count for all possible k-mers sequence of lenght half_k
    """
    gc_dict = {}
    for index in range(4**half_k):
        v1_c, v1_g = calculate_gc_index(index, half_k)
        gc_dict[index] = (v1_c, v1_g)
    return gc_dict

def generate_cg_lists(half_k):
    c_list = []
    g_list = []
    for i in range(4 ** half_k):
        c = g = 0
        for j in range(half_k):
            base = (i >> (2 * j)) & 3
            c += (base == 1)
            g += (base == 3)
        c_list.append(c)
        g_list.append(g)
    return c_list, g_list
