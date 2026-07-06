"""Functions aimed at motifs analysis of nullomeric sequences from bit files"""


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
        complement_index_dict(dict): dictionary pairing indexes that represent complementary k-mer sequences
    """
    m = 4**half_k
    mask = int(2 * half_k, 2)
    complement_index_dict = {}

    for number in range(m):
        comp = number ^ mask
        rev = 0
        temp = comp
        for _ in range(half_k):
            rev = (rev << 2) | (temp & 3)
            temp >>= 2
        complement_index_dict[number] = rev

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
