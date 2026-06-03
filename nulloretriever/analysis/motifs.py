"""Functions aimed at motifs analysis of nullomeric sequences from bit files"""

import struct

def calculate_gpc_index(index, half_k):
    """Calculate cpg occurrence, last and first base of each index sequence
    Args:
        index(int): index of sequence to be analyzed
        half_k(int): original size of the k-mer sequence
    Returns:
        cpg(int): total count of tuples consisting in C nucleotides directly followed by G nucleotides in the index original sequence
        c(bool): marks if given sequence ends in a C nucleotide
        g(bool): marks if given sequence starts in a G nucleotide
    """
    cpg = 0
    c = False
    g = False
    for i in range(half_k):
        base = (index // (4 ** (half_k - i - 1))) % 4
        if i == 0 and base == 3:
            g = True
        if base == 2:
            c = True
        elif base == 3 and c:
            cpg += 1
            c = False
        else:
            c = False
    return cpg,c,g

def generate_cpg_dict(half_k):
    """Generates a dict with cpg informations by sequence indexes
    Args:
        half_k(int): original k-mer size
    Returns:
        cpg_dict(dict): dict contaning for each possible index for all k-mers with size half_k it's total cpg count, if it starts with G or ends in C 
    """
    cpg_dict = {}
    for index in range(4**half_k):
        cpg, c, g = calculate_gpc_index(index, half_k)
        cpg_dict[index] = (cpg,c,g)
    return cpg_dict

def retrieve_nullomers_cpg_stats(filename):
    """Retrieve multiple CpG stats for nullomers in a bit file
    Args:
        filename(str): location of bit file containing the nullomeric sequences
    Returns:
        cpg_stats(arr): array containing the following statistics:
            cpg_tot(int): total occurences of CpG dinucleotides
            null_with_cpg(int): total number of nullomers with at least one CpG dinucleotide occurrence
            cpg_global_mean(float): relation between null_with_cpg and total nullomer count
            cpg_count_mean(float): relation between cpg_tot and null_with_cpg
    Counts the percentage of nullomers that have at least 1 CpG dinucleotide.
    """
    count = 0
    with open(filename, 'rb') as f:
        # Skip header
        f.seek(6)  # Skip magic(4) + version(2)
        l_bytes = f.read(2)
        half_k = struct.unpack('<H', l_bytes)[0]
        byte_to_format = {1: 'B', 2: 'H', 4: 'I', 8: 'Q'}
        byte_size = struct.unpack('<B', f.read(1))[0]
        byte_format = byte_to_format[byte_size]
        cpg_dict = generate_cpg_dict(half_k)
        cpg_tot = 0
        null_with_cpg = 0
        try:
            while True:
                end_c = 0
                cpg_exists = 0
                index_bytes = f.read(byte_size)
                if len(index_bytes) < byte_size:
                    break
                v1 = struct.unpack(f'<{byte_format}', index_bytes)[0]
                if cpg_dict[v1][1] == 1:
                    end_c = 1
                if cpg_dict[v1][0] != 0:
                    cpg_exists = 1
                nullomer_count_bytes = f.read(byte_size)
                if len(nullomer_count_bytes) < byte_size:
                    break
                nullomer_count = struct.unpack(f'<{byte_format}', nullomer_count_bytes)[0]
                i = 0
                while i < nullomer_count:
                    end_g = 0
                    nullomer_byte = f.read(byte_size)
                    nullomer_index = struct.unpack(f'<{byte_format}', nullomer_byte)[0]
                    try:
                        cpg_tot += cpg_dict[v1][0] + cpg_dict[nullomer_index][0]
                        if cpg_dict[nullomer_index][2] == 1:
                            end_g = 1
                        if end_c and end_g:
                                cpg_tot +=1
                        if cpg_exists:
                            null_with_cpg += 1
                        elif end_c and end_g:
                            null_with_cpg += 1
                        else:
                            if cpg_dict[nullomer_index][0] > 0:
                                null_with_cpg += 1
                        i += 1
                    except:
                        print("DEBUG: NAO ACHOU CHAVE")
                        print(v1,nullomer_index,byte_format, nullomer_byte, half_k)
                count += nullomer_count     
        except (struct.error, OSError):
            pass
    try:
        cpg_count_mean = cpg_tot/null_with_cpg
    except ZeroDivisionError:
        cpg_count_mean = 0
    try:
        cpg_global_mean = (null_with_cpg/count)*100
    except ZeroDivisionError:
        cpg_global_mean = 0
    cpg_stats = {
        "total": cpg_tot,
        "nullomers_with_cpg": null_with_cpg,
        "global_mean": cpg_global_mean,
        "mean_nullomers_with_cpg": cpg_count_mean
    }
    return cpg_stats

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

def retrieve_palindrome_stats(filename):
    """Retrieve palindromic sequences statistics from a nullomer bit file
    Args:
        filename(str): path of nullomer bit file
    Returns:
        palindrome_count(int): total of found nullomers that are palindromic
        palindrome_relative(float): relation between palindromic nullomers and non palindromic nullomers found
    """
    palindrome_count = 0
    with open(filename, 'rb') as f:
        # Skip header
        f.seek(6)  # Skip magic(4) + version(2)
        l_bytes = f.read(2)
        half_k = struct.unpack('<H', l_bytes)[0]
        byte_to_format = {1: 'B', 2: 'H', 4: 'I', 8: 'Q'}
        byte_size = struct.unpack('<B', f.read(1))[0]
        byte_format = byte_to_format[byte_size]
        total_null = 0 
        complement_index_dict = generate_complement_index_dict(half_k)
        try:
            while True:
                index_bytes = f.read(byte_size)
                if len(index_bytes) < byte_size:
                    break
                v1 = struct.unpack(f'<{byte_format}', index_bytes)[0]
                v1_comp = complement_index_dict[v1]
                nullomer_count_bytes = f.read(byte_size)
                if len(nullomer_count_bytes) < byte_size:
                    break
                nullomer_count = struct.unpack(f'<{byte_format}', nullomer_count_bytes)[0]
                total_null += nullomer_count
                # Skip v2 indices (nullomer IDs)
                total_bytes = nullomer_count * byte_size
                i = 0
                while i < nullomer_count:
                    nullomer_byte = f.read(byte_size)
                    total_bytes -= byte_size
                    v2_index = struct.unpack(f'<{byte_format}', nullomer_byte)[0]
                    if v2_index == v1_comp:
                        palindrome_count += 1
                        f.seek(total_bytes,1)
                        break
                    i += 1                
        except (struct.error, OSError):
            pass
        try:
            palindrome_relative = (palindrome_count/total_null) * 100
        except ZeroDivisionError:
            palindrome_relative = 0
    palindrome_stats = {
        "count": palindrome_count,
        "relative_fraction": palindrome_relative
    }
    return palindrome_stats

def is_homopolymer(index,half_k):
    """Checks if given index is a homopolymer
    Args:
        index(int): index relative to a k-mer sequence
        half_k(int): original length of the k-mer sequence
    Returns:
        same(bool): boolean indicating if the index represents a homopolymeric sequence or not
    """
    same = False
    count = 0
    current = 0
    for i in range(half_k):
        base = (index // (4 ** (half_k - i - 1))) % 4
        if i == 0:
            current = base
        if base == current:
            count += 1
            continue
        else:
            same = False
            break
    if count == half_k:
        same = True
    return same

def generate_homopolymer_array(half_k):
    """Generates an array containing all homopolymer indexes for given half_k
    Args:
        half_k(int): original size of k-mer sequences
    Returns:
        homopolymer_array(arr): all possible indexes that represent homopolymeric sequences
    """
    homopolymer_array = []
    for i in range(4**half_k):
        homopolymer = is_homopolymer(i, half_k)
        if homopolymer:
            homopolymer_array.append(i)
    return homopolymer_array

def retrieve_homopolymer_stats(filename):
    """Retrieves which homopolymers can be found as nullomers of given organism
    Args:
        filename(str): path to nullomer data in binary
    Returns:
        found_homopolymers(arr): which homopolymers are present in the nullomer collection of given organism
    Obs:
        In this initial approach, a homopolymer is considered only if all bases of the sequence with length k are the same
    """
    byte_to_format = {1: 'B', 2: 'H', 4: 'I', 8: 'Q'}
    found_homopolymers = []
    with open(filename, 'rb') as f:
        # Skip header
        f.seek(6)  # Skip magic(4) + version(2)
        half_k = struct.unpack('<H', f.read(2))[0]
        half_k = int(half_k)
        byte_size = struct.unpack('<B', f.read(1))[0]
        byte_format = byte_to_format[byte_size]
        homopolymers = generate_homopolymer_array(half_k)
        try:
            while True:
                index_bytes = f.read(byte_size)
                if len(index_bytes) < byte_size:
                    break
                v1 = struct.unpack(f'<{byte_format}', index_bytes)[0]
                nullomer_count_bytes = f.read(byte_size)
                if len(nullomer_count_bytes) < byte_size:
                    break
                nullomer_count = struct.unpack(f'<{byte_format}', nullomer_count_bytes)[0]
                total_bytes = nullomer_count * byte_size
                i = 0
                if v1 in homopolymers:
                    while i < nullomer_count:
                        nullomer_byte = f.read(byte_size)
                        total_bytes -= byte_size
                        v2_index = struct.unpack(f'<{byte_format}', nullomer_byte)[0]
                        if v2_index == v1:
                            found_homopolymers.append(v1)
                            f.seek(total_bytes,1)
                            break
                        i += 1
                else:
                    f.seek(total_bytes,1)
        except (struct.error, OSError):
            pass
    return found_homopolymers
